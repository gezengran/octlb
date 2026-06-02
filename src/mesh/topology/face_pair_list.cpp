#include "src/mesh/topology/face_pair_list.h"

#include <stdexcept>

#include <p8est_bits.h>
#include <p8est_iterate.h>

#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/forest/octree_forest_access.h"

namespace octlb {
namespace {

struct FaceBuildContext {
  p8est_t* forest = nullptr;
  p8est_ghost_t* ghost = nullptr;
  int my_rank = 0;
  std::vector<SameLevelFace>* same_level = nullptr;
  std::vector<CoarseFineFace>* coarse_fine = nullptr;
};

int GhostOwnerRank(const p8est_ghost_t* ghost, p4est_locidx_t ghost_index) {
  for (int rank = 0; rank < ghost->mpisize; ++rank) {
    const p4est_locidx_t begin = ghost->proc_offsets[rank];
    const p4est_locidx_t end = ghost->proc_offsets[rank + 1];
    if (ghost_index >= begin && ghost_index < end) {
      return rank;
    }
  }
  return -1;
}

OctantId LocalOctantIndex(const p8est_t* forest, p4est_topidx_t treeid,
                          const p8est_quadrant_t* quadrant) {
  p8est_tree_t* tree = p8est_tree_array_index(forest->trees, treeid);
  p8est_quadrant_t* first =
      p8est_quadrant_array_index(&tree->quadrants, 0);
  return static_cast<OctantId>(quadrant - first);
}

bool ResolveFullSide(const FaceBuildContext& ctx,
                     const p8est_iter_face_side_t& side, OctantId* out_id,
                     int* out_rank) {
  if (side.is_hanging != 0) {
    return false;
  }
  if (side.is.full.quad == nullptr) {
    return false;
  }
  if (side.is.full.is_ghost) {
    *out_id = static_cast<OctantId>(side.is.full.quadid);
    *out_rank = GhostOwnerRank(ctx.ghost, side.is.full.quadid);
    return *out_rank >= 0;
  }
  *out_id = LocalOctantIndex(ctx.forest, side.treeid, side.is.full.quad);
  *out_rank = ctx.my_rank;
  return true;
}

bool ResolveHangingSide(const FaceBuildContext& ctx,
                        const p8est_iter_face_side_t& side,
                        OctantId fine_ids[4], int remote_ranks[4]) {
  if (side.is_hanging == 0) {
    return false;
  }
  for (int i = 0; i < 4; ++i) {
    if (side.is.hanging.quad[i] == nullptr) {
      return false;
    }
    if (side.is.hanging.is_ghost[i]) {
      fine_ids[i] = static_cast<OctantId>(side.is.hanging.quadid[i]);
      remote_ranks[i] = GhostOwnerRank(ctx.ghost, side.is.hanging.quadid[i]);
      if (remote_ranks[i] < 0) {
        return false;
      }
    } else {
      fine_ids[i] =
          LocalOctantIndex(ctx.forest, side.treeid, side.is.hanging.quad[i]);
      remote_ranks[i] = ctx.my_rank;
    }
  }
  return true;
}

bool SideHasLocalQuadrant(const FaceBuildContext& ctx,
                          const p8est_iter_face_side_t& side) {
  if (side.is_hanging != 0) {
    for (int i = 0; i < 4; ++i) {
      if (side.is.hanging.quad[i] != nullptr &&
          side.is.hanging.is_ghost[i] == 0) {
        return true;
      }
    }
    return false;
  }
  return side.is.full.quad != nullptr && side.is.full.is_ghost == 0;
}

void FaceCallback(p8est_iter_face_info_t* info, void* user_data) {
  auto* ctx = static_cast<FaceBuildContext*>(user_data);
  if (info->tree_boundary != 0) {
    return;
  }
  if (info->sides.elem_count != 2) {
    return;
  }

  p8est_iter_face_side_t* side0 =
      p8est_iter_fside_array_index_int(&info->sides, 0);
  p8est_iter_face_side_t* side1 =
      p8est_iter_fside_array_index_int(&info->sides, 1);

  const bool local0 = SideHasLocalQuadrant(*ctx, *side0);
  const bool local1 = SideHasLocalQuadrant(*ctx, *side1);
  if (!local0 && !local1) {
    return;
  }

  if (side0->is_hanging == side1->is_hanging) {
    if (side0->is_hanging != 0) {
      return;
    }
    OctantId id0 = 0;
    OctantId id1 = 0;
    int rank0 = -1;
    int rank1 = -1;
    if (!ResolveFullSide(*ctx, *side0, &id0, &rank0) ||
        !ResolveFullSide(*ctx, *side1, &id1, &rank1)) {
      return;
    }
    const p8est_iter_face_side_t* local_side = local0 ? side0 : side1;
    const p8est_iter_face_side_t* remote_side = local0 ? side1 : side0;
    OctantId local_id = local0 ? id0 : id1;
    OctantId remote_id = local0 ? id1 : id0;
    int remote_rank = local0 ? rank1 : rank0;
    if (!local0) {
      local_id = id1;
      remote_id = id0;
      remote_rank = rank0;
    }
    (void)remote_side;
    ctx->same_level->push_back(
        {local_id, static_cast<FaceDir>(local_side->face), remote_id,
         remote_rank});
    return;
  }

  const p8est_iter_face_side_t* coarse_side =
      side0->is_hanging == 0 ? side0 : side1;
  const p8est_iter_face_side_t* fine_side =
      side0->is_hanging != 0 ? side0 : side1;

  OctantId coarse_id = 0;
  int coarse_rank = -1;
  OctantId fine_ids[4] = {};
  int remote_ranks[4] = {};
  if (!ResolveFullSide(*ctx, *coarse_side, &coarse_id, &coarse_rank) ||
      !ResolveHangingSide(*ctx, *fine_side, fine_ids, remote_ranks)) {
    return;
  }
  if (!SideHasLocalQuadrant(*ctx, *coarse_side) &&
      !SideHasLocalQuadrant(*ctx, *fine_side)) {
    return;
  }

  CoarseFineFace entry{};
  entry.coarse_id = coarse_id;
  entry.normal = static_cast<FaceDir>(coarse_side->face);
  for (int i = 0; i < 4; ++i) {
    entry.fine_ids[i] = fine_ids[i];
    entry.remote_ranks[i] = remote_ranks[i];
  }
  (void)coarse_rank;
  ctx->coarse_fine->push_back(entry);
}

void BuildFromForest(const OctreeForest& forest, FaceBuildContext* ctx) {
  ctx->forest = MeshForestAccess::Forest(forest);
  ctx->ghost = MeshForestAccess::Ghost(forest);
  if (ctx->forest == nullptr) {
    throw std::runtime_error("FacePairList: null forest");
  }
  if (ctx->ghost == nullptr) {
    throw std::runtime_error(
        "FacePairList: ghost layer missing; call partition() first");
  }
  MPI_Comm_rank(ctx->forest->mpicomm, &ctx->my_rank);
  p8est_iterate(ctx->forest, ctx->ghost, ctx, nullptr, FaceCallback, nullptr,
                nullptr);
}

}  // namespace

FacePairList::FacePairList(const OctreeForest& forest) {
  FaceBuildContext ctx;
  ctx.same_level = &same_level_;
  ctx.coarse_fine = &coarse_fine_;
  BuildFromForest(forest, &ctx);
}

const std::vector<SameLevelFace>& FacePairList::same_level_faces() const {
  return same_level_;
}

const std::vector<CoarseFineFace>& FacePairList::coarse_fine_faces() const {
  return coarse_fine_;
}

void FacePairList::rebuild(const OctreeForest& forest) {
  same_level_.clear();
  coarse_fine_.clear();
  FaceBuildContext ctx;
  ctx.same_level = &same_level_;
  ctx.coarse_fine = &coarse_fine_;
  BuildFromForest(forest, &ctx);
}

}  // namespace octlb
