#include "src/mesh/topology/face_pair_list.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

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
  std::unordered_map<uint64_t, int> comm_tag_by_face_key;
  std::unordered_map<int, uint64_t> tag_owner;
  std::unordered_set<uint64_t> cross_rank_face_keys;
  std::vector<SameLevelFace>* same_level = nullptr;
  std::vector<CoarseFineFace>* coarse_fine = nullptr;
  std::vector<TreeBoundaryFace>* tree_boundary = nullptr;
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

uint64_t QuadrantKey(p4est_topidx_t treeid, const p8est_quadrant_t* q) {
  return (static_cast<uint64_t>(treeid) << 48) |
         (static_cast<uint64_t>(q->level) << 40) |
         (static_cast<uint64_t>(q->x) << 30) |
         (static_cast<uint64_t>(q->y) << 20) |
         (static_cast<uint64_t>(q->z) << 10);
}

uint64_t Mix64(uint64_t h) {
  h ^= h >> 33;
  h *= 0xff51afd7ed558ccdULL;
  h ^= h >> 33;
  h *= 0xc4ceb9fe1a85ec53ULL;
  h ^= h >> 33;
  return h;
}

uint64_t CanonicalFaceKey(uint64_t ka, uint64_t kb, int face) {
  if (ka > kb) {
    std::swap(ka, kb);
    face ^= 1;
  }
  uint64_t h = ka;
  h = Mix64(h ^ kb);
  h = Mix64(h ^ static_cast<uint64_t>(face));
  return h;
}

int AssignCommTag(FaceBuildContext* ctx, uint64_t face_key) {
  const auto existing = ctx->comm_tag_by_face_key.find(face_key);
  if (existing != ctx->comm_tag_by_face_key.end()) {
    return existing->second;
  }

  int tag = static_cast<int>(Mix64(face_key) & 0x7FFFFFFF);
  if (tag <= 0) {
    tag = 1;
  }
  uint64_t salt = face_key;
  while (ctx->tag_owner.count(tag) != 0 &&
         ctx->tag_owner[tag] != face_key) {
    salt = Mix64(salt ^ static_cast<uint64_t>(tag));
    tag = static_cast<int>(salt & 0x7FFFFFFF);
    if (tag <= 0) {
      tag = 1;
    }
  }
  ctx->tag_owner[tag] = face_key;
  ctx->comm_tag_by_face_key.emplace(face_key, tag);
  return tag;
}

uint64_t SymmetricFaceKey(const p8est_quadrant_t* qa,
                          const p8est_quadrant_t* qb, p4est_topidx_t tree_a,
                          p4est_topidx_t tree_b, int face) {
  return CanonicalFaceKey(QuadrantKey(tree_a, qa), QuadrantKey(tree_b, qb),
                          face);
}

int MakeSymmetricCommTag(FaceBuildContext* ctx, const p8est_quadrant_t* qa,
                         const p8est_quadrant_t* qb, p4est_topidx_t tree_a,
                         p4est_topidx_t tree_b, int face) {
  if (qa == nullptr || qb == nullptr) {
    return -1;
  }
  return AssignCommTag(ctx,
                       SymmetricFaceKey(qa, qb, tree_a, tree_b, face));
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
  auto* build = static_cast<FaceBuildContext*>(user_data);

  if (info->tree_boundary != 0) {
    if (build->tree_boundary == nullptr || info->sides.elem_count == 0) {
      return;
    }
    for (p4est_locidx_t si = 0; si < info->sides.elem_count; ++si) {
      p8est_iter_face_side_t* side =
          p8est_iter_fside_array_index_int(&info->sides, si);
      if (!SideHasLocalQuadrant(*build, *side)) {
        continue;
      }
      OctantId local_id = 0;
      int rank = -1;
      if (!ResolveFullSide(*build, *side, &local_id, &rank)) {
        continue;
      }
      if (rank != build->my_rank) {
        continue;
      }
      build->tree_boundary->push_back(
          {local_id, static_cast<FaceDir>(side->face)});
    }
    return;
  }

  if (info->sides.elem_count != 2) {
    return;
  }

  p8est_iter_face_side_t* side0 =
      p8est_iter_fside_array_index_int(&info->sides, 0);
  p8est_iter_face_side_t* side1 =
      p8est_iter_fside_array_index_int(&info->sides, 1);

  const bool local0 = SideHasLocalQuadrant(*build, *side0);
  const bool local1 = SideHasLocalQuadrant(*build, *side1);
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
    if (!ResolveFullSide(*build, *side0, &id0, &rank0) ||
        !ResolveFullSide(*build, *side1, &id1, &rank1)) {
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
    const int comm_tag = MakeSymmetricCommTag(
        build, side0->is.full.quad, side1->is.full.quad, side0->treeid,
        side1->treeid, side0->face);
    if (comm_tag < 0) {
      return;
    }
    const uint64_t face_key = SymmetricFaceKey(
        side0->is.full.quad, side1->is.full.quad, side0->treeid,
        side1->treeid, side0->face);
    if (remote_rank != build->my_rank &&
        !build->cross_rank_face_keys.insert(face_key).second) {
      return;
    }
    const FaceDir local_dir = static_cast<FaceDir>(local_side->face);
    const FaceDir remote_dir = static_cast<FaceDir>(remote_side->face);
    build->same_level->push_back(
        {local_id, local_dir, remote_id, remote_rank, comm_tag});
    if (remote_rank == build->my_rank) {
      build->same_level->push_back(
          {remote_id, remote_dir, local_id, remote_rank, comm_tag});
    }
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
  if (!ResolveFullSide(*build, *coarse_side, &coarse_id, &coarse_rank) ||
      !ResolveHangingSide(*build, *fine_side, fine_ids, remote_ranks)) {
    return;
  }
  if (!SideHasLocalQuadrant(*build, *coarse_side) &&
      !SideHasLocalQuadrant(*build, *fine_side)) {
    return;
  }

  CoarseFineFace entry{};
  entry.coarse_id = coarse_id;
  entry.normal = static_cast<FaceDir>(coarse_side->face);
  entry.coarse_remote_rank = coarse_rank;
  for (int i = 0; i < 4; ++i) {
    entry.fine_ids[i] = fine_ids[i];
    entry.remote_ranks[i] = remote_ranks[i];
    entry.comm_tags[i] = MakeSymmetricCommTag(
        build, coarse_side->is.full.quad, fine_side->is.hanging.quad[i],
        coarse_side->treeid, fine_side->treeid, coarse_side->face);
  }
  build->coarse_fine->push_back(entry);
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
  ctx.tree_boundary = &tree_boundary_;
  BuildFromForest(forest, &ctx);
}

const std::vector<SameLevelFace>& FacePairList::same_level_faces() const {
  return same_level_;
}

const std::vector<CoarseFineFace>& FacePairList::coarse_fine_faces() const {
  return coarse_fine_;
}

const std::vector<TreeBoundaryFace>& FacePairList::tree_boundary_faces() const {
  return tree_boundary_;
}

void FacePairList::rebuild(const OctreeForest& forest) {
  same_level_.clear();
  coarse_fine_.clear();
  tree_boundary_.clear();
  FaceBuildContext ctx;
  ctx.same_level = &same_level_;
  ctx.coarse_fine = &coarse_fine_;
  ctx.tree_boundary = &tree_boundary_;
  BuildFromForest(forest, &ctx);
}

}  // namespace octlb
