#include "src/mesh/topology/face_pair_list.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include <mpi.h>
#include <p8est_bits.h>
#include <p8est_iterate.h>

#include "src/common/comm_tag_assigner.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/forest/octree_forest_access.h"

namespace octlb {
namespace {

// A pending tag fill: a face whose comm_tag must be written after the assigner
// runs. Tags are assigned in a second phase (after p8est_iterate) so that, per
// rank-pair, all faces are known and a symmetric sorted-index tag can be given.
struct PendingTagFill {
  uint64_t face_key;
  enum Kind { kSameLevel, kCoarseFine } kind;
  std::size_t record_index;  // index into same_level_ or coarse_fine_
  int slot;                  // 0 for same-level; 0..3 for coarse-fine
};

struct FaceBuildContext {
  p8est_t* forest = nullptr;
  p8est_ghost_t* ghost = nullptr;
  int my_rank = 0;
  int tag_ub = 0x3FFFFFFF;
  CommTagAssigner* tagger = nullptr;
  std::vector<PendingTagFill>* pending = nullptr;  // set by BuildFromForest
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

// Canonical, symmetric key for a face shared by two quadrants. Both endpoints
// of a cross-rank face compute the same value (operands are sorted), so the
// key is usable as the per-rank-pair face identity for tag assignment.
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

uint64_t SymmetricFaceKey(const p8est_quadrant_t* qa,
                          const p8est_quadrant_t* qb, p4est_topidx_t tree_a,
                          p4est_topidx_t tree_b, int face) {
  return CanonicalFaceKey(QuadrantKey(tree_a, qa), QuadrantKey(tree_b, qb),
                          face);
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

// Register a face (key, peer) for tag assignment and record where its tag must
// be written back in phase two.
void RegisterFace(FaceBuildContext* ctx, uint64_t face_key, int peer_rank,
                  PendingTagFill::Kind kind, std::size_t record_index,
                  int slot) {
  ctx->tagger->add_face(face_key, peer_rank);
  ctx->pending->push_back({face_key, kind, record_index, slot});
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
    const uint64_t face_key = SymmetricFaceKey(
        side0->is.full.quad, side1->is.full.quad, side0->treeid,
        side1->treeid, side0->face);
    if (remote_rank != build->my_rank &&
        !build->cross_rank_face_keys.insert(face_key).second) {
      return;
    }
    const FaceDir local_dir = static_cast<FaceDir>(local_side->face);
    const FaceDir remote_dir = static_cast<FaceDir>(remote_side->face);
    const int peer_rank =
        remote_rank != build->my_rank ? remote_rank : build->my_rank;

    build->same_level->push_back(
        {local_id, local_dir, remote_id, remote_rank, /*comm_tag=*/0});
    RegisterFace(build, face_key, peer_rank, PendingTagFill::kSameLevel,
                 build->same_level->size() - 1, /*slot=*/0);
    if (remote_rank == build->my_rank) {
      build->same_level->push_back(
          {remote_id, remote_dir, local_id, remote_rank, /*comm_tag=*/0});
      RegisterFace(build, face_key, peer_rank, PendingTagFill::kSameLevel,
                   build->same_level->size() - 1, /*slot=*/0);
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
    entry.comm_tags[i] = 0;  // filled in phase two
  }
  build->coarse_fine->push_back(entry);
  const std::size_t entry_index = build->coarse_fine->size() - 1;

  const bool owns_coarse = SideHasLocalQuadrant(*build, *coarse_side);
  for (int i = 0; i < 4; ++i) {
    const uint64_t face_key_i = SymmetricFaceKey(
        coarse_side->is.full.quad, fine_side->is.hanging.quad[i],
        coarse_side->treeid, fine_side->treeid, coarse_side->face);
    int peer_rank = -1;
    if (owns_coarse) {
      // Coarse side participates with every fine slot (peer = fine_i's rank).
      peer_rank = remote_ranks[i];
    } else if (fine_side->is.hanging.is_ghost[i] == 0) {
      // Fine side: only the slots this rank owns (peer = coarse's rank).
      peer_rank = coarse_rank;
    }
    if (peer_rank < 0) {
      continue;  // not owned by this rank -> tag unused here, leave 0
    }
    RegisterFace(build, face_key_i, peer_rank, PendingTagFill::kCoarseFine,
                 entry_index, i);
  }
}

int QueryTagUb(MPI_Comm comm) {
  int tag_ub = 0x3FFFFFFF;
  int exists = 0;
  int* value = nullptr;
  if (MPI_Comm_get_attr(comm, MPI_TAG_UB, &value, &exists) == MPI_SUCCESS &&
      exists && value != nullptr && *value > 0) {
    tag_ub = *value;
  }
  return tag_ub;
}

void BuildFromForest(const OctreeForest& forest, FaceBuildContext* ctx,
                     std::vector<PendingTagFill>* pending) {
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
  ctx->tag_ub = QueryTagUb(ctx->forest->mpicomm);

  CommTagAssigner tagger(ctx->tag_ub);
  ctx->tagger = &tagger;
  ctx->pending = pending;
  p8est_iterate(ctx->forest, ctx->ghost, ctx, nullptr, FaceCallback, nullptr,
                nullptr);
  ctx->tagger = nullptr;
  ctx->pending = nullptr;

  tagger.assign();
  for (const PendingTagFill& p : *pending) {
    const int tag = tagger.tag_for(p.face_key);
    if (p.kind == PendingTagFill::kSameLevel) {
      (*ctx->same_level)[p.record_index].comm_tag = tag;
    } else {
      (*ctx->coarse_fine)[p.record_index].comm_tags[p.slot] = tag;
    }
  }
}

}  // namespace

FacePairList::FacePairList(const OctreeForest& forest) {
  FaceBuildContext ctx;
  ctx.same_level = &same_level_;
  ctx.coarse_fine = &coarse_fine_;
  ctx.tree_boundary = &tree_boundary_;
  std::vector<PendingTagFill> pending;
  BuildFromForest(forest, &ctx, &pending);
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
  std::vector<PendingTagFill> pending;
  BuildFromForest(forest, &ctx, &pending);
}

}  // namespace octlb