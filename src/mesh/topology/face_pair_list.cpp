#include "src/mesh/topology/face_pair_list.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
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
  enum Kind { kSameLevel, kCoarseFine, kCrossRankEdge } kind;
  std::size_t record_index;  // index into same_level_ / coarse_fine_ / cross_rank_edges_
  int slot;                  // 0 for same-level; 0..3 for coarse-fine; 0 for edge
};

struct FaceBuildContext {
  p8est_t* forest = nullptr;
  const OctreeForest* forest_obj = nullptr;  // for MeshForestAccess::QuadrantBounds
  p8est_ghost_t* ghost = nullptr;
  int my_rank = 0;
  int tag_ub = 0x3FFFFFFF;
  CommTagAssigner* tagger = nullptr;
  std::vector<PendingTagFill>* pending = nullptr;  // set by BuildFromForest
  std::unordered_set<uint64_t> cross_rank_face_keys;
  std::unordered_set<uint64_t> cross_rank_edge_keys;
  std::vector<SameLevelFace>* same_level = nullptr;
  std::vector<CoarseFineFace>* coarse_fine = nullptr;
  std::vector<TreeBoundaryFace>* tree_boundary = nullptr;
  std::vector<CrossRankEdge>* cross_rank_edges = nullptr;
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
    if (remote_rank != build->my_rank &&
        std::getenv("OCTLB_FACE_DEBUG") != nullptr) {
      const int loc_lvl = static_cast<int>(
          (local0 ? side0 : side1)->is.full.quad->level);
      const int rem_lvl = static_cast<int>(
          (local0 ? side1 : side0)->is.full.quad->level);
      std::fprintf(stderr,
                   "[sface r%d] peer=%d key=%llu locL=%d remL=%d loc=%d\n",
                   build->my_rank, remote_rank,
                   static_cast<unsigned long long>(face_key), loc_lvl, rem_lvl,
                   local0 ? 0 : 1);
      std::fflush(stderr);
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
  // Defect 5: snapshot the physical geometry of each side now, while the
  // iterate ghost layer holding the remote side is alive. The coarse side is
  // full; the 4 fine slots hang off it. coarse_side/fine_side quadrant pointers
  // are non-null here (ResolveFullSide/ResolveHangingSide returned true).
  entry.coarse_bounds = MeshForestAccess::QuadrantBounds(
      *build->forest_obj, coarse_side->treeid, coarse_side->is.full.quad);
  entry.coarse_level = static_cast<int>(coarse_side->is.full.quad->level);
  entry.fine_level = static_cast<int>(fine_side->is.hanging.quad[0]->level);
  for (int i = 0; i < 4; ++i) {
    entry.fine_ids[i] = fine_ids[i];
    entry.remote_ranks[i] = remote_ranks[i];
    entry.comm_tags[i] = 0;  // filled in phase two
    entry.fine_bounds[i] = MeshForestAccess::QuadrantBounds(
        *build->forest_obj, fine_side->treeid, fine_side->is.hanging.quad[i]);
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

// ② Resolve an edge side to (local_id | ghost quadid, rank). Same field
// layout as a face side (treeid + is.full.quad/quadid/is_ghost). A null quad
// means the edge-diagonal neighbour is absent from the ghost layer -- skip.
bool ResolveEdgeSide(const FaceBuildContext& ctx,
                     const p8est_iter_edge_side_t& side, OctantId* out_id,
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

// Canonical, symmetric key for an edge shared by two quadrants. The edge is
// identified by the two FACE AXES that meet at it (d/2), which are the same
// from both sides (an x-y edge is x-y regardless of Min/Max orientation), so
// both ranks compute the same key. The two quadrant keys are sorted so the
// key is symmetric -> both ranks of a cross-rank edge register the same key
// with CommTagAssigner and receive the same tag.
uint64_t CanonicalEdgeKey(uint64_t ka, uint64_t kb, int axis1, int axis2) {
  if (ka > kb) {
    std::swap(ka, kb);
  }
  if (axis1 > axis2) {
    std::swap(axis1, axis2);
  }
  uint64_t h = Mix64(ka);
  h = Mix64(h ^ kb);
  h = Mix64(h ^ (static_cast<uint64_t>(axis1) | (static_cast<uint64_t>(axis2)
                                                 << 8)));
  return h;
}

uint64_t SymmetricEdgeKey(const p8est_quadrant_t* qa, const p8est_quadrant_t* qb,
                          p4est_topidx_t tree_a, p4est_topidx_t tree_b,
                          int face_a0, int face_a1) {
  return CanonicalEdgeKey(QuadrantKey(tree_a, qa), QuadrantKey(tree_b, qb),
                          face_a0 / 2, face_a1 / 2);
}

// ② Stage B: enumerate cross-rank edge-diagonal neighbours. p8est fires this
// once per edge; a non-hanging same-level edge has exactly two sides (the two
// octants sharing the edge). We record only cross-rank edges (same-rank edges
// are filled by GhostSchedule Stage A via face-neighbour composition). Hanging
// edges (2:1 across the edge) are skipped -- cross-rank hanging-edge
// prolongation is a LevelCoupler concern, not same-level halo exchange.
void EdgeCallback(p8est_iter_edge_info_t* info, void* user_data) {
  auto* build = static_cast<FaceBuildContext*>(user_data);

  if (info->tree_boundary != 0) {
    return;  // tree-boundary edge: no in-domain diagonal neighbour
  }
  // An interior edge of a same-level forest is shared by 4 octants (one in each
  // quadrant around the edge). The 4 pair into 2 DIAGONAL pairs (octants whose
  // faces at the edge are opposite, i.e. they share ONLY the edge -- not a
  // face); the other 4 pairings share a face and are handled by face exchange.
  // For each diagonal pair with one local and one cross-rank octant, record a
  // directed cross-rank edge. Hanging edges (2:1) are deferred (LevelCoupler).
  const int n = static_cast<int>(info->sides.elem_count);
  if (n != 4) {
    return;  // only the uniform 4-octant same-level edge for now
  }

  struct SideInfo {
    p8est_iter_edge_side_t* side;
    OctantId id;
    int rank;
    int f0;
    int f1;
    bool valid;
  };
  SideInfo sides[4];
  for (int s = 0; s < 4; ++s) {
    sides[s].side = p8est_iter_eside_array_index_int(&info->sides, s);
    if (sides[s].side->is_hanging != 0) {
      return;  // hanging edge (2:1) deferred
    }
    sides[s].f0 = sides[s].side->faces[0];
    sides[s].f1 = sides[s].side->faces[1];
    sides[s].valid = ResolveEdgeSide(*build, *sides[s].side, &sides[s].id,
                                     &sides[s].rank);
  }

  // Pair diagonals: side j is the diagonal partner of side i iff its two faces
  // are the opposites (XOR 1) of side i's faces (in either order).
  for (int i = 0; i < 4; ++i) {
    if (!sides[i].valid) continue;
    const int of0 = sides[i].f0 ^ 1;
    const int of1 = sides[i].f1 ^ 1;
    for (int j = i + 1; j < 4; ++j) {
      if (!sides[j].valid) continue;
      const bool opposite =
          (sides[j].f0 == of0 && sides[j].f1 == of1) ||
          (sides[j].f0 == of1 && sides[j].f1 == of0);
      if (!opposite) continue;

      const bool i_local = sides[i].rank == build->my_rank;
      const bool j_local = sides[j].rank == build->my_rank;
      if (!i_local && !j_local) break;  // neither owned here
      if (i_local == j_local) break;    // same-rank pair: Stage A handles it

      const int li = i_local ? i : j;
      const int ri = i_local ? j : i;
      const int remote_rank = sides[ri].rank;
      const uint64_t edge_key = SymmetricEdgeKey(
          sides[i].side->is.full.quad, sides[j].side->is.full.quad,
          sides[i].side->treeid, sides[j].side->treeid, sides[li].f0,
          sides[li].f1);
      if (!build->cross_rank_edge_keys.insert(edge_key).second) break;

      const FaceDir d1 = static_cast<FaceDir>(sides[li].f0);
      const FaceDir d2 = static_cast<FaceDir>(sides[li].f1);
      build->cross_rank_edges->push_back(
          {sides[li].id, d1, d2, sides[ri].id, remote_rank, /*comm_tag=*/0});
      RegisterFace(build, edge_key, remote_rank, PendingTagFill::kCrossRankEdge,
                   build->cross_rank_edges->size() - 1, /*slot=*/0);
      break;
    }
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
  ctx->forest_obj = &forest;
  if (ctx->forest == nullptr) {
    throw std::runtime_error("FacePairList: null forest");
  }
  // The forest's face-connected ghost (built by partition()) is the contract
  // that partition() has run; the face callback uses it so that SameLevelFace /
  // CoarseFineFace remote_id remains an index into THIS face-ghost (consumers
  // like GhostSchedule and test_ghost_topology look it up there).
  p8est_ghost_t* face_ghost = MeshForestAccess::Ghost(forest);
  if (face_ghost == nullptr) {
    throw std::runtime_error(
        "FacePairList: ghost layer missing; call partition() first");
  }
  MPI_Comm_rank(ctx->forest->mpicomm, &ctx->my_rank);
  ctx->tag_ub = QueryTagUb(ctx->forest->mpicomm);

  CommTagAssigner tagger(ctx->tag_ub);
  ctx->tagger = &tagger;
  ctx->pending = pending;

  // Pass 1: face pairs with the face-connected ghost (remote_id indexes it).
  ctx->ghost = face_ghost;
  p8est_iterate(ctx->forest, face_ghost, ctx, nullptr, FaceCallback, nullptr,
                nullptr);

  // ② Stage B pass 2: edge pairs with a local EDGE-connected ghost. The
  // face-connected ghost lacks edge-diagonal neighbours; the edge ghost
  // contains them. CrossRankEdge.remote_id indexes the edge ghost, but no
  // consumer looks it up in the face ghost (GhostSchedule Stage B exchanges by
  // rank+tag, not by ghost index), so the two ghost index spaces are kept
  // separate and consistent. Face + edge keys share one CommTagAssigner so a
  // face and an edge to the same peer can't collide.
  p8est_ghost_t* edge_ghost =
      p8est_ghost_new(ctx->forest, P8EST_CONNECT_EDGE);
  if (edge_ghost == nullptr) {
    throw std::runtime_error("FacePairList: edge ghost build failed");
  }
  ctx->ghost = edge_ghost;
  p8est_iterate(ctx->forest, edge_ghost, ctx, nullptr, nullptr, EdgeCallback,
                nullptr);
  p8est_ghost_destroy(edge_ghost);
  ctx->ghost = nullptr;

  ctx->tagger = nullptr;
  ctx->pending = nullptr;
  ctx->forest_obj = nullptr;

  tagger.assign();
  for (const PendingTagFill& p : *pending) {
    const int tag = tagger.tag_for(p.face_key);
    if (p.kind == PendingTagFill::kSameLevel) {
      (*ctx->same_level)[p.record_index].comm_tag = tag;
    } else if (p.kind == PendingTagFill::kCoarseFine) {
      (*ctx->coarse_fine)[p.record_index].comm_tags[p.slot] = tag;
    } else {  // kCrossRankEdge
      (*ctx->cross_rank_edges)[p.record_index].comm_tag = tag;
    }
  }

  ctx->forest = nullptr;
}

}  // namespace

FacePairList::FacePairList(const OctreeForest& forest) {
  FaceBuildContext ctx;
  ctx.same_level = &same_level_;
  ctx.coarse_fine = &coarse_fine_;
  ctx.tree_boundary = &tree_boundary_;
  ctx.cross_rank_edges = &cross_rank_edges_;
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

const std::vector<CrossRankEdge>& FacePairList::cross_rank_edges() const {
  return cross_rank_edges_;
}

void FacePairList::rebuild(const OctreeForest& forest) {
  same_level_.clear();
  coarse_fine_.clear();
  tree_boundary_.clear();
  cross_rank_edges_.clear();
  FaceBuildContext ctx;
  ctx.same_level = &same_level_;
  ctx.coarse_fine = &coarse_fine_;
  ctx.tree_boundary = &tree_boundary_;
  ctx.cross_rank_edges = &cross_rank_edges_;
  std::vector<PendingTagFill> pending;
  BuildFromForest(forest, &ctx, &pending);
}

}  // namespace octlb