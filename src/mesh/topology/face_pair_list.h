#ifndef OCTLB_SRC_MESH_TOPOLOGY_FACE_PAIR_LIST_H_
#define OCTLB_SRC_MESH_TOPOLOGY_FACE_PAIR_LIST_H_

#include <vector>

#include "src/common/bounding_box.h"
#include "src/common/types.h"

namespace octlb {

class OctreeForest;

struct SameLevelFace {
  OctantId local_id;
  FaceDir dir;
  OctantId remote_id;
  int remote_rank;
  int comm_tag;  // Symmetric MPI tag for cross-rank halo exchange (T05).
};

struct CoarseFineFace {
  OctantId coarse_id;
  OctantId fine_ids[4];
  FaceDir normal;
  int coarse_remote_rank;
  int remote_ranks[4];
  int comm_tags[4];  // Symmetric MPI tag per fine slot (T06).
  // Defect 5 (cross-rank coarse-fine): physical bounds + level of each side,
  // captured at FacePairList build time (during p8est_iterate, when the ghost
  // layer holding the remote side is alive). The remote side's quadid is a
  // transient ghost-array index invalid after construction, so the geometry
  // is snapshotted here for LevelCoupler -- it never re-resolves a quadid.
  BoundingBox coarse_bounds;
  int coarse_level;
  BoundingBox fine_bounds[4];
  int fine_level;  // all 4 fine slots share one level (2:1 balance)
};

/** Domain outer face: octant touches the physical tree boundary (T09-W1). */
struct TreeBoundaryFace {
  OctantId octant_id;
  FaceDir face_dir;
};

// ② cross-rank edge-diagonal neighbour (T11 W3 Stage B): octant `local_id`
// shares an edge (intersection of faces d1, d2) with `remote_id` on
// `remote_rank`. Same-rank edge pairs are NOT here -- GhostSchedule Stage A
// fills those by composing same-rank face neighbours; this list is cross-rank
// only (remote_rank != owner rank). comm_tag is symmetric across the two ranks
// (CommTagAssigner over a canonical edge key), shared with the face tags so a
// face and an edge to the same peer can't collide. Pure topology: no physics.
struct CrossRankEdge {
  OctantId local_id;
  FaceDir d1;
  FaceDir d2;
  OctantId remote_id;
  int remote_rank;
  int comm_tag;
};

class FacePairList {
 public:
  explicit FacePairList(const OctreeForest& forest);

  const std::vector<SameLevelFace>& same_level_faces() const;
  const std::vector<CoarseFineFace>& coarse_fine_faces() const;
  const std::vector<TreeBoundaryFace>& tree_boundary_faces() const;
  const std::vector<CrossRankEdge>& cross_rank_edges() const;

  void rebuild(const OctreeForest& forest);

 private:
  std::vector<SameLevelFace> same_level_;
  std::vector<CoarseFineFace> coarse_fine_;
  std::vector<TreeBoundaryFace> tree_boundary_;
  std::vector<CrossRankEdge> cross_rank_edges_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_TOPOLOGY_FACE_PAIR_LIST_H_
