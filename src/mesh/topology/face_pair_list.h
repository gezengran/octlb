#ifndef OCTLB_SRC_MESH_TOPOLOGY_FACE_PAIR_LIST_H_
#define OCTLB_SRC_MESH_TOPOLOGY_FACE_PAIR_LIST_H_

#include <vector>

#include "src/common/types.h"

namespace octlb {

class OctreeForest;

struct SameLevelFace {
  OctantId local_id;
  FaceDir dir;
  OctantId remote_id;
  int remote_rank;
};

struct CoarseFineFace {
  OctantId coarse_id;
  OctantId fine_ids[4];
  FaceDir normal;
  int remote_ranks[4];
};

class FacePairList {
 public:
  explicit FacePairList(const OctreeForest& forest);

  const std::vector<SameLevelFace>& same_level_faces() const;
  const std::vector<CoarseFineFace>& coarse_fine_faces() const;

  void rebuild(const OctreeForest& forest);

 private:
  std::vector<SameLevelFace> same_level_;
  std::vector<CoarseFineFace> coarse_fine_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_TOPOLOGY_FACE_PAIR_LIST_H_
