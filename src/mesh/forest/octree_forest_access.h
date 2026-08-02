#ifndef OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_ACCESS_H_
#define OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_ACCESS_H_

/** Mesh-module-only accessors; not for solver/. */

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest_internal.h"

namespace octlb {

class OctreeForest;

struct MeshForestAccess {
  static p8est_t* Forest(const OctreeForest& forest);
  static p8est_ghost_t* Ghost(const OctreeForest& forest);
  // Physical bounds of any quadrant (local OR ghost) identified by its treeid
  // + quadrant pointer. Works on ghost quadrants captured during p8est_iterate,
  // whose quadid is a transient ghost-array index and cannot be resolved later
  // via OctreeForest::quadrant_bounds(OctantId). Used by FacePairList to capture
  // coarse/fine side bounds at iterate time (defect 5: cross-rank coarse-fine).
  static BoundingBox QuadrantBounds(const OctreeForest& forest,
                                    p4est_topidx_t treeid,
                                    const p8est_quadrant_t* q);
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_ACCESS_H_
