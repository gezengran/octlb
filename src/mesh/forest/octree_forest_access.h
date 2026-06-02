#ifndef OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_ACCESS_H_
#define OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_ACCESS_H_

/** Mesh-module-only accessors; not for solver/. */

#include "src/mesh/forest/octree_forest_internal.h"

namespace octlb {

class OctreeForest;

struct MeshForestAccess {
  static p8est_t* Forest(const OctreeForest& forest);
  static p8est_ghost_t* Ghost(const OctreeForest& forest);
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_ACCESS_H_
