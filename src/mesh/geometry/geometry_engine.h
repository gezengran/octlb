#ifndef OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_ENGINE_H_
#define OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_ENGINE_H_

#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/geometry/geometry_assembly.h"
#include "src/mesh/geometry/geometry_config.h"
#include "src/mesh/geometry/material_field.h"

namespace octlb {

/**
 * Geometry-driven static AMR and voxel material tagging.
 *
 * Modifies \p forest (refine/balance/partition) and returns a local
 * MaterialField. Does not construct FacePairList; after build() the caller
 * must rebuild FacePairList, GhostSchedule, LevelCoupler, and TimeLoop caches.
 */
class GeometryEngine {
 public:
  MaterialField build(OctreeForest& forest, const GeometryAssembly& assembly,
                      const GeometryConfig& config) const;
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_ENGINE_H_
