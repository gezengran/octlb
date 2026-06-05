#ifndef OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_ADAPTIVE_REFINEMENT_H_
#define OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_ADAPTIVE_REFINEMENT_H_

#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/geometry/geometry_config.h"
#include "src/mesh/geometry/triangle_soup.h"

namespace octlb {

BoundingBox extended_geometry_bbox(const TriangleSoup& soup,
                                   const GeometryConfig& config);

void resolve_bounding(OctreeForest& forest, const BoundingBox& geom_bbox,
                      const GeometryConfig& config);

void resolve_surface(OctreeForest& forest, const TriangleSoup& surface_soup,
                     const GeometryConfig& config);

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_ADAPTIVE_REFINEMENT_H_
