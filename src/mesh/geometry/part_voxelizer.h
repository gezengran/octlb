#ifndef OCTLB_SRC_MESH_GEOMETRY_PART_VOXELIZER_H_
#define OCTLB_SRC_MESH_GEOMETRY_PART_VOXELIZER_H_

#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/geometry/cgal_surface_mesh.h"
#include "src/mesh/geometry/geometry_config.h"
#include "src/mesh/geometry/geometry_types.h"
#include "src/mesh/geometry/material_field.h"

namespace octlb {

void validate_geometry_part(const GeometryPart& part);

MaterialField voxelize_part(const OctreeForest& forest, const GeometryPart& part,
                            const GeometryConfig& config);

/** Merges \p overlay into \p base (used after sorting parts by ascending priority).
 *
 * The merge rule is role-dependent:
 *  - kInternalChannel DEFINES the fluid domain: its kFluid (lumen) carves fluid
 *    out of whatever sits below (e.g. a lower-priority obstacle), so it
 *    overwrites UNCONDITIONALLY -- every cell takes the channel's material.
 *  - kExternalObstacle ADDS an obstacle into the existing fluid: its kSolid
 *    (interior) and kBoundary (surface) claim those cells, but its kFluid means
 *    "outside my obstacle" -- it must NOT erase a lower-priority channel's wall
 *    (kBoundary) or earth (kSolid). So it overwrites ONLY where it claims
 *    kSolid/kBoundary (claim-based). The unconditional rule would let the
 *    obstacle's kFluid erase the channel walls/earth, turning a ducted flow
 *    into a full-cube fluid and corrupting the drag. */
void merge_material_field(MaterialField* base, const MaterialField& overlay,
                          GeometryPartRole role);

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_PART_VOXELIZER_H_
