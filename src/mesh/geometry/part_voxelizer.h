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

/** Overwrites \p base with \p overlay (used after sorting parts by ascending priority). */
void merge_material_field(MaterialField* base, const MaterialField& overlay);

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_PART_VOXELIZER_H_
