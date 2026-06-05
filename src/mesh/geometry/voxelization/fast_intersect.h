#ifndef OCTLB_SRC_MESH_GEOMETRY_VOXELIZATION_FAST_INTERSECT_H_
#define OCTLB_SRC_MESH_GEOMETRY_VOXELIZATION_FAST_INTERSECT_H_

#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/geometry/cgal_surface_mesh.h"  // CgalTriangleSurface

namespace octlb {

/** Returns true per octant index when the octant box intersects the surface. */
std::vector<bool> fast_surface_intersect_octants(
    const std::vector<BoundingBox>& octant_boxes,
    const CgalTriangleSurface& surface, scalar bound_width);

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_VOXELIZATION_FAST_INTERSECT_H_
