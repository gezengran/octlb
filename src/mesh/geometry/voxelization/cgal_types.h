#ifndef OCTLB_SRC_MESH_GEOMETRY_VOXELIZATION_CGAL_TYPES_H_
#define OCTLB_SRC_MESH_GEOMETRY_VOXELIZATION_CGAL_TYPES_H_

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/box_intersection_d.h>

#include <vector>

namespace octlb {
namespace voxelization {

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using CgalPoint = Kernel::Point_3;
using CgalTriangle = Kernel::Triangle_3;
using CgalBbox = CGAL::Bbox_3;
using Polyhedron = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>;
using SideOfTriangleMesh = CGAL::Side_of_triangle_mesh<Polyhedron, Kernel>;

using BoxIterator = std::vector<CgalBbox>::iterator;
using CgalBoxHandle =
    CGAL::Box_intersection_d::Box_with_handle_d<double, 3, BoxIterator>;

}  // namespace voxelization
}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_VOXELIZATION_CGAL_TYPES_H_
