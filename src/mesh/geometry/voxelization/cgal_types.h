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

// AABB tree over the triangle vector for accelerated box/ray intersection
// queries (voxelization hot path). Built once per surface in from_soup; the
// owning vector must outlive the tree (it does -- it is a member, populated in
// from_soup before the tree and never reallocated after).
using TriangleConstIterator = std::vector<CgalTriangle>::const_iterator;
using AABBPrimitive =
    CGAL::AABB_triangle_primitive<Kernel, TriangleConstIterator>;
using AABBTraits = CGAL::AABB_traits<Kernel, AABBPrimitive>;
using AABBTree = CGAL::AABB_tree<AABBTraits>;

}  // namespace voxelization
}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_VOXELIZATION_CGAL_TYPES_H_
