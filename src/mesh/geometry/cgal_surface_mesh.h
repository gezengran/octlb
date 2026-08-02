#ifndef OCTLB_SRC_MESH_GEOMETRY_CGAL_SURFACE_MESH_H_
#define OCTLB_SRC_MESH_GEOMETRY_CGAL_SURFACE_MESH_H_

#include <memory>
#include <optional>
#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/geometry/geometry_types.h"
#include "src/mesh/geometry/triangle_soup.h"
#include "src/mesh/geometry/voxelization/cgal_types.h"

namespace octlb {

/** Triangle batch for fast surface intersection (no inside/outside queries). */
class CgalTriangleSurface {
 public:
  static CgalTriangleSurface from_soup(const TriangleSoup& soup);

  const BoundingBox& bounding_box() const { return bbox_; }
  const std::vector<voxelization::CgalTriangle>& cgal_triangles() const {
    return cgal_triangles_;
  }
  bool intersects_box(const BoundingBox& box) const;

 private:
  BoundingBox bbox_{};
  std::vector<voxelization::CgalTriangle> cgal_triangles_;
  // AABB tree over cgal_triangles_ for O(log N) box intersection (accelerates
  // the per-cell voxelization intersects_box query, the brute-force hot path).
  voxelization::AABBTree aabb_tree_;
  // Degenerate (zero-area) triangles excluded from the AABB tree (CGAL AABB
  // requires non-degenerate primitives). intersects_box tests them directly to
  // preserve the original brute-force semantics (which tested every triangle).
  std::vector<voxelization::CgalTriangle> degenerate_triangles_;
};

/** Closed CGAL surface for inside/outside classification (voxelization). */
class CgalSurfaceMesh {
 public:
  static CgalSurfaceMesh from_soup(const TriangleSoup& soup);

  const BoundingBox& bounding_box() const { return bbox_; }
  const std::vector<voxelization::CgalTriangle>& cgal_triangles() const {
    return cgal_triangles_;
  }

  bool is_closed() const;
  bool is_inside(scalar x, scalar y, scalar z) const;
  bool intersects_box(const BoundingBox& box) const;
  scalar outside_sample_fraction(int samples_per_axis) const;

 private:
  CgalSurfaceMesh() = default;

  BoundingBox bbox_{};
  std::vector<voxelization::CgalTriangle> cgal_triangles_;
  voxelization::Polyhedron polyhedron_;
  std::optional<voxelization::SideOfTriangleMesh> side_;
  // AABB tree over cgal_triangles_ accelerating intersects_box (O(log N) vs the
  // brute-force O(triangles) hot path) and the ray-triangle count inside
  // RayParityInside (bbox-cull + exact do_intersect, preserving the exact
  // degenerate-skip parity count). Built in from_soup after cgal_triangles_.
  voxelization::AABBTree aabb_tree_;
  // Degenerate (zero-area) triangles excluded from the AABB tree (CGAL AABB
  // precondition). intersects_box tests them directly to preserve the original
  // brute-force semantics; RayParityInside always skipped them (the tree holds
  // only non-degenerate triangles, so the AABB parity count matches exactly).
  std::vector<voxelization::CgalTriangle> degenerate_triangles_;
};

TriangleSoup merge_assembly_soup(const GeometryAssembly& assembly);

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_CGAL_SURFACE_MESH_H_
