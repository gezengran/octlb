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
};

TriangleSoup merge_assembly_soup(const GeometryAssembly& assembly);

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_CGAL_SURFACE_MESH_H_
