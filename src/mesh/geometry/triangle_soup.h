#ifndef OCTLB_SRC_MESH_GEOMETRY_TRIANGLE_SOUP_H_
#define OCTLB_SRC_MESH_GEOMETRY_TRIANGLE_SOUP_H_

#include <array>
#include <cstddef>
#include <vector>

#include "src/common/bounding_box.h"

namespace octlb {

/** One STL triangle: 3 vertices (x,y,z) and a unit normal. */
struct Triangle {
  std::array<scalar, 3> normal{};
  std::array<scalar, 3> v0{};
  std::array<scalar, 3> v1{};
  std::array<scalar, 3> v2{};
};

/** Triangle mesh in physical coordinates with cached axis-aligned bounds. */
class TriangleSoup {
 public:
  const std::vector<Triangle>& triangles() const { return triangles_; }
  const BoundingBox& bounding_box() const { return bbox_; }
  bool empty() const { return triangles_.empty(); }

  void clear();
  void add_triangle(const Triangle& tri);
  void merge(const TriangleSoup& other);

  /** Recompute \c bbox_ from all vertices. */
  void update_bounding_box();

 private:
  std::vector<Triangle> triangles_;
  BoundingBox bbox_{};
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_TRIANGLE_SOUP_H_
