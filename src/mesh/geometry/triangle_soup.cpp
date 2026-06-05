#include "src/mesh/geometry/triangle_soup.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace octlb {

void TriangleSoup::clear() {
  triangles_.clear();
  bbox_ = BoundingBox{};
}

void TriangleSoup::add_triangle(const Triangle& tri) {
  triangles_.push_back(tri);
  update_bounding_box();
}

void TriangleSoup::merge(const TriangleSoup& other) {
  triangles_.insert(triangles_.end(), other.triangles().begin(),
                    other.triangles().end());
  update_bounding_box();
}

void TriangleSoup::update_bounding_box() {
  if (triangles_.empty()) {
    bbox_ = BoundingBox{};
    return;
  }
  scalar xmin = std::numeric_limits<scalar>::max();
  scalar ymin = std::numeric_limits<scalar>::max();
  scalar zmin = std::numeric_limits<scalar>::max();
  scalar xmax = std::numeric_limits<scalar>::lowest();
  scalar ymax = std::numeric_limits<scalar>::lowest();
  scalar zmax = std::numeric_limits<scalar>::lowest();
  for (const Triangle& tri : triangles_) {
    for (const auto* v : {&tri.v0, &tri.v1, &tri.v2}) {
      xmin = std::min(xmin, (*v)[0]);
      ymin = std::min(ymin, (*v)[1]);
      zmin = std::min(zmin, (*v)[2]);
      xmax = std::max(xmax, (*v)[0]);
      ymax = std::max(ymax, (*v)[1]);
      zmax = std::max(zmax, (*v)[2]);
    }
  }
  bbox_.x_min = xmin;
  bbox_.y_min = ymin;
  bbox_.z_min = zmin;
  bbox_.x_max = xmax;
  bbox_.y_max = ymax;
  bbox_.z_max = zmax;
}

}  // namespace octlb
