#include "src/mesh/geometry/voxelization/fast_intersect.h"

#include <CGAL/box_intersection_d.h>

namespace octlb {
namespace {

class BoxIntersectionCallback {
 public:
  BoxIntersectionCallback(std::vector<bool>* flags,
                          const std::vector<voxelization::CgalBbox>& oct_boxes,
                          const std::vector<voxelization::CgalBbox>& tri_boxes,
                          const std::vector<voxelization::CgalTriangle>& tris)
      : flags_(flags),
        oct_boxes_(oct_boxes),
        tri_boxes_(tri_boxes),
        tris_(tris) {}

  void operator()(const voxelization::CgalBoxHandle& tri_box,
                  const voxelization::CgalBoxHandle& oct_box) const {
    const std::size_t oct_id =
        static_cast<std::size_t>(oct_box.handle() - oct_boxes_.begin());
    const std::size_t tri_id =
        static_cast<std::size_t>(tri_box.handle() - tri_boxes_.begin());
    if ((*flags_)[oct_id]) {
      return;
    }
    if (CGAL::do_intersect(tris_[tri_id], *(oct_box.handle()))) {
      (*flags_)[oct_id] = true;
    }
  }

 private:
  std::vector<bool>* flags_;
  const std::vector<voxelization::CgalBbox>& oct_boxes_;
  const std::vector<voxelization::CgalBbox>& tri_boxes_;
  const std::vector<voxelization::CgalTriangle>& tris_;
};

}  // namespace

std::vector<bool> fast_surface_intersect_octants(
    const std::vector<BoundingBox>& octant_boxes,
    const CgalTriangleSurface& surface,
    scalar bound_width) {
  std::vector<bool> crossed(octant_boxes.size(), false);
  const auto& tris = surface.cgal_triangles();
  if (tris.empty()) {
    return crossed;
  }

  std::vector<voxelization::CgalBbox> tri_boxes;
  tri_boxes.reserve(tris.size());
  for (const auto& tri : tris) {
    tri_boxes.push_back(tri.bbox());
  }

  std::vector<voxelization::CgalBbox> oct_boxes;
  oct_boxes.reserve(octant_boxes.size());
  for (const BoundingBox& oct : octant_boxes) {
    oct_boxes.emplace_back(oct.x_min - bound_width, oct.y_min - bound_width,
                           oct.z_min - bound_width, oct.x_max + bound_width,
                           oct.y_max + bound_width, oct.z_max + bound_width);
  }

  std::vector<voxelization::CgalBoxHandle> tri_handles;
  for (auto it = tri_boxes.begin(); it != tri_boxes.end(); ++it) {
    tri_handles.emplace_back(*it, it);
  }
  std::vector<voxelization::CgalBoxHandle> oct_handles;
  for (auto it = oct_boxes.begin(); it != oct_boxes.end(); ++it) {
    oct_handles.emplace_back(*it, it);
  }

  BoxIntersectionCallback callback(&crossed, oct_boxes, tri_boxes, tris);
  CGAL::box_intersection_all_pairs_d(
      tri_handles.begin(), tri_handles.end(), oct_handles.begin(),
      oct_handles.end(), callback);
  return crossed;
}

}  // namespace octlb
