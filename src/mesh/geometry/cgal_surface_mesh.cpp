#include "src/mesh/geometry/cgal_surface_mesh.h"

#include <CGAL/intersections.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <map>
#include <stdexcept>
#include <tuple>

namespace octlb {
namespace {

using voxelization::CgalPoint;
using voxelization::Kernel;
using Ray = Kernel::Ray_3;
using CgalTriangle = voxelization::CgalTriangle;

void BuildTrianglesFromSoup(const TriangleSoup& soup, BoundingBox* bbox,
                            std::vector<CgalTriangle>* tris) {
  *bbox = soup.bounding_box();
  tris->reserve(soup.triangles().size());
  for (const Triangle& tri : soup.triangles()) {
    const CgalPoint p0(tri.v0[0], tri.v0[1], tri.v0[2]);
    const CgalPoint p1(tri.v1[0], tri.v1[1], tri.v1[2]);
    const CgalPoint p2(tri.v2[0], tri.v2[1], tri.v2[2]);
    tris->emplace_back(p0, p1, p2);
  }
}

std::size_t WeldVertex(std::vector<CgalPoint>* points,
                       std::map<std::tuple<long long, long long, long long>,
                                std::size_t>* index,
                       const CgalPoint& p) {
  const long long scale = 1000000;
  const auto key = std::make_tuple(
      static_cast<long long>(std::llround(CGAL::to_double(p.x()) * scale)),
      static_cast<long long>(std::llround(CGAL::to_double(p.y()) * scale)),
      static_cast<long long>(std::llround(CGAL::to_double(p.z()) * scale)));
  const auto it = index->find(key);
  if (it != index->end()) {
    return it->second;
  }
  const std::size_t id = points->size();
  points->push_back(p);
  index->emplace(key, id);
  return id;
}

bool IsDegenerate(const CgalTriangle& tri) {
  return CGAL::collinear(tri.vertex(0), tri.vertex(1), tri.vertex(2));
}

bool RayParityInsideDirection(const std::vector<CgalTriangle>& tris,
                            const CgalPoint& query,
                            const Kernel::Direction_3& dir) {
  const Ray ray(query, dir);
  int hits = 0;
  for (const CgalTriangle& tri : tris) {
    if (IsDegenerate(tri)) {
      continue;
    }
    if (CGAL::do_intersect(ray, tri)) {
      ++hits;
    }
  }
  return (hits % 2) == 1;
}

bool RayParityInside(const std::vector<CgalTriangle>& tris,
                     const CgalPoint& query) {
  int inside_votes = 0;
  const Kernel::Direction_3 dirs[3] = {Kernel::Direction_3(1, 0, 0),
                                       Kernel::Direction_3(0, 1, 0),
                                       Kernel::Direction_3(0, 0, 1)};
  for (const auto& dir : dirs) {
    if (RayParityInsideDirection(tris, query, dir)) {
      ++inside_votes;
    }
  }
  return inside_votes >= 2;
}

}  // namespace

CgalTriangleSurface CgalTriangleSurface::from_soup(const TriangleSoup& soup) {
  if (soup.empty()) {
    throw std::runtime_error("CgalTriangleSurface: empty triangle soup");
  }
  CgalTriangleSurface mesh;
  BuildTrianglesFromSoup(soup, &mesh.bbox_, &mesh.cgal_triangles_);
  return mesh;
}

bool CgalTriangleSurface::intersects_box(const BoundingBox& box) const {
  const voxelization::CgalBbox query(box.x_min, box.y_min, box.z_min, box.x_max,
                                    box.y_max, box.z_max);
  for (const auto& tri : cgal_triangles_) {
    if (CGAL::do_intersect(tri, query)) {
      return true;
    }
  }
  return false;
}

CgalSurfaceMesh CgalSurfaceMesh::from_soup(const TriangleSoup& soup) {
  if (soup.empty()) {
    throw std::runtime_error("CgalSurfaceMesh: empty triangle soup");
  }
  CgalSurfaceMesh mesh;
  BuildTrianglesFromSoup(soup, &mesh.bbox_, &mesh.cgal_triangles_);

  std::vector<CgalPoint> points;
  std::vector<std::vector<std::size_t>> faces;
  std::map<std::tuple<long long, long long, long long>, std::size_t> index;
  points.reserve(soup.triangles().size() * 3);
  faces.reserve(soup.triangles().size());

  for (const Triangle& tri : soup.triangles()) {
    const CgalPoint p0(tri.v0[0], tri.v0[1], tri.v0[2]);
    const CgalPoint p1(tri.v1[0], tri.v1[1], tri.v1[2]);
    const CgalPoint p2(tri.v2[0], tri.v2[1], tri.v2[2]);
    faces.push_back({WeldVertex(&points, &index, p0), WeldVertex(&points, &index, p1),
                     WeldVertex(&points, &index, p2)});
  }

  namespace PMP = CGAL::Polygon_mesh_processing;
  PMP::repair_polygon_soup(points, faces);
  PMP::orient_polygon_soup(points, faces);
  if (PMP::is_polygon_soup_a_polygon_mesh(faces)) {
    PMP::polygon_soup_to_polygon_mesh(points, faces, mesh.polyhedron_);
    PMP::stitch_borders(mesh.polyhedron_);
    if (mesh.polyhedron_.is_closed()) {
      mesh.side_.emplace(mesh.polyhedron_);
    }
  }
  return mesh;
}

bool CgalSurfaceMesh::is_closed() const {
  return static_cast<bool>(side_);
}

bool CgalSurfaceMesh::is_inside(scalar x, scalar y, scalar z) const {
  const CgalPoint query(x, y, z);
  const bool ray = RayParityInside(cgal_triangles_, query);
  if (side_) {
    const bool cgal = (*side_)(query) == CGAL::ON_BOUNDED_SIDE;
    return cgal == ray ? cgal : ray;
  }
  return ray;
}

bool CgalSurfaceMesh::intersects_box(const BoundingBox& box) const {
  const voxelization::CgalBbox query(box.x_min, box.y_min, box.z_min, box.x_max,
                                    box.y_max, box.z_max);
  for (const auto& tri : cgal_triangles_) {
    if (CGAL::do_intersect(tri, query)) {
      return true;
    }
  }
  return false;
}

scalar CgalSurfaceMesh::outside_sample_fraction(int samples_per_axis) const {
  if (samples_per_axis < 2) {
    return 0.0;
  }
  int outside_count = 0;
  int total = 0;
  const scalar dx =
      (bbox_.x_max - bbox_.x_min) / static_cast<scalar>(samples_per_axis - 1);
  const scalar dy =
      (bbox_.y_max - bbox_.y_min) / static_cast<scalar>(samples_per_axis - 1);
  const scalar dz =
      (bbox_.z_max - bbox_.z_min) / static_cast<scalar>(samples_per_axis - 1);
  for (int iz = 0; iz < samples_per_axis; ++iz) {
    for (int iy = 0; iy < samples_per_axis; ++iy) {
      for (int ix = 0; ix < samples_per_axis; ++ix) {
        const scalar x = bbox_.x_min + static_cast<scalar>(ix) * dx;
        const scalar y = bbox_.y_min + static_cast<scalar>(iy) * dy;
        const scalar z = bbox_.z_min + static_cast<scalar>(iz) * dz;
        if (!is_inside(x, y, z)) {
          ++outside_count;
        }
        ++total;
      }
    }
  }
  return total > 0 ? static_cast<scalar>(outside_count) / static_cast<scalar>(total)
                   : 0.0;
}

TriangleSoup merge_assembly_soup(const GeometryAssembly& assembly) {
  TriangleSoup merged;
  for (const GeometryPart& part : assembly.parts) {
    merged.merge(part.soup);
  }
  return merged;
}

}  // namespace octlb
