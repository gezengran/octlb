#ifndef OCTLB_TESTS_UNIT_MESH_GEOMETRY_FIXTURES_H_
#define OCTLB_TESTS_UNIT_MESH_GEOMETRY_FIXTURES_H_

#include <cmath>

#include "src/mesh/geometry/geometry_types.h"
#include "src/mesh/geometry/triangle_soup.h"

namespace octlb {
namespace test_geom {

inline void AddTriangle(TriangleSoup* soup, const std::array<scalar, 3>& a,
                        const std::array<scalar, 3>& b,
                        const std::array<scalar, 3>& c) {
  Triangle tri;
  tri.v0 = a;
  tri.v1 = b;
  tri.v2 = c;
  const scalar ux = b[0] - a[0];
  const scalar uy = b[1] - a[1];
  const scalar uz = b[2] - a[2];
  const scalar vx = c[0] - a[0];
  const scalar vy = c[1] - a[1];
  const scalar vz = c[2] - a[2];
  tri.normal[0] = uy * vz - uz * vy;
  tri.normal[1] = uz * vx - ux * vz;
  tri.normal[2] = ux * vy - uy * vx;
  const scalar len = std::sqrt(tri.normal[0] * tri.normal[0] +
                               tri.normal[1] * tri.normal[1] +
                               tri.normal[2] * tri.normal[2]);
  if (len > 0.0) {
    tri.normal[0] /= len;
    tri.normal[1] /= len;
    tri.normal[2] /= len;
  }
  soup->add_triangle(tri);
}

inline void AddAxisAlignedBox(TriangleSoup* soup, scalar xmin, scalar ymin,
                              scalar zmin, scalar xmax, scalar ymax,
                              scalar zmax, bool inward_normals) {
  const auto p000 = std::array<scalar, 3>{xmin, ymin, zmin};
  const auto p100 = std::array<scalar, 3>{xmax, ymin, zmin};
  const auto p010 = std::array<scalar, 3>{xmin, ymax, zmin};
  const auto p110 = std::array<scalar, 3>{xmax, ymax, zmin};
  const auto p001 = std::array<scalar, 3>{xmin, ymin, zmax};
  const auto p101 = std::array<scalar, 3>{xmax, ymin, zmax};
  const auto p011 = std::array<scalar, 3>{xmin, ymax, zmax};
  const auto p111 = std::array<scalar, 3>{xmax, ymax, zmax};
  auto tri = [&](const std::array<scalar, 3>& a, const std::array<scalar, 3>& b,
                 const std::array<scalar, 3>& c) {
    if (inward_normals) {
      AddTriangle(soup, a, c, b);
    } else {
      AddTriangle(soup, a, b, c);
    }
  };
  tri(p000, p100, p010);
  tri(p100, p110, p010);
  tri(p100, p110, p101);
  tri(p100, p101, p110);
  tri(p110, p011, p101);
  tri(p110, p011, p111);
  tri(p010, p110, p011);
  tri(p010, p011, p001);
  tri(p001, p011, p101);
  tri(p001, p101, p100);
  tri(p000, p010, p001);
  tri(p000, p001, p100);
}

inline TriangleSoup SolidCube(scalar xmin, scalar ymin, scalar zmin, scalar xmax,
                              scalar ymax, scalar zmax) {
  TriangleSoup soup;
  AddAxisAlignedBox(&soup, xmin, ymin, zmin, xmax, ymax, zmax, false);
  return soup;
}

inline void FillHollowChannelPart(GeometryPart* part, scalar outer_min,
                                  scalar outer_max, scalar inner_min,
                                  scalar inner_max) {
  AddAxisAlignedBox(&part->soup, outer_min, outer_min, outer_min, outer_max,
                    outer_max, outer_max, false);
  AddAxisAlignedBox(&part->inner_cavity_soup, inner_min, inner_min, inner_min,
                    inner_max, inner_max, inner_max, true);
}

inline TriangleSoup CenteredSphere(scalar cx, scalar cy, scalar cz, scalar r,
                                   int segments) {
  TriangleSoup soup;
  for (int lat = 0; lat < segments; ++lat) {
    const scalar theta0 =
        static_cast<scalar>(lat) * static_cast<scalar>(M_PI) / segments;
    const scalar theta1 =
        static_cast<scalar>(lat + 1) * static_cast<scalar>(M_PI) / segments;
    for (int lon = 0; lon < segments; ++lon) {
      const scalar phi0 =
          static_cast<scalar>(lon) * 2.0 * static_cast<scalar>(M_PI) / segments;
      const scalar phi1 = static_cast<scalar>(lon + 1) * 2.0 *
                          static_cast<scalar>(M_PI) / segments;
      auto sph = [&](scalar theta, scalar phi) {
        const scalar st = std::sin(theta);
        return std::array<scalar, 3>{
            cx + r * st * std::cos(phi), cy + r * st * std::sin(phi),
            cz + r * std::cos(theta)};
      };
      const auto p00 = sph(theta0, phi0);
      const auto p01 = sph(theta0, phi1);
      const auto p10 = sph(theta1, phi0);
      const auto p11 = sph(theta1, phi1);
      AddTriangle(&soup, p00, p10, p01);
      AddTriangle(&soup, p01, p10, p11);
    }
  }
  return soup;
}

inline GeometryConfig DefaultConfig() {
  GeometryConfig cfg;
  cfg.max_level = 6;
  cfg.resolve_surface_times = 4;
  cfg.bound_width = 0.05;
  cfg.cell_width = 4;
  return cfg;
}

}  // namespace test_geom
}  // namespace octlb

#endif  // OCTLB_TESTS_UNIT_MESH_GEOMETRY_FIXTURES_H_
