#include "src/solver/lbm/bouzidi_link_data.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "core/platform/platform.h"
#include "descriptor/descriptor.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/geometry/triangle_soup.h"

namespace octlb {
namespace {

using Descriptor = olb::descriptors::D3Q19<>;

bool IsSolidLike(MaterialKind kind) {
  return kind == MaterialKind::kSolid || kind == MaterialKind::kBoundary;
}

void ClosestPointOnSegment(scalar px, scalar py, scalar pz, scalar ax, scalar ay,
                           scalar az, scalar bx, scalar by, scalar bz,
                           scalar* out_x, scalar* out_y, scalar* out_z) {
  const scalar abx = bx - ax;
  const scalar aby = by - ay;
  const scalar abz = bz - az;
  const scalar apx = px - ax;
  const scalar apy = py - ay;
  const scalar apz = pz - az;
  const scalar ab2 = abx * abx + aby * aby + abz * abz;
  scalar t = 0.0;
  if (ab2 > 0.0) {
    t = (apx * abx + apy * aby + apz * abz) / ab2;
    t = std::clamp(t, scalar{0}, scalar{1});
  }
  *out_x = ax + t * abx;
  *out_y = ay + t * aby;
  *out_z = az + t * abz;
}

void ClosestPointOnTriangle(const Triangle& tri, scalar px, scalar py, scalar pz,
                            scalar* out_x, scalar* out_y, scalar* out_z) {
  scalar best_dist2 = std::numeric_limits<scalar>::max();
  scalar bx = tri.v0[0];
  scalar by = tri.v0[1];
  scalar bz = tri.v0[2];

  const auto try_edge = [&](scalar ax, scalar ay, scalar az, scalar ex,
                              scalar ey, scalar ez) {
    scalar cx = 0.0;
    scalar cy = 0.0;
    scalar cz = 0.0;
    ClosestPointOnSegment(px, py, pz, ax, ay, az, ex, ey, ez, &cx, &cy, &cz);
    const scalar dx = px - cx;
    const scalar dy = py - cy;
    const scalar dz = pz - cz;
    const scalar d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < best_dist2) {
      best_dist2 = d2;
      bx = cx;
      by = cy;
      bz = cz;
    }
  };

  try_edge(tri.v0[0], tri.v0[1], tri.v0[2], tri.v1[0], tri.v1[1], tri.v1[2]);
  try_edge(tri.v1[0], tri.v1[1], tri.v1[2], tri.v2[0], tri.v2[1], tri.v2[2]);
  try_edge(tri.v2[0], tri.v2[1], tri.v2[2], tri.v0[0], tri.v0[1], tri.v0[2]);

  *out_x = bx;
  *out_y = by;
  *out_z = bz;
}

TriangleSoup MergeAssemblySoup(const GeometryAssembly& assembly) {
  TriangleSoup merged;
  for (const GeometryPart& part : assembly.parts) {
    merged.merge(part.soup);
  }
  return merged;
}

void CellCenter(const BoundingBox& bounds, int nx, int ny, int nz, int i, int j,
                int k, scalar* cx, scalar* cy, scalar* cz) {
  const scalar dx = (bounds.x_max - bounds.x_min) / static_cast<scalar>(nx);
  const scalar dy = (bounds.y_max - bounds.y_min) / static_cast<scalar>(ny);
  const scalar dz = (bounds.z_max - bounds.z_min) / static_cast<scalar>(nz);
  *cx = bounds.x_min + (static_cast<scalar>(i) + scalar{0.5}) * dx;
  *cy = bounds.y_min + (static_cast<scalar>(j) + scalar{0.5}) * dy;
  *cz = bounds.z_min + (static_cast<scalar>(k) + scalar{0.5}) * dz;
}

MaterialKind NeighborMaterial(const MaterialField& material, OctantId id, int i,
                              int j, int k, int di, int dj, int dk) {
  const int ni = i + di;
  const int nj = j + dj;
  const int nk = k + dk;
  if (ni < 0 || ni >= material.nx() || nj < 0 || nj >= material.ny() ||
      nk < 0 || nk >= material.nz()) {
    return MaterialKind::kFluid;
  }
  return material.at(id, ni, nj, nk);
}

}  // namespace

BouzidiLinkData::BouzidiLinkData(label num_octants, int nx, int ny, int nz,
                                   int q)
    : num_octants_(num_octants),
      nx_(nx),
      ny_(ny),
      nz_(nz),
      q_(q),
      q_frac_(static_cast<std::size_t>(num_octants) *
                  static_cast<std::size_t>(nx * ny * nz * q),
              1.0) {}

std::size_t BouzidiLinkData::index(OctantId id, int ix, int iy, int iz,
                                   int iPop) const {
  const std::size_t cells =
      static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_) *
      static_cast<std::size_t>(nz_);
  return static_cast<std::size_t>(id) * cells * static_cast<std::size_t>(q_) +
         static_cast<std::size_t>((ix * ny_ + iy) * nz_ + iz) *
             static_cast<std::size_t>(q_) +
         static_cast<std::size_t>(iPop);
}

double BouzidiLinkData::q_frac(OctantId id, int ix, int iy, int iz,
                               int iPop) const {
  return q_frac_[index(id, ix, iy, iz, iPop)];
}

void BouzidiLinkData::set_q_frac(OctantId id, int ix, int iy, int iz, int iPop,
                                 double q) {
  q_frac_[index(id, ix, iy, iz, iPop)] = q;
}

BouzidiLinkData BouzidiLinkData::Build(const OctreeForest& forest,
                                       const MaterialField& material,
                                       const GeometryAssembly& assembly) {
  const int nx = material.nx();
  const int ny = material.ny();
  const int nz = material.nz();
  BouzidiLinkData data(material.num_octants(), nx, ny, nz, Descriptor::q);
  const TriangleSoup soup = MergeAssemblySoup(assembly);
  if (soup.empty()) {
    return data;
  }

  for (label oid = 0; oid < material.num_octants(); ++oid) {
    const OctantId id = static_cast<OctantId>(oid);
    const BoundingBox bounds = forest.quadrant_bounds(id);
    const scalar dx =
        (bounds.x_max - bounds.x_min) / static_cast<scalar>(nx);

    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          if (material.at(id, i, j, k) != MaterialKind::kFluid) {
            continue;
          }

          scalar cx = 0.0;
          scalar cy = 0.0;
          scalar cz = 0.0;
          CellCenter(bounds, nx, ny, nz, i, j, k, &cx, &cy, &cz);

          scalar wx = cx;
          scalar wy = cy;
          scalar wz = cz;
          scalar best_dist2 = std::numeric_limits<scalar>::max();
          for (const Triangle& tri : soup.triangles()) {
            scalar px = 0.0;
            scalar py = 0.0;
            scalar pz = 0.0;
            ClosestPointOnTriangle(tri, cx, cy, cz, &px, &py, &pz);
            const scalar d2 = (cx - px) * (cx - px) + (cy - py) * (cy - py) +
                              (cz - pz) * (cz - pz);
            if (d2 < best_dist2) {
              best_dist2 = d2;
              wx = px;
              wy = py;
              wz = pz;
            }
          }

          for (int iPop = 1; iPop < Descriptor::q; ++iPop) {
            const int qOpp = olb::descriptors::opposite<Descriptor>(iPop);
            int ex = 0;
            int ey = 0;
            int ez = 0;
            for (int d = 0; d < Descriptor::d; ++d) {
              const int c = olb::descriptors::c<Descriptor>(qOpp, d);
              if (d == 0) {
                ex = c;
              } else if (d == 1) {
                ey = c;
              } else {
                ez = c;
              }
            }
            const MaterialKind nb = NeighborMaterial(material, id, i, j, k, ex, ey,
                                                     ez);
            if (!IsSolidLike(nb)) {
              continue;
            }

            const scalar dir_x = static_cast<scalar>(ex) * dx;
            const scalar dir_y = static_cast<scalar>(ey) * dx;
            const scalar dir_z = static_cast<scalar>(ez) * dx;
            const scalar dn2 = dir_x * dir_x + dir_y * dir_y + dir_z * dir_z;
            if (dn2 <= 1.0e-30) {
              continue;
            }

            const scalar to_wall_x = wx - cx;
            const scalar to_wall_y = wy - cy;
            const scalar to_wall_z = wz - cz;
            const scalar proj =
                (to_wall_x * dir_x + to_wall_y * dir_y + to_wall_z * dir_z) /
                dn2;
            if (proj > 1.0e-8 && proj < 1.0 - 1.0e-8) {
              data.set_q_frac(id, i, j, k, iPop, proj);
            }
          }
        }
      }
    }
  }

  return data;
}

}  // namespace octlb
