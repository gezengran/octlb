#include "src/mesh/geometry/part_voxelizer.h"

#include <cmath>

#include "src/mesh/geometry/geometry_build_error.h"

namespace octlb {
namespace {

constexpr scalar kMinLumenOutsideFraction = 0.05;
constexpr int kWatertightSampleResolution = 5;

MaterialKind MapExternal(bool intersects, bool inside) {
  if (intersects) {
    return MaterialKind::kBoundary;
  }
  return inside ? MaterialKind::kSolid : MaterialKind::kFluid;
}

MaterialKind MapInternalChannel(bool intersects, bool in_wall, bool in_lumen) {
  if (intersects) {
    return MaterialKind::kBoundary;
  }
  if (in_lumen) {
    return MaterialKind::kFluid;
  }
  if (in_wall) {
    return MaterialKind::kSolid;
  }
  return MaterialKind::kFluid;
}

BoundingBox CellBounds(const BoundingBox& oct, int nx, int ny, int nz, int i,
                       int j, int k) {
  const scalar dx = (oct.x_max - oct.x_min) / static_cast<scalar>(nx);
  const scalar dy = (oct.y_max - oct.y_min) / static_cast<scalar>(ny);
  const scalar dz = (oct.z_max - oct.z_min) / static_cast<scalar>(nz);
  BoundingBox cell;
  cell.x_min = oct.x_min + static_cast<scalar>(i) * dx;
  cell.y_min = oct.y_min + static_cast<scalar>(j) * dy;
  cell.z_min = oct.z_min + static_cast<scalar>(k) * dz;
  cell.x_max = cell.x_min + dx;
  cell.y_max = cell.y_min + dy;
  cell.z_max = cell.z_min + dz;
  return cell;
}

}  // namespace

void validate_geometry_part(const GeometryPart& part) {
  if (part.soup.empty()) {
    throw GeometryBuildError("geometry part has no triangles: " + part.name);
  }
  if (part.role == GeometryPartRole::kInternalChannel) {
    if (part.inner_cavity_soup.empty()) {
      throw GeometryBuildError(
          "kInternalChannel requires inner_cavity_soup: " + part.name);
    }
    const CgalSurfaceMesh outer = CgalSurfaceMesh::from_soup(part.soup);
    const scalar outside_frac =
        outer.outside_sample_fraction(kWatertightSampleResolution);
    if (outside_frac < kMinLumenOutsideFraction) {
      throw GeometryBuildError(
          "kInternalChannel mesh has insufficient exterior/lumen volume "
          "(solid-wall channel expected): " +
          part.name);
    }
  }
}

MaterialField voxelize_part(const OctreeForest& forest, const GeometryPart& part,
                            const GeometryConfig& config) {
  validate_geometry_part(part);

  const CgalSurfaceMesh outer_mesh = CgalSurfaceMesh::from_soup(part.soup);
  std::optional<CgalSurfaceMesh> inner_mesh;
  if (part.role == GeometryPartRole::kInternalChannel) {
    inner_mesh.emplace(CgalSurfaceMesh::from_soup(part.inner_cavity_soup));
  }

  const int nx = config.cell_width;
  const int ny = config.cell_width;
  const int nz = config.cell_width;
  MaterialField field(forest.local_num_octants(), nx, ny, nz);

  for (label id = 0; id < forest.local_num_octants(); ++id) {
    const BoundingBox oct = forest.quadrant_bounds(id);
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          const BoundingBox cell = CellBounds(oct, nx, ny, nz, i, j, k);
          const scalar cx = 0.5 * (cell.x_min + cell.x_max);
          const scalar cy = 0.5 * (cell.y_min + cell.y_max);
          const scalar cz = 0.5 * (cell.z_min + cell.z_max);
          const bool intersects = outer_mesh.intersects_box(cell) ||
                                  (inner_mesh && inner_mesh->intersects_box(cell));

          MaterialKind kind = MaterialKind::kFluid;
          if (part.role == GeometryPartRole::kExternalObstacle) {
            kind = MapExternal(intersects,
                               outer_mesh.is_inside(cx, cy, cz));
          } else {
            const bool in_lumen = inner_mesh->is_inside(cx, cy, cz);
            const bool in_wall =
                outer_mesh.is_inside(cx, cy, cz) && !in_lumen;
            kind = MapInternalChannel(intersects, in_wall, in_lumen);
          }
          field.set(id, i, j, k, kind);
        }
      }
    }
  }
  return field;
}

void merge_material_field(MaterialField* base, const MaterialField& overlay,
                          GeometryPartRole role) {
  if (base->num_octants() != overlay.num_octants() ||
      base->nx() != overlay.nx() || base->ny() != overlay.ny() ||
      base->nz() != overlay.nz()) {
    throw GeometryBuildError("MaterialField merge size mismatch");
  }
  // kInternalChannel defines the fluid domain -> unconditional overwrite (its
  // lumen kFluid carves fluid out of a lower-priority obstacle; its wall
  // kSolid/surface kBoundary likewise take precedence). kExternalObstacle only
  // ADDS itself -> claim-based: overwrite only where it classifies kSolid/kBoundary
  // (its own interior/surface), never where it says kFluid ("outside my obstacle"),
  // so a lower-priority channel's wall/earth survives underneath. See the header
  // doc for the full rationale.
  const bool claim_based = (role == GeometryPartRole::kExternalObstacle);
  for (label id = 0; id < base->num_octants(); ++id) {
    for (int k = 0; k < base->nz(); ++k) {
      for (int j = 0; j < base->ny(); ++j) {
        for (int i = 0; i < base->nx(); ++i) {
          const MaterialKind ov = overlay.at(id, i, j, k);
          if (claim_based && ov == MaterialKind::kFluid) {
            continue;  // obstacle "outside" -- defer to the lower-priority part
          }
          base->set(id, i, j, k, ov);
        }
      }
    }
  }
}

}  // namespace octlb
