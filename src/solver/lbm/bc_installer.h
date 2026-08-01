#ifndef OCTLB_SRC_SOLVER_LBM_BC_INSTALLER_H_
#define OCTLB_SRC_SOLVER_LBM_BC_INSTALLER_H_

#include <vector>

#include "src/common/types.h"
#include "src/mesh/geometry/geometry_types.h"
#include "src/mesh/geometry/material_field.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/bc_kind.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {

// BcInstaller: setup-time per-cell BcKind stamping. OctLB equivalent of
// OpenLB's setVelocityBoundary(mat)/setPressureBoundary(mat) -- it resolves
// (geometry MaterialField + TreeBoundaryFace face roles + BouzidiLinkData
// markers) into a per-cell BcKind. The Mesh side still exposes only geometry
// (MaterialKind); the BC role is a solver-side decision, so the seam holds.
//
// Header-only inline: octlb_lbm's only .cpp remains block_lattice.cpp.
namespace bc {

using InstallerLattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

// Geometry material -> BcKind. MaterialKind is geometry-only; the cylinder
// surface cut cells (kBoundary) become kBouzidi; interior solid -> kSolid;
// fluid -> kBulk.
inline BcKind BcKindFromMaterial(MaterialKind mk) {
  switch (mk) {
    case MaterialKind::kFluid:
      return BcKind::kBulk;
    case MaterialKind::kSolid:
      return BcKind::kSolid;
    case MaterialKind::kBoundary:
      return BcKind::kBouzidi;
  }
  return BcKind::kBulk;
}

// Domain face spec type -> BcKind (the per-face role assigned at setup).
inline BcKind BcKindFromSpecType(DomainBcType t) {
  switch (t) {
    case DomainBcType::kNoSlip:
      return BcKind::kBounceBack;
    case DomainBcType::kMovingLid:
      return BcKind::kMovingBounceBack;
    case DomainBcType::kInterpolatedVelocity:
      return BcKind::kVelocityDirichlet;
    case DomainBcType::kOutflow:
      return BcKind::kOutflow;
    case DomainBcType::kInterpolatedPressure:
      return BcKind::kPressureDirichlet;
  }
  return BcKind::kBulk;
}

inline DomainBcSpec FindInstallerSpec(const std::vector<DomainBcSpec>& specs,
                                      FaceDir dir) {
  for (const DomainBcSpec& spec : specs) {
    if (spec.face == dir) {
      return spec;
    }
  }
  return DomainBcSpec{};
}

// Stamp per-cell BcKind from the geometry material field (per octant).
// Replaces the legacy MaterialKind->CellKind direct map.
inline void StampFromMaterial(OctantId id, InstallerLattice& lat,
                              const MaterialField& material) {
  for (int k = 0; k < material.nz(); ++k) {
    for (int j = 0; j < material.ny(); ++j) {
      for (int i = 0; i < material.nx(); ++i) {
        lat.set_bc_kind(i, j, k, BcKindFromMaterial(material.at(id, i, j, k)));
      }
    }
  }
}

inline void StampFromMaterial(BlockCollection<InstallerLattice>& blocks,
                              const MaterialField& material) {
  for (label oid = 0; oid < material.num_octants(); ++oid) {
    StampFromMaterial(static_cast<OctantId>(oid),
                      blocks[static_cast<OctantId>(oid)], material);
  }
}

// Stamp the cells sitting on each domain-outer (tree-boundary) face with the
// BcKind implied by that face's spec. Replaces MarkDomainBoundaryCellKinds.
// Geometry-aware: a face cell already stamped kSolid (e.g. the solid exterior of
// a cubic domain whose channel is a carved sub-region) is left alone -- the
// spec only applies to fluid face cells (kBulk / kBouzidi), not to solid cells.
inline void StampTreeBoundaryCells(BlockCollection<InstallerLattice>& blocks,
                                   const std::vector<TreeBoundaryFace>& faces,
                                   const std::vector<DomainBcSpec>& specs,
                                   int nx, int ny, int nz) {
  for (const TreeBoundaryFace& face : faces) {
    InstallerLattice& lat = blocks[face.octant_id];
    const BcKind kind = BcKindFromSpecType(
        FindInstallerSpec(specs, face.face_dir).type);
    // Stamp the face slab with the spec kind, skipping cells the geometry
    // already marked solid (cubic-domain carved-channel inlet/outlet faces).
    auto stamp = [&](int ix, int iy, int iz) {
      if (lat.bc_kind(ix, iy, iz) != BcKind::kSolid) {
        lat.set_bc_kind(ix, iy, iz, kind);
      }
    };
    switch (face.face_dir) {
      case FaceDir::kXMin:
        for (int iy = 0; iy < ny; ++iy)
          for (int iz = 0; iz < nz; ++iz) stamp(0, iy, iz);
        break;
      case FaceDir::kXMax:
        for (int iy = 0; iy < ny; ++iy)
          for (int iz = 0; iz < nz; ++iz) stamp(nx - 1, iy, iz);
        break;
      case FaceDir::kYMin:
        for (int ix = 0; ix < nx; ++ix)
          for (int iz = 0; iz < nz; ++iz) stamp(ix, 0, iz);
        break;
      case FaceDir::kYMax:
        for (int ix = 0; ix < nx; ++ix)
          for (int iz = 0; iz < nz; ++iz) stamp(ix, ny - 1, iz);
        break;
      case FaceDir::kZMin:
        for (int ix = 0; ix < nx; ++ix)
          for (int iy = 0; iy < ny; ++iy) stamp(ix, iy, 0);
        break;
      case FaceDir::kZMax:
        for (int ix = 0; ix < nx; ++ix)
          for (int iy = 0; iy < ny; ++iy) stamp(ix, iy, nz - 1);
        break;
    }
  }
}

}  // namespace bc

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BC_INSTALLER_H_