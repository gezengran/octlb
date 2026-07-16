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
//
// boundary_lattice_mode is a transient migration hint (removed in R4): legacy
// mode (cylinder3d) leaves tree-boundary cells as kBulk -- walls/inlet/outlet
// are still enforced by the legacy ghost-fill path. boundary_lattice mode
// (cavity3d) stamps the face slabs so per-cell dispatch takes over.
inline void StampTreeBoundaryCells(BlockCollection<InstallerLattice>& blocks,
                                   const std::vector<TreeBoundaryFace>& faces,
                                   const std::vector<DomainBcSpec>& specs,
                                   int nx, int ny, int nz,
                                   bool boundary_lattice_mode) {
  if (!boundary_lattice_mode) {
    return;
  }
  for (const TreeBoundaryFace& face : faces) {
    InstallerLattice& lat = blocks[face.octant_id];
    const BcKind kind = BcKindFromSpecType(
        FindInstallerSpec(specs, face.face_dir).type);
    switch (face.face_dir) {
      case FaceDir::kXMin:
        for (int iy = 0; iy < ny; ++iy)
          for (int iz = 0; iz < nz; ++iz) lat.set_bc_kind(0, iy, iz, kind);
        break;
      case FaceDir::kXMax:
        for (int iy = 0; iy < ny; ++iy)
          for (int iz = 0; iz < nz; ++iz) lat.set_bc_kind(nx - 1, iy, iz, kind);
        break;
      case FaceDir::kYMin:
        for (int ix = 0; ix < nx; ++ix)
          for (int iz = 0; iz < nz; ++iz) lat.set_bc_kind(ix, 0, iz, kind);
        break;
      case FaceDir::kYMax:
        for (int ix = 0; ix < nx; ++ix)
          for (int iz = 0; iz < nz; ++iz) lat.set_bc_kind(ix, ny - 1, iz, kind);
        break;
      case FaceDir::kZMin:
        for (int ix = 0; ix < nx; ++ix)
          for (int iy = 0; iy < ny; ++iy) lat.set_bc_kind(ix, iy, 0, kind);
        break;
      case FaceDir::kZMax:
        for (int ix = 0; ix < nx; ++ix)
          for (int iy = 0; iy < ny; ++iy) lat.set_bc_kind(ix, iy, nz - 1, kind);
        break;
    }
  }
}

}  // namespace bc

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BC_INSTALLER_H_