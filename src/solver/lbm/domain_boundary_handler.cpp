#include "src/solver/lbm/domain_boundary_handler.h"

#include <vector>

#include "dynamics/lbm.h"
#include "src/solver/lbm/bc_dispatcher.h"
#include "src/solver/lbm/boundary/interpolated_velocity.h"

namespace octlb {

using Descriptor = olb::descriptors::D3Q19<>;

ConcreteDomainBoundaryHandler::ConcreteDomainBoundaryHandler(
    BlockCollection<DomainBoundaryLattice>& blocks,
    const std::vector<TreeBoundaryFace>& faces,
    const std::vector<DomainBcSpec>& specs, int nx, int ny, int nz,
    double omega, OverlapPaddingMode padding_mode)
    : blocks_(blocks),
      faces_(faces),
      specs_(specs),
      nx_(nx),
      ny_(ny),
      nz_(nz),
      omega_(omega),
      padding_mode_(padding_mode) {}

// Per-cell dispatch is the only path now: each boundary cell collides per its
// own BcKind. Called by the flat advance path (the recursive advance uses
// collide_interleaved_with); kBulk cells are left to BlockLattice::collide so
// the flat path does not double-collide them.
void ConcreteDomainBoundaryHandler::apply(CollideRhoStats* rho_stats,
                                           double average_rho,
                                           bool use_const_rho_bgk) {
  std::vector<bool> seen(static_cast<std::size_t>(blocks_.size()), false);
  for (const TreeBoundaryFace& face : faces_) {
    if (seen[static_cast<std::size_t>(face.octant_id)]) {
      continue;
    }
    seen[static_cast<std::size_t>(face.octant_id)] = true;
    DomainBoundaryLattice& lat = blocks_[face.octant_id];
    for (int ix = 0; ix < nx_; ++ix) {
      for (int iy = 0; iy < ny_; ++iy) {
        for (int iz = 0; iz < nz_; ++iz) {
          if (lat.bc_kind(ix, iy, iz) == BcKind::kBulk) {
            continue;  // bulk handled by BlockLattice::collide in the flat path
          }
          BcDispatcher::collide(lat, ix, iy, iz, nx_, ny_, nz_, omega_,
                                 specs_, rho_stats, average_rho,
                                 use_const_rho_bgk);
        }
      }
    }
  }
}

bool ConcreteDomainBoundaryHandler::collide_boundary_before_bulk() const {
  return true;
}

void ConcreteDomainBoundaryHandler::collide_interleaved_with(
    BlockLattice<double, olb::descriptors::D3Q19<>>& lat,
    CollideRhoStats* rho_stats, double average_rho, bool use_const_rho_bgk) {
  // OpenLB Dominant spatial order: per-cell dispatch walks iX,iY,iZ so each
  // boundary cell collides before its inward fluid neighbors. Centralized in
  // BcDispatcher so mixed BCs coexist.
  for (int ix = 0; ix < nx_; ++ix) {
    for (int iy = 0; iy < ny_; ++iy) {
      for (int iz = 0; iz < nz_; ++iz) {
        BcDispatcher::collide(lat, ix, iy, iz, nx_, ny_, nz_, omega_, specs_,
                              rho_stats, average_rho, use_const_rho_bgk);
      }
    }
  }
  lat.collide_overlap_padding_bgk(omega_);
}

void ConcreteDomainBoundaryHandler::apply_post_stream() {
  std::vector<bool> seen(static_cast<std::size_t>(blocks_.size()), false);
  for (const TreeBoundaryFace& face : faces_) {
    if (seen[static_cast<std::size_t>(face.octant_id)]) {
      continue;
    }
    seen[static_cast<std::size_t>(face.octant_id)] = true;
    DomainBoundaryLattice& lat = blocks_[face.octant_id];
    boundary::ApplyInterpolatedVelocityBoundaryCells<double, Descriptor>(
        lat, nx_, ny_, nz_, omega_, specs_, padding_mode_);
  }
}

}  // namespace octlb