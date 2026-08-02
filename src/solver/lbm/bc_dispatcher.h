#ifndef OCTLB_SRC_SOLVER_LBM_BC_DISPATCHER_H_
#define OCTLB_SRC_SOLVER_LBM_BC_DISPATCHER_H_

#include <vector>

#include "src/solver/lbm/bc_kind.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/interpolated_velocity.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {

// Centralized per-cell BC dispatch. Each entry switches on lat.bc_kind(ix,iy,iz)
// and delegates to the existing stateless kernels (algorithm reuse, not rewrite).
// Header-only inline: octlb_lbm's only .cpp remains block_lattice.cpp.
//
// Concrete over the production lattice type (DomainBoundaryLattice), matching
// DomainBoundaryHandler. Arms:
//   kBulk              -> BGK
//   kBounceBack        -> full-way reflection
//   kVelocityDirichlet -> Dirichlet collide (delegates to CollideDirichletBoundaryCellAt)
// kMovingBounceBack / kBouzidi / kPressureDirichlet / kOutflow / kSolid are no-op
// here (Bouzidi/bounce-back act at stream; pressure arm lands in R2; moving-wall
// and outflow in R3+).
struct BcDispatcher {
  using Descriptor = olb::descriptors::D3Q19<>;

  // Per-cell collide. Mirrors OpenLB per-cell dynamics collide: the cell's
  // BcKind selects the operator. nx/ny/nz are the block extent (Dirichlet arm
  // needs face/edge topology). Specs supply prescribed values for Dirichlet
  // arms.
  static inline void collide(DomainBoundaryLattice& lat, int ix, int iy, int iz,
                             int nx, int ny, int nz, double omega,
                             const std::vector<DomainBcSpec>& specs,
                             CollideRhoStats* rho_stats = nullptr,
                             double average_rho = 1.0,
                             bool use_const_rho_bgk = false,
                             double t = 0.0) {
    const BcKind kind = lat.bc_kind(ix, iy, iz);
    switch (kind) {
      case BcKind::kBulk: {
        if (use_const_rho_bgk) {
          // ConstRhoBGK: global average-rho correction (cavity3d bulk). The
          // BlockLattice helper self-gates on bc_kind==kBulk and feeds stats.
          lat.collide_const_rho_at(ix, iy, iz, omega, average_rho, rho_stats);
        } else {
          auto cell = lat.get(ix, iy, iz);
          double rho = 0.0;
          double u[Descriptor::d]{};
          cell.computeRhoU(rho, u);
          olb::lbm<Descriptor>::bgkCollision(cell, rho, u, omega);
        }
        return;
      }
      case BcKind::kBounceBack: {
        auto cell = lat.get(ix, iy, iz);
        // Full-way bounce-back collision: reflect every population to its
        // opposite (rest population maps to itself).
        double tmp[Descriptor::q];
        for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
          tmp[iPop] = cell[iPop];
        }
        for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
          const int opp = olb::descriptors::opposite<Descriptor>(iPop);
          cell[iPop] = tmp[opp];
        }
        return;
      }
      case BcKind::kVelocityDirichlet: {
        boundary::CollideDirichletBoundaryCellAt<double, Descriptor,
                                                  DomainBoundaryLattice>(
            lat, ix, iy, iz, nx, ny, nz, omega, specs, t, rho_stats);
        return;
      }
      default:
        // kMovingBounceBack / kBouzidi / kPressureDirichlet / kOutflow / kSolid:
        // no collide-time action (Bouzidi/bounce-back at stream; pressure in R2;
        // moving-wall/outflow in R3+).
        return;
    }
  }
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BC_DISPATCHER_H_