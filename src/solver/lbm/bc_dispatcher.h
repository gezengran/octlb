#ifndef OCTLB_SRC_SOLVER_LBM_BC_DISPATCHER_H_
#define OCTLB_SRC_SOLVER_LBM_BC_DISPATCHER_H_

#include <vector>

#include "src/solver/lbm/bc_kind.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {

// Centralized per-cell BC dispatch. Each entry switches on lat.bc_kind(ix,iy,iz)
// and delegates to the existing stateless kernels (algorithm reuse, not rewrite).
// Header-only inline: octlb_lbm's only .cpp remains block_lattice.cpp.
//
// Concrete over the production lattice type (DomainBoundaryLattice), matching
// DomainBoundaryHandler. R0 arms: kBulk (BGK), kBounceBack (full-way
// reflection). Further arms are added by later waves (R1 wires this into the
// handler).
struct BcDispatcher {
  using Descriptor = olb::descriptors::D3Q19<>;

  // Per-cell collide. Mirrors OpenLB per-cell dynamics collide: the cell's
  // BcKind selects the operator. Specs supply prescribed values for Dirichlet
  // arms (unused by kBulk/kBounceBack).
  static inline void collide(DomainBoundaryLattice& lat, int ix, int iy, int iz,
                             double omega,
                             const std::vector<DomainBcSpec>& specs,
                             CollideRhoStats* rho_stats = nullptr,
                             double /*average_rho*/ = 1.0,
                             bool /*use_const_rho_bgk*/ = false) {
    (void)specs;
    (void)rho_stats;
    const BcKind kind = lat.bc_kind(ix, iy, iz);
    auto cell = lat.get(ix, iy, iz);
    switch (kind) {
      case BcKind::kBulk: {
        double rho = 0.0;
        double u[Descriptor::d]{};
        cell.computeRhoU(rho, u);
        olb::lbm<Descriptor>::bgkCollision(cell, rho, u, omega);
        return;
      }
      case BcKind::kBounceBack: {
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
      default:
        // kMovingBounceBack/kBouzidi/kVelocityDirichlet/kPressureDirichlet/
        // kOutflow/kSolid arms are added in R1+; no-op until wired.
        return;
    }
  }
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BC_DISPATCHER_H_