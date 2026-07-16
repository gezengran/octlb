#ifndef OCTLB_SRC_SOLVER_LBM_DRAG_POSTPROCESSOR_H_
#define OCTLB_SRC_SOLVER_LBM_DRAG_POSTPROCESSOR_H_

#include <array>

#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/bc_kind.h"

namespace octlb {

// Momentum-exchange drag post-processor for a single BlockLattice.
//
// Fullway bounce-back MEM: for every fluid cell and every direction iPop whose
// neighbour (cell + c_i) is a solid/boundary cell, the obstacle gains momentum
//   F_obstacle += 2 * f_i(cell) * c_i
// along that boundary link (the factor 2 accounts for the incident population
// and its bounced-back opposite; fullway BB stores f_opp = f_i at the solid).
// The post-processor returns the force on the FLUID (the reaction by which the
// obstacle extracts momentum from the flow): force_on_fluid = -F_obstacle. For
// a flow in +x past an obstacle, force_on_fluid.x < 0.
//
// Reads only fluid-cell populations (no solid/ghost storage dependency), so it
// is unit-testable with a hand-built lattice. Bouzidi-corrected MEM (using the
// partial-boundary q_frac) is deferred to a later wave; the fullway form is
// sufficient for the W1 sanity oracle.
class MomentumExchangeDrag {
 public:
  explicit MomentumExchangeDrag(
      const BlockLattice<double, olb::descriptors::D3Q19<>>& lat)
      : lat_(lat) {}

  // Force on the fluid, summed over all boundary links, in lattice units.
  std::array<double, 3> force_on_fluid() const {
    using Descriptor = olb::descriptors::D3Q19<>;
    std::array<double, 3> f_obstacle{0.0, 0.0, 0.0};
    const int nx = lat_.nx();
    const int ny = lat_.ny();
    const int nz = lat_.nz();
    for (int iz = 0; iz < nz; ++iz) {
      for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
          if (lat_.bc_kind(ix, iy, iz) != BcKind::kBulk) {
            continue;
          }
          const double* cell = lat_.populations_at_halo(ix + 1, iy + 1, iz + 1);
          for (int iPop = 1; iPop < Descriptor::q; ++iPop) {
            const int cx = olb::descriptors::c<Descriptor>(iPop, 0);
            const int cy = olb::descriptors::c<Descriptor>(iPop, 1);
            const int cz = olb::descriptors::c<Descriptor>(iPop, 2);
            const int nxn = ix + cx;
            const int nyn = iy + cy;
            const int nzn = iz + cz;
            if (nxn < 0 || nxn >= nx || nyn < 0 || nyn >= ny || nzn < 0 ||
                nzn >= nz) {
              continue;  // halo neighbour, not a local obstacle link
            }
            const BcKind nk = lat_.bc_kind(nxn, nyn, nzn);
            if (nk != BcKind::kSolid && nk != BcKind::kBouzidi) {
              continue;
            }
            const double fi = cell[iPop];
            f_obstacle[0] += 2.0 * fi * static_cast<double>(cx);
            f_obstacle[1] += 2.0 * fi * static_cast<double>(cy);
            f_obstacle[2] += 2.0 * fi * static_cast<double>(cz);
          }
        }
      }
    }
    return {-f_obstacle[0], -f_obstacle[1], -f_obstacle[2]};
  }

  // Drag coefficient Cd = 2 * (-force_on_fluid . flow_dir) / (rho * u^2 * area),
  // where flow_dir is the unit direction of the free stream. Positive for a
  // body in a flow (force_on_fluid opposes the stream).
  double drag_coefficient(double rho, double u, double area,
                          std::array<double, 3> flow_dir) const {
    const std::array<double, 3> f = force_on_fluid();
    const double fdot = f[0] * flow_dir[0] + f[1] * flow_dir[1] +
                        f[2] * flow_dir[2];
    return 2.0 * (-fdot) / (rho * u * u * area);
  }

 private:
  const BlockLattice<double, olb::descriptors::D3Q19<>>& lat_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_DRAG_POSTPROCESSOR_H_