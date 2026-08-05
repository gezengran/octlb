#ifndef OCTLB_SRC_SOLVER_LBM_DRAG_POSTPROCESSOR_H_
#define OCTLB_SRC_SOLVER_LBM_DRAG_POSTPROCESSOR_H_

#include <array>

#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/bc_kind.h"

namespace octlb {

// Momentum-exchange drag post-processor for a single BlockLattice.
//
// Correct MEM (Latt et al. 2008) for a boundary link whose neighbour (cell+c_i)
// is a solid/Bouzidi cell: the obstacle gains the momentum the fluid sent toward
// the wall MINUS what bounced back, so along direction c_i
//   F_obstacle += (f_i^outgoing + f_bar^bounced) * c_i
// where, under this code's pull scheme:
//   f_i^outgoing  = the fluid cell's PRE-stream (post-collide) population in
//                   direction iPop (toward the wall) -- read from the post-collide
//                   snapshot (BlockLattice::post_collide_populations_at_halo).
//   f_bar^bounced = the fluid cell's POST-stream population in direction
//                   opposite(iPop) (away from the wall) -- stream() pulled it from
//                   the wall cell, where the bounce-back / Bouzidi interpolation
//                   already encoded the wall location (fullway stores f_opp = f_i
//                   for a stationary wall; Bouzidi stores the q-weighted value).
// So the formula needs BOTH the snapshot (for f_i) and the live post-stream
// array (for f_bar). Reading only one -- as the legacy 2*f_i form did -- is wrong:
//   - 2*live[iPop]      : reads the INTERIOR cell's streamed-in outgoing, not F's
//                        own -> large timing error (~6.6x on the cylinder).
//   - 2*snapshot[iPop] : F's own outgoing only; for Bouzidi links the bounced
//                        f_bar != f_i, so it ignores the q-weighting (~2.7x).
// The (snapshot + live[opp]) form is exact for both fullway and Bouzidi links
// without an explicit q_frac in the formula -- stream() already applied it.
//
// If no post-collide snapshot is available (unit tests that set populations by
// hand and never call stream()), f_i falls back to the live array; at the
// equilibrium rest state used by the unit tests live == post-collide, so the
// known-force checks still hold.
//
// The post-processor returns the force on the FLUID (the reaction by which the
// obstacle extracts momentum from the flow): force_on_fluid = -F_obstacle. For
// a flow in +x past an obstacle, force_on_fluid.x < 0.
class MomentumExchangeDrag {
 public:
  explicit MomentumExchangeDrag(
      const BlockLattice<double, olb::descriptors::D3Q19<>>& lat)
      : lat_(lat) {}

  // Force on the fluid, summed over all boundary links, in lattice units.
  // Counts every kSolid/kBouzidi neighbour as part of the obstacle (whole-body
  // drag). Use force_on_fluid_if when only a subset of the boundary -- e.g. just
  // the cylinder surface in a ducted cylinder3d, excluding the duct walls --
  // should contribute to the drag.
  std::array<double, 3> force_on_fluid() const {
    return force_on_fluid_if(
        [](int /*ix*/, int /*iy*/, int /*iz*/) { return true; });
  }

  // Filtered MEM: same sum as force_on_fluid, but a boundary neighbour only
  // contributes when is_obstacle_cell(nxn, nyn, nzn) is true (indices are the
  // neighbour's cell indices in this block). This lets a caller restrict the drag
  // to one body when several kBouzidi/kSolid boundaries share the lattice -- e.g.
  // a cylinder3d where both the curved cylinder surface and the flat duct walls
  // are stamped kBouzidi; without the filter the duct-wall shear (a long strip
  // of Bouzidi cells along the whole duct) dominates the sum and inflates Cd ~100x.
  template <class ObstaclePred>
  std::array<double, 3> force_on_fluid_if(ObstaclePred is_obstacle_cell) const {
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
          // f_i (outgoing, toward wall): post-collide snapshot if available,
          // else the live array (unit-test rest-state fallback).
          const double* snap =
              lat_.post_collide_populations_at_halo(ix + 1, iy + 1, iz + 1);
          const double* live =
              lat_.populations_at_halo(ix + 1, iy + 1, iz + 1);
          const double* f_out_src = snap != nullptr ? snap : live;
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
            if (!is_obstacle_cell(nxn, nyn, nzn)) {
              continue;  // boundary link belongs to a different body -- skip
            }
            const int opp = olb::descriptors::opposite<Descriptor>(iPop);
            const double f_out = f_out_src[iPop];  // toward-wall outgoing
            const double f_bar = live[opp];        // bounced-back incoming
            f_obstacle[0] += (f_out + f_bar) * static_cast<double>(cx);
            f_obstacle[1] += (f_out + f_bar) * static_cast<double>(cy);
            f_obstacle[2] += (f_out + f_bar) * static_cast<double>(cz);
          }
        }
      }
    }
    return {-f_obstacle[0], -f_obstacle[1], -f_obstacle[2]};
  }

  // Legacy fullway MEM (2 * live[iPop] * c_i), kept ONLY for the T11 timing
  // diagnostic: it reads the post-stream live array, which on a boundary link
  // holds the interior cell's streamed-in outgoing -- the wrong timing. The
  // difference between this and force_on_fluid_if quantifies the timing + q
  // correction. Do NOT use for the real Cd.
  template <class ObstaclePred>
  std::array<double, 3> force_on_fluid_if_legacy(
      ObstaclePred is_obstacle_cell) const {
    using Descriptor = olb::descriptors::D3Q19<>;
    std::array<double, 3> f_obstacle{0.0, 0.0, 0.0};
    const int nx = lat_.nx();
    const int ny = lat_.ny();
    const int nz = lat_.nz();
    for (int iz = 0; iz < nz; ++iz) {
      for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
          if (lat_.bc_kind(ix, iy, iz) != BcKind::kBulk) continue;
          const double* cell = lat_.populations_at_halo(ix + 1, iy + 1, iz + 1);
          for (int iPop = 1; iPop < Descriptor::q; ++iPop) {
            const int cx = olb::descriptors::c<Descriptor>(iPop, 0);
            const int cy = olb::descriptors::c<Descriptor>(iPop, 1);
            const int cz = olb::descriptors::c<Descriptor>(iPop, 2);
            const int nxn = ix + cx, nyn = iy + cy, nzn = iz + cz;
            if (nxn < 0 || nxn >= nx || nyn < 0 || nyn >= ny || nzn < 0 ||
                nzn >= nz) continue;
            const BcKind nk = lat_.bc_kind(nxn, nyn, nzn);
            if (nk != BcKind::kSolid && nk != BcKind::kBouzidi) continue;
            if (!is_obstacle_cell(nxn, nyn, nzn)) continue;
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