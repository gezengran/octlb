#include <gtest/gtest.h>

#include <array>

#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/drag_postprocessor.h"

namespace octlb {
namespace {

using Lattice = BlockLattice<double, olb::descriptors::D3Q19<>>;
using Descriptor = olb::descriptors::D3Q19<>;

// D3Q19 equilibrium f_i = t_i * (rho + 3*(c.u) + 4.5*(c.u)^2 - 1.5*u^2)
// for u = (ux, 0, 0) (invCs2 = 3). Fixture helper, not under test.
double Feq(int iPop, double rho, double ux) {
  const int cx = olb::descriptors::c<Descriptor>(iPop, 0);
  const double cu = static_cast<double>(cx) * ux;
  const double usqr = ux * ux;
  return olb::descriptors::t<double, Descriptor>(iPop) *
         (rho + 3.0 * cu + 4.5 * cu * cu - 1.5 * usqr);
}

// One isolated boundary link: fluid cell (0,0,0) at equilibrium rho=1, u=0
// (f_i = t_i), with its only solid neighbour at +x. Diagonal/out-of-bounds
// neighbours are not boundary links, so the single +x face link (iPop 9,
// c=(1,0,0), w=1/18) gives F_obstacle = 2*f_eq*c = 2*(1/18)*(1,0,0) = 1/9 in x.
// force_on_fluid = -F_obstacle -> fx = -1/9.
TEST(MomentumExchangeDrag, OneBoundaryLink_KnownForce) {
  Lattice lat(3, 1, 1);
  lat.set_bc_kind(1, 0, 0, BcKind::kSolid);
  lat.set_bc_kind(2, 0, 0, BcKind::kSolid);
  // cell (0,0,0) at equilibrium rho=1, u=0: f_i = weight t_i
  double* p = lat.populations_at_halo(1, 1, 1);
  for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    p[iPop] = olb::descriptors::t<double, Descriptor>(iPop);
  }

  MomentumExchangeDrag drag(lat);
  const std::array<double, 3> f = drag.force_on_fluid();
  EXPECT_NEAR(f[0], -1.0 / 9.0, 1e-12);
  EXPECT_NEAR(f[1], 0.0, 1e-12);
  EXPECT_NEAR(f[2], 0.0, 1e-12);
}

// A free stream in +x past an obstacle: the obstacle extracts +x momentum from
// the fluid, so force_on_fluid.x < 0 and Cd > 0.
TEST(MomentumExchangeDrag, Sign_OpposesFlow) {
  Lattice lat(3, 1, 1);
  lat.set_bc_kind(1, 0, 0, BcKind::kSolid);
  lat.set_bc_kind(2, 0, 0, BcKind::kSolid);
  const double rho = 1.0;
  const double U = 0.1;  // free stream in +x
  double* p = lat.populations_at_halo(1, 1, 1);
  for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    p[iPop] = Feq(iPop, rho, U);
  }

  MomentumExchangeDrag drag(lat);
  const std::array<double, 3> f = drag.force_on_fluid();
  EXPECT_LT(f[0], 0.0) << "drag on fluid must oppose the +x free stream";
  EXPECT_NEAR(f[1], 0.0, 1e-12);
  EXPECT_NEAR(f[2], 0.0, 1e-12);

  const double area = 1.0;
  const double Cd = drag.drag_coefficient(rho, U, area, {1.0, 0.0, 0.0});
  EXPECT_GT(Cd, 0.0);
}

}  // namespace
}  // namespace octlb