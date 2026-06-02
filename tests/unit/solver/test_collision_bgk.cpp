#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <numeric>

#include "block_lattice.h"

namespace octlb {
namespace {

using T          = double;
using Descriptor = olb::descriptors::D3Q19<>;

constexpr int kN     = 8;
constexpr int kSteps = 100;
constexpr T   kOmega = 1.4;
constexpr T   kRho0  = 1.0;

// ── helpers ──────────────────────────────────────────────────────────────────

BlockLattice<T, Descriptor> MakeUniform() {
  BlockLattice<T, Descriptor> lat(kN, kN, kN);
  const std::array<T, 3> u0{T{0}, T{0}, T{0}};
  lat.initialize(kRho0, u0.data());
  return lat;
}

struct Conserved {
  T mass{};
  T momentum[3]{};
};

Conserved Integrate(BlockLattice<T, Descriptor>& lat) {
  Conserved q;
  for (int i = 0; i < lat.nx(); ++i)
    for (int j = 0; j < lat.ny(); ++j)
      for (int k = 0; k < lat.nz(); ++k) {
        auto cell = lat.get(i, j, k);
        T rho{}, u[3]{};
        cell.computeRhoU(rho, u);
        q.mass        += rho;
        q.momentum[0] += rho * u[0];
        q.momentum[1] += rho * u[1];
        q.momentum[2] += rho * u[2];
      }
  return q;
}

// ── Test 1: BGK conserves total mass ─────────────────────────────────────────

TEST(CollisionBgk, ConservesMassForUniformState) {
  auto lat = MakeUniform();
  const T mass_before = Integrate(lat).mass;

  for (int t = 0; t < kSteps; ++t) {
    lat.collide(kOmega);
    lat.fill_periodic_halo();  // ghost ← post-collision values; stream() reads ghost
    lat.stream();
  }

  EXPECT_NEAR(Integrate(lat).mass, mass_before, 1e-10);
}

// ── RED test 2: BGK conserves total momentum in x, y, z ──────────────────────

TEST(CollisionBgk, ConservesMomentumXYZForUniformState) {
  auto lat = MakeUniform();
  const auto before = Integrate(lat);

  for (int t = 0; t < kSteps; ++t) {
    lat.collide(kOmega);
    lat.fill_periodic_halo();
    lat.stream();
  }

  const auto after = Integrate(lat);
  EXPECT_NEAR(after.momentum[0], before.momentum[0], 1e-10);
  EXPECT_NEAR(after.momentum[1], before.momentum[1], 1e-10);
  EXPECT_NEAR(after.momentum[2], before.momentum[2], 1e-10);
}

// ── RED test 3: sinusoidal velocity perturbation relaxes to equilibrium ───────
//
// Initial condition: u_x = A * sin(2*pi*i/N) — a pure k=(2pi/N,0,0) mode.
// This mode is physically damped by BGK (not a spurious lattice mode).
//
// Diffusion timescale tau = N^2 / (4*pi^2 * nu) ≈ 64 / (4*pi^2 * 0.071) ≈ 23 steps.
// After 400 steps ≈ 17 tau, the amplitude decays by exp(-17) ≈ 4e-8.
// Initial amplitude = 2*A = 0.1, expected residual << 1e-5.

TEST(CollisionBgk, RelaxesSinusoidalPerturbationToEquilibrium) {
  BlockLattice<T, Descriptor> lat(kN, kN, kN);

  // Set u_x = A*sin(2*pi*i/N) at every cell; all other cells at equilibrium.
  const T kA = T{0.05};
  const T k_wave = T{2} * std::acos(T{-1}) / static_cast<T>(kN);
  for (int i = 0; i < kN; ++i)
    for (int j = 0; j < kN; ++j)
      for (int k = 0; k < kN; ++k) {
        const T ux = kA * std::sin(k_wave * static_cast<T>(i));
        const T u_p[3] = {ux, T{0}, T{0}};
        const T uSqr = ux * ux;
        auto cell = lat.get(i, j, k);
        for (int iPop = 0; iPop < Descriptor::q; ++iPop)
          cell[iPop] = olb::equilibrium<Descriptor>::secondOrder(
              iPop, kRho0, u_p, uSqr);
      }

  // Measure spread = max(ux) - min(ux) over all cells.
  auto ux_spread = [&]() {
    T mn = 1e30, mx = -1e30;
    for (int i = 0; i < lat.nx(); ++i)
      for (int j = 0; j < lat.ny(); ++j)
        for (int k = 0; k < lat.nz(); ++k) {
          T rho{}, u[3]{};
          lat.get(i, j, k).computeRhoU(rho, u);
          mn = std::min(mn, u[0]);
          mx = std::max(mx, u[0]);
        }
    return mx - mn;
  };

  const T initial_spread = ux_spread();  // ≈ 2*kA = 0.1

  for (int t = 0; t < 400; ++t) {
    lat.collide(kOmega);
    lat.fill_periodic_halo();
    lat.stream();
  }

  const T final_spread = ux_spread();
  EXPECT_GT(initial_spread, final_spread) << "spread should decrease";
  EXPECT_LT(final_spread, 1e-5)          << "sinusoidal mode should damp to near-zero";
}

}  // namespace
}  // namespace octlb
