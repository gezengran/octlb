#include <gtest/gtest.h>

#include <array>
#include <memory>

#include "src/solver/lbm/boundary/inlet_velocity_field.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {
namespace {

using namespace octlb::boundary;

// Poiseuille duct profile on a square cross-section: at full ramp the
// centerline cell (in-plane centre) reaches u_peak * flow_dir, wall-adjacent
// cells are slower per the parabolic profile, and the off-flow components are
// zero. Verified against the analytic parab(s) = 4*c*(1-c), c=(i+0.5)/N.
TEST(InletVelocityField, PoiseuilleProfile_MatchesAnalytic_AtFullRamp) {
  // Odd in-plane extent N=7 so a cell centre sits exactly at c=0.5 (peak).
  constexpr int N = 7;
  constexpr double kPeak = 0.045;            // lattice centreline velocity
  const std::array<double, 3> flow_dir{1.0, 0.0, 0.0};
  PoiseuilleInletProfile profile(FaceDir::kXMin, N, N, kPeak, flow_dir,
                                 /*ramp_end_t=*/0.0);

  auto parab = [](int i, int n) {
    const double c = (i + 0.5) / static_cast<double>(n);
    return 4.0 * c * (1.0 - c);
  };

  // Centre cell (iy=3, iz=3): both parabolas = 1 -> u = peak * flow_dir.
  {
    double u[3] = {-1.0, -1.0, -1.0};
    profile.velocity(/*ix=*/0, /*iy=*/3, /*iz=*/3, /*t=*/0.0, u);
    EXPECT_NEAR(u[0], kPeak, 1e-12);
    EXPECT_NEAR(u[1], 0.0, 1e-12);
    EXPECT_NEAR(u[2], 0.0, 1e-12);
  }

  // Edge cell (centre in y, wall-adjacent in z): parab_y = 1, parab_z < 1.
  {
    double u[3] = {-1.0, -1.0, -1.0};
    profile.velocity(/*ix=*/0, /*iy=*/3, /*iz=*/0, /*t=*/0.0, u);
    const double expected = kPeak * parab(3, N) * parab(0, N);
    EXPECT_NEAR(u[0], expected, 1e-12);
    EXPECT_NEAR(u[1], 0.0, 1e-12);
    EXPECT_NEAR(u[2], 0.0, 1e-12);
  }

  // Corner cell (wall-adjacent in both y and z): parab_y * parab_z.
  {
    double u[3] = {-1.0, -1.0, -1.0};
    profile.velocity(/*ix=*/0, /*iy=*/0, /*iz=*/0, /*t=*/0.0, u);
    const double expected = kPeak * parab(0, N) * parab(0, N);
    EXPECT_NEAR(u[0], expected, 1e-12);
    EXPECT_NEAR(u[1], 0.0, 1e-12);
    EXPECT_NEAR(u[2], 0.0, 1e-12);
  }
}

// Time ramp: with ramp_end_t > 0 the velocity scales by OpenLB's
// PolynomialStartScale (smootherstep 10x^3 - 15x^4 + 6x^5, x = t/ramp_end).
// At t=0 -> 0, t=ramp_end/2 -> 1/2 (smootherstep is symmetric), t>=ramp_end ->
// full velocity (matching the no-ramp profile). Sampled at the centre cell so
// the spatial parabolas are both 1 and only the ramp factor is exercised.
TEST(InletVelocityField, PoiseuilleProfile_RampScalesVelocityOverTime) {
  constexpr int N = 7;
  constexpr double kPeak = 0.045;
  const std::array<double, 3> flow_dir{1.0, 0.0, 0.0};
  constexpr double kRampEnd = 10.0;
  PoiseuilleInletProfile profile(FaceDir::kXMin, N, N, kPeak, flow_dir,
                                 kRampEnd);

  auto u_at = [&](double t) {
    double u[3] = {-1.0, -1.0, -1.0};
    profile.velocity(/*ix=*/0, /*iy=*/3, /*iz=*/3, t, u);
    return u[0];
  };

  EXPECT_NEAR(u_at(0.0), 0.0, 1e-12) << "ramp starts from zero";
  EXPECT_NEAR(u_at(kRampEnd / 2.0), kPeak / 2.0, 1e-12)
      << "smootherstep at x=0.5 is exactly 0.5";
  EXPECT_NEAR(u_at(kRampEnd), kPeak, 1e-12)
      << "ramp reaches full velocity at ramp_end_t";
  EXPECT_NEAR(u_at(2.0 * kRampEnd), kPeak, 1e-12)
      << "velocity holds after ramp_end_t (clamp)";

  // At ramp end the profile must match the no-ramp (full) profile everywhere,
  // not just at the centre.
  PoiseuilleInletProfile full(FaceDir::kXMin, N, N, kPeak, flow_dir, 0.0);
  for (int iy : {0, 1, 3, 5, 6}) {
    for (int iz : {0, 3, 6}) {
      double ur[3] = {-1, -1, -1};
      double uf[3] = {-1, -1, -1};
      profile.velocity(0, iy, iz, kRampEnd, ur);
      full.velocity(0, iy, iz, 0.0, uf);
      EXPECT_NEAR(ur[0], uf[0], 1e-12) << "ramped matches full at iy=" << iy;
    }
  }
}

// The BC velocity lookup: when a spec carries an inlet_field, the prescribed
// per-cell velocity comes from the field (spatially varying); otherwise it
// falls back to the constant u_wall. Backward compatible -- existing specs
// without an inlet_field keep their constant-u behaviour.
TEST(InletVelocityField, PrescribedVelocity_UsesInletFieldWhenSet) {
  constexpr int N = 7;
  constexpr double kPeak = 0.045;
  const std::array<double, 3> flow_dir{1.0, 0.0, 0.0};

  DomainBcSpec spec_field;
  spec_field.face = FaceDir::kXMin;
  spec_field.type = DomainBcType::kInterpolatedVelocity;
  spec_field.inlet_field = std::make_shared<PoiseuilleInletProfile>(
      FaceDir::kXMin, N, N, kPeak, flow_dir, 0.0);

  // Centre cell: both parabolas = 1, no ramp -> full peak along flow_dir.
  {
    double u[3] = {-1.0, -1.0, -1.0};
    PrescribedVelocity(spec_field, /*ix=*/0, /*iy=*/3, /*iz=*/3, /*t=*/0.0, u);
    EXPECT_NEAR(u[0], kPeak, 1e-12);
    EXPECT_NEAR(u[1], 0.0, 1e-12);
    EXPECT_NEAR(u[2], 0.0, 1e-12);
  }
  // Corner cell: parabolic value (spatial variation reaches the lookup).
  {
    auto parab = [](int i, int n) {
      const double c = (i + 0.5) / static_cast<double>(n);
      return 4.0 * c * (1.0 - c);
    };
    double u[3] = {-1.0, -1.0, -1.0};
    PrescribedVelocity(spec_field, 0, 0, 0, 0.0, u);
    EXPECT_NEAR(u[0], kPeak * parab(0, N) * parab(0, N), 1e-12);
  }

  // Without inlet_field: falls back to constant u_wall.
  DomainBcSpec spec_const;
  spec_const.face = FaceDir::kXMin;
  spec_const.type = DomainBcType::kInterpolatedVelocity;
  spec_const.u_wall = {0.1, 0.0, 0.0};
  {
    double u[3] = {-1.0, -1.0, -1.0};
    PrescribedVelocity(spec_const, 0, 3, 3, 0.0, u);
    EXPECT_NEAR(u[0], 0.1, 1e-12);
    EXPECT_NEAR(u[1], 0.0, 1e-12);
    EXPECT_NEAR(u[2], 0.0, 1e-12);
  }
}

}  // namespace
}  // namespace octlb