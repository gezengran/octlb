#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <string>

#include "examples/cylinder3d/cylinder3d_case.h"

namespace octlb {
namespace {

int CylinderSteps(int default_steps) {
  const char* env = std::getenv("OCTLB_CYL_STEPS");
  return env != nullptr ? std::atoi(env) : default_steps;
}

std::string CylinderStlPath() {
  return std::string(OCTLB_TEST_DATA_DIR) + "/mesh/cylinder_st.stl";
}

// Schäfer-Turek sanity config (uniform, no AMR). omega/u_inlet are relaxed
// from the strict Re=20 values so the transient is well-behaved on the coarse
// uniform grid; the sanity oracle only checks runs / mass bounded / drag sign.
Cylinder3dConfig SanityConfig(const std::string& stl_path) {
  Cylinder3dConfig cfg;
  cfg.mode = Cylinder3dMode::kUniform;
  cfg.stl_path = stl_path;
  cfg.n = 8;
  cfg.max_level = 2;          // 4^3 = 64 same-level octants
  cfg.omega = 1.0 / 0.6;
  cfg.u_inlet = 0.05;
  cfg.rho0 = 1.0;
  cfg.ramp_end_t = 0.0;
  return cfg;
}

// W3 (T11) Schäfer-Turek cylinder3d uniform sanity. The oracle ladder is
// sanity -> order-of-magnitude -> Cd<1% (W3-e/W4); this wave only asserts the
// pipeline runs on the migrated cubic-domain + carved-channel geometry.
TEST(Cylinder3dUniform, MultiRank_RunsNoCrash) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "two-rank uniform test";
  }

  Cylinder3dCase cas(MPI_COMM_WORLD, SanityConfig(CylinderStlPath()));
  cas.advance_steps(CylinderSteps(20));

  EXPECT_FALSE(cas.has_non_finite_velocity())
      << "uniform cylinder3d produced non-finite velocity";
}

// Open flow-through domain (velocity inlet + pressure outlet): strict mass
// conservation does not hold. The sanity gate is a BOUNDED relative drift.
TEST(Cylinder3dUniform, MassBoundedDrift) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "two-rank uniform test";
  }

  Cylinder3dCase cas(MPI_COMM_WORLD, SanityConfig(CylinderStlPath()));
  const double m0 = cas.total_mass(MPI_COMM_WORLD);
  cas.advance_steps(50);
  const double m1 = cas.total_mass(MPI_COMM_WORLD);

  EXPECT_TRUE(std::isfinite(m1));
  ASSERT_GT(m0, 0.0);
  EXPECT_LT(std::abs(m1 - m0) / m0, 0.2) << "fluid mass drift blowup";
}

// Drag on the cylinder opposes the +x free stream: force_on_fluid.x < 0, so
// Cd = 2*(-force.flow_dir)/(rho*u^2*A) > 0.
TEST(Cylinder3dUniform, DragFinitePositiveSign) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "two-rank uniform test";
  }

  Cylinder3dCase cas(MPI_COMM_WORLD, SanityConfig(CylinderStlPath()));
  cas.advance_steps(100);  // let the +x stream reach and pass the cylinder

  const double Cd = cas.drag_coefficient(MPI_COMM_WORLD, /*area=*/1.0,
                                         /*flow_dir=*/{1.0, 0.0, 0.0});
  EXPECT_TRUE(std::isfinite(Cd));
  EXPECT_GT(Cd, 0.0) << "drag must oppose the +x free stream";
}

}  // namespace
}  // namespace octlb