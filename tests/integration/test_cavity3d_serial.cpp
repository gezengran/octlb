#include <gtest/gtest.h>
#include <mpi.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

#include "examples/cavity3d/cavity3d_case.h"

namespace octlb {
namespace {

constexpr int kSmokeSteps = 50;

bool IsSingleRank() {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size == 1;
}

int CavityStepsFromEnv(int default_steps) {
  if (const char* env = std::getenv("OCTLB_CAVITY_STEPS")) {
    return std::atoi(env);
  }
  return default_steps;
}

double MidBulkFluidL2VsOpenLbInterp(const CenterlineColumn& col) {
  return BandFluidL2VsOpenLbInterp(col, 5, 21);
}

}  // namespace

TEST(Cavity3dSerial, Smoke_RunsWithoutCrash) {
  if (!IsSingleRank()) {
    GTEST_SKIP() << "single-rank test";
  }

  Cavity3dCase cavity(MPI_COMM_WORLD, UnitConverter::OpenLbCavity3dDefaults());
  EXPECT_NO_THROW(cavity.advance_steps(kSmokeSteps));
  EXPECT_FALSE(cavity.has_non_finite_velocity());
}

TEST(Cavity3dSerial, Smoke_MovingLid_NonZeroUx) {
  if (!IsSingleRank()) {
    GTEST_SKIP() << "single-rank test";
  }

  Cavity3dCase cavity(MPI_COMM_WORLD, UnitConverter::OpenLbCavity3dDefaults());
  cavity.advance_steps(kSmokeSteps);
  EXPECT_TRUE(cavity.any_ux_above(0.0));
}

// PRD #12: OctLB vs OpenLB cavity3d centerline L2 < 2% at convergence.
TEST(Cavity3dSerial, OpenLbCenterline_RelativeL2_Below2Pct) {
  if (!IsSingleRank()) {
    GTEST_SKIP() << "single-rank test";
  }

  const int steps = CavityStepsFromEnv(kOpenLbCavity3dConvergedSteps);
  const UnitConverter converter = UnitConverter::OpenLbCavity3dDefaults();
  Cavity3dCase cavity(MPI_COMM_WORLD, converter);
  cavity.advance_steps(steps);

  const auto ghia_diag =
      AnalyzeGhiaCenterline(cavity.sample_ghia_centerline());
  const CenterlineColumn column = AnalyzeCenterlineColumn(cavity);
  const double l2_openlb = RelativeL2VsOpenLb(ghia_diag.simulated);
  const double mid_bulk_l2 = MidBulkFluidL2VsOpenLbInterp(column);

  if (l2_openlb >= 0.02) {
    PrintGhiaDiagnostics(converter, steps, ghia_diag);
    std::cout << "  mid_bulk_L2_openlb=" << mid_bulk_l2
              << " loop.average_rho=" << cavity.loop.average_rho() << '\n';
    PrintCenterlinePointComparison(converter, steps,
                                   AnalyzeCenterlinePoints(cavity));
  }

  EXPECT_LT(l2_openlb, 0.02)
      << "OpenLB-ref centerline L2=" << l2_openlb
      << " mid_bulk_L2=" << mid_bulk_l2;

  if (CavityWriteVtkEnabled()) {
    const std::string dir =
        (std::filesystem::temp_directory_path() / "octlb_cavity3d_vtk")
            .string();
    std::filesystem::create_directories(dir);
    cavity.write_vtk_timestep(MPI_COMM_WORLD, steps, dir);
  }
}

}  // namespace octlb
