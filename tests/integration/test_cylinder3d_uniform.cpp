#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "examples/cylinder3d/cylinder3d_case.h"

namespace octlb {
namespace {

int CylinderSteps(int default_steps) {
  const char* env = std::getenv("OCTLB_CYL_STEPS");
  return env != nullptr ? std::atoi(env) : default_steps;
}

std::string CylinderStlPath() {
  return std::string(OCTLB_TEST_DATA_DIR) + "/mesh/cylinder.stl";
}

// W2 (T11) uniform cylinder3d sanity. The oracle ladder is sanity -> order-of-
// magnitude -> Cd<1% (W3/W4); this wave only asserts the pipeline runs.
TEST(Cylinder3dUniform, MultiRank_RunsNoCrash) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "two-rank uniform test";
  }

  const std::string stl_path = CylinderStlPath();
  Cylinder3dCase cas(MPI_COMM_WORLD, stl_path, /*n=*/4, /*omega=*/1.0 / 0.6,
                     /*u_inlet=*/0.05, /*rho0=*/1.0);
  cas.advance_steps(CylinderSteps(20));

  EXPECT_FALSE(cas.has_non_finite_velocity())
      << "uniform cylinder3d produced non-finite velocity";
}

// Open flow-through domain (Zou-He inlet + outflow outlet): strict mass
// conservation does not hold. The sanity gate is a BOUNDED relative drift --
// the initial transient (inlet driving flow from u=0, pressure field forming)
// redistributes density, but a real instability would drift far past this.
TEST(Cylinder3dUniform, MassBoundedDrift) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "two-rank uniform test";
  }

  Cylinder3dCase cas(MPI_COMM_WORLD, CylinderStlPath(), /*n=*/4,
                     /*omega=*/1.0 / 0.6, /*u_inlet=*/0.05, /*rho0=*/1.0);
  const double m0 = cas.total_mass(MPI_COMM_WORLD);
  cas.advance_steps(50);
  const double m1 = cas.total_mass(MPI_COMM_WORLD);

  EXPECT_TRUE(std::isfinite(m1));
  ASSERT_GT(m0, 0.0);
  EXPECT_LT(std::abs(m1 - m0) / m0, 0.2) << "fluid mass drift blowup";
}

// Drag on the cylinder opposes the +x free stream: force_on_fluid.x < 0, so
// Cd = 2*(-force.flow_dir)/(rho*u^2*A) > 0. Area is the projected area (its
// exact value matters for W4 accuracy, not the sign).
TEST(Cylinder3dUniform, DragFinitePositiveSign) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "two-rank uniform test";
  }

  Cylinder3dCase cas(MPI_COMM_WORLD, CylinderStlPath(), /*n=*/4,
                     /*omega=*/1.0 / 0.6, /*u_inlet=*/0.05, /*rho0=*/1.0);
  cas.advance_steps(100);  // let the +x stream reach and pass the cylinder

  const double Cd = cas.drag_coefficient(MPI_COMM_WORLD, /*area=*/1.0,
                                         /*flow_dir=*/{1.0, 0.0, 0.0});
  EXPECT_TRUE(std::isfinite(Cd));
  EXPECT_GT(Cd, 0.0) << "drag must oppose the +x free stream";
}

// ② edge-ghost probe -- DISABLED: this multi-vs-single comparison is confounded
// by per-octant vs whole-grid voxelization differences at the cylinder boundary
// (invisible at step 0 because solid and fluid-at-rest both store f=0 -> rho=1,
// u=0 under the f=f_eq-t convention, then diverge at step 1). It does NOT
// cleanly isolate ②. Kept for reference; ② needs a clean re-probe (no-cylinder
// multi-vs-single, or a with/without-edge-exchange toggle).
TEST(Cylinder3dUniform, DISABLED_NoEdgeArtifact_MatchesSingleBlock) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank edge-ghost probe";
  }

  Cylinder3dCase multi(MPI_COMM_WORLD, CylinderStlPath(), /*n=*/4, 1.0 / 0.6,
                       0.05, 1.0, /*max_level=*/1);
  Cylinder3dCase single(MPI_COMM_WORLD, CylinderStlPath(), /*n=*/8, 1.0 / 0.6,
                        0.05, 1.0, /*max_level=*/0);

  const BoundingBox d = CylinderDomain();
  const double dx = (d.x_max - d.x_min) / 8.0;
  const double dy = (d.y_max - d.y_min) / 8.0;
  const double dz = (d.z_max - d.z_min) / 8.0;

  // Compare multi vs single over all 8^3 cell centres; return mismatch count
  // and max errors. Called at three checkpoints to localize the divergence.
  auto compare = [&](const char* label) {
    double max_rho = 0.0;
    double max_ux = 0.0;
    int bad = 0;
    for (int iz = 0; iz < 8; ++iz) {
      for (int iy = 0; iy < 8; ++iy) {
        for (int ix = 0; ix < 8; ++ix) {
          const double x = d.x_min + (ix + 0.5) * dx;
          const double y = d.y_min + (iy + 0.5) * dy;
          const double z = d.z_min + (iz + 0.5) * dz;
          const auto [rm, um] = multi.sample_rho_ux(x, y, z);
          const auto [rs, us] = single.sample_rho_ux(x, y, z);
          max_rho = std::max(max_rho, std::abs(rm - rs));
          max_ux = std::max(max_ux, std::abs(um - us));
          if (std::abs(rm - rs) > 1e-6 || std::abs(um - us) > 1e-6) {
            ++bad;
          }
        }
      }
    }
    std::cout << "[NoEdgeArtifact] " << label
              << ": bad=" << bad << " max|dρ|=" << max_rho
              << " max|dux|=" << max_ux << "\n";
    return bad;
  };

  // Checkpoint 0: before any advance. If this mismatches, multi and single
  // initialize differently (voxelization/BC/Bouzidi) -- an approach flaw, NOT
  // ②. If it matches, the init is identical and any later divergence is purely
  // from time-stepping (② edge-ghost).
  const int bad0 = compare("step0");
  EXPECT_EQ(bad0, 0) << "multi vs single differ at init -- approach flaw";

  multi.advance_steps(1);
  single.advance_steps(1);
  const int bad1 = compare("step1");

  // Categorize the step-1 mismatched cells to localize the cause:
  //   kind[fluid,boundary,solid]   -- boundary-heavy => Bouzidi, fluid => ghost
  //   edge_count[interior,face,edge,corner] -- edge/corner => ②, face => face xchg
  {
    int by_kind[3] = {0, 0, 0};  // kFluid=0, kSolid=1, kBoundary=2
    int by_ec[4] = {0, 0, 0, 0};  // interior, face(1), edge(2), corner(3)
    for (int iz = 0; iz < 8; ++iz) {
      for (int iy = 0; iy < 8; ++iy) {
        for (int ix = 0; ix < 8; ++ix) {
          const double x = d.x_min + (ix + 0.5) * dx;
          const double y = d.y_min + (iy + 0.5) * dy;
          const double z = d.z_min + (iz + 0.5) * dz;
          const auto dm = multi.sample_diag(x, y, z);
          const auto ds = single.sample_diag(x, y, z);
          if (std::abs(dm.rho - ds.rho) > 1e-6 ||
              std::abs(dm.ux - ds.ux) > 1e-6) {
            by_kind[static_cast<int>(dm.kind) % 3]++;
            const int ec = dm.edge_count < 0 ? 0
                           : dm.edge_count > 3 ? 3
                                               : dm.edge_count;
            by_ec[ec]++;
          }
        }
      }
    }
    std::cout << "[NoEdgeArtifact] step1 breakdown: kind[fluid,boundary,solid]=["
              << by_kind[0] << "," << by_kind[1] << "," << by_kind[2]
              << "] edge_count[interior,face,edge,corner]=[" << by_ec[0] << ","
              << by_ec[1] << "," << by_ec[2] << "," << by_ec[3] << "]\n";
  }

  const int remaining = CylinderSteps(50) - 1;
  multi.advance_steps(remaining);
  single.advance_steps(remaining);
  const int badN = compare("stepN");

  EXPECT_EQ(badN, 0)
      << "multi-block differs from single-block at " << badN
      << " cells (② edge-ghost); max|dρ|/|dux| see checkpoint log above";
}

}  // namespace
}  // namespace octlb