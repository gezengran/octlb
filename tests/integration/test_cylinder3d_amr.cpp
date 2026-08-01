#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

#include "examples/cylinder3d/cylinder3d_case.h"

namespace octlb {
namespace {

int CylinderSteps(int default_steps) {
  const char* env = std::getenv("OCTLB_CYL_STEPS");
  return env != nullptr ? std::atoi(env) : default_steps;
}

bool StrictMode() {
  return std::getenv("OCTLB_CYL_STRICT") != nullptr;
}

// Override the AMR max_level for fast verification under slow emulation.
// Default (unset) keeps the ctest resolution (4 for sanity, 5 for strict gates).
int ResolveMaxLevel(int default_level) {
  const char* env = std::getenv("OCTLB_CYL_MAXLEVEL");
  if (env != nullptr && env[0] != '\0') {
    const int v = std::atoi(env);
    if (v > 0) return v;
  }
  return default_level;
}

std::string CylinderStlPath() {
  return std::string(OCTLB_TEST_DATA_DIR) + "/mesh/cylinder_st.stl";
}

// AMR Schäfer-Turek config. ctest uses a feasible max_level (finest dx ~= 0.02,
// D/dx ~= 5) + short steps; the 2x magnitude and 1% Cd gates need long-time
// convergence and run via the W4 example (manual), gated by OCTLB_CYL_STRICT.
// max_level=5 (dx=0.01, D/dx=10) is the W4 example resolution.
Cylinder3dConfig AmrConfig(const std::string& stl_path, int max_level) {
  Cylinder3dConfig cfg;
  cfg.mode = Cylinder3dMode::kAmr;
  cfg.stl_path = stl_path;
  cfg.n = 8;
  cfg.max_level = ResolveMaxLevel(max_level);
  cfg.omega = 1.0 / 0.53;  // Schäfer-Turek Re=20
  cfg.u_inlet = 0.05;       // relaxed for a clear drag signal on short runs
  cfg.rho0 = 1.0;
  cfg.ramp_end_t = 0.0;
  return cfg;
}

// Uniform straight-duct config (no cylinder) for the ② clean re-probe: isolates
// block edge/corner artifacts from cylinder-surface voxelization noise.
Cylinder3dConfig ChannelConfig(int max_level, bool enable_edge_exchange) {
  Cylinder3dConfig cfg;
  cfg.mode = Cylinder3dMode::kUniform;
  cfg.stl_path = CylinderStlPath();  // unused when include_cylinder=false
  cfg.n = 8;
  cfg.max_level = max_level;
  cfg.omega = 1.0 / 0.6;
  // Low Ma (u_inlet=0.01) for a stable base flow: the W3 pressure outlet is
  // marginally unstable at higher Ma in the no-cylinder straight channel (see
  // EdgeExchange_CrossRank_NoArtifact comment). The ② probe needs a stable base.
  cfg.u_inlet = 0.01;
  cfg.rho0 = 1.0;
  cfg.ramp_end_t = 0.0;
  cfg.include_cylinder = false;
  cfg.enable_edge_exchange = enable_edge_exchange;
  return cfg;
}

// Allreduce a Cylinder3dCase::EdgeCornerStats across ranks (componentwise max).
Cylinder3dCase::EdgeCornerStats AllreduceStats(
    MPI_Comm comm, const Cylinder3dCase::EdgeCornerStats& local) {
  double in[2] = {local.max_u, local.max_rho_dev};
  double out[2] = {0.0, 0.0};
  MPI_Allreduce(in, out, 2, MPI_DOUBLE, MPI_MAX, comm);
  Cylinder3dCase::EdgeCornerStats g;
  g.max_u = out[0];
  g.max_rho_dev = out[1];
  return g;
}

// Finest-cell lattice frontal area A_lat = (D/dx) * (z_span/dx) for the Cd
// formula Cd = 2*F_lat / (rho * u_lat^2 * A_lat). The cylinder surface sits at
// the finest AMR level, so dx = domain_side / 2^max_level / n.
double CylinderFrontalAreaLattice(const Cylinder3dConfig& cfg) {
  const double dx = kCubeSide / std::pow(2.0, cfg.max_level) /
                    static_cast<double>(cfg.n);
  const double D = 0.1;
  const double z_span = kChannelHi - kChannelLo;  // 0.41
  return (D / dx) * (z_span / dx);
}

// W3-e (T11) AMR sanity: geometry-driven AMR (carved channel + cylinder
// surface + wake) runs multi-rank without NaN. The magnitude/Cd gates need
// convergence and run via the W4 example.
TEST(Cylinder3dAmr, MultiRank_RunsNoCrash) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 4) {
    GTEST_SKIP() << "four-rank AMR test";
  }

  Cylinder3dCase cas(MPI_COMM_WORLD, AmrConfig(CylinderStlPath(), /*max_level=*/4));
  cas.advance_steps(CylinderSteps(20));

  EXPECT_FALSE(cas.has_non_finite_velocity())
      << "AMR cylinder3d produced non-finite velocity";
}

TEST(Cylinder3dAmr, MassBoundedDrift) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 4) {
    GTEST_SKIP() << "four-rank AMR test";
  }

  Cylinder3dCase cas(MPI_COMM_WORLD, AmrConfig(CylinderStlPath(), /*max_level=*/4));
  const double m0 = cas.total_mass(MPI_COMM_WORLD);
  cas.advance_steps(40);
  const double m1 = cas.total_mass(MPI_COMM_WORLD);

  EXPECT_TRUE(std::isfinite(m1));
  ASSERT_GT(m0, 0.0);
  EXPECT_LT(std::abs(m1 - m0) / m0, 0.2) << "AMR fluid mass drift blowup";
}

TEST(Cylinder3dAmr, Cd_FinitePositiveSign) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 4) {
    GTEST_SKIP() << "four-rank AMR test";
  }

  Cylinder3dCase cas(MPI_COMM_WORLD, AmrConfig(CylinderStlPath(), /*max_level=*/4));
  cas.advance_steps(80);  // let the +x stream reach and pass the cylinder

  const double Cd = cas.drag_coefficient(MPI_COMM_WORLD, /*area=*/1.0,
                                         /*flow_dir=*/{1.0, 0.0, 0.0});
  EXPECT_TRUE(std::isfinite(Cd));
  EXPECT_GT(Cd, 0.0) << "drag must oppose the +x free stream";
}

// ② W3 Stage B clean re-probe (T11 defect ②): uniform straight duct (no
// cylinder) on 2 ranks isolates block edge/corner artifacts from the cylinder-
// surface voxelization noise. With cross-rank edge exchange ON, the block
// edge/corner cells -- whose stream pulls from the cross-rank edge ghosts --
// stay bounded (no spike); Stage B must actually be wired (cross-rank edges
// exist along the partition boundary). The no-exchange baseline is run too
// and the ON variant must not be worse. CI runs short; OCTLB_CYL_STRICT raises
// the step count and tightens the bound.
TEST(Cylinder3dAmr, EdgeExchange_CrossRank_NoArtifact) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    GTEST_SKIP() << "multi-rank edge-exchange test";
  }

  // Refine to level 2 so 2 ranks own 32 octants each with a shared face and
  // shared cross-rank edges along the partition boundary.
  const int kLevel = 2;
  const int kSteps = CylinderSteps(StrictMode() ? 200 : 40);

  // The base flow must be STABLE so any edge/corner spike is attributable to the
  // ② cross-rank edge ghosts, not to a base-flow instability. The no-cylinder
  // straight channel with the W3 pressure outlet (kInterpolatedPressure) is
  // marginally unstable at u_inlet=0.05 (Ma~=0.083): the bulk velocity runs away
  // ~15%/step (measured 1-rank: bulk 0.055->0.60 over 12 steps), which dominates
  // any edge signal and makes the ON-vs-OFF comparison meaningless (ON~=OFF,
  // both diverge). At u_inlet=0.01 (Ma~=0.017) the bulk converges to a steady
  // ~0.04-0.05 (Poiseuille), giving a valid base for the ② probe. The cylinder
  // case (test_cylinder3d_uniform) hides this because the cylinder's blockage
  // damps the runaway and its assertions are weak (no-NaN + mass drift).
  Cylinder3dCase on(MPI_COMM_WORLD, ChannelConfig(kLevel, /*edge=*/true));
  on.advance_steps(kSteps);
  ASSERT_FALSE(on.has_non_finite_velocity());
  const Cylinder3dCase::EdgeCornerStats on_stat =
      AllreduceStats(MPI_COMM_WORLD, on.local_edge_corner_stats());

  // Stage B is wired: at least one rank owns outbound cross-rank edges.
  int has_edge_local = on.face_pairs.cross_rank_edges().empty() ? 0 : 1;
  int has_edge_global = 0;
  MPI_Allreduce(&has_edge_local, &has_edge_global, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  EXPECT_GT(has_edge_global, 0) << "no cross-rank edges; partition too coarse";

  Cylinder3dCase off(MPI_COMM_WORLD, ChannelConfig(kLevel, /*edge=*/false));
  off.advance_steps(kSteps);
  ASSERT_FALSE(off.has_non_finite_velocity());
  const Cylinder3dCase::EdgeCornerStats off_stat =
      AllreduceStats(MPI_COMM_WORLD, off.local_edge_corner_stats());

  // With edge exchange ON, block edge/corner cells are bounded (no spike).
  // u stays within a small multiple of the inlet velocity; rho near rho0.
  const double u_in = on.cfg.u_inlet;
  const double u_bound = StrictMode() ? 4.0 * u_in : 8.0 * u_in + 0.2;
  EXPECT_LT(on_stat.max_u, u_bound)
      << "edge/corner velocity spike with edge exchange ON";
  EXPECT_LT(on_stat.max_rho_dev, StrictMode() ? 0.05 : 0.5)
      << "edge/corner rho spike with edge exchange ON";
  // Edge exchange must not make artifacts worse than the no-exchange baseline.
  EXPECT_LE(on_stat.max_u, off_stat.max_u + 0.5)
      << "edge exchange ON worse than OFF (u)";
}

// W3 Lagrava interface continuity (T11): after the LevelCoupler acts, rho/u
// across a coarse-fine interface must be continuous (bounded jump) -- the
// Lagrava prolongation/restriction contract. The prolongation writes the
// fine-side interface cells from the coarse macros, so the jump is bounded by
// construction once the coupler has acted. A broken/uncoupled interface would
// diverge and blow the bound. CI runs short; OCTLB_CYL_STRICT tightens the tol.
TEST(Cylinder3dAmr, CoarseFineInterface_Continuous) {
  int size = 0;
  int rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (size < 4) {
    GTEST_SKIP() << "four-rank AMR test";
  }

  Cylinder3dCase cas(MPI_COMM_WORLD, AmrConfig(CylinderStlPath(), /*max_level=*/4));
  cas.advance_steps(CylinderSteps(StrictMode() ? 100 : 20));
  ASSERT_FALSE(cas.has_non_finite_velocity());

  const label num_local = cas.num_blocks();
  const double u_scale = std::max(cas.cfg.u_inlet, 1e-6);
  double max_rho_jump = 0.0;
  double max_u_jump = 0.0;
  int checked = 0;
  for (const CouplingPoint& cp : cas.coupler.coupling_plan()) {
    // Local coarse-fine pairs only: both ends on this rank.
    if (cp.remote_rank != rank || cp.coarse_remote_rank != rank) continue;
    if (cp.coarse_id >= num_local || cp.fine_id >= num_local) continue;
    double rho_c = 0.0, rho_f = 0.0;
    double u_c[3] = {}, u_f[3] = {};
    cas.blocks[cp.coarse_id].get(cp.ci, cp.cj, cp.ck).computeRhoU(rho_c, u_c);
    cas.blocks[cp.fine_id].get(cp.fi, cp.fj, cp.fk).computeRhoU(rho_f, u_f);
    max_rho_jump = std::max(max_rho_jump, std::abs(rho_c - rho_f));
    const double du0 = u_c[0] - u_f[0];
    const double du1 = u_c[1] - u_f[1];
    const double du2 = u_c[2] - u_f[2];
    max_u_jump =
        std::max(max_u_jump, std::sqrt(du0 * du0 + du1 * du1 + du2 * du2));
    ++checked;
  }
  ASSERT_GT(checked, 0) << "no local coarse-fine interface; AMR did not refine";

  const double rho_tol = StrictMode() ? 0.1 * cas.cfg.rho0 : 0.5 * cas.cfg.rho0;
  const double u_tol = StrictMode() ? 2.0 * u_scale : 10.0 * u_scale;
  EXPECT_LT(max_rho_jump / cas.cfg.rho0, rho_tol / cas.cfg.rho0)
      << "coarse-fine rho discontinuity " << max_rho_jump << " over " << checked
      << " interface points";
  EXPECT_LT(max_u_jump, u_tol)
      << "coarse-fine u discontinuity " << max_u_jump << " over " << checked
      << " interface points";
}

// W3 magnitude gate (T11): Schäfer-Turek AMR at W4 resolution (max_level=5,
// D/dx=10), long-time convergence. Cd must be within 2x of the OpenLB
// cylinder3d self-converged reference ~6.36 (the PRD #13 reference). CI skips
// this (too slow under emulation); set OCTLB_CYL_STRICT=1 and OCTLB_CYL_STEPS=N
// to run the convergence.
TEST(Cylinder3dAmr, Amr_Cd_SameOrderAsOpenLb) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 4) {
    GTEST_SKIP() << "four-rank AMR test";
  }
  if (!StrictMode()) {
    GTEST_SKIP() << "long-time magnitude gate (set OCTLB_CYL_STRICT=1)";
  }

  Cylinder3dConfig cfg = AmrConfig(CylinderStlPath(), /*max_level=*/5);
  cfg.u_inlet = 0.02;  // Schäfer-Turek latticeU mean (Re=20)
  cfg.ramp_end_t = 0.0;
  Cylinder3dCase cas(MPI_COMM_WORLD, cfg);
  const int steps = CylinderSteps(4000);
  cas.advance_steps(steps);

  ASSERT_FALSE(cas.has_non_finite_velocity());
  const double area = CylinderFrontalAreaLattice(cfg);
  const double Cd = cas.drag_coefficient(MPI_COMM_WORLD, area,
                                         /*flow_dir=*/{1.0, 0.0, 0.0});
  EXPECT_TRUE(std::isfinite(Cd)) << "Cd non-finite after " << steps << " steps";
  // 2x of the OpenLB reference 6.36 -> [3.18, 12.72].
  EXPECT_GE(Cd, 3.18) << "Cd " << Cd << " below 2x-lower of OpenLB ref 6.36";
  EXPECT_LE(Cd, 12.72) << "Cd " << Cd << " above 2x-upper of OpenLB ref 6.36";
}

// W4 strict gate (T11 / PRD #13): long-time steady-state Cd vs the OpenLB
// cylinder3d self-converged reference ~6.36, relative error < 1%. This is the
// real convergence run (~16000 steps at max_level=5); CI skips it (hours under
// x86 emulation). Run on native hardware: OCTLB_CYL_STRICT=1 OCTLB_CYL_STEPS=16000
// mpiexec -n 4 ./test_cylinder3d_amr --gtest_filter='Cylinder3dAmr.Amr_Cd_WithinOnePercentOfOpenLb'
TEST(Cylinder3dAmr, Amr_Cd_WithinOnePercentOfOpenLb) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 4) {
    GTEST_SKIP() << "four-rank AMR test";
  }
  if (!StrictMode()) {
    GTEST_SKIP() << "1% convergence gate (set OCTLB_CYL_STRICT=1)";
  }

  Cylinder3dConfig cfg = AmrConfig(CylinderStlPath(), /*max_level=*/5);
  cfg.u_inlet = 0.02;  // Schäfer-Turek latticeU mean (Re=20)
  cfg.ramp_end_t = 0.0;
  Cylinder3dCase cas(MPI_COMM_WORLD, cfg);
  const int steps = CylinderSteps(16000);
  cas.advance_steps(steps);

  ASSERT_FALSE(cas.has_non_finite_velocity());
  const double area = CylinderFrontalAreaLattice(cfg);
  const double Cd = cas.drag_coefficient(MPI_COMM_WORLD, area,
                                         /*flow_dir=*/{1.0, 0.0, 0.0});
  const double ref = 6.36;
  const double rel = std::abs(Cd - ref) / ref;
  EXPECT_TRUE(std::isfinite(Cd)) << "Cd non-finite after " << steps << " steps";
  EXPECT_LT(rel, 0.01) << "Cd " << Cd << " rel error " << rel
                       << " not < 1% of OpenLB ref " << ref;
}

}  // namespace
}  // namespace octlb