// W4 (T11) cylinder3d example: Schäfer-Turek flow past a cylinder, geometry-
// driven AMR, Bouzidi surface BC, momentum-exchange drag. Runs to steady state
// and writes VTK + a Cd time series CSV for manual OpenLB cross-check (PRD
// integration test #13: Cd vs the OpenLB cylinder3d self-converged reference
// ~6.36, relative error < 1% at convergence -- a long-time native-CPU run;
// under x86 emulation this host is too slow for the full convergence).
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <mpi.h>

#include "examples/cylinder3d/cylinder3d_case.h"

namespace {

std::string Opt(int argc, char** argv, const std::string& flag,
                const std::string& dflt) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (std::string(argv[i]) == flag) return argv[i + 1];
  }
  return dflt;
}

int ResolveSteps(int argc, char** argv) {
  const std::string s = Opt(argc, argv, "--steps", "");
  if (!s.empty()) return std::atoi(s.c_str());
  if (const char* env = std::getenv("OCTLB_CYL_STEPS")) return std::atoi(env);
  // Schäfer-Turek MAX_PHYS_T=16, dt=0.001 -> ~16000 lattice steps. Long-time
  // steady state; the OpenLB reference Cd converges around physT=16.
  return 16000;
}

// Finest-cell lattice frontal area A_lat = (D/dx)*(z_span/dx) for Cd =
// 2*F_lat/(rho*u_lat^2*A_lat). The cylinder surface sits at the finest AMR
// level, so dx = domain_side / 2^max_level / n. See test_cylinder3d_amr.
double FrontalAreaLattice(const octlb::Cylinder3dConfig& cfg) {
  const double dx = octlb::kCubeSide / std::pow(2.0, cfg.max_level) /
                    static_cast<double>(cfg.n);
  const double D = 0.1;
  const double z_span = octlb::kChannelHi - octlb::kChannelLo;  // 0.41
  return (D / dx) * (z_span / dx);
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const std::string stl_path = Opt(
      argc, argv, "--stl",
      std::string(OCTLB_TEST_DATA_DIR) + "/mesh/cylinder_st.stl");
  const int max_level = std::atoi(Opt(argc, argv, "--max-level", "5").c_str());
  const int vtk_interval =
      std::atoi(Opt(argc, argv, "--vtk-interval", "2000").c_str());
  const int cd_interval =
      std::atoi(Opt(argc, argv, "--cd-interval", "100").c_str());
  const std::string output_dir = Opt(argc, argv, "--output", "cylinder3d_vtk");
  const int steps = ResolveSteps(argc, argv);

  octlb::Cylinder3dConfig cfg;
  cfg.mode = octlb::Cylinder3dMode::kAmr;
  if (std::getenv("OCTLB_CYL_UNIFORM")) cfg.mode = octlb::Cylinder3dMode::kUniform;
  cfg.stl_path = stl_path;
  cfg.n = 8;
  cfg.max_level = max_level;
  cfg.omega = 1.0 / 0.53;   // Schäfer-Turek Re=20, tau=0.53
  cfg.u_inlet = 0.02;       // Schäfer-Turek U_max peak (Re=20; duct mean = 4/9 * this)
  cfg.rho0 = 1.0;
  cfg.ramp_end_t = 0.0;
  cfg.charU = 0.2;
  cfg.poiseuille_inlet = true;  // W4: Schäfer-Turek Poiseuille duct inlet
  cfg.refine_channel_to_finest = true;  // T11: resolve inlet+channel at finest so Poiseuille develops
  if (std::getenv("OCTLB_CYL_NO_EDGE")) cfg.enable_edge_exchange = false;

  if (rank == 0) {
    std::cout << "cylinder3d: max_level=" << max_level << " n=" << cfg.n
              << " omega=" << cfg.omega << " u_inlet=" << cfg.u_inlet
              << " ranks=" << size << " steps=" << steps
              << " OpenLB_ref_Cd~6.36" << std::endl;
  }

  octlb::Cylinder3dCase cas(MPI_COMM_WORLD, cfg);

  // T11 MEM timing: the correct MEM reads the pre-stream post-collide snapshot
  // for the outgoing f_i, so the time loop must snapshot every block right before
  // stream() on any step where the drag will be computed. Snapshots on drag
  // steps only (cheap: one extra memcpy per block per drag step).
  const bool mem_timing =
      std::getenv("OCTLB_CYL_MEM_TIMING") != nullptr;

  if (std::getenv("OCTLB_CYL_AUDIT")) {
    cas.audit_spurious_fluid(MPI_COMM_WORLD, std::cout);
  }

  const double area = FrontalAreaLattice(cfg);
  std::ofstream cd_csv;
  if (rank == 0) {
    cd_csv.open(output_dir + "/cd.csv");
    cd_csv << "iT,Cd\n";
  }

  for (int step = 1; step <= steps; ++step) {
    // advance_steps() enables the post-collide snapshot itself, so the MEM in
    // drag_coefficient() reads the correct (pre-stream outgoing + post-stream
    // bounced) timing on every drag step.
    cas.advance_steps(1);
    if (step % cd_interval == 0) {
      const double Cd = cas.drag_coefficient(MPI_COMM_WORLD, area,
                                             /*flow_dir=*/{1.0, 0.0, 0.0});
      const auto fs = cas.flow_stats(MPI_COMM_WORLD);
      if (rank == 0) {
        cd_csv << step << ',' << Cd << '\n' << std::flush;
        std::cout << "iT=" << step << " Cd=" << Cd
                  << " mean_ux=" << fs.mean_ux
                  << " mean_rho=" << fs.mean_rho
                  << " max_umag=" << fs.max_umag
                  << " n=" << fs.n << std::endl;
      }
      if (mem_timing) {
        // Legacy 2*live MEM (wrong timing) for before/after comparison.
        const double Cd_legacy = cas.drag_coefficient(MPI_COMM_WORLD, area,
                                                       /*flow_dir=*/{1.0, 0.0, 0.0},
                                                       /*legacy=*/true);
        if (rank == 0) {
          std::cout << "[mem-timing] Cd_correct=" << Cd
                    << " Cd_legacy(2*live)=" << Cd_legacy << '\n';
        }
      }
      if (std::getenv("OCTLB_CYL_MAXU")) {
        cas.report_max_umag_cell(MPI_COMM_WORLD, std::cout);
      }
      if (std::getenv("OCTLB_CYL_DRAG_DEBUG")) {
        cas.drag_breakdown(MPI_COMM_WORLD, std::cout);
      }
      if (std::getenv("OCTLB_CYL_FLUX")) {
        const double flux = cas.inlet_adjacent_mass_flux(MPI_COMM_WORLD);
        if (rank == 0) {
          std::cout << "[flux] inlet-adjacent mass flux=" << flux << '\n';
        }
      }
    }
    if (step % vtk_interval == 0) {
      cas.write_vtk_timestep(MPI_COMM_WORLD, step, output_dir);
    }
  }

  cas.write_vtk_timestep(MPI_COMM_WORLD, steps, output_dir);
  if (rank == 0) {
    if (cas.has_non_finite_velocity()) {
      std::cerr << "warning: non-finite velocity after " << steps << " steps\n";
    }
    cd_csv.close();
    std::cout << "Wrote VTK + cd.csv to " << output_dir << std::endl;
  }

  MPI_Finalize();
  return 0;
}