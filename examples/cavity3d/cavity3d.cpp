#include <cstdlib>
#include <iostream>
#include <string>

#include <mpi.h>

#include "examples/cavity3d/cavity3d_case.h"

namespace {

std::string OutputDir(int argc, char** argv) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (std::string(argv[i]) == "--output") {
      return argv[i + 1];
    }
  }
  return "cavity3d_vtk";
}

// OpenLB cavity3d sets iT_max=getLatticeTime(100) but ValueTracer stops at
// iT=5269; running the full 30000 steps without convergence check diverges.
int ResolveSteps(int argc, char** argv) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (std::string(argv[i]) == "--steps") {
      return std::atoi(argv[i + 1]);
    }
  }
  if (const char* env = std::getenv("OCTLB_CAVITY_STEPS")) {
    return std::atoi(env);
  }
  return octlb::kOpenLbCavity3dConvergedSteps;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const octlb::UnitConverter converter = octlb::UnitConverter::OpenLbCavity3dDefaults();
  octlb::Cavity3dCase cavity(MPI_COMM_WORLD, converter);
  const int steps = ResolveSteps(argc, argv);

  if (rank == 0) {
    std::cout << "cavity3d: N=" << converter.resolution()
              << " tau=" << converter.lattice_relaxation_time()
              << " Re=" << converter.reynolds()
              << " iT=" << steps << std::endl;
  }

  cavity.advance_steps(steps);

  const std::string output_dir = OutputDir(argc, argv);
  cavity.write_vtk_timestep(MPI_COMM_WORLD, steps, output_dir);
  if (rank == 0) {
    if (cavity.has_non_finite_velocity()) {
      std::cerr << "warning: non-finite velocity in lattice after " << steps
                << " steps\n";
    }
    cavity.write_centerline_csv(output_dir + "/centerline.csv");
    std::cout << "Wrote VTK to " << output_dir << std::endl;
  }

  MPI_Finalize();
  return 0;
}
