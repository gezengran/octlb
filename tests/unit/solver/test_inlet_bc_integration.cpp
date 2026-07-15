#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <memory>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/ghost_schedule.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/inlet_velocity_field.h"
#include "src/solver/lbm/domain_boundary_handler.h"
#include "src/solver/lbm/level_coupler.h"
#include "src/solver/lbm/time_loop/time_loop.h"

namespace octlb {
namespace {

using Lattice = BlockLattice<double, olb::descriptors::D3Q19<>>;
using Descriptor = olb::descriptors::D3Q19<>;
using namespace octlb::boundary;

// In-plane extent is odd so a cell centre sits at the cross-section peak.
constexpr int kN = 7;
constexpr double kPeak = 0.045;

BoundingBox UnitCubeDomain() { return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}; }

BlockCollection<Lattice> MakeSingleBlockLattice() {
  return BlockCollection<Lattice>(1, [](OctantId) {
    Lattice lat(kN, kN, kN, 1);
    const double u0[3] = {0.0, 0.0, 0.0};
    lat.initialize(1.0, u0);
    return lat;
  });
}

double Parab(int i, int n) {
  const double c = (i + 0.5) / static_cast<double>(n);
  return 4.0 * c * (1.0 - c);
}

// Legacy (non-boundary-lattice) velocity path: the handler fills the inlet
// ghost so that the ghost cell's velocity matches the prescribed per-cell
// Poiseuille profile. Verifies inlet_field reaches the BC, not just the lookup
// helper.
TEST(InletBcIntegration, LegacyAppliesPoiseuilleToGhost) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();
  const FacePairList pairs(forest);
  auto blocks = MakeSingleBlockLattice();
  Lattice& lat = blocks[0];
  lat.collide(1.0);

  DomainBcSpec spec;
  spec.face = FaceDir::kXMin;
  spec.type = DomainBcType::kInterpolatedVelocity;
  spec.inlet_field = std::make_shared<PoiseuilleInletProfile>(
      FaceDir::kXMin, kN, kN, kPeak, std::array<double, 3>{1.0, 0.0, 0.0},
      /*ramp_end_t=*/0.0);

  ConcreteDomainBoundaryHandler handler(blocks, pairs.tree_boundary_faces(),
                                       {spec}, kN, kN, kN,
                                       /*omega=*/1.0, /*boundary_lattice_mode=*/false);
  handler.set_time(0.0);
  handler.apply();

  for (int iy = 0; iy < kN; ++iy) {
    for (int iz = 0; iz < kN; ++iz) {
      const double* ghost = lat.populations_at_halo(0, iy + 1, iz + 1);
      CellProxy<double, Descriptor> cell(const_cast<double*>(ghost));
      double rho = 0.0;
      double u[3] = {};
      cell.computeRhoU(rho, u);
      const double expected = kPeak * Parab(iy, kN) * Parab(iz, kN);
      EXPECT_NEAR(u[0], expected, 1e-9)
          << "inlet ghost ux at iy=" << iy << " iz=" << iz;
      EXPECT_NEAR(u[1], 0.0, 1e-9);
      EXPECT_NEAR(u[2], 0.0, 1e-9);
    }
  }
}

// TimeLoop threads time to the boundary handler so an inlet_field's ramp can
// respond to the advancing simulation. Verifies the plumbing (set_time is
// called with an incrementing lattice step per advance_one); the ramp shape
// itself is covered by the profile tests and the per-cell application by the
// legacy ghost test above.
TEST(InletBcIntegration, TimeLoopThreadsTimeToBoundaryHandler) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();
  const FacePairList pairs(forest);
  BlockCollection<Lattice> blocks(forest.local_num_octants(), [](OctantId) {
    Lattice lat(kN, kN, kN, 1);
    const double u0[3] = {0.0, 0.0, 0.0};
    lat.initialize(1.0, u0);
    return lat;
  });
  GhostSchedule<Lattice> ghosts(MPI_COMM_WORLD, pairs, blocks, kN, kN, kN);
  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       /*omega=*/1.0);
  NoOpDomainBoundaryHandler domain_bc;
  TimeLoop loop(forest, blocks, ghosts, coupler, domain_bc, /*omega=*/1.0);

  // current_time() is the 0-based step index used in the most recent advance
  // (OpenLB iT semantics: the first step uses t=0 so the ramp starts at zero).
  EXPECT_NEAR(domain_bc.current_time(), 0.0, 1e-12);
  loop.advance_one();
  EXPECT_NEAR(domain_bc.current_time(), 0.0, 1e-12) << "step 0 used t=0";
  loop.advance_one();
  EXPECT_NEAR(domain_bc.current_time(), 1.0, 1e-12) << "step 1 used t=1";
  loop.advance_one();
  EXPECT_NEAR(domain_bc.current_time(), 2.0, 1e-12) << "step 2 used t=2";
}

}  // namespace
}  // namespace octlb