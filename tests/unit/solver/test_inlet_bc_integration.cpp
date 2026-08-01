#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <memory>
#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/ghost_schedule.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/bc_installer.h"
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

// Per-cell dispatch (R4): the inlet is a velocity-Dirichlet FD cell whose
// prescribed u comes from the inlet_field (Poiseuille profile). After collide +
// PostStream FD, the inlet boundary cell's velocity matches the profile.
TEST(InletBcIntegration, FdAppliesPoiseuilleToInletCell) {
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

  DomainBcSpec spec;
  spec.face = FaceDir::kXMin;
  spec.type = DomainBcType::kInterpolatedVelocity;
  spec.inlet_field = std::make_shared<PoiseuilleInletProfile>(
      FaceDir::kXMin, kN, kN, kPeak, std::array<double, 3>{1.0, 0.0, 0.0},
      /*ramp_end_t=*/0.0);
  const std::vector<DomainBcSpec> specs = {spec};
  bc::StampTreeBoundaryCells(blocks, pairs.tree_boundary_faces(), specs, kN,
                              kN, kN);

  ConcreteDomainBoundaryHandler handler(blocks, pairs.tree_boundary_faces(),
                                        specs, kN, kN, kN, /*omega=*/1.0);
  handler.set_time(0.0);
  handler.collide_interleaved_with(lat, nullptr, 1.0, false);
  handler.apply_post_stream();

  // Flat inlet-face cells (interior of the face, not the wall edges/corners)
  // are kVelocityDirichlet and prescribe the Poiseuille profile. Edge/corner
  // cells belong to the walls and are not asserted here.
  for (int iy = 1; iy < kN - 1; ++iy) {
    for (int iz = 1; iz < kN - 1; ++iz) {
      double rho = 0.0;
      double u[3] = {};
      lat.get(0, iy, iz).computeRhoU(rho, u);
      const double expected = kPeak * Parab(iy, kN) * Parab(iz, kN);
      EXPECT_NEAR(u[0], expected, 1e-9)
          << "inlet cell ux at iy=" << iy << " iz=" << iz;
      EXPECT_NEAR(u[1], 0.0, 1e-9);
      EXPECT_NEAR(u[2], 0.0, 1e-9);
    }
  }
}

// W3-b (遗留 3c): the FD InterpolatedVelocity path must thread the handler's
// current_time into the inlet_field -- PrescribedBoundaryU previously hardcoded
// t=0.0, so a ramped Poiseuille inlet was stuck at u=0 (Ramp(0)=0). With ramp_end_t
// saturated (t >= ramp_end_t -> ramp=1), the inlet cell must prescribe the full
// peak*parab*parab, proving the FD path honours the real time.
TEST(InletBcIntegration, FdHonorsPoiseuilleRampTime) {
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

  constexpr double kRampEnd = 10.0;
  DomainBcSpec spec;
  spec.face = FaceDir::kXMin;
  spec.type = DomainBcType::kInterpolatedVelocity;
  spec.inlet_field = std::make_shared<PoiseuilleInletProfile>(
      FaceDir::kXMin, kN, kN, kPeak, std::array<double, 3>{1.0, 0.0, 0.0},
      /*ramp_end_t=*/kRampEnd);
  const std::vector<DomainBcSpec> specs = {spec};
  bc::StampTreeBoundaryCells(blocks, pairs.tree_boundary_faces(), specs, kN,
                              kN, kN);

  ConcreteDomainBoundaryHandler handler(blocks, pairs.tree_boundary_faces(),
                                        specs, kN, kN, kN, /*omega=*/1.0);
  // Saturate the ramp: t == ramp_end_t -> smootherstep(1) = 1 (full velocity).
  handler.set_time(kRampEnd);
  handler.collide_interleaved_with(lat, nullptr, 1.0, false);
  handler.apply_post_stream();

  for (int iy = 1; iy < kN - 1; ++iy) {
    for (int iz = 1; iz < kN - 1; ++iz) {
      double rho = 0.0;
      double u[3] = {};
      lat.get(0, iy, iz).computeRhoU(rho, u);
      const double expected = kPeak * Parab(iy, kN) * Parab(iz, kN);
      EXPECT_NEAR(u[0], expected, 1e-9)
          << "ramped inlet cell ux at iy=" << iy << " iz=" << iz;
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