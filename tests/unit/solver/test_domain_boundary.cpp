#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/bc_installer.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {
namespace {

using Lattice = BlockLattice<double, olb::descriptors::D3Q19<>>;
using Descriptor = olb::descriptors::D3Q19<>;

constexpr int kN = 8;

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

BlockCollection<Lattice> MakeSingleBlockLattice() {
  return BlockCollection<Lattice>(1, [](OctantId) {
    Lattice lat(kN, kN, kN, 1);
    const double u0[3] = {0.0, 0.0, 0.0};
    lat.initialize(1.0, u0);
    return lat;
  });
}

// Per-cell dispatch (R4): a no-slip face is a velocity-Dirichlet FD cell with
// prescribed u=0. After collide + PostStream FD, the boundary cell's velocity
// is zero (no-slip), with no ghost-fill involved.
TEST(DomainBoundary, NoSlip_FdBoundaryCellPrescribesZeroU) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();
  const FacePairList pairs(forest);
  auto blocks = MakeSingleBlockLattice();

  DomainBcSpec spec;
  spec.face = FaceDir::kXMin;
  spec.type = DomainBcType::kInterpolatedVelocity;  // FD no-slip (u=0)
  const std::vector<DomainBcSpec> specs = {spec};
  bc::StampTreeBoundaryCells(blocks, pairs.tree_boundary_faces(), specs, kN,
                              kN, kN);

  ConcreteDomainBoundaryHandler handler(blocks, pairs.tree_boundary_faces(),
                                        specs, kN, kN, kN);
  handler.collide_interleaved_with(blocks[0], nullptr, 1.0, false);
  handler.apply_post_stream();

  const int iy = kN / 2;
  const int iz = kN / 2;
  double rho = 0.0;
  double u[3] = {};
  blocks[0].get(0, iy, iz).computeRhoU(rho, u);
  EXPECT_NEAR(u[0], 0.0, 1e-12) << "no-slip wall cell must prescribe u=0";
  EXPECT_NEAR(u[1], 0.0, 1e-12);
  EXPECT_NEAR(u[2], 0.0, 1e-12);
}

// Per-cell dispatch (R4): the moving lid is a velocity-Dirichlet FD cell with
// prescribed u_lid. After collide + PostStream FD, the interior top-lid cell's
// velocity matches u_lid.
TEST(DomainBoundary, MovingLid_FdBoundaryCellPrescribesLidU) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();
  const FacePairList pairs(forest);
  auto blocks = MakeSingleBlockLattice();

  constexpr double kUwall = 0.1;
  DomainBcSpec spec;
  spec.face = FaceDir::kYMax;
  spec.type = DomainBcType::kInterpolatedVelocity;
  spec.u_wall = {kUwall, 0.0, 0.0};
  const std::vector<DomainBcSpec> specs = {spec};
  bc::StampTreeBoundaryCells(blocks, pairs.tree_boundary_faces(), specs, kN,
                              kN, kN);

  ConcreteDomainBoundaryHandler handler(blocks, pairs.tree_boundary_faces(),
                                        specs, kN, kN, kN);
  handler.collide_interleaved_with(blocks[0], nullptr, 1.0, false);
  handler.apply_post_stream();

  const int ix = kN / 2;
  const int iz = kN / 2;
  double rho = 0.0;
  double u[3] = {};
  blocks[0].get(ix, kN - 1, iz).computeRhoU(rho, u);
  EXPECT_NEAR(u[0], kUwall, 1e-12) << "moving-lid cell must prescribe u_lid";
  EXPECT_NEAR(u[1], 0.0, 1e-12);
  EXPECT_NEAR(u[2], 0.0, 1e-12);
}

}  // namespace
}  // namespace octlb