#include <gtest/gtest.h>
#include <mpi.h>

#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/bc_installer.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/outflow_bc.h"
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

// The zero-gradient outflow kernel (still a stateless helper) copies the
// interior boundary populations into the outer ghost verbatim.
TEST(OutflowBc, ZeroGradientGhostEqualsInterior) {
  constexpr int Q = Descriptor::q;
  std::vector<double> interior(Q);
  std::vector<double> ghost(Q, 0.0);
  for (int i = 0; i < Q; ++i) {
    interior[i] = static_cast<double>(i + 1) * 0.1;
  }
  boundary::ApplyOutflowGhost<double, Descriptor>(ghost.data(), interior.data());
  for (int i = 0; i < Q; ++i) {
    EXPECT_DOUBLE_EQ(ghost[i], interior[i]);
  }
}

// Per-cell dispatch (R4): the outlet is a pressure-Dirichlet FD cell with
// prescribed rho=1.0 (p=0). After collide + PostStream FD, the outlet boundary
// cell reconstructs to rho==1.0.
TEST(DomainBoundary, Outlet_PressureZero_ReconstructsRhoOne) {
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
  spec.face = FaceDir::kXMax;
  spec.type = DomainBcType::kInterpolatedPressure;
  spec.rho_target = 1.0;
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
  blocks[0].get(kN - 1, iy, iz).computeRhoU(rho, u);
  EXPECT_NEAR(rho, 1.0, 1e-12) << "pressure outlet must prescribe rho=1 (p=0)";
}

}  // namespace
}  // namespace octlb