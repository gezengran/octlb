#include <gtest/gtest.h>
#include <mpi.h>

#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/outflow_bc.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {
namespace {

using Lattice = BlockLattice<double, olb::descriptors::D3Q19<>>;
using Descriptor = olb::descriptors::D3Q19<>;

constexpr int kN = 4;

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

// Zero-gradient outflow copies the interior boundary populations into the outer
// ghost verbatim, so a uniform outflow is not reflected.
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

// The handler dispatches a kOutflow spec to the zero-gradient kernel on the
// matching tree face: the +x outer ghost equals the adjacent interior.
TEST(DomainBoundary, OutflowHandler_ZeroGradientNoReflection) {
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
  spec.face = FaceDir::kXMax;
  spec.type = DomainBcType::kOutflow;
  ConcreteDomainBoundaryHandler handler(blocks, pairs.tree_boundary_faces(),
                                       {spec}, kN, kN, kN);
  handler.apply();

  for (int iy = 0; iy < kN; ++iy) {
    for (int iz = 0; iz < kN; ++iz) {
      const double* ghost = lat.populations_at_halo(kN + 1, iy + 1, iz + 1);
      const double* interior = lat.populations_at_halo(kN, iy + 1, iz + 1);
      for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
        EXPECT_DOUBLE_EQ(ghost[iPop], interior[iPop])
            << "outflow ghost must mirror interior at iy=" << iy << " iz=" << iz;
      }
    }
  }
}

}  // namespace
}  // namespace octlb