#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {
namespace {

using Lattice = BlockLattice<double, olb::descriptors::D3Q19<>>;
using Descriptor = olb::descriptors::D3Q19<>;

constexpr int kN = 4;

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

bool PopPointsIntoDomainFromGhost(int iPop, FaceDir dir) {
  const int cx = olb::descriptors::c<Descriptor>(iPop, 0);
  const int cy = olb::descriptors::c<Descriptor>(iPop, 1);
  const int cz = olb::descriptors::c<Descriptor>(iPop, 2);
  switch (dir) {
    case FaceDir::kXMin:
      return cx > 0;
    case FaceDir::kXMax:
      return cx < 0;
    case FaceDir::kYMin:
      return cy > 0;
    case FaceDir::kYMax:
      return cy < 0;
    case FaceDir::kZMin:
      return cz > 0;
    case FaceDir::kZMax:
      return cz < 0;
  }
  return false;
}

void ComputeRhoU(const double* f, double* rho, double* u) {
  *rho = 1.0;
  u[0] = u[1] = u[2] = 0.0;
  for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    *rho += f[iPop];
    u[0] += f[iPop] * olb::descriptors::c<Descriptor>(iPop, 0);
    u[1] += f[iPop] * olb::descriptors::c<Descriptor>(iPop, 1);
    u[2] += f[iPop] * olb::descriptors::c<Descriptor>(iPop, 2);
  }
  u[0] /= *rho;
  u[1] /= *rho;
  u[2] /= *rho;
}

BlockCollection<Lattice> MakeSingleBlockLattice() {
  return BlockCollection<Lattice>(1, [](OctantId) {
    Lattice lat(kN, kN, kN, 1);
    const double u0[3] = {0.0, 0.0, 0.0};
    lat.initialize(1.0, u0);
    return lat;
  });
}

}  // namespace

TEST(DomainBoundary, NoSlip_FillsGhostAfterCollide) {
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
  spec.type = DomainBcType::kNoSlip;
  ConcreteDomainBoundaryHandler handler(blocks, pairs.tree_boundary_faces(), {spec},
                                      kN, kN, kN);
  handler.apply();

  for (int iy = 0; iy < kN; ++iy) {
    for (int iz = 0; iz < kN; ++iz) {
      const double* ghost = lat.populations_at_halo(0, iy + 1, iz + 1);
      const double* interior = lat.populations_at_halo(1, iy + 1, iz + 1);
      for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
        if (!PopPointsIntoDomainFromGhost(iPop, FaceDir::kXMin)) {
          continue;
        }
        const int opp = olb::descriptors::opposite<Descriptor>(iPop);
        EXPECT_DOUBLE_EQ(ghost[iPop], interior[opp]);
      }
    }
  }
}

TEST(DomainBoundary, MovingLid_ZouHe) {
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

  constexpr double kUwall = 0.1;
  DomainBcSpec spec;
  spec.face = FaceDir::kYMax;
  spec.type = DomainBcType::kMovingLid;
  spec.u_wall = {kUwall, 0.0, 0.0};
  ConcreteDomainBoundaryHandler handler(blocks, pairs.tree_boundary_faces(), {spec},
                                      kN, kN, kN);
  handler.apply();

  for (int ix = 0; ix < kN; ++ix) {
    for (int iz = 0; iz < kN; ++iz) {
      const double* ghost = lat.populations_at_halo(ix + 1, kN + 1, iz + 1);
      double rho = 0.0;
      double u[3] = {};
      ComputeRhoU(ghost, &rho, u);
      EXPECT_NEAR(u[0], kUwall, 1e-10);
      EXPECT_NEAR(u[1], 0.0, 1e-10);
      EXPECT_NEAR(u[2], 0.0, 1e-10);
    }
  }
}

}  // namespace octlb
