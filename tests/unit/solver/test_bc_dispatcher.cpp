#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/bc_dispatcher.h"
#include "src/solver/lbm/bc_kind.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/interpolated_velocity.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {
namespace {

using T = double;
using Descriptor = olb::descriptors::D3Q19<>;
using Lattice = BlockLattice<T, Descriptor>;

constexpr int kQ = Descriptor::q;
constexpr int kN = 8;

BoundingBox UnitCubeDomain() { return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}; }

BlockCollection<Lattice> MakeSingleBlockLattice() {
  return BlockCollection<Lattice>(1, [](OctantId) {
    Lattice lat(kN, kN, kN, 1);
    const double u0[3] = {0.0, 0.0, 0.0};
    lat.initialize(1.0, u0);
    return lat;
  });
}

// Per-cell BC dispatch (R0): a cell stamped kBounceBack must collide as a
// full-way bounce-back reflection -- f[i] <- f[opposite(i)] -- independent of
// any global mode flag. Tests behavior (the reflected populations), not the
// internal switch structure.
TEST(BcDispatcher, BounceBack_CollidesAsBounceBack) {
  Lattice lat(4, 4, 4, 1);
  const std::array<T, 3> u0{T{0}, T{0}, T{0}};
  lat.initialize(T{1}, u0.data());

  lat.set_bc_kind(1, 1, 1, BcKind::kBounceBack);
  ASSERT_EQ(lat.bc_kind(1, 1, 1), BcKind::kBounceBack);

  auto cell = lat.get(1, 1, 1);
  std::vector<T> before(static_cast<std::size_t>(kQ));
  for (int iPop = 0; iPop < kQ; ++iPop) {
    cell[iPop] = T{0.01} * static_cast<T>(iPop + 1) - T{0.05};
    before[static_cast<std::size_t>(iPop)] = cell[iPop];
  }

  const std::vector<DomainBcSpec> specs;  // unused by the bounce-back arm
  BcDispatcher::collide(lat, 1, 1, 1, 4, 4, 4, /*omega=*/1.0, specs);

  cell = lat.get(1, 1, 1);
  for (int iPop = 0; iPop < kQ; ++iPop) {
    const int opp = olb::descriptors::opposite<Descriptor>(iPop);
    EXPECT_DOUBLE_EQ(cell[iPop], before[static_cast<std::size_t>(opp)])
        << "iPop=" << iPop << " must equal pre-collision opposite";
  }
}

// R1: the handler collides per cell via BcDispatcher, so kBounceBack and
// kVelocityDirichlet cells in the same block each collide per their own kind
// without crosstalk. In R0 the handler's uniform Dirichlet path skipped
// kBounceBack cells (no reflection), so this is red until dispatch is wired.
TEST(BcDispatcher, MixedKinds_PerCellDispatch) {
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

  // Mixed kinds in one block: an interior bounce-back cell and a face
  // velocity-Dirichlet cell.
  const int bb_ix = 1, bb_iy = 1, bb_iz = 1;
  const int vd_ix = 0, vd_iy = 4, vd_iz = 4;
  lat.set_bc_kind(bb_ix, bb_iy, bb_iz, BcKind::kBounceBack);
  lat.set_bc_kind(vd_ix, vd_iy, vd_iz, BcKind::kVelocityDirichlet);

  auto plant = [](Lattice& l, int ix, int iy, int iz, std::vector<T>& snap) {
    auto c = l.get(ix, iy, iz);
    snap.resize(static_cast<std::size_t>(kQ));
    for (int iPop = 0; iPop < kQ; ++iPop) {
      c[iPop] = T{0.01} * static_cast<T>((ix * 64 + iy * 8 + iz) + iPop + 1) -
                T{0.03};
      snap[static_cast<std::size_t>(iPop)] = c[iPop];
    }
  };
  std::vector<T> bb_before;
  std::vector<T> vd_before;
  plant(lat, bb_ix, bb_iy, bb_iz, bb_before);
  plant(lat, vd_ix, vd_iy, vd_iz, vd_before);

  // Reference: what the velocity-Dirichlet cell should collide to (the arm
  // delegates to CollideDirichletBoundaryCellAt). Computed on a pristine copy
  // so the in-place handler dispatch can be compared against it.
  Lattice vd_ref(kN, kN, kN, 1);
  {
    const double u0[3] = {0.0, 0.0, 0.0};
    vd_ref.initialize(1.0, u0);
    vd_ref.set_bc_kind(vd_ix, vd_iy, vd_iz, BcKind::kVelocityDirichlet);
    auto c = vd_ref.get(vd_ix, vd_iy, vd_iz);
    for (int iPop = 0; iPop < kQ; ++iPop) {
      c[iPop] = vd_before[static_cast<std::size_t>(iPop)];
    }
  }
  DomainBcSpec spec;
  spec.face = FaceDir::kXMin;
  spec.type = DomainBcType::kInterpolatedVelocity;
  const std::vector<DomainBcSpec> specs = {spec};
  boundary::CollideDirichletBoundaryCellAt<T, Descriptor, Lattice>(
      vd_ref, vd_ix, vd_iy, vd_iz, kN, kN, kN, /*omega=*/1.0, specs,
      /*t=*/0.0, nullptr);
  std::vector<T> vd_expected(static_cast<std::size_t>(kQ));
  {
    auto c = vd_ref.get(vd_ix, vd_iy, vd_iz);
    for (int iPop = 0; iPop < kQ; ++iPop) {
      vd_expected[static_cast<std::size_t>(iPop)] = c[iPop];
    }
  }

  ConcreteDomainBoundaryHandler handler(blocks, pairs.tree_boundary_faces(),
                                        specs, kN, kN, kN, /*omega=*/1.0);
  handler.collide_interleaved_with(lat, /*rho_stats=*/nullptr,
                                   /*average_rho=*/1.0,
                                   /*use_const_rho_bgk=*/false);

  // kBounceBack cell: full-way reflection f[i] <- before[opp].
  auto bb_cell = lat.get(bb_ix, bb_iy, bb_iz);
  for (int iPop = 0; iPop < kQ; ++iPop) {
    const int opp = olb::descriptors::opposite<Descriptor>(iPop);
    EXPECT_DOUBLE_EQ(bb_cell[iPop], bb_before[static_cast<std::size_t>(opp)])
        << "kBounceBack cell must reflect under per-cell dispatch; iPop="
        << iPop;
  }

  // kVelocityDirichlet cell: collides to the same Dirichlet result the arm
  // delegates to, independent of the neighboring kBounceBack cell.
  auto vd_cell = lat.get(vd_ix, vd_iy, vd_iz);
  for (int iPop = 0; iPop < kQ; ++iPop) {
    EXPECT_NEAR(vd_cell[iPop], vd_expected[static_cast<std::size_t>(iPop)],
                1e-12)
        << "kVelocityDirichlet cell must Dirichlet-collide; iPop=" << iPop;
  }
}

}  // namespace
}  // namespace octlb