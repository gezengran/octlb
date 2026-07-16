#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "src/solver/lbm/bc_dispatcher.h"
#include "src/solver/lbm/bc_kind.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {
namespace {

using T = double;
using Descriptor = olb::descriptors::D3Q19<>;
using Lattice = BlockLattice<T, Descriptor>;

constexpr int kQ = Descriptor::q;

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

  // Plant known, distinct incident populations on the stamped cell.
  auto cell = lat.get(1, 1, 1);
  std::vector<T> before(static_cast<std::size_t>(kQ));
  for (int iPop = 0; iPop < kQ; ++iPop) {
    cell[iPop] = T{0.01} * static_cast<T>(iPop + 1) - T{0.05};
    before[static_cast<std::size_t>(iPop)] = cell[iPop];
  }

  const std::vector<DomainBcSpec> specs;  // unused by the bounce-back arm
  BcDispatcher::collide(lat, 1, 1, 1, /*omega=*/1.0, specs);

  cell = lat.get(1, 1, 1);
  for (int iPop = 0; iPop < kQ; ++iPop) {
    const int opp = olb::descriptors::opposite<Descriptor>(iPop);
    EXPECT_DOUBLE_EQ(cell[iPop], before[static_cast<std::size_t>(opp)])
        << "iPop=" << iPop << " must equal pre-collision opposite";
  }
}

}  // namespace
}  // namespace octlb