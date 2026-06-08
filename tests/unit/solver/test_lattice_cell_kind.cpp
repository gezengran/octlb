#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "block_lattice.h"
#include "src/solver/lbm/cell_kind.h"

namespace octlb {
namespace {

using T = double;
using Descriptor = olb::descriptors::D3Q19<>;

std::vector<T> SnapshotPopulations(BlockLattice<T, Descriptor>& lat) {
  const int n = lat.nx() * lat.ny() * lat.nz() * Descriptor::q;
  std::vector<T> out(static_cast<std::size_t>(n));
  for (int ix = 0; ix < lat.nx(); ++ix) {
    for (int iy = 0; iy < lat.ny(); ++iy) {
      for (int iz = 0; iz < lat.nz(); ++iz) {
        auto cell = lat.get(ix, iy, iz);
        for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
          const std::size_t flat =
              static_cast<std::size_t>((ix * lat.ny() + iy) * lat.nz() + iz) *
                  Descriptor::q +
              static_cast<std::size_t>(iPop);
          out[flat] = cell[iPop];
        }
      }
    }
  }
  return out;
}

}  // namespace

TEST(LatticeCellKind, DefaultAllFluid) {
  BlockLattice<T, Descriptor> lat(4, 4, 4);
  const std::array<T, 3> u0{T{0}, T{0}, T{0}};
  lat.initialize(T{1}, u0.data());

  auto cell = lat.get(1, 1, 1);
  cell[5] += T{0.05};
  const T perturbed = cell[5];

  lat.collide(T{1.0});
  cell = lat.get(1, 1, 1);
  EXPECT_NE(cell[5], perturbed);
}

TEST(LatticeCellKind, SolidCell_SkipsCollide) {
  BlockLattice<T, Descriptor> lat(4, 4, 4);
  const std::array<T, 3> u0{T{0}, T{0}, T{0}};
  lat.initialize(T{1}, u0.data());
  lat.set_cell_kind(1, 1, 1, CellKind::kSolid);
  lat.set_cell_kind(2, 2, 2, CellKind::kSolid);

  const auto before = SnapshotPopulations(lat);
  lat.collide(T{1.0});
  const auto after = SnapshotPopulations(lat);

  EXPECT_EQ(before, after);
}

}  // namespace octlb
