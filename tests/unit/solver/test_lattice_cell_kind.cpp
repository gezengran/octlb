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

TEST(YminYmaxStreamOnlyPadding, IncludesIzFaceCornersOnYmaxSlab) {
  constexpr int nx = 31;
  constexpr int ny = 31;
  constexpr int nz = 31;
  EXPECT_TRUE(IsYminYmaxStreamOnlyPadding(15, ny, 15, nx, ny, nz));
  EXPECT_TRUE(IsYminYmaxStreamOnlyPadding(15, ny, 0, nx, ny, nz));
  EXPECT_TRUE(IsYminYmaxStreamOnlyPadding(15, ny, nz - 1, nx, ny, nz));
  EXPECT_TRUE(IsYminYmaxStreamOnlyPadding(0, ny, 0, nx, ny, nz));
  EXPECT_FALSE(OverlapPaddingMaterialNonZero(15, ny, 0, nx, ny, nz));
  EXPECT_FALSE(OverlapPaddingMaterialNonZero(15, -1, 15, nx, ny, nz));
  EXPECT_FALSE(OverlapPaddingMaterialNonZero(31, 30, 31, nx, ny, nz));
  EXPECT_FALSE(OverlapPaddingMaterialNonZero(-3, 5, 33, nx, ny, nz));
}

TEST(YminYmaxStreamOnlyPadding, YmaxIzCornerSkipsMirrorFill) {
  constexpr int n = 8;
  BlockLattice<T, Descriptor> lat(n, n, n, 1);
  const std::array<T, 3> u0{T{0}, T{0}, T{0}};
  lat.initialize(T{0}, u0.data());

  const int ix = 4;
  const int iz = 0;
  const int iy_core = n - 2;
  const int iy_pad = n;
  T* core = lat.populations_at_halo(ix + 1, iy_core + 1, iz + 1);
  T* ymax_corner =
      lat.populations_at_halo(ix + 1, iy_pad + 1, iz + 1);
  ASSERT_NE(core, nullptr);
  ASSERT_NE(ymax_corner, nullptr);
  for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    core[iPop] = T{1} + static_cast<T>(iPop) * T{0.01};
    ymax_corner[iPop] = T{0};
  }

  lat.fill_overlap_padding_from_core();

  for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    EXPECT_EQ(ymax_corner[iPop], T{0})
        << "iPop=" << iPop << " iy=" << iy_pad << " iz=" << iz;
  }
}

}  // namespace octlb
