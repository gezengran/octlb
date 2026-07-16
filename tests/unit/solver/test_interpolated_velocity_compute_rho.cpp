#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/interpolated_velocity.h"
#include "src/solver/lbm/bc_kind.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {
namespace {

using T = double;
using Descriptor = olb::descriptors::D3Q19<>;
using Lattice = BlockLattice<T, Descriptor>;

std::vector<DomainBcSpec> AllInterpolatedVelocitySpecs() {
  std::vector<DomainBcSpec> specs;
  for (const FaceDir face :
       {FaceDir::kXMin, FaceDir::kXMax, FaceDir::kYMin, FaceDir::kYMax,
        FaceDir::kZMin, FaceDir::kZMax}) {
    DomainBcSpec spec;
    spec.face = face;
    spec.type = DomainBcType::kInterpolatedVelocity;
    specs.push_back(spec);
  }
  return specs;
}

// OpenLB lbm::computeRho on f-t shifted storage: rho = 1 + sum_i f_i.
T OpenLbBulkRhoFromPops(const T* f) {
  T rho = T{1};
  for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    rho += f[iPop];
  }
  return rho;
}

Lattice MakeInitializedBoundaryLattice(int n, int halo = 1) {
  Lattice lat(n, n, n, halo);
  const std::array<T, 3> u0{T{0}, T{0}, T{0}};
  lat.initialize(T{1}, u0.data());
  boundary::MarkDomainBoundaryBcKinds(lat, n, n, n);
  return lat;
}

TEST(InterpolatedVelocityComputeRho, CornerEdge_Equilibrium_ReturnsOneNotRawSum) {
  constexpr int kN = 8;
  Lattice lat = MakeInitializedBoundaryLattice(kN);
  const std::vector<DomainBcSpec> specs = AllInterpolatedVelocitySpecs();
  boundary::detail::BoundaryLatticeView<T, Descriptor, Lattice> view(
      lat, kN, kN, kN, specs);

  const std::array<std::array<int, 3>, 4> corner_edge_cells{{
      {0, 0, 0},
      {kN - 1, 0, kN - 1},
      {0, kN / 2, kN - 1},
      {0, kN - 1, kN - 1},
  }};

  for (const auto& xyz : corner_edge_cells) {
    const int ix = xyz[0];
    const int iy = xyz[1];
    const int iz = xyz[2];
    ASSERT_EQ(lat.bc_kind(ix, iy, iz), BcKind::kVelocityDirichlet);
    EXPECT_GE(boundary::detail::BoundaryFaceCount(ix, iy, iz, kN, kN, kN), 2)
        << "cell (" << ix << ',' << iy << ',' << iz << ")";

    T raw_sum = T{0};
    for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
      raw_sum += lat.get(ix, iy, iz)[iPop];
    }
    EXPECT_NEAR(raw_sum, T{0}, 1e-14)
        << "cell (" << ix << ',' << iy << ',' << iz << ")";

    const T rho = view.ComputeRho(ix, iy, iz);
    EXPECT_NEAR(rho, T{1}, 1e-14)
        << "cell (" << ix << ',' << iy << ',' << iz << ")";
    EXPECT_NE(rho, raw_sum)
        << "ComputeRho must use OpenLB sum+1, not raw sum";
  }
}

TEST(InterpolatedVelocityComputeRho, CornerEdge_MatchesCellProxyAndOpenLbFormula) {
  constexpr int kN = 8;
  Lattice lat = MakeInitializedBoundaryLattice(kN);
  const std::vector<DomainBcSpec> specs = AllInterpolatedVelocitySpecs();
  boundary::detail::BoundaryLatticeView<T, Descriptor, Lattice> view(
      lat, kN, kN, kN, specs);

  constexpr int kIx = 0;
  constexpr int kIy = kN / 2;
  constexpr int kIz = kN - 1;

  auto cell = lat.get(kIx, kIy, kIz);
  cell[2] = T{0.04};
  cell[7] = T{-0.02};
  cell[12] = T{0.01};

  const T* f = &cell[0];
  T rho_cell = T{0};
  T u[3]{};
  cell.computeRhoU(rho_cell, u);

  EXPECT_NEAR(view.ComputeRho(kIx, kIy, kIz), rho_cell, 1e-14);
  EXPECT_NEAR(view.ComputeRho(kIx, kIy, kIz), OpenLbBulkRhoFromPops(f), 1e-14);
}

TEST(InterpolatedVelocityComputeRho, FlatFace_StillUsesVelocityBoundaryRho) {
  constexpr int kN = 8;
  Lattice lat = MakeInitializedBoundaryLattice(kN);
  const std::vector<DomainBcSpec> specs = AllInterpolatedVelocitySpecs();
  boundary::detail::BoundaryLatticeView<T, Descriptor, Lattice> view(
      lat, kN, kN, kN, specs);

  constexpr int kIx = kN / 2;
  constexpr int kIy = 0;
  constexpr int kIz = kN / 2;

  int direction = 0;
  int orientation = 0;
  ASSERT_TRUE(boundary::detail::FlatBoundaryFaceInfo(kIx, kIy, kIz, kN, kN, kN,
                                                     direction, orientation));
  ASSERT_EQ(direction, 1);
  ASSERT_EQ(orientation, -1);

  auto cell = lat.get(kIx, kIy, kIz);
  cell[3] = T{0.02};
  cell[9] = T{-0.01};

  T u_wall[3]{};
  boundary::detail::PrescribedBoundaryU(kIx, kIy, kIz, kN, kN, kN, specs,
                                        u_wall);
  const T rho_expected = boundary::detail::VelocityBoundaryRhoFromPop<
      T, Descriptor>(direction, orientation, &cell[0], u_wall);

  EXPECT_NEAR(view.ComputeRho(kIx, kIy, kIz), rho_expected, 1e-14);
}

TEST(InterpolatedVelocityComputeRho, CavityEdgeCell_02930_EquilibriumRhoOne) {
  constexpr int kN = 31;
  Lattice lat = MakeInitializedBoundaryLattice(kN, 3);
  const std::vector<DomainBcSpec> specs = AllInterpolatedVelocitySpecs();
  boundary::detail::BoundaryLatticeView<T, Descriptor, Lattice> view(
      lat, kN, kN, kN, specs);

  constexpr int kIx = 0;
  constexpr int kIy = kN - 2;
  constexpr int kIz = kN - 1;

  T raw_sum = T{0};
  for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    raw_sum += lat.get(kIx, kIy, kIz)[iPop];
  }
  EXPECT_NEAR(raw_sum, T{0}, 1e-14);

  EXPECT_NEAR(view.ComputeRho(kIx, kIy, kIz), T{1}, 1e-14);
}

// R2: pressure outlet (kPressureDirichlet, p=0). With a known uniform interior
// field (rho=1, u=(u0,0,0)) and the outlet face stamped kPressureDirichlet,
// PostStream reconstruction must reset the outlet cell to rho==1.0 (p=0) and
// u extrapolated from the interior (== interior u for a uniform field). The
// cell is perturbed first so the reconstruction is observable.
TEST(InterpolatedVelocityComputeRho, Outlet_PressureZero_ReconstructsRhoOneAndInteriorU) {
  constexpr int kN = 8;
  constexpr double kU0 = 0.05;
  Lattice lat(kN, kN, kN, 1);
  {
    const double u0[3] = {kU0, 0.0, 0.0};
    lat.initialize(1.0, u0);
  }
  // Stamp the whole +x (outlet) face as kPressureDirichlet so tangential
  // neighbors are also pressure cells and extrapolate from the interior.
  for (int iy = 0; iy < kN; ++iy) {
    for (int iz = 0; iz < kN; ++iz) {
      lat.set_bc_kind(kN - 1, iy, iz, BcKind::kPressureDirichlet);
    }
  }

  DomainBcSpec spec;
  spec.face = FaceDir::kXMax;
  spec.type = DomainBcType::kInterpolatedPressure;
  spec.rho_target = 1.0;  // p = 0
  const std::vector<DomainBcSpec> specs = {spec};
  boundary::detail::BoundaryLatticeView<T, Descriptor, Lattice> view(
      lat, kN, kN, kN, specs);

  // Perturb the outlet cell off the target so reconstruction is observable.
  constexpr int kIx = kN - 1;
  constexpr int kIy = kN / 2;
  constexpr int kIz = kN / 2;
  auto cell = lat.get(kIx, kIy, kIz);
  for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    cell[iPop] = T{0.02} * static_cast<T>(iPop) - T{0.18};  // arbitrary non-eq
  }

  boundary::detail::ApplyPlaneFdBoundary<T, Descriptor>(
      view, &lat.get(kIx, kIy, kIz)[0], kIx, kIy, kIz,
      /*direction=*/0, /*orientation=*/1, /*omega=*/1.0);

  cell = lat.get(kIx, kIy, kIz);
  T rho = T{0};
  T u[Descriptor::d]{};
  cell.computeRhoU(rho, u);
  EXPECT_NEAR(rho, T{1}, 1e-12) << "pressure outlet must prescribe rho=1 (p=0)";
  EXPECT_NEAR(u[0], kU0, 1e-12) << "outlet u must extrapolate the interior u";
  EXPECT_NEAR(u[1], T{0}, 1e-12);
  EXPECT_NEAR(u[2], T{0}, 1e-12);
}

}  // namespace
}  // namespace octlb
