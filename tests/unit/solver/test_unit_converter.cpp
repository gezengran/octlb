#include <gtest/gtest.h>

#include <cmath>

#include "src/solver/lbm/unit_converter/unit_converter.h"

namespace octlb {

// D3Q19: cs^2 = 1/3 => invCs2 = 3 (same as olb::descriptors::invCs2 for D3Q19).
constexpr scalar kD3Q19InvCs2 = 3.0;

TEST(UnitConverter, OpenLbDefaults_Re1000) {
  const UnitConverter conv = UnitConverter::OpenLbCavity3dDefaults();
  // OpenLB cavity3d defaults: U=L=1, nu=0.001 => Re=U*L/nu=1000 (Ghia table is Re=100 label).
  EXPECT_NEAR(conv.reynolds(), 1000.0, 1e-10);
}

TEST(UnitConverter, TauOmegaConsistent) {
  const UnitConverter conv = UnitConverter::OpenLbCavity3dDefaults();
  EXPECT_NEAR(conv.omega(), 1.0 / conv.lattice_relaxation_time(), 1e-14);
}

TEST(UnitConverter, LatticeLidVelocity) {
  const UnitConverter conv = UnitConverter::OpenLbCavity3dDefaults();
  constexpr int kN = 30;
  constexpr scalar kTau = 0.509;
  constexpr scalar kL = 1.0;
  constexpr scalar kU = 1.0;
  constexpr scalar kNu = 0.001;
  const scalar dx = kL / kN;
  const scalar dt = (kTau - 0.5) / kD3Q19InvCs2 * dx * dx / kNu;
  const scalar u_lid = kU / (dx / dt);
  EXPECT_NEAR(conv.char_lattice_velocity(), u_lid, 1e-14);
  EXPECT_NEAR(conv.char_lattice_velocity(), 0.1, 1e-14);
}

TEST(UnitConverter, LatticeTime_MaxPhysT) {
  const UnitConverter conv = UnitConverter::OpenLbCavity3dDefaults();
  const scalar dt = conv.phys_delta_t();
  const std::size_t expected =
      static_cast<std::size_t>(100.0 / dt + 0.5);
  EXPECT_EQ(conv.get_lattice_time(100.0), expected);
  EXPECT_EQ(expected, 30000u);
}

}  // namespace octlb
