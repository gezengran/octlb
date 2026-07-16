#ifndef OCTLB_SRC_SOLVER_LBM_BC_KIND_H_
#define OCTLB_SRC_SOLVER_LBM_BC_KIND_H_

#include <cstdint>

namespace octlb {

// Per-cell BC category (OctLB-flavored: per-cell uint8, no virtual dispatch).
// Stamped once at setup by BcInstaller; the core loops dispatch on this kind.
//
// This supersedes the overloaded CellKind::kBoundary, which conflated cavity3d
// tree-face Dirichlet cells and cylinder3d Bouzidi surface cells. Each BcKind
// names one BC role, so mixed BCs coexist within a single block.
enum class BcKind : std::uint8_t {
  kBulk,               // Fluid BGK.
  kSolid,              // Inert solid.
  kBounceBack,         // Static wall (full-way bounce-back collide).
  kMovingBounceBack,   // Moving wall (bounce-back + wall velocity, cavity3d lid).
  kBouzidi,            // Curved surface (cylinder), stream-time Bouzidi pull.
  kVelocityDirichlet,  // Inlet FD (prescribed u); stream bounce-back pull + PostStream FD.
  kPressureDirichlet,  // Outlet FD (prescribed rho); stream bounce-back pull + PostStream FD.
  kOutflow,            // Zero-gradient do-nothing outflow.
};

// Reflecting-boundary pull set: cells whose own stream reflects wall-pointing
// links (half-way bounce-back pull), matching the legacy CellKind::kBoundary
// stream behavior. PostStream FD then overwrites the Dirichlet variants.
// kBouzidi is included because a Bouzidi surface cell reflects its own
// wall-pointing links (Bouzidi interpolation fires on the *source* side).
inline bool BcKindReflectsOnPull(BcKind k) {
  return k == BcKind::kVelocityDirichlet || k == BcKind::kBouzidi;
}

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BC_KIND_H_