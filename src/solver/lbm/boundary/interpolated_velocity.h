#ifndef OCTLB_SRC_SOLVER_LBM_BOUNDARY_INTERPOLATED_VELOCITY_H_
#define OCTLB_SRC_SOLVER_LBM_BOUNDARY_INTERPOLATED_VELOCITY_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "core/platform/platform.h"
#include "core/util.h"
#include "descriptor/descriptor.h"
#include "dynamics/lbm.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/bc_kind.h"
#include "src/solver/lbm/domain_boundary_handler.h"

namespace octlb {
namespace boundary {

// OpenLB InterpolatedVelocity::getNeighborhoodRadius() (PostStream overlap).
inline constexpr int kInterpolatedVelocityOverlap = 2;

namespace detail {

using Descriptor = olb::descriptors::D3Q19<>;
namespace TensorIdx = olb::util::tensorIndices3D;

template <typename T>
inline T BoundaryGradient(T u0, T u1, T u2) {
  return (-T{3} * u0 + T{4} * u1 - T{1} * u2) / T{2};
}

template <typename T>
inline T CentralGradient(T u_p1, T u_m1) {
  return (u_p1 - u_m1) / T{2};
}

inline DomainBcSpec FindSpec(const std::vector<DomainBcSpec>& specs, FaceDir dir) {
  for (const DomainBcSpec& spec : specs) {
    if (spec.face == dir) {
      return spec;
    }
  }
  return DomainBcSpec{};
}

inline bool IsBoundaryLatticeCell(int ix, int iy, int iz, int nx, int ny,
                                  int nz) {
  return ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 ||
         iz == nz - 1;
}

// Interior lattice indices (ix,iy,iz); not halo coordinates (see pack_face).
template <typename Lattice>
inline auto PopulationsAt(Lattice& lat, int ix, int iy, int iz) {
  return &lat.get(ix, iy, iz)[0];
}

template <typename Lattice>
inline auto PopulationsAt(const Lattice& lat, int ix, int iy, int iz) {
  return &lat.get(ix, iy, iz)[0];
}

inline bool IsInteriorTopLidCell(int ix, int iy, int iz, int nx, int ny,
                                 int nz) {
  return iy == ny - 1 && ix > 0 && ix < nx - 1 && iz > 0 && iz < nz - 1;
}

inline void PrescribedBoundaryU(int ix, int iy, int iz, int nx, int ny, int nz,
                                const std::vector<DomainBcSpec>& specs,
                                double u[3]) {
  u[0] = u[1] = u[2] = 0.0;
  // A velocity-Dirichlet face with an inlet_field prescribes u per cell (e.g.
  // cylinder3d's uniform/ Poiseuille inlet). Pick the field on any face the
  // cell touches (flat inlet cell -> its face; inlet/wall corner -> the inlet).
  const struct {
    bool on;
    FaceDir dir;
  } faces[6] = {
      {ix == 0, FaceDir::kXMin},    {ix == nx - 1, FaceDir::kXMax},
      {iy == 0, FaceDir::kYMin},    {iy == ny - 1, FaceDir::kYMax},
      {iz == 0, FaceDir::kZMin},    {iz == nz - 1, FaceDir::kZMax}};
  for (const auto& f : faces) {
    if (!f.on) {
      continue;
    }
    const DomainBcSpec spec = FindSpec(specs, f.dir);
    if (spec.inlet_field) {
      double u_d[3]{};
      spec.inlet_field->velocity(ix, iy, iz, 0.0, u_d);
      for (int d = 0; d < 3; ++d) u[d] = u_d[d];
      return;
    }
  }
  // OpenLB cavity3d fallback (no inlet_field): material 3 (moving lid) covers
  // the interior top face only; top edges/corners stay material 2 with u=0.
  if (IsInteriorTopLidCell(ix, iy, iz, nx, ny, nz)) {
    const DomainBcSpec spec = FindSpec(specs, FaceDir::kYMax);
    for (int d = 0; d < 3; ++d) {
      u[d] = spec.u_wall[d];
    }
  }
}

inline int BoundaryFaceCount(int ix, int iy, int iz, int nx, int ny, int nz) {
  int count = 0;
  if (ix == 0 || ix == nx - 1) {
    ++count;
  }
  if (iy == 0 || iy == ny - 1) {
    ++count;
  }
  if (iz == 0 || iz == nz - 1) {
    ++count;
  }
  return count;
}

// Returns true for flat-face boundary cells (exactly one domain face).
inline bool FlatBoundaryFaceInfo(int ix, int iy, int iz, int nx, int ny, int nz,
                                 int& direction, int& orientation) {
  int count = 0;
  direction = 0;
  orientation = 0;
  if (ix == 0) {
    ++count;
    direction = 0;
    orientation = -1;
  } else if (ix == nx - 1) {
    ++count;
    direction = 0;
    orientation = 1;
  }
  if (iy == 0) {
    ++count;
    direction = 1;
    orientation = -1;
  } else if (iy == ny - 1) {
    ++count;
    direction = 1;
    orientation = 1;
  }
  if (iz == 0) {
    ++count;
    direction = 2;
    orientation = -1;
  } else if (iz == nz - 1) {
    ++count;
    direction = 2;
    orientation = 1;
  }
  return count == 1;
}

// OpenLB velocityBMRho for a fixed face direction/orientation.
template <int direction, int orientation, typename T, typename DESCRIPTOR>
inline T VelocityBoundaryRhoFromPopImpl(const T* f, const T u[DESCRIPTOR::d]) {
  constexpr auto on_wall =
      olb::util::populationsContributingToVelocity<DESCRIPTOR, direction, 0>();
  constexpr auto normal =
      olb::util::populationsContributingToVelocity<DESCRIPTOR, direction,
                                                   orientation>();

  T rho_on_wall = T{0};
  for (auto iPop : on_wall) {
    rho_on_wall += f[iPop];
  }

  T rho_normal = T{0};
  for (auto iPop : normal) {
    rho_normal += f[iPop];
  }

  return (T{2} * rho_normal + rho_on_wall + T{1}) /
         (T{1} + static_cast<T>(orientation) * u[direction]);
}

template <typename T, typename DESCRIPTOR>
inline T VelocityBoundaryRhoFromPop(int direction, int orientation, const T* f,
                                    const T u[DESCRIPTOR::d]) {
  if (direction == 0) {
    if (orientation < 0) {
      return VelocityBoundaryRhoFromPopImpl<0, -1, T, DESCRIPTOR>(f, u);
    }
    return VelocityBoundaryRhoFromPopImpl<0, 1, T, DESCRIPTOR>(f, u);
  }
  if (direction == 1) {
    if (orientation < 0) {
      return VelocityBoundaryRhoFromPopImpl<1, -1, T, DESCRIPTOR>(f, u);
    }
    return VelocityBoundaryRhoFromPopImpl<1, 1, T, DESCRIPTOR>(f, u);
  }
  if (orientation < 0) {
    return VelocityBoundaryRhoFromPopImpl<2, -1, T, DESCRIPTOR>(f, u);
  }
  return VelocityBoundaryRhoFromPopImpl<2, 1, T, DESCRIPTOR>(f, u);
}

template <typename T, typename DESCRIPTOR, typename Lattice>
class BoundaryLatticeView {
 public:
  BoundaryLatticeView(Lattice& lat, int nx, int ny, int nz,
                      const std::vector<DomainBcSpec>& specs)
      : lat_(lat),
        nx_(nx),
        ny_(ny),
        nz_(nz),
        h_(lat.halo_width()),
        specs_(specs) {}

  bool InDomain(int ix, int iy, int iz) const {
    return ix >= 0 && ix < nx_ && iy >= 0 && iy < ny_ && iz >= 0 && iz < nz_;
  }

  bool InOverlapStencil(int ix, int iy, int iz) const {
    const int stencil =
        std::min(h_, boundary::kInterpolatedVelocityOverlap);
    return ix >= -stencil && ix < nx_ + stencil && iy >= -stencil &&
           iy < ny_ + stencil && iz >= -stencil && iz < nz_ + stencil;
  }

  // OpenLB NoDynamics (material 0) on overlap padding: rho=1, u=0 for FD stencils.
  bool IsPaddingOverlapCell(int ix, int iy, int iz) const {
    return !InDomain(ix, iy, iz) && h_ > 0 && InOverlapStencil(ix, iy, iz);
  }

  static void NoDynamicsRhoU(T& rho, T u[DESCRIPTOR::d]) {
    rho = T{1};
    for (int d = 0; d < DESCRIPTOR::d; ++d) {
      u[d] = T{0};
    }
  }

  CellProxy<T, DESCRIPTOR> CellAtLatticeR(int ix, int iy, int iz) const {
    return CellProxy<T, DESCRIPTOR>(
        const_cast<T*>(lat_.populations_at_halo(ix + h_, iy + h_, iz + h_)));
  }

  // Map exterior lattice indices to interior mirror points for FD stencils.
  // Fallback when halo_width()==0 (no overlap padding).
  static int MirrorExteriorIndex(int i, int n) {
    if (i < 0) {
      return -i;
    }
    if (i >= n) {
      return 2 * (n - 1) - i;
    }
    return i;
  }

  void MirrorExteriorVelocity(int ix, int iy, int iz, T u[DESCRIPTOR::d]) const {
    const int mx = MirrorExteriorIndex(ix, nx_);
    const int my = MirrorExteriorIndex(iy, ny_);
    const int mz = MirrorExteriorIndex(iz, nz_);
    T u_in[DESCRIPTOR::d]{};
    ComputeInteriorU(mx, my, mz, u_in);

    const bool flip_x = ix < 0 || ix >= nx_;
    const bool flip_y = iy < 0 || iy >= ny_;
    const bool flip_z = iz < 0 || iz >= nz_;
    u[0] = flip_x ? -u_in[0] : u_in[0];
    u[1] = flip_y ? -u_in[1] : u_in[1];
    u[2] = flip_z ? -u_in[2] : u_in[2];
  }

  void ComputeInteriorU(int ix, int iy, int iz, T u[DESCRIPTOR::d]) const {
    const BcKind kind = lat_.bc_kind(ix, iy, iz);
    if (kind == BcKind::kVelocityDirichlet) {
      PrescribedU(ix, iy, iz, u);
      return;
    }
    if (kind == BcKind::kPressureDirichlet) {
      ExtrapolatedU(ix, iy, iz, u);
      return;
    }
    T rho = T{0};
    lat_.get(ix, iy, iz).computeRhoU(rho, u);
  }

  // FD-extrapolated u for a pressure cell: the inward (interior) neighbor's u.
  // For a uniform interior field this equals the interior u exactly.
  void ExtrapolatedU(int ix, int iy, int iz, T u[DESCRIPTOR::d]) const {
    int direction = 0;
    int orientation = 0;
    if (!FlatBoundaryFaceInfo(ix, iy, iz, nx_, ny_, nz_, direction,
                              orientation)) {
      for (int d = 0; d < DESCRIPTOR::d; ++d) u[d] = T{0};
      return;
    }
    int nix = ix;
    int niy = iy;
    int niz = iz;
    if (direction == 0) {
      nix = ix - orientation;
    } else if (direction == 1) {
      niy = iy - orientation;
    } else {
      niz = iz - orientation;
    }
    ComputeU(nix, niy, niz, u);
  }

  // Prescribed outlet density for a pressure cell (p = cs^2 * (rho - 1)).
  // Reads the face spec's rho_target (default 1.0 -> p=0).
  T PrescribedRho(int ix, int iy, int iz) const {
    int direction = 0;
    int orientation = 0;
    if (FlatBoundaryFaceInfo(ix, iy, iz, nx_, ny_, nz_, direction,
                             orientation)) {
      const FaceDir face = FaceDirFromDirectionOrientation(direction,
                                                            orientation);
      return static_cast<T>(FindSpec(specs_, face).rho_target);
    }
    return T{1};
  }

  static FaceDir FaceDirFromDirectionOrientation(int direction,
                                                 int orientation) {
    if (direction == 0) {
      return orientation < 0 ? FaceDir::kXMin : FaceDir::kXMax;
    }
    if (direction == 1) {
      return orientation < 0 ? FaceDir::kYMin : FaceDir::kYMax;
    }
    return orientation < 0 ? FaceDir::kZMin : FaceDir::kZMax;
  }

  void PrescribedU(int ix, int iy, int iz, T u[DESCRIPTOR::d]) const {
    double u_d[3]{};
    PrescribedBoundaryU(ix, iy, iz, nx_, ny_, nz_, specs_, u_d);
    for (int d = 0; d < DESCRIPTOR::d; ++d) {
      u[d] = static_cast<T>(u_d[d]);
    }
  }

  void ComputeRhoU(int ix, int iy, int iz, T& rho, T u[DESCRIPTOR::d]) const {
    if (IsPaddingOverlapCell(ix, iy, iz)) {
      NoDynamicsRhoU(rho, u);
      return;
    }
    if (!InDomain(ix, iy, iz)) {
      MirrorExteriorVelocity(ix, iy, iz, u);
      rho = T{1};
      return;
    }
    const BcKind kind = lat_.bc_kind(ix, iy, iz);
    if (kind == BcKind::kVelocityDirichlet) {
      PrescribedU(ix, iy, iz, u);
      rho = ComputeRho(ix, iy, iz);
      return;
    }
    if (kind == BcKind::kPressureDirichlet) {
      rho = PrescribedRho(ix, iy, iz);
      ExtrapolatedU(ix, iy, iz, u);
      return;
    }
    lat_.get(ix, iy, iz).computeRhoU(rho, u);
  }

  void ComputeU(int ix, int iy, int iz, T u[DESCRIPTOR::d]) const {
    if (IsPaddingOverlapCell(ix, iy, iz)) {
      T rho = T{0};
      NoDynamicsRhoU(rho, u);
      return;
    }
    if (!InDomain(ix, iy, iz)) {
      MirrorExteriorVelocity(ix, iy, iz, u);
      return;
    }
    // On boundary lattice cells, use prescribed Dirichlet velocity instead of
    // reconstructing from populations. This matches OpenLB, where the boundary
    // dynamics stores the wall velocity in an external VELOCITY field and uses
    // it for computeU().
    const BcKind kind = lat_.bc_kind(ix, iy, iz);
    if (kind == BcKind::kVelocityDirichlet) {
      PrescribedU(ix, iy, iz, u);
      return;
    }
    if (kind == BcKind::kPressureDirichlet) {
      ExtrapolatedU(ix, iy, iz, u);
      return;
    }
    T rho = T{0};
    lat_.get(ix, iy, iz).computeRhoU(rho, u);
  }

  T ComputeRho(int ix, int iy, int iz) const {
    if (IsPaddingOverlapCell(ix, iy, iz)) {
      return T{1};
    }
    if (!InDomain(ix, iy, iz)) {
      const int mx = MirrorExteriorIndex(ix, nx_);
      const int my = MirrorExteriorIndex(iy, ny_);
      const int mz = MirrorExteriorIndex(iz, nz_);
      return ComputeRho(mx, my, mz);
    }
    const BcKind kind = lat_.bc_kind(ix, iy, iz);
    if (kind == BcKind::kPressureDirichlet) {
      return PrescribedRho(ix, iy, iz);
    }
    const T* f = detail::PopulationsAt(lat_, ix, iy, iz);
    if (kind == BcKind::kVelocityDirichlet) {
      int direction = 0;
      int orientation = 0;
      if (FlatBoundaryFaceInfo(ix, iy, iz, nx_, ny_, nz_, direction,
                             orientation)) {
        T u[DESCRIPTOR::d]{};
        PrescribedU(ix, iy, iz, u);
        return VelocityBoundaryRhoFromPop<T, DESCRIPTOR>(direction, orientation,
                                                       f, u);
      }
      // OpenLB BulkDensity on external corner/edge: lbm::computeRho = sum(f)+1.
    }
    T rho = T{0};
    T u[DESCRIPTOR::d]{};
    lat_.get(ix, iy, iz).computeRhoU(rho, u);
    return rho;
  }

 private:
  Lattice& lat_;
  int nx_;
  int ny_;
  int nz_;
  int h_;
  const std::vector<DomainBcSpec>& specs_;
};

template <typename T, typename DESCRIPTOR>
inline void WriteFromPi(T* f, T rho, const T* u, const T* pi) {
  const T u_sqr = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    f[iPop] = olb::equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, u_sqr) +
              olb::equilibrium<DESCRIPTOR>::template fromPiToFneq<T>(iPop, pi);
  }
}

// OpenLB InnerCornerDensity3D / InnerEdgeDensity3D (collide-time momenta).
template <typename T, typename DESCRIPTOR>
inline T InnerCornerDensityFromPop(const T* f, const T u[DESCRIPTOR::d],
                                   int x_normal, int y_normal, int z_normal) {
  const T rho_x =
      VelocityBoundaryRhoFromPop<T, DESCRIPTOR>(0, x_normal, f, u);
  const T rho_y =
      VelocityBoundaryRhoFromPop<T, DESCRIPTOR>(1, y_normal, f, u);
  const T rho_z =
      VelocityBoundaryRhoFromPop<T, DESCRIPTOR>(2, z_normal, f, u);
  return (rho_x + rho_y + rho_z) / T{3};
}

template <typename T, typename DESCRIPTOR>
inline T InnerEdgeDensityFromPop(const T* f, const T u[DESCRIPTOR::d], int plane,
                                 int normal1, int normal2) {
  const int direction1 = (plane + 1) % 3;
  const int direction2 = (plane + 2) % 3;
  const T rho1 =
      VelocityBoundaryRhoFromPop<T, DESCRIPTOR>(direction1, normal1, f, u);
  const T rho2 =
      VelocityBoundaryRhoFromPop<T, DESCRIPTOR>(direction2, normal2, f, u);
  return (rho1 + rho2) / T{2};
}

// OpenLB InnerCornerStress3D / InnerEdgeStress3D (BulkStress on regularized pops).
template <typename T, typename DESCRIPTOR>
inline void InnerCornerStressFromPop(const T* f, T rho, const T u[DESCRIPTOR::d],
                                     int x_normal, int y_normal, int z_normal,
                                     T pi[olb::util::TensorVal<DESCRIPTOR>::n]) {
  T f_mod[DESCRIPTOR::q];
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    f_mod[iPop] = f[iPop];
  }
  const int v[3] = {-x_normal, -y_normal, -z_normal};
  const int unknown_f = olb::util::findVelocity<DESCRIPTOR>(v);
  if (unknown_f != DESCRIPTOR::q) {
    const int opposite_f = olb::descriptors::opposite<DESCRIPTOR>(unknown_f);
    const T u_sqr = olb::util::normSqr<T, DESCRIPTOR::d>(u);
    f_mod[unknown_f] =
        f_mod[opposite_f] -
        olb::equilibrium<DESCRIPTOR>::secondOrder(opposite_f, rho, u, u_sqr) +
        olb::equilibrium<DESCRIPTOR>::secondOrder(unknown_f, rho, u, u_sqr);
  }
  octlb::CellProxy<T, DESCRIPTOR> cell(f_mod);
  olb::lbm<DESCRIPTOR>::computeStress(cell, rho, u, pi);
}

template <typename T, typename DESCRIPTOR>
inline void InnerEdgeStressFromPop(const T* f, T rho, const T u[DESCRIPTOR::d],
                                   int plane, int normal1, int normal2,
                                   T pi[olb::util::TensorVal<DESCRIPTOR>::n]) {
  T f_mod[DESCRIPTOR::q];
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    f_mod[iPop] = f[iPop];
  }
  const T u_sqr = olb::util::normSqr<T, DESCRIPTOR::d>(u);
  const int direction1 = (plane + 1) % 3;
  const int direction2 = (plane + 2) % 3;
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    if (olb::descriptors::c<DESCRIPTOR>(iPop, direction1) == -normal1 &&
        olb::descriptors::c<DESCRIPTOR>(iPop, direction2) == -normal2) {
      const int opp = olb::descriptors::opposite<DESCRIPTOR>(iPop);
      f_mod[iPop] = f_mod[opp] -
                    olb::equilibrium<DESCRIPTOR>::secondOrder(opp, rho, u,
                                                                u_sqr) +
                    olb::equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u,
                                                              u_sqr);
    }
  }
  octlb::CellProxy<T, DESCRIPTOR> cell(f_mod);
  olb::lbm<DESCRIPTOR>::computeStress(cell, rho, u, pi);
}

// OpenLB CombinedRLBdynamics collide: f <- f_eq + pi_neq, then BGK/ConstRhoBGK.
template <typename T, typename DESCRIPTOR>
inline void CombinedRlbThenBgkCollide(octlb::CellProxy<T, DESCRIPTOR>& cell, T rho,
                                      const T u[DESCRIPTOR::d],
                                      const T pi[olb::util::TensorVal<DESCRIPTOR>::n],
                                      T omega, T average_rho = T{1},
                                      bool use_const_rho_bgk = false) {
  const T u_sqr = olb::util::normSqr<T, DESCRIPTOR::d>(u);
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] =
        olb::equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, u_sqr) +
        olb::equilibrium<DESCRIPTOR>::template fromPiToFneq<T>(iPop, pi);
  }
  if (use_const_rho_bgk) {
    CollideConstRhoBgkWithMacroscopic(cell, omega, average_rho, rho, u);
  } else {
    olb::lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  }
}

inline void EdgeTopologyFromFaces(bool on_xmin, bool on_xmax, bool on_ymin,
                                  bool on_ymax, bool on_zmin, bool on_zmax,
                                  int& plane, int& normal1, int& normal2) {
  const int x_normal = on_xmin ? -1 : on_xmax ? 1 : 0;
  const int y_normal = on_ymin ? -1 : on_ymax ? 1 : 0;
  const int z_normal = on_zmin ? -1 : on_zmax ? 1 : 0;
  if (x_normal != 0 && y_normal != 0) {
    plane = 2;
    normal1 = x_normal;
    normal2 = y_normal;
  } else if (x_normal != 0 && z_normal != 0) {
    plane = 1;
    normal1 = x_normal;
    normal2 = z_normal;
  } else {
    plane = 0;
    normal1 = y_normal;
    normal2 = z_normal;
  }
}

template <typename T, typename DESCRIPTOR, typename View>
inline void ApplyPlaneFdBoundary(View& view, T* f, int ix, int iy, int iz,
                                 int direction, int orientation, T omega) {
  T rho = T{0};
  T u[DESCRIPTOR::d]{};
  // Match OpenLB PlaneFdBoundaryProcessor3D: single computeRhoU before gradients.
  view.ComputeRhoU(ix, iy, iz, rho, u);

  T dx_u[DESCRIPTOR::d]{};
  T dy_u[DESCRIPTOR::d]{};
  T dz_u[DESCRIPTOR::d]{};

  const int inward_x = (direction == 0) ? -orientation : 0;
  const int inward_y = (direction == 1) ? -orientation : 0;
  const int inward_z = (direction == 2) ? -orientation : 0;

  auto normal_grad = [&](T vel_deriv[DESCRIPTOR::d]) {
    T u0[DESCRIPTOR::d]{};
    T u1[DESCRIPTOR::d]{};
    T u2[DESCRIPTOR::d]{};
    view.ComputeU(ix, iy, iz, u0);
    view.ComputeU(ix + inward_x, iy + inward_y, iz + inward_z, u1);
    view.ComputeU(ix + 2 * inward_x, iy + 2 * inward_y, iz + 2 * inward_z,
                  u2);
    for (int d = 0; d < DESCRIPTOR::d; ++d) {
      vel_deriv[d] = -orientation * BoundaryGradient(u0[d], u1[d], u2[d]);
    }
  };

  auto tangential_grad = [&](int tdx, int tdy, int tdz,
                             T vel_deriv[DESCRIPTOR::d]) {
    T u_p1[DESCRIPTOR::d]{};
    T u_m1[DESCRIPTOR::d]{};
    // Match OpenLB cell.neighbor(±1): out-of-domain neighbors use prescribed u.
    view.ComputeU(ix + tdx, iy + tdy, iz + tdz, u_p1);
    view.ComputeU(ix - tdx, iy - tdy, iz - tdz, u_m1);
    for (int d = 0; d < DESCRIPTOR::d; ++d) {
      vel_deriv[d] = CentralGradient(u_p1[d], u_m1[d]);
    }
  };

  if (direction == 0) {
    normal_grad(dx_u);
    tangential_grad(0, 1, 0, dy_u);
    tangential_grad(0, 0, 1, dz_u);
  } else if (direction == 1) {
    tangential_grad(1, 0, 0, dx_u);
    normal_grad(dy_u);
    tangential_grad(0, 0, 1, dz_u);
  } else {
    tangential_grad(1, 0, 0, dx_u);
    tangential_grad(0, 1, 0, dy_u);
    normal_grad(dz_u);
  }

  const T inv_cs2 = olb::descriptors::invCs2<T, DESCRIPTOR>();
  const T s_to_pi = -rho / inv_cs2 / omega;
  T pi[olb::util::TensorVal<DESCRIPTOR>::n]{};
  pi[TensorIdx::xx] = T{2} * dx_u[0] * s_to_pi;
  pi[TensorIdx::yy] = T{2} * dy_u[1] * s_to_pi;
  pi[TensorIdx::zz] = T{2} * dz_u[2] * s_to_pi;
  pi[TensorIdx::xy] = (dx_u[1] + dy_u[0]) * s_to_pi;
  pi[TensorIdx::xz] = (dx_u[2] + dz_u[0]) * s_to_pi;
  pi[TensorIdx::yz] = (dy_u[2] + dz_u[1]) * s_to_pi;

  WriteFromPi<T, DESCRIPTOR>(f, rho, u, pi);
}

template <typename T, typename DESCRIPTOR, typename View>
inline void ApplyCornerBoundary(View& view, T* f, int ix, int iy, int iz,
                                int x_normal, int y_normal, int z_normal,
                                T omega) {
  const T rho100 = view.ComputeRho(ix - x_normal, iy, iz);
  const T rho010 = view.ComputeRho(ix, iy - y_normal, iz);
  const T rho001 = view.ComputeRho(ix, iy, iz - z_normal);
  const T rho200 = view.ComputeRho(ix - 2 * x_normal, iy, iz);
  const T rho020 = view.ComputeRho(ix, iy - 2 * y_normal, iz);
  const T rho002 = view.ComputeRho(ix, iy, iz - 2 * z_normal);
  const T rho = T{4} / T{9} * (rho001 + rho010 + rho100) -
                T{1} / T{9} * (rho002 + rho020 + rho200);

  auto normal_grad = [&](int axis, int normal, T vel_deriv[DESCRIPTOR::d]) {
    const int ox = axis == 0 ? -normal : 0;
    const int oy = axis == 1 ? -normal : 0;
    const int oz = axis == 2 ? -normal : 0;
    T u0[DESCRIPTOR::d]{};
    T u1[DESCRIPTOR::d]{};
    T u2[DESCRIPTOR::d]{};
    view.ComputeU(ix, iy, iz, u0);
    view.ComputeU(ix + ox, iy + oy, iz + oz, u1);
    view.ComputeU(ix + 2 * ox, iy + 2 * oy, iz + 2 * oz, u2);
    for (int d = 0; d < DESCRIPTOR::d; ++d) {
      vel_deriv[d] = -normal * BoundaryGradient(u0[d], u1[d], u2[d]);
    }
  };

  T dx_u[DESCRIPTOR::d]{};
  T dy_u[DESCRIPTOR::d]{};
  T dz_u[DESCRIPTOR::d]{};
  T u[DESCRIPTOR::d]{};
  normal_grad(0, x_normal, dx_u);
  normal_grad(1, y_normal, dy_u);
  normal_grad(2, z_normal, dz_u);
  view.ComputeU(ix, iy, iz, u);

  const T inv_cs2 = olb::descriptors::invCs2<T, DESCRIPTOR>();
  const T s_to_pi = -rho / inv_cs2 / omega;
  T pi[olb::util::TensorVal<DESCRIPTOR>::n]{};
  pi[TensorIdx::xx] = T{2} * dx_u[0] * s_to_pi;
  pi[TensorIdx::yy] = T{2} * dy_u[1] * s_to_pi;
  pi[TensorIdx::zz] = T{2} * dz_u[2] * s_to_pi;
  pi[TensorIdx::xy] = (dx_u[1] + dy_u[0]) * s_to_pi;
  pi[TensorIdx::xz] = (dx_u[2] + dz_u[0]) * s_to_pi;
  pi[TensorIdx::yz] = (dy_u[2] + dz_u[1]) * s_to_pi;
  WriteFromPi<T, DESCRIPTOR>(f, rho, u, pi);
}

template <typename T, typename DESCRIPTOR, typename View>
inline void ApplyEdgeBoundary(View& view, T* f, int ix, int iy, int iz,
                              int plane, int normal1, int normal2, T omega) {
  const int direction1 = (plane + 1) % 3;
  const int direction2 = (plane + 2) % 3;

  auto step_rho = [&](int s1, int s2) {
    int coords[3] = {0, 0, 0};
    coords[direction1] = -normal1 * s1;
    coords[direction2] = -normal2 * s2;
    return view.ComputeRho(ix + coords[0], iy + coords[1], iz + coords[2]);
  };

  const T rho10 = step_rho(1, 0);
  const T rho01 = step_rho(0, 1);
  const T rho20 = step_rho(2, 0);
  const T rho02 = step_rho(0, 2);
  const T rho = T{2} / T{3} * (rho01 + rho10) - T{1} / T{6} * (rho02 + rho20);

  T dA_uB[3][3]{};
  auto central_along = [&](int axis, T vel_deriv[DESCRIPTOR::d]) {
    const int tdx = axis == 0 ? 1 : 0;
    const int tdy = axis == 1 ? 1 : 0;
    const int tdz = axis == 2 ? 1 : 0;
    T u_p1[DESCRIPTOR::d]{};
    T u_m1[DESCRIPTOR::d]{};
    view.ComputeU(ix + tdx, iy + tdy, iz + tdz, u_p1);
    view.ComputeU(ix - tdx, iy - tdy, iz - tdz, u_m1);
    for (int d = 0; d < DESCRIPTOR::d; ++d) {
      vel_deriv[d] = CentralGradient(u_p1[d], u_m1[d]);
    }
  };
  auto normal_along = [&](int axis, int normal, T vel_deriv[DESCRIPTOR::d]) {
    const int ox = axis == 0 ? -normal : 0;
    const int oy = axis == 1 ? -normal : 0;
    const int oz = axis == 2 ? -normal : 0;
    T u0[DESCRIPTOR::d]{};
    T u1[DESCRIPTOR::d]{};
    T u2[DESCRIPTOR::d]{};
    view.ComputeU(ix, iy, iz, u0);
    view.ComputeU(ix + ox, iy + oy, iz + oz, u1);
    view.ComputeU(ix + 2 * ox, iy + 2 * oy, iz + 2 * oz, u2);
    for (int d = 0; d < DESCRIPTOR::d; ++d) {
      vel_deriv[d] = -normal * BoundaryGradient(u0[d], u1[d], u2[d]);
    }
  };

  central_along(plane, dA_uB[plane]);
  normal_along(direction1, normal1, dA_uB[direction1]);
  normal_along(direction2, normal2, dA_uB[direction2]);

  T u[DESCRIPTOR::d]{};
  view.ComputeU(ix, iy, iz, u);

  const T inv_cs2 = olb::descriptors::invCs2<T, DESCRIPTOR>();
  const T s_to_pi = -rho / inv_cs2 / omega;
  T pi[olb::util::TensorVal<DESCRIPTOR>::n]{};
  pi[TensorIdx::xx] = T{2} * dA_uB[0][0] * s_to_pi;
  pi[TensorIdx::yy] = T{2} * dA_uB[1][1] * s_to_pi;
  pi[TensorIdx::zz] = T{2} * dA_uB[2][2] * s_to_pi;
  pi[TensorIdx::xy] = (dA_uB[0][1] + dA_uB[1][0]) * s_to_pi;
  pi[TensorIdx::xz] = (dA_uB[0][2] + dA_uB[2][0]) * s_to_pi;
  pi[TensorIdx::yz] = (dA_uB[1][2] + dA_uB[2][1]) * s_to_pi;

  WriteFromPi<T, DESCRIPTOR>(f, rho, u, pi);
}

inline int OutwardOrientation(FaceDir face) {
  switch (face) {
    case FaceDir::kXMin:
      return -1;
    case FaceDir::kXMax:
      return 1;
    case FaceDir::kYMin:
      return -1;
    case FaceDir::kYMax:
      return 1;
    case FaceDir::kZMin:
      return -1;
    case FaceDir::kZMax:
      return 1;
  }
  return 0;
}

inline int FaceDirection(FaceDir face) {
  switch (face) {
    case FaceDir::kXMin:
    case FaceDir::kXMax:
      return 0;
    case FaceDir::kYMin:
    case FaceDir::kYMax:
      return 1;
    case FaceDir::kZMin:
    case FaceDir::kZMax:
      return 2;
  }
  return 0;
}

}  // namespace detail

inline void MarkDomainBoundaryBcKinds(
    BlockLattice<double, olb::descriptors::D3Q19<>>& lat, int nx, int ny,
    int nz) {
  // Transitional geometric marker (R0): stamps the closed-cavity outer shell
  // as kVelocityDirichlet and the interior as kBulk. Replaced by
  // bc::StampTreeBoundaryCells (face + spec) in R4 once cavity3d migrates to
  // per-cell dispatch.
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        if (detail::IsBoundaryLatticeCell(ix, iy, iz, nx, ny, nz)) {
          lat.set_bc_kind(ix, iy, iz, BcKind::kVelocityDirichlet);
        } else {
          lat.set_bc_kind(ix, iy, iz, BcKind::kBulk);
        }
      }
    }
  }
}

// OpenLB InterpolatedVelocity (default MixinDynamics=BGKdynamics): flat shells use
// exchange_momenta<BasicDirichletVelocityBoundaryTuple> -> plain BGK collide;
// external edge/corner use exchange_momenta<FixedVelocityBoundaryTuple> -> plain BGK.
// CombinedRLBdynamics is only for InternalEdge/InternalCorner (concave nodes), not
// on the closed cavity outer shell.
// Pre-collide rho statistic for ConstRho average (matches BGK collide input).
template <typename T, typename DESCRIPTOR, typename Lattice>
inline T BoundaryCellStatisticRho(Lattice& lat, int ix, int iy, int iz, int nx,
                                  int ny, int nz,
                                  const std::vector<DomainBcSpec>& specs) {
  double u_wall_d[3]{};
  detail::PrescribedBoundaryU(ix, iy, iz, nx, ny, nz, specs, u_wall_d);
  T u[DESCRIPTOR::d]{};
  for (int d = 0; d < DESCRIPTOR::d; ++d) {
    u[d] = static_cast<T>(u_wall_d[d]);
  }

  const bool on_xmin = ix == 0;
  const bool on_xmax = ix == nx - 1;
  const bool on_ymin = iy == 0;
  const bool on_ymax = iy == ny - 1;
  const bool on_zmin = iz == 0;
  const bool on_zmax = iz == nz - 1;
  const int face_count = static_cast<int>(on_xmin) +
                         static_cast<int>(on_xmax) +
                         static_cast<int>(on_ymin) +
                         static_cast<int>(on_ymax) +
                         static_cast<int>(on_zmin) +
                         static_cast<int>(on_zmax);

  const T* f = detail::PopulationsAt(lat, ix, iy, iz);
  T rho = T{0};
  if (face_count >= 2) {
    T u_tmp[DESCRIPTOR::d]{};
    lat.get(ix, iy, iz).computeRhoU(rho, u_tmp);
  } else {
    int direction = 0;
    int orientation = 0;
    detail::FlatBoundaryFaceInfo(ix, iy, iz, nx, ny, nz, direction, orientation);
    rho = detail::VelocityBoundaryRhoFromPop<T, DESCRIPTOR>(
        direction, orientation, f, u);
  }
  return rho;
}

// Single boundary lattice cell collide (OpenLB Dominant spatial order).
template <typename T, typename DESCRIPTOR, typename Lattice>
inline void CollideDirichletBoundaryCellAt(Lattice& lat, int ix, int iy, int iz,
                                           int nx, int ny, int nz, T omega,
                                           const std::vector<DomainBcSpec>& specs,
                                           CollideRhoStats* rho_stats = nullptr) {
  if (lat.bc_kind(ix, iy, iz) != BcKind::kVelocityDirichlet) {
    return;
  }
  double u_wall_d[3]{};
  detail::PrescribedBoundaryU(ix, iy, iz, nx, ny, nz, specs, u_wall_d);
  T u[DESCRIPTOR::d]{};
  for (int d = 0; d < DESCRIPTOR::d; ++d) {
    u[d] = static_cast<T>(u_wall_d[d]);
  }

  T rho = BoundaryCellStatisticRho<T, DESCRIPTOR, Lattice>(
      lat, ix, iy, iz, nx, ny, nz, specs);
  auto cell = lat.get(ix, iy, iz);
  olb::lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);

  if (rho_stats != nullptr) {
    rho_stats->add(static_cast<double>(rho));
  }
}

template <typename T, typename DESCRIPTOR, typename Lattice>
inline void CollideDirichletBoundaryCells(Lattice& lat, int nx, int ny, int nz,
                                          T omega,
                                          const std::vector<DomainBcSpec>& specs,
                                          CollideRhoStats* rho_stats = nullptr,
                                          T /*average_rho*/ = T{1},
                                          bool /*use_const_rho_bgk*/ = false) {
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        CollideDirichletBoundaryCellAt<T, DESCRIPTOR, Lattice>(
            lat, ix, iy, iz, nx, ny, nz, omega, specs, rho_stats);
      }
    }
  }
}

// OpenLB geometry: padding cells outside the closed cavity surface have material 0.
inline bool PaddingMaterialNonZero(int ix, int iy, int iz, int nx, int ny,
                                   int nz) {
  return OverlapPaddingMaterialNonZero(ix, iy, iz, nx, ny, nz);
}

// PostStream: addPoints2CommBC communicate (padding already streamed in stream()).
template <typename Lattice>
inline void FillInterpolatedVelocityOverlapPadding(Lattice& lat, int nx, int ny,
                                                   int nz) {
  lat.fill_overlap_padding_bc_post_stream(
      [&](int ix, int iy, int iz) {
        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz) {
          return false;
        }
        return BcKindIsFdBoundary(lat.bc_kind(ix, iy, iz));
      },
      [&](int ix, int iy, int iz) {
        return PaddingMaterialNonZero(ix, iy, iz, nx, ny, nz);
      });
}

template <typename Lattice>
inline void FillOverlapPaddingForMode(Lattice& lat, int nx, int ny, int nz,
                                      OverlapPaddingMode mode) {
  if (mode == OverlapPaddingMode::kMirror) {
    lat.fill_overlap_padding_from_core();
  } else if (mode == OverlapPaddingMode::kHybrid) {
    // OpenLB PostStream: addPoints2CommBC partner copy on mat!=0 padding only.
    // Stream-only mat-0 cells (incl. ymin/ymax slabs) keep stream() rotate values.
    lat.fill_overlap_padding_bc_post_stream(
        [&](int ix, int iy, int iz) {
          if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz) {
            return false;
          }
          return BcKindIsFdBoundary(lat.bc_kind(ix, iy, iz));
        },
        [&](int ix, int iy, int iz) {
          return PaddingMaterialNonZero(ix, iy, iz, nx, ny, nz) &&
                 !IsYminYmaxStreamOnlyPadding(ix, iy, iz, nx, ny, nz);
        });
  } else {
    FillInterpolatedVelocityOverlapPadding(lat, nx, ny, nz);
  }
}

// OpenLB InterpolatedVelocity (Skordos FD) on explicit boundary lattice cells.
template <typename T, typename DESCRIPTOR, typename Lattice>
inline void ApplyInterpolatedVelocityBoundaryCells(
    Lattice& lat, int nx, int ny, int nz, T omega,
    const std::vector<DomainBcSpec>& specs,
    OverlapPaddingMode padding_mode = OverlapPaddingMode::kHybrid) {
  FillOverlapPaddingForMode(lat, nx, ny, nz, padding_mode);
  detail::BoundaryLatticeView<T, DESCRIPTOR, Lattice> view(lat, nx, ny, nz,
                                                           specs);

  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        if (!BcKindIsFdBoundary(lat.bc_kind(ix, iy, iz))) {
          continue;
        }

        const bool on_xmin = ix == 0;
        const bool on_xmax = ix == nx - 1;
        const bool on_ymin = iy == 0;
        const bool on_ymax = iy == ny - 1;
        const bool on_zmin = iz == 0;
        const bool on_zmax = iz == nz - 1;
        const int face_count = static_cast<int>(on_xmin) +
                               static_cast<int>(on_xmax) +
                               static_cast<int>(on_ymin) +
                               static_cast<int>(on_ymax) +
                               static_cast<int>(on_zmin) +
                               static_cast<int>(on_zmax);

        T* f = detail::PopulationsAt(lat, ix, iy, iz);

        if (face_count == 3) {
          const int x_normal = on_xmin ? -1 : on_xmax ? 1 : 0;
          const int y_normal = on_ymin ? -1 : on_ymax ? 1 : 0;
          const int z_normal = on_zmin ? -1 : on_zmax ? 1 : 0;
          detail::ApplyCornerBoundary<T, DESCRIPTOR>(view, f, ix, iy, iz,
                                                     x_normal, y_normal,
                                                     z_normal, omega);
          continue;
        }

        if (face_count == 2) {
          const int x_normal = on_xmin ? -1 : on_xmax ? 1 : 0;
          const int y_normal = on_ymin ? -1 : on_ymax ? 1 : 0;
          const int z_normal = on_zmin ? -1 : on_zmax ? 1 : 0;
          int plane = 0;
          int normal1 = 0;
          int normal2 = 0;
          if (x_normal != 0 && y_normal != 0) {
            plane = 2;
            normal1 = x_normal;
            normal2 = y_normal;
          } else if (x_normal != 0 && z_normal != 0) {
            plane = 1;
            normal1 = x_normal;
            normal2 = z_normal;
          } else {
            plane = 0;
            normal1 = y_normal;
            normal2 = z_normal;
          }
          detail::ApplyEdgeBoundary<T, DESCRIPTOR>(view, f, ix, iy, iz, plane,
                                                   normal1, normal2, omega);
          continue;
        }

        FaceDir face = FaceDir::kXMin;
        if (on_xmax) {
          face = FaceDir::kXMax;
        } else if (on_ymax) {
          face = FaceDir::kYMax;
        } else if (on_ymin) {
          face = FaceDir::kYMin;
        } else if (on_zmax) {
          face = FaceDir::kZMax;
        } else if (on_zmin) {
          face = FaceDir::kZMin;
        }

        detail::ApplyPlaneFdBoundary<T, DESCRIPTOR>(
            view, f, ix, iy, iz, detail::FaceDirection(face),
            detail::OutwardOrientation(face), omega);
      }
    }
  }
}

}  // namespace boundary
}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BOUNDARY_INTERPOLATED_VELOCITY_H_
