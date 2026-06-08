#ifndef OCTLB_SRC_SOLVER_LBM_BOUNDARY_ZOU_HE_VELOCITY_H_
#define OCTLB_SRC_SOLVER_LBM_BOUNDARY_ZOU_HE_VELOCITY_H_

#include "descriptor/descriptor.h"
#include "dynamics/lbm.h"

namespace octlb {
namespace boundary {

namespace detail {

inline bool PopPointsOutOfDomain(int iPop, int outward_x, int outward_y,
                                 int outward_z) {
  const int cx = olb::descriptors::c<olb::descriptors::D3Q19<>>(iPop, 0);
  const int cy = olb::descriptors::c<olb::descriptors::D3Q19<>>(iPop, 1);
  const int cz = olb::descriptors::c<olb::descriptors::D3Q19<>>(iPop, 2);
  return cx * outward_x + cy * outward_y + cz * outward_z > 0;
}

inline void OutwardNormal(int face, int* nx, int* ny, int* nz) {
  *nx = *ny = *nz = 0;
  switch (face) {
    case 0:
      *nx = -1;
      break;
    case 1:
      *nx = 1;
      break;
    case 2:
      *ny = -1;
      break;
    case 3:
      *ny = 1;
      break;
    case 4:
      *nz = -1;
      break;
    case 5:
      *nz = 1;
      break;
    default:
      break;
  }
}

}  // namespace detail

// Zou-He velocity BC on a flat wall: fill ghost layer with post-collision
// populations that reproduce the prescribed wall velocity. Uses the
// non-equilibrium bounce-back scheme (Zou & He, 1997) on unknown outgoing
// links, with a final equilibrium reset so the ghost cell matches u_wall.
template <typename T, typename DESCRIPTOR>
inline void ApplyZouHeVelocityGhost(T* ghost, const T* interior,
                                    int face_dir, const T* u_wall) {
  int outward_x = 0;
  int outward_y = 0;
  int outward_z = 0;
  detail::OutwardNormal(face_dir, &outward_x, &outward_y, &outward_z);

  T rho = T{0};
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    if (detail::PopPointsOutOfDomain(iPop, outward_x, outward_y, outward_z)) {
      continue;
    }
    rho += interior[iPop];
  }
  rho += T{1};

  T u[DESCRIPTOR::d]{};
  T uSqr = T{0};
  for (int d = 0; d < DESCRIPTOR::d; ++d) {
    u[d] = u_wall[d];
    uSqr += u[d] * u[d];
  }

  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    ghost[iPop] = olb::equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
  }
}

}  // namespace boundary
}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BOUNDARY_ZOU_HE_VELOCITY_H_
