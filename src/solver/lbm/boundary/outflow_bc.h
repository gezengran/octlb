#ifndef OCTLB_SRC_SOLVER_LBM_BOUNDARY_OUTFLOW_BC_H_
#define OCTLB_SRC_SOLVER_LBM_BOUNDARY_OUTFLOW_BC_H_

#include "descriptor/descriptor.h"

namespace octlb {
namespace boundary {

// Zero-gradient (do-nothing) outflow: copy the adjacent interior boundary
// populations into the outer ghost layer. A pull-stream then reads back the
// same populations that left the domain, so a uniform outflow produces no
// reflection. This is the W1 outflow choice; a convective outflow
// (df/dt + U df/dx = 0) is deferred to a later wave.
template <typename T, typename DESCRIPTOR>
inline void ApplyOutflowGhost(T* ghost, const T* interior) {
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    ghost[iPop] = interior[iPop];
  }
}

}  // namespace boundary
}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BOUNDARY_OUTFLOW_BC_H_