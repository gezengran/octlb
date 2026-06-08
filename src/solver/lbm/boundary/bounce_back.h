#ifndef OCTLB_SRC_SOLVER_LBM_BOUNDARY_BOUNCE_BACK_H_
#define OCTLB_SRC_SOLVER_LBM_BOUNDARY_BOUNCE_BACK_H_

#include "descriptor/descriptor.h"

namespace octlb {
namespace boundary {

// Post-collision no-slip: fill ghost populations that stream into the domain
// from the adjacent interior boundary cell (half-way bounce-back).
template <typename T, typename DESCRIPTOR>
inline void ApplyNoSlipGhost(T* ghost, const T* interior) {
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    int cx = 0;
    int cy = 0;
    int cz = 0;
    for (int d = 0; d < DESCRIPTOR::d; ++d) {
      const int c = olb::descriptors::c<DESCRIPTOR>(iPop, d);
      if (d == 0) {
        cx = c;
      } else if (d == 1) {
        cy = c;
      } else {
        cz = c;
      }
    }
    if (cx == 0 && cy == 0 && cz == 0) {
      continue;
    }
    const int opp = olb::descriptors::opposite<DESCRIPTOR>(iPop);
    ghost[iPop] = interior[opp];
  }
}

}  // namespace boundary
}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BOUNDARY_BOUNCE_BACK_H_
