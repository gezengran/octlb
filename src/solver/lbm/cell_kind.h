#ifndef OCTLB_SRC_SOLVER_LBM_CELL_KIND_H_
#define OCTLB_SRC_SOLVER_LBM_CELL_KIND_H_

#include <cstdint>

namespace octlb {

enum class CellKind : std::uint8_t {
  kFluid,
  kSolid,
  kBoundary,
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_CELL_KIND_H_
