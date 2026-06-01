#ifndef OCTLB_SRC_COMMON_TYPES_H_
#define OCTLB_SRC_COMMON_TYPES_H_

#include <cstdint>

namespace octlb {

using label = std::int32_t;
using scalar = double;

/** Local linear index of an owned quadrant (0 … local_num_octants()-1). */
using OctantId = label;

}  // namespace octlb

#endif  // OCTLB_SRC_COMMON_TYPES_H_
