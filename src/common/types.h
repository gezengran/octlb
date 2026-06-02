#ifndef OCTLB_SRC_COMMON_TYPES_H_
#define OCTLB_SRC_COMMON_TYPES_H_

#include <cstdint>

namespace octlb {

using label = std::int32_t;
using scalar = double;

/** Local linear index of an owned quadrant (0 … local_num_octants()-1). */
using OctantId = label;

/** Face direction; indices match p4est face numbering (x-, x+, y-, y+, z-, z+). */
enum class FaceDir : int {
  kXMin = 0,
  kXMax = 1,
  kYMin = 2,
  kYMax = 3,
  kZMin = 4,
  kZMax = 5,
};

}  // namespace octlb

#endif  // OCTLB_SRC_COMMON_TYPES_H_
