#ifndef OCTLB_SRC_COMMON_BOUNDING_BOX_H_
#define OCTLB_SRC_COMMON_BOUNDING_BOX_H_

#include "src/common/types.h"

namespace octlb {

struct BoundingBox {
  scalar x_min = 0;
  scalar y_min = 0;
  scalar z_min = 0;
  scalar x_max = 0;
  scalar y_max = 0;
  scalar z_max = 0;
};

}  // namespace octlb

#endif  // OCTLB_SRC_COMMON_BOUNDING_BOX_H_
