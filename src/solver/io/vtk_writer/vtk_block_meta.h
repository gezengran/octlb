#ifndef OCTLB_SRC_SOLVER_IO_VTK_WRITER_VTK_BLOCK_META_H_
#define OCTLB_SRC_SOLVER_IO_VTK_WRITER_VTK_BLOCK_META_H_

#include "src/common/bounding_box.h"
#include "src/common/types.h"

namespace octlb {

struct VtkBlockMeta {
  OctantId id = 0;
  int nx = 0;
  int ny = 0;
  int nz = 0;
  BoundingBox bounds{};
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_IO_VTK_WRITER_VTK_BLOCK_META_H_
