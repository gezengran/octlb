#ifndef OCTLB_SRC_SOLVER_IO_VTK_WRITER_STRUCTURED_GRID_WRITER_H_
#define OCTLB_SRC_SOLVER_IO_VTK_WRITER_STRUCTURED_GRID_WRITER_H_

#include <span>
#include <string>

#include "src/solver/io/vtk_writer/vtk_block_meta.h"
#include "src/solver/io/vtk_writer/vtk_cell_field.h"

namespace octlb {

/** Writes one StructuredGrid `.vts` file (CellData, binary Float64). */
void WriteStructuredGridVts(const std::string& filename, const VtkBlockMeta& meta,
                            std::span<const VtkCellFieldView> fields);

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_IO_VTK_WRITER_STRUCTURED_GRID_WRITER_H_
