#ifndef OCTLB_SRC_SOLVER_IO_VTK_WRITER_AMR_VTK_WRITER_H_
#define OCTLB_SRC_SOLVER_IO_VTK_WRITER_AMR_VTK_WRITER_H_

#include <functional>
#include <string>
#include <vector>

#include <mpi.h>

#include "src/mesh/forest/octree_forest.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/io/vtk_writer/structured_grid_writer.h"
#include "src/solver/io/vtk_writer/vtk_cell_field.h"

namespace octlb {

/** Writes one `.vts` per local octant; optional W3 index files via WriteVtmAndPvd. */
class AmrVtkWriter {
 public:
  AmrVtkWriter(MPI_Comm comm, const OctreeForest& forest, int nx, int ny, int nz,
               std::string output_dir, std::string base_name);

  template <typename BlockT>
  std::vector<std::string> WriteTimestep(
      int iT, const BlockCollection<BlockT>& blocks,
      std::span<const VtkCellFieldView> fields);

  /** Gathers paths from all ranks and writes `.vtm` + appends `.pvd` on rank 0. */
  void WriteVtmAndPvd(int iT, const std::vector<std::string>& local_vts_paths);

  const std::string& output_dir() const { return output_dir_; }
  const std::string& base_name() const { return base_name_; }

 private:
  std::string BlockVtsFilename(OctantId id, int iT) const;

  MPI_Comm comm_;
  const OctreeForest& forest_;
  int nx_ = 0;
  int ny_ = 0;
  int nz_ = 0;
  std::string output_dir_;
  std::string base_name_;
};

template <typename BlockT>
std::vector<std::string> AmrVtkWriter::WriteTimestep(
    int iT, const BlockCollection<BlockT>& blocks,
    std::span<const VtkCellFieldView> fields) {
  std::vector<std::string> written;
  const label n = forest_.local_num_octants();
  for (label o = 0; o < n; ++o) {
    const OctantId id = static_cast<OctantId>(o);
    (void)blocks[id];
    VtkBlockMeta meta;
    meta.id = id;
    meta.nx = nx_;
    meta.ny = ny_;
    meta.nz = nz_;
    meta.bounds = forest_.quadrant_bounds(id);
    const std::string path = BlockVtsFilename(id, iT);
    WriteStructuredGridVts(path, meta, fields);
    written.push_back(path);
  }
  return written;
}

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_IO_VTK_WRITER_AMR_VTK_WRITER_H_
