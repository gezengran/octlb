#ifndef OCTLB_SRC_SOLVER_IO_VTK_WRITER_PARALLEL_VTK_INDEX_H_
#define OCTLB_SRC_SOLVER_IO_VTK_WRITER_PARALLEL_VTK_INDEX_H_

#include <mpi.h>

#include <string>
#include <vector>

namespace octlb {

int GetMpiRank(MPI_Comm comm);
int GetMpiSize(MPI_Comm comm);

/** Makes paths relative to \a output_dir for `.vtm` references. */
std::vector<std::string> CollectVtsRelativePaths(
    MPI_Comm comm, const std::vector<std::string>& local_absolute_paths,
    const std::string& output_dir);

/** Rank 0 writes `{base_name}_T{iT:05d}.vtm` listing all blocks. */
void WriteVtmIndex(MPI_Comm comm, int iT, const std::string& output_dir,
                   const std::string& base_name,
                   const std::vector<std::string>& vts_relative_paths);

/** Appends one timestep entry to an existing or new `.pvd` (rank 0 only). */
void AppendPvdTimestep(int iT, const std::string& pvd_path,
                       const std::string& vtm_relative_path);

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_IO_VTK_WRITER_PARALLEL_VTK_INDEX_H_
