#include "src/solver/io/vtk_writer/amr_vtk_writer.h"

#include <cstdio>
#include <filesystem>

#include "src/solver/io/vtk_writer/parallel_vtk_index.h"

namespace octlb {
namespace {

std::string TimestepTag(int iT) {
  char buf[16];
  std::snprintf(buf, sizeof(buf), "%05d", iT);
  return buf;
}

}  // namespace

AmrVtkWriter::AmrVtkWriter(MPI_Comm comm, const OctreeForest& forest, int nx,
                           int ny, int nz, std::string output_dir,
                           std::string base_name)
    : comm_(comm),
      forest_(forest),
      nx_(nx),
      ny_(ny),
      nz_(nz),
      output_dir_(std::move(output_dir)),
      base_name_(std::move(base_name)) {
  std::filesystem::create_directories(output_dir_);
}

std::string AmrVtkWriter::BlockVtsFilename(OctantId id, int iT) const {
  // Local OctantId repeats across ranks after partition; prefix MPI rank for uniqueness.
  return output_dir_ + "/" + base_name_ + "_r" + std::to_string(GetMpiRank(comm_)) +
         "_oct" + std::to_string(id) + "_T" + TimestepTag(iT) + ".vts";
}

void AmrVtkWriter::WriteVtmAndPvd(int iT,
                                  const std::vector<std::string>& local_vts_paths) {
  const std::vector<std::string> relative =
      CollectVtsRelativePaths(comm_, local_vts_paths, output_dir_);
  WriteVtmIndex(comm_, iT, output_dir_, base_name_, relative);
  if (GetMpiRank(comm_) == 0) {
    const std::string vtm_rel = base_name_ + "_T" + TimestepTag(iT) + ".vtm";
    AppendPvdTimestep(iT, output_dir_ + "/" + base_name_ + ".pvd", vtm_rel);
  }
}

}  // namespace octlb
