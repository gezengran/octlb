#include "src/solver/io/vtk_writer/parallel_vtk_index.h"

#include <cstring>

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace octlb {
namespace {

std::string TimestepTag(int iT) {
  char buf[16];
  std::snprintf(buf, sizeof(buf), "%05d", iT);
  return buf;
}

std::filesystem::path NormalizeDir(const std::string& dir) {
  std::filesystem::path p(dir);
  if (p.empty()) {
    return std::filesystem::path(".");
  }
  return p;
}

std::vector<char> PackStrings(const std::vector<std::string>& items) {
  std::vector<char> buf;
  auto append_i32 = [&buf](std::int32_t v) {
    const char* p = reinterpret_cast<const char*>(&v);
    buf.insert(buf.end(), p, p + sizeof(v));
  };
  append_i32(static_cast<std::int32_t>(items.size()));
  for (const std::string& s : items) {
    append_i32(static_cast<std::int32_t>(s.size()));
    buf.insert(buf.end(), s.begin(), s.end());
  }
  return buf;
}

std::vector<std::string> UnpackStrings(const std::vector<char>& buf) {
  std::vector<std::string> items;
  std::size_t pos = 0;
  auto read_i32 = [&]() -> std::int32_t {
    if (pos + sizeof(std::int32_t) > buf.size()) {
      throw std::runtime_error("CollectVtsRelativePaths: corrupt packed buffer");
    }
    std::int32_t v = 0;
    std::memcpy(&v, buf.data() + pos, sizeof(v));
    pos += sizeof(std::int32_t);
    return v;
  };
  const std::int32_t count = read_i32();
  items.reserve(static_cast<std::size_t>(count));
  for (std::int32_t i = 0; i < count; ++i) {
    const std::int32_t len = read_i32();
    if (pos + static_cast<std::size_t>(len) > buf.size()) {
      throw std::runtime_error("CollectVtsRelativePaths: corrupt string length");
    }
    items.emplace_back(buf.data() + pos, buf.data() + pos + len);
    pos += static_cast<std::size_t>(len);
  }
  return items;
}

}  // namespace

int GetMpiRank(MPI_Comm comm) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  return rank;
}

int GetMpiSize(MPI_Comm comm) {
  int size = 1;
  MPI_Comm_size(comm, &size);
  return size;
}

std::vector<std::string> CollectVtsRelativePaths(
    MPI_Comm comm, const std::vector<std::string>& local_absolute_paths,
    const std::string& output_dir) {
  const int rank = GetMpiRank(comm);
  const int size = GetMpiSize(comm);
  const auto out_root = NormalizeDir(output_dir);

  std::vector<std::string> local_rel;
  local_rel.reserve(local_absolute_paths.size());
  for (const std::string& abs : local_absolute_paths) {
    std::filesystem::path rel =
        std::filesystem::relative(std::filesystem::path(abs), out_root);
    local_rel.push_back(rel.generic_string());
  }

  const std::vector<char> local_packed = PackStrings(local_rel);
  const int local_bytes = static_cast<int>(local_packed.size());

  std::vector<int> byte_counts(size, 0);
  MPI_Allgather(&local_bytes, 1, MPI_INT, byte_counts.data(), 1, MPI_INT, comm);

  std::vector<int> displs(size, 0);
  int total_bytes = 0;
  for (int r = 0; r < size; ++r) {
    displs[r] = total_bytes;
    total_bytes += byte_counts[r];
  }

  std::vector<char> gathered(static_cast<std::size_t>(total_bytes));
  MPI_Allgatherv(local_packed.data(), local_bytes, MPI_CHAR, gathered.data(),
                 byte_counts.data(), displs.data(), MPI_CHAR, comm);

  std::vector<std::string> all;
  for (int r = 0; r < size; ++r) {
    const std::vector<char> chunk(
        gathered.begin() + displs[r],
        gathered.begin() + displs[r] + byte_counts[r]);
    std::vector<std::string> part = UnpackStrings(chunk);
    all.insert(all.end(), part.begin(), part.end());
  }
  (void)rank;
  return all;
}

void WriteVtmIndex(MPI_Comm comm, int iT, const std::string& output_dir,
                   const std::string& base_name,
                   const std::vector<std::string>& vts_relative_paths) {
  if (GetMpiRank(comm) != 0) {
    return;
  }

  std::filesystem::create_directories(output_dir);
  const std::string vtm_path =
      output_dir + "/" + base_name + "_T" + TimestepTag(iT) + ".vtm";

  std::ofstream out(vtm_path, std::ios::trunc);
  if (!out) {
    throw std::runtime_error("WriteVtmIndex: failed to open " + vtm_path);
  }

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
         "byte_order=\"LittleEndian\">\n";
  out << "<vtkMultiBlockDataSet>\n";
  for (std::size_t i = 0; i < vts_relative_paths.size(); ++i) {
    out << "<Block index=\"" << i << "\">\n";
    out << "<DataSet index=\"0\" file=\"" << vts_relative_paths[i] << "\"/>\n";
    out << "</Block>\n";
  }
  out << "</vtkMultiBlockDataSet>\n";
  out << "</VTKFile>\n";
}

void AppendPvdTimestep(int iT, const std::string& pvd_path,
                       const std::string& vtm_relative_path) {
  const bool exists = std::filesystem::exists(pvd_path);
  if (!exists) {
    std::ofstream create(pvd_path, std::ios::trunc);
    if (!create) {
      throw std::runtime_error("AppendPvdTimestep: failed to create " + pvd_path);
    }
    create << "<?xml version=\"1.0\"?>\n";
    create << "<VTKFile type=\"Collection\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">\n";
    create << "<Collection>\n";
    create << "</Collection>\n";
    create << "</VTKFile>\n";
  }

  std::fstream out(pvd_path, std::ios::in | std::ios::out);
  if (!out) {
    throw std::runtime_error("AppendPvdTimestep: failed to open " + pvd_path);
  }
  out.seekp(-25, std::ios::end);
  out << "<DataSet timestep=\"" << iT << "\" group=\"\" part=\"0\" "
         "file=\"" << vtm_relative_path << "\"/>\n";
  out << "</Collection>\n";
  out << "</VTKFile>\n";
}

}  // namespace octlb
