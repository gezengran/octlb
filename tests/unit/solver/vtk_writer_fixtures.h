#ifndef OCTLB_TESTS_UNIT_SOLVER_VTK_WRITER_FIXTURES_H_
#define OCTLB_TESTS_UNIT_SOLVER_VTK_WRITER_FIXTURES_H_

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "src/solver/io/vtk_writer/vtk_binary_codec.h"

namespace octlb {
namespace vtk_test {

/** Shared output directory for multi-rank VTK tests (same path on every rank). */
inline std::string MpiSharedTempDir(const std::string& tag) {
  return (std::filesystem::temp_directory_path() / ("octlb_vtk_" + tag))
      .string();
}

inline std::string ReadFile(const std::string& path) {
  std::ifstream in(path);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

inline std::string ExtractFirstBase64Payload(const std::string& xml,
                                             const std::string& section_tag) {
  const std::string open = "<" + section_tag;
  const std::size_t sec = xml.find(open);
  if (sec == std::string::npos) {
    return {};
  }
  const std::size_t data_start = xml.find("<DataArray", sec);
  const std::size_t gt = xml.find('>', data_start);
  const std::size_t content_start = xml.find_first_not_of(" \n\r\t", gt + 1);
  const std::size_t close = xml.find("</DataArray>", content_start);
  if (content_start == std::string::npos || close == std::string::npos) {
    return {};
  }
  std::string payload = xml.substr(content_start, close - content_start);
  payload.erase(std::remove_if(payload.begin(), payload.end(),
                               [](unsigned char c) {
                                 return c == '\n' || c == '\r' || c == ' ' ||
                                        c == '\t';
                               }),
                payload.end());
  return payload;
}

inline std::vector<double> DecodePointsArray(const std::string& xml) {
  const std::string payload = ExtractFirstBase64Payload(xml, "Points");
  return vtk_binary::DecodeFloat64Array(payload);
}

inline std::vector<double> DecodeCellArrayNamed(const std::string& xml,
                                                const std::string& name) {
  const std::regex re(
      "<DataArray[^>]*Name=\"" + name +
      "\"[^>]*format=\"binary\"[^>]*>\\s*([^<]+)</DataArray>");
  std::smatch m;
  if (!std::regex_search(xml, m, re)) {
    return {};
  }
  std::string payload = m[1].str();
  payload.erase(std::remove_if(payload.begin(), payload.end(),
                               [](unsigned char c) {
                                 return c == '\n' || c == '\r' || c == ' ' ||
                                        c == '\t';
                               }),
                payload.end());
  return vtk_binary::DecodeFloat64Array(payload);
}

inline bool ApproxEq(double a, double b, double eps = 1e-10) {
  return std::abs(a - b) <= eps;
}

/** Constant scalar field for W1/W2 behavior tests. */
class DummyScalarField {
 public:
  explicit DummyScalarField(double value = 1.0) : value_(value) {}

  std::string_view vtk_name() const { return "dummy_scalar"; }
  int vtk_components() const { return 1; }

  void sample_cell(int i, int j, int k, double* out) const {
    (void)i;
    (void)j;
    (void)k;
    out[0] = value_;
  }

 private:
  double value_ = 1.0;
};

/** Vector field encoding cell index for multi-field tests. */
class DummyVectorField {
 public:
  std::string_view vtk_name() const { return "dummy_vector"; }
  int vtk_components() const { return 3; }

  void sample_cell(int i, int j, int k, double* out) const {
    out[0] = static_cast<double>(i);
    out[1] = static_cast<double>(j);
    out[2] = static_cast<double>(k);
  }
};

/** Interior-only field; halo marker value must never appear in VTK output. */
class InteriorMarkedField {
 public:
  explicit InteriorMarkedField(int nx, int ny, int nz, int halo = 1)
      : nx_(nx), ny_(ny), nz_(nz), halo_(halo) {
    const int sx = nx + 2 * halo;
    const int sy = ny + 2 * halo;
    const int sz = nz + 2 * halo;
    storage_.assign(static_cast<std::size_t>(sx * sy * sz), kHaloMarker);
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          At(i, j, k) = static_cast<double>(i + 100 * j + 10000 * k);
        }
      }
    }
  }

  std::string_view vtk_name() const { return "interior_marker"; }
  int vtk_components() const { return 1; }

  void sample_cell(int i, int j, int k, double* out) const {
    out[0] = At(i, j, k);
  }

  double HaloMarker() const { return kHaloMarker; }

 private:
  static constexpr double kHaloMarker = -999.0;

  std::size_t Index(int i, int j, int k) const {
    const int hi = i + halo_;
    const int hj = j + halo_;
    const int hk = k + halo_;
    const int sx = nx_ + 2 * halo_;
    const int sy = ny_ + 2 * halo_;
    return static_cast<std::size_t>(hi + sx * (hj + sy * hk));
  }

  double& At(int i, int j, int k) { return storage_[Index(i, j, k)]; }
  double At(int i, int j, int k) const { return storage_[Index(i, j, k)]; }

  int nx_ = 0;
  int ny_ = 0;
  int nz_ = 0;
  int halo_ = 0;
  std::vector<double> storage_;
};

}  // namespace vtk_test
}  // namespace octlb

#endif  // OCTLB_TESTS_UNIT_SOLVER_VTK_WRITER_FIXTURES_H_
