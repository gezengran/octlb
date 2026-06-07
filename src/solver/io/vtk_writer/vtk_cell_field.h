#ifndef OCTLB_SRC_SOLVER_IO_VTK_WRITER_VTK_CELL_FIELD_H_
#define OCTLB_SRC_SOLVER_IO_VTK_WRITER_VTK_CELL_FIELD_H_

#include <concepts>
#include <functional>
#include <string_view>

namespace octlb {

template <typename T>
concept VtkCellField3D = requires(const T& f, int i, int j, int k, double* out) {
  { f.vtk_name() } -> std::convertible_to<std::string_view>;
  { f.vtk_components() } -> std::same_as<int>;
  { f.sample_cell(i, j, k, out) } -> std::same_as<void>;
};

/** Type-erased view for writing multiple fields in one .vts file. */
struct VtkCellFieldView {
  std::string_view name;
  int components = 1;
  std::function<void(int i, int j, int k, double* out)> sample;

  template <VtkCellField3D F>
  static VtkCellFieldView From(const F& field) {
    return {field.vtk_name(), field.vtk_components(),
            [&field](int i, int j, int k, double* out) {
              field.sample_cell(i, j, k, out);
            }};
  }
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_IO_VTK_WRITER_VTK_CELL_FIELD_H_
