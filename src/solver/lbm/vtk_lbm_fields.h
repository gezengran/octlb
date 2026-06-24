#ifndef OCTLB_SRC_SOLVER_LBM_VTK_LBM_FIELDS_H_
#define OCTLB_SRC_SOLVER_LBM_VTK_LBM_FIELDS_H_

#include <array>
#include <string_view>

#include "block_lattice.h"
#include "src/solver/io/vtk_writer/vtk_cell_field.h"

namespace octlb {

/** Samples LBM velocity (u) in physical units on interior cells. */
template <typename T, typename DESCRIPTOR>
class VtkPhysVelocityField {
 public:
  VtkPhysVelocityField(const BlockLattice<T, DESCRIPTOR>& lattice,
                       const UnitConverter& converter)
      : lattice_(lattice),
        scale_(converter.char_phys_velocity() /
               converter.char_lattice_velocity()) {}

  std::string_view vtk_name() const { return "physVelocity"; }
  int vtk_components() const { return 3; }

  void sample_cell(int i, int j, int k, double* out) const {
    T rho{};
    T u[3]{};
    lattice_.get(i, j, k).computeRhoU(rho, u);
    for (int d = 0; d < 3; ++d) {
      out[d] = static_cast<double>(u[d]) * scale_;
    }
  }

 private:
  const BlockLattice<T, DESCRIPTOR>& lattice_;
  double scale_ = 1.0;
};

/** Samples LBM velocity (u) on interior cells for VTK CellData output. */
template <typename T, typename DESCRIPTOR>
class VtkVelocityField {
 public:
  explicit VtkVelocityField(const BlockLattice<T, DESCRIPTOR>& lattice)
      : lattice_(lattice) {}

  std::string_view vtk_name() const { return "velocity"; }
  int vtk_components() const { return 3; }

  void sample_cell(int i, int j, int k, double* out) const {
    T rho{};
    T u[3]{};
    lattice_.get(i, j, k).computeRhoU(rho, u);
    for (int d = 0; d < 3; ++d) {
      out[d] = static_cast<double>(u[d]);
    }
  }

 private:
  const BlockLattice<T, DESCRIPTOR>& lattice_;
};

/** Samples LBM density (rho) on interior cells. */
template <typename T, typename DESCRIPTOR>
class VtkDensityField {
 public:
  explicit VtkDensityField(const BlockLattice<T, DESCRIPTOR>& lattice)
      : lattice_(lattice) {}

  std::string_view vtk_name() const { return "density"; }
  int vtk_components() const { return 1; }

  void sample_cell(int i, int j, int k, double* out) const {
    T rho{};
    T u[3]{};
    lattice_.get(i, j, k).computeRhoU(rho, u);
    out[0] = static_cast<double>(rho);
  }

 private:
  const BlockLattice<T, DESCRIPTOR>& lattice_;
};

/** Samples LBM pressure p = rho * cs^2 on interior cells. */
template <typename T, typename DESCRIPTOR>
class VtkPressureField {
 public:
  explicit VtkPressureField(const BlockLattice<T, DESCRIPTOR>& lattice)
      : lattice_(lattice) {}

  std::string_view vtk_name() const { return "pressure"; }
  int vtk_components() const { return 1; }

  void sample_cell(int i, int j, int k, double* out) const {
    T rho{};
    T u[3]{};
    lattice_.get(i, j, k).computeRhoU(rho, u);
    const T cs2 = T{1} / static_cast<T>(DESCRIPTOR::invCs2);
    out[0] = static_cast<double>(rho * cs2);
    (void)u;
  }

 private:
  const BlockLattice<T, DESCRIPTOR>& lattice_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_VTK_LBM_FIELDS_H_
