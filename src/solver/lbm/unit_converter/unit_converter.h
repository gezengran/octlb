#ifndef OCTLB_SRC_SOLVER_LBM_UNIT_CONVERTER_UNIT_CONVERTER_H_
#define OCTLB_SRC_SOLVER_LBM_UNIT_CONVERTER_UNIT_CONVERTER_H_

#include <cstddef>

#include "src/common/types.h"

namespace octlb {

/** OpenLB `UnitConverterFromResolutionAndRelaxationTime` formulas (D3Q19). */
class UnitConverter {
 public:
  UnitConverter(int resolution, scalar lattice_relaxation_time,
                scalar char_phys_length, scalar char_phys_velocity,
                scalar phys_viscosity, scalar phys_density);

  int resolution() const { return resolution_; }
  scalar lattice_relaxation_time() const { return lattice_relaxation_time_; }
  scalar char_phys_length() const { return char_phys_length_; }
  scalar char_phys_velocity() const { return char_phys_velocity_; }
  scalar phys_viscosity() const { return phys_viscosity_; }
  scalar phys_density() const { return phys_density_; }

  scalar phys_delta_x() const { return phys_delta_x_; }
  scalar phys_delta_t() const { return phys_delta_t_; }
  scalar conversion_velocity() const { return conversion_velocity_; }
  scalar char_lattice_velocity() const { return char_lattice_velocity_; }
  scalar omega() const { return omega_; }

  scalar reynolds() const {
    return char_phys_velocity_ * char_phys_length_ / phys_viscosity_;
  }

  std::size_t get_lattice_time(scalar phys_t) const {
    return static_cast<std::size_t>(phys_t / phys_delta_t_ + 0.5);
  }

  /** OpenLB `examples/laminar/cavity3d` default parameters. */
  static UnitConverter OpenLbCavity3dDefaults();

 private:
  int resolution_ = 0;
  scalar lattice_relaxation_time_ = 0.0;
  scalar char_phys_length_ = 0.0;
  scalar char_phys_velocity_ = 0.0;
  scalar phys_viscosity_ = 0.0;
  scalar phys_density_ = 0.0;

  scalar phys_delta_x_ = 0.0;
  scalar phys_delta_t_ = 0.0;
  scalar conversion_velocity_ = 0.0;
  scalar char_lattice_velocity_ = 0.0;
  scalar omega_ = 0.0;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_UNIT_CONVERTER_UNIT_CONVERTER_H_
