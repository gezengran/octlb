#include "core/platform/platform.h"
#include "descriptor/descriptor.h"

#include "src/solver/lbm/unit_converter/unit_converter.h"

namespace octlb {

UnitConverter::UnitConverter(int resolution, scalar lattice_relaxation_time,
                             scalar char_phys_length, scalar char_phys_velocity,
                             scalar phys_viscosity, scalar phys_density)
    : resolution_(resolution),
      lattice_relaxation_time_(lattice_relaxation_time),
      char_phys_length_(char_phys_length),
      char_phys_velocity_(char_phys_velocity),
      phys_viscosity_(phys_viscosity),
      phys_density_(phys_density) {
  using Descriptor = olb::descriptors::D3Q19<>;
  const scalar inv_cs2 = olb::descriptors::invCs2<scalar, Descriptor>();

  phys_delta_x_ = char_phys_length_ / static_cast<scalar>(resolution_);
  phys_delta_t_ = (lattice_relaxation_time_ - 0.5) / inv_cs2 * phys_delta_x_ *
                  phys_delta_x_ / phys_viscosity_;
  conversion_velocity_ = phys_delta_x_ / phys_delta_t_;
  char_lattice_velocity_ = char_phys_velocity_ / conversion_velocity_;
  omega_ = 1.0 / lattice_relaxation_time_;
}

UnitConverter UnitConverter::OpenLbCavity3dDefaults() {
  return UnitConverter(/*resolution=*/30, /*tau=*/0.509,
                       /*L=*/1.0, /*U=*/1.0, /*nu=*/0.001, /*rho=*/1.0);
}

}  // namespace octlb
