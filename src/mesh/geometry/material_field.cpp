#include "src/mesh/geometry/material_field.h"

#include <stdexcept>

namespace octlb {

MaterialField::MaterialField(label num_octants, int nx, int ny, int nz)
    : num_octants_(num_octants),
      nx_(nx),
      ny_(ny),
      nz_(nz),
      data_(static_cast<std::size_t>(num_octants) * static_cast<std::size_t>(nx) *
                static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz),
            MaterialKind::kFluid) {}

std::size_t MaterialField::index(OctantId id, int i, int j, int k) const {
  if (id < 0 || id >= num_octants_) {
    throw std::out_of_range("MaterialField OctantId out of range");
  }
  if (i < 0 || i >= nx_ || j < 0 || j >= ny_ || k < 0 || k >= nz_) {
    throw std::out_of_range("MaterialField cell index out of range");
  }
  const std::size_t cells_per_oct =
      static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_) *
      static_cast<std::size_t>(nz_);
  return static_cast<std::size_t>(id) * cells_per_oct +
         static_cast<std::size_t>(k) * static_cast<std::size_t>(nx_) *
             static_cast<std::size_t>(ny_) +
         static_cast<std::size_t>(j) * static_cast<std::size_t>(nx_) +
         static_cast<std::size_t>(i);
}

MaterialKind MaterialField::at(OctantId id, int i, int j, int k) const {
  return data_[index(id, i, j, k)];
}

void MaterialField::set(OctantId id, int i, int j, int k, MaterialKind kind) {
  data_[index(id, i, j, k)] = kind;
}

}  // namespace octlb
