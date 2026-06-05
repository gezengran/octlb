#ifndef OCTLB_SRC_MESH_GEOMETRY_MATERIAL_FIELD_H_
#define OCTLB_SRC_MESH_GEOMETRY_MATERIAL_FIELD_H_

#include <vector>

#include "src/common/types.h"
#include "src/mesh/geometry/geometry_types.h"

namespace octlb {

/** Per-octant Nx×Ny×Nz material tags (local octants only). */
class MaterialField {
 public:
  MaterialField() = default;
  MaterialField(label num_octants, int nx, int ny, int nz);

  MaterialKind at(OctantId id, int i, int j, int k) const;
  void set(OctantId id, int i, int j, int k, MaterialKind kind);

  label num_octants() const { return num_octants_; }
  int nx() const { return nx_; }
  int ny() const { return ny_; }
  int nz() const { return nz_; }

 private:
  label num_octants_ = 0;
  int nx_ = 0;
  int ny_ = 0;
  int nz_ = 0;
  std::vector<MaterialKind> data_;
  std::size_t index(OctantId id, int i, int j, int k) const;
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_MATERIAL_FIELD_H_
