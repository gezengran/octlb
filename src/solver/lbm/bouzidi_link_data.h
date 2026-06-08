#ifndef OCTLB_SRC_SOLVER_LBM_BOUZIDI_LINK_DATA_H_
#define OCTLB_SRC_SOLVER_LBM_BOUZIDI_LINK_DATA_H_

#include <vector>

#include "src/common/types.h"
#include "src/mesh/geometry/geometry_types.h"
#include "src/mesh/geometry/material_field.h"

namespace octlb {

class OctreeForest;

/** Init-time cache of Bouzidi fractional link lengths per fluid cell/direction. */
class BouzidiLinkData {
 public:
  BouzidiLinkData() = default;
  BouzidiLinkData(label num_octants, int nx, int ny, int nz, int q);

  double q_frac(OctantId id, int ix, int iy, int iz, int iPop) const;
  void set_q_frac(OctantId id, int ix, int iy, int iz, int iPop, double q);

  static BouzidiLinkData Build(const OctreeForest& forest,
                               const MaterialField& material,
                               const GeometryAssembly& assembly);

 private:
  label num_octants_ = 0;
  int nx_ = 0;
  int ny_ = 0;
  int nz_ = 0;
  int q_ = 0;
  std::vector<double> q_frac_;
  std::size_t index(OctantId id, int ix, int iy, int iz, int iPop) const;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BOUZIDI_LINK_DATA_H_
