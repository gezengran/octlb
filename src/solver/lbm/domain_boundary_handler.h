#ifndef OCTLB_SRC_SOLVER_LBM_DOMAIN_BOUNDARY_HANDLER_H_
#define OCTLB_SRC_SOLVER_LBM_DOMAIN_BOUNDARY_HANDLER_H_

#include <array>
#include <vector>

#include "src/common/types.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/block_lattice.h"

namespace octlb {

enum class DomainBcType {
  kNoSlip,
  kMovingLid,
};

struct DomainBcSpec {
  FaceDir face = FaceDir::kXMin;
  DomainBcType type = DomainBcType::kNoSlip;
  std::array<double, 3> u_wall{{0.0, 0.0, 0.0}};
};

class DomainBoundaryHandler {
 public:
  virtual ~DomainBoundaryHandler() = default;
  virtual void apply() = 0;
};

class NoOpDomainBoundaryHandler : public DomainBoundaryHandler {
 public:
  void apply() override {}
};

using DomainBoundaryLattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

class ConcreteDomainBoundaryHandler : public DomainBoundaryHandler {
 public:
  ConcreteDomainBoundaryHandler(
      BlockCollection<DomainBoundaryLattice>& blocks,
      const std::vector<TreeBoundaryFace>& faces,
      const std::vector<DomainBcSpec>& specs, int nx, int ny, int nz);

  void apply() override;

 private:
  BlockCollection<DomainBoundaryLattice>& blocks_;
  std::vector<TreeBoundaryFace> faces_;
  std::vector<DomainBcSpec> specs_;
  int nx_;
  int ny_;
  int nz_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_DOMAIN_BOUNDARY_HANDLER_H_
