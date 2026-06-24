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
  kInterpolatedVelocity,
};

// InterpolatedVelocity overlap padding fill before PlaneFd (PostStream stage).
enum class OverlapPaddingMode {
  kPostStream,  // stream snapshot + addPoints2CommBC communicate (all faces)
  kMirror,      // mirror-from-core everywhere (diagnostic baseline)
  kHybrid,      // ymin/ymax exterior padding slabs stream-only; lateral mirror
};

struct DomainBcSpec {
  FaceDir face = FaceDir::kXMin;
  DomainBcType type = DomainBcType::kNoSlip;
  std::array<double, 3> u_wall{{0.0, 0.0, 0.0}};
};

class DomainBoundaryHandler {
 public:
  virtual ~DomainBoundaryHandler() = default;
  // Pre-stream: ghost fill (legacy) or Dirichlet collide on boundary lattice cells.
  virtual void apply(CollideRhoStats* rho_stats = nullptr,
                     double average_rho = 1.0,
                     bool use_const_rho_bgk = false) = 0;
  // Post-stream: InterpolatedVelocity on boundary lattice cells (OpenLB PostStream).
  virtual void apply_post_stream() {}
  // OpenLB Dominant collide walks iX,iY,iZ with boundary before inward fluid neighbors.
  virtual bool collide_boundary_before_bulk() const { return false; }
  // Spatial iX,iY,iZ collide matching OpenLB ConcreteBlockCollisionO::applyDominant.
  virtual void collide_interleaved_with(
      BlockLattice<double, olb::descriptors::D3Q19<>>& lat,
      CollideRhoStats* rho_stats, double average_rho, bool use_const_rho_bgk) {}
};

class NoOpDomainBoundaryHandler : public DomainBoundaryHandler {
 public:
  void apply(CollideRhoStats* /*rho_stats*/ = nullptr,
             double /*average_rho*/ = 1.0,
             bool /*use_const_rho_bgk*/ = false) override {}
};

using DomainBoundaryLattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

class ConcreteDomainBoundaryHandler : public DomainBoundaryHandler {
 public:
  ConcreteDomainBoundaryHandler(
      BlockCollection<DomainBoundaryLattice>& blocks,
      const std::vector<TreeBoundaryFace>& faces,
      const std::vector<DomainBcSpec>& specs, int nx, int ny, int nz,
      double omega = 1.0, bool boundary_lattice_mode = false,
      OverlapPaddingMode padding_mode = OverlapPaddingMode::kHybrid);

  void apply(CollideRhoStats* rho_stats = nullptr, double average_rho = 1.0,
             bool use_const_rho_bgk = false) override;
  void apply_post_stream() override;
  bool collide_boundary_before_bulk() const override;
  void collide_interleaved_with(
      BlockLattice<double, olb::descriptors::D3Q19<>>& lat,
      CollideRhoStats* rho_stats, double average_rho,
      bool use_const_rho_bgk) override;

 private:
  void ApplyLegacyFaceBc(DomainBoundaryLattice& lat, FaceDir dir,
                         const DomainBcSpec& spec);
  bool UsesInterpolatedVelocity() const;

  BlockCollection<DomainBoundaryLattice>& blocks_;
  std::vector<TreeBoundaryFace> faces_;
  std::vector<DomainBcSpec> specs_;
  int nx_;
  int ny_;
  int nz_;
  double omega_;
  bool boundary_lattice_mode_;
  OverlapPaddingMode padding_mode_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_DOMAIN_BOUNDARY_HANDLER_H_
