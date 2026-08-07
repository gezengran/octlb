#ifndef OCTLB_SRC_SOLVER_LBM_DOMAIN_BOUNDARY_HANDLER_H_
#define OCTLB_SRC_SOLVER_LBM_DOMAIN_BOUNDARY_HANDLER_H_

#include <array>
#include <memory>
#include <vector>

#include "src/common/types.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/inlet_velocity_field.h"

namespace octlb {

enum class DomainBcType {
  kNoSlip,
  kMovingLid,
  kInterpolatedVelocity,
  kOutflow,  // Zero-gradient (do-nothing) outflow BC (W1, T11).
  kInterpolatedPressure,  // Pressure-outlet FD (prescribed rho) BC (T11 W3).
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
  // Optional per-cell, time-dependent inlet velocity. When set, velocity-type
  // faces (kInterpolatedVelocity / kMovingLid) prescribe u per cell from the
  // field instead of the constant u_wall. Empty -> backward compatible.
  std::shared_ptr<boundary::InletVelocityField> inlet_field;
  // Prescribed outlet density for kInterpolatedPressure faces (p = cs^2 * (
  // rho - 1 ), so rho_target=1.0 is p=0). Default 1.0.
  double rho_target = 1.0;
};

// Per-cell prescribed velocity for a BC spec: the inlet_field's value at
// (ix, iy, iz, t) when set, otherwise the constant u_wall. Single source for
// both the legacy Zou-He and the InterpolatedVelocity BC paths. When the field
// is_physical() and a block phys_origin is supplied, evaluate from the cell's
// domain position so a physical-coordinate profile (e.g. channel Poiseuille)
// is correct on a channel patch carved into a larger face.
inline void PrescribedVelocity(const DomainBcSpec& spec, int ix, int iy,
                               int iz, double t, double u[3],
                               const double origin[3] = nullptr,
                               double cell_width = 0.0) {
  if (spec.inlet_field) {
    if (origin != nullptr && spec.inlet_field->is_physical()) {
      const double px = origin[0] + (ix + 0.5) * cell_width;
      const double py = origin[1] + (iy + 0.5) * cell_width;
      const double pz = origin[2] + (iz + 0.5) * cell_width;
      spec.inlet_field->velocity_phys(px, py, pz, t, u);
    } else {
      spec.inlet_field->velocity(ix, iy, iz, t, u);
    }
    return;
  }
  for (int d = 0; d < 3; ++d) {
    u[d] = spec.u_wall[d];
  }
}

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
  // Current simulation time (lattice steps when threaded from TimeLoop). Used
  // by inlet velocity fields for ramp-up. Default 0 -> no ramp effect.
  void set_time(double t) { current_time_ = t; }
  double current_time() const { return current_time_; }

 protected:
  double current_time_ = 0.0;
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
      double omega = 1.0,
      OverlapPaddingMode padding_mode = OverlapPaddingMode::kHybrid);

  // Per-cell dispatch is now the only path; apply() is retained for the flat
  // advance path and collides non-bulk boundary cells per their BcKind.
  void apply(CollideRhoStats* rho_stats = nullptr, double average_rho = 1.0,
             bool use_const_rho_bgk = false) override;
  void apply_post_stream() override;
  bool collide_boundary_before_bulk() const override;
  void collide_interleaved_with(
      BlockLattice<double, olb::descriptors::D3Q19<>>& lat,
      CollideRhoStats* rho_stats, double average_rho,
      bool use_const_rho_bgk) override;

 private:
  BlockCollection<DomainBoundaryLattice>& blocks_;
  std::vector<TreeBoundaryFace> faces_;
  std::vector<DomainBcSpec> specs_;
  int nx_;
  int ny_;
  int nz_;
  double omega_;
  OverlapPaddingMode padding_mode_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_DOMAIN_BOUNDARY_HANDLER_H_
