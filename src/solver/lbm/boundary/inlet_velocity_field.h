#ifndef OCTLB_SRC_SOLVER_LBM_BOUNDARY_INLET_VELOCITY_FIELD_H_
#define OCTLB_SRC_SOLVER_LBM_BOUNDARY_INLET_VELOCITY_FIELD_H_

#include <array>

#include "src/common/types.h"

namespace octlb {
namespace boundary {

// Per-cell, time-dependent inlet velocity field. Pure functional: no LBM
// state, no side effects. Consulted by DomainBoundaryHandler on velocity-type
// inlet faces (kInterpolatedVelocity / kMovingLid) to prescribe a spatially
// varying, possibly ramped velocity -- e.g. a Poiseuille duct profile.
//
// (ix, iy, iz) are interior lattice indices within the block that owns the
// inlet face; t is the current time in the same units as ramp_end_t (lattice
// steps when threaded from TimeLoop). Writes u[3].
class InletVelocityField {
 public:
  virtual ~InletVelocityField() = default;
  virtual void velocity(int ix, int iy, int iz, double t,
                        double u[3]) const = 0;
};

// Parabolic (Poiseuille) duct profile on a rectangular cross-section of the
// inlet face: u = ramp(t) * u_peak * parab(in0) * parab(in1) * flow_dir, where
// parab(i, N) = 4*c*(1-c) with c = (i + 0.5)/N the cell-centre normalised
// coordinate (0 at the wall, 1 at the cross-section centre). The two in-plane
// axes are determined by face_dir; flow_dir is taken as a unit vector.
//
// ramp_end_t <= 0 disables the ramp (full velocity from t = 0). The ramp shape
// for ramp_end_t > 0 is added with the ramp-time test (behaviour 2).
class PoiseuilleInletProfile : public InletVelocityField {
 public:
  PoiseuilleInletProfile(FaceDir face_dir, int n_inplane0, int n_inplane1,
                         double u_peak, std::array<double, 3> flow_dir,
                         double ramp_end_t = 0.0)
      : face_dir_(face_dir),
        n_inplane0_(n_inplane0),
        n_inplane1_(n_inplane1),
        u_peak_(u_peak),
        flow_dir_(flow_dir),
        ramp_end_t_(ramp_end_t) {}

  void velocity(int ix, int iy, int iz, double t, double u[3]) const override {
    const double r = Ramp(t);
    const double p0 = Parab(InPlane0(ix, iy, iz), n_inplane0_);
    const double p1 = Parab(InPlane1(ix, iy, iz), n_inplane1_);
    const double mag = r * u_peak_ * p0 * p1;
    for (int d = 0; d < 3; ++d) {
      u[d] = mag * flow_dir_[d];
    }
  }

 private:
  static double Parab(int i, int n) {
    const double c = (static_cast<double>(i) + 0.5) / static_cast<double>(n);
    return 4.0 * c * (1.0 - c);
  }

  // OpenLB PolynomialStartScale: smootherstep 10x^3 - 15x^4 + 6x^5 with
  // x = t / ramp_end_t, clamped to [0, 1]. ramp_end_t <= 0 disables the ramp.
  double Ramp(double t) const {
    if (ramp_end_t_ <= 0.0) {
      return 1.0;
    }
    double x = t / ramp_end_t_;
    if (x < 0.0) {
      x = 0.0;
    } else if (x > 1.0) {
      x = 1.0;
    }
    const double x2 = x * x;
    const double x3 = x2 * x;
    return x3 * (10.0 + x * (-15.0 + 6.0 * x));
  }

  // The two in-plane lattice indices for the inlet face (the face-normal index
  // is irrelevant to the profile and is ignored).
  int InPlane0(int ix, int iy, int /*iz*/) const {
    switch (face_dir_) {
      case FaceDir::kXMin:
      case FaceDir::kXMax:
        return iy;
      case FaceDir::kYMin:
      case FaceDir::kYMax:
        return ix;
      default:  // kZMin / kZMax
        return ix;
    }
  }

  int InPlane1(int /*ix*/, int iy, int iz) const {
    switch (face_dir_) {
      case FaceDir::kXMin:
      case FaceDir::kXMax:
        return iz;
      case FaceDir::kYMin:
      case FaceDir::kYMax:
        return iz;
      default:  // kZMin / kZMax
        return iy;
    }
  }

  FaceDir face_dir_;
  int n_inplane0_;
  int n_inplane1_;
  double u_peak_;
  std::array<double, 3> flow_dir_;
  double ramp_end_t_;
};

// Spatially uniform inlet velocity: u = ramp(t) * u_const. Used by cases that
// prescribe a constant free-stream inlet (e.g. cylinder3d Schäfer-Turek sanity)
// on a velocity-Dirichlet face, without a duct profile.
class UniformInletProfile : public InletVelocityField {
 public:
  UniformInletProfile(std::array<double, 3> u_const, double ramp_end_t = 0.0)
      : u_const_(u_const), ramp_end_t_(ramp_end_t) {}

  void velocity(int /*ix*/, int /*iy*/, int /*iz*/, double t,
                double u[3]) const override {
    const double r = Ramp(t);
    for (int d = 0; d < 3; ++d) {
      u[d] = r * u_const_[d];
    }
  }

 private:
  double Ramp(double t) const {
    if (ramp_end_t_ <= 0.0) {
      return 1.0;
    }
    double x = t / ramp_end_t_;
    if (x < 0.0) {
      x = 0.0;
    } else if (x > 1.0) {
      x = 1.0;
    }
    const double x2 = x * x;
    const double x3 = x2 * x;
    return x3 * (10.0 + x * (-15.0 + 6.0 * x));
  }

  std::array<double, 3> u_const_;
  double ramp_end_t_;
};

}  // namespace boundary
}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BOUNDARY_INLET_VELOCITY_FIELD_H_