#include "src/solver/lbm/domain_boundary_handler.h"

#include <vector>

#include "dynamics/lbm.h"
#include "src/solver/lbm/boundary/bounce_back.h"
#include "src/solver/lbm/boundary/interpolated_velocity.h"
#include "src/solver/lbm/boundary/outflow_bc.h"
#include "src/solver/lbm/boundary/zou_he_velocity.h"

namespace octlb {
namespace {

using Descriptor = olb::descriptors::D3Q19<>;

DomainBcSpec FindSpec(const std::vector<DomainBcSpec>& specs, FaceDir dir) {
  for (const DomainBcSpec& spec : specs) {
    if (spec.face == dir) {
      return spec;
    }
  }
  return DomainBcSpec{};
}

}  // namespace

ConcreteDomainBoundaryHandler::ConcreteDomainBoundaryHandler(
    BlockCollection<DomainBoundaryLattice>& blocks,
    const std::vector<TreeBoundaryFace>& faces,
    const std::vector<DomainBcSpec>& specs, int nx, int ny, int nz,
    double omega, bool boundary_lattice_mode, OverlapPaddingMode padding_mode)
    : blocks_(blocks),
      faces_(faces),
      specs_(specs),
      nx_(nx),
      ny_(ny),
      nz_(nz),
      omega_(omega),
      boundary_lattice_mode_(boundary_lattice_mode),
      padding_mode_(padding_mode) {}

bool ConcreteDomainBoundaryHandler::UsesInterpolatedVelocity() const {
  for (const DomainBcSpec& spec : specs_) {
    if (spec.type == DomainBcType::kInterpolatedVelocity) {
      return true;
    }
  }
  return false;
}

void ConcreteDomainBoundaryHandler::ApplyLegacyFaceBc(
    DomainBoundaryLattice& lat, FaceDir dir, const DomainBcSpec& spec) {
  const int face = static_cast<int>(dir);
  switch (dir) {
    case FaceDir::kXMin:
      for (int iy = 0; iy < ny_; ++iy) {
        for (int iz = 0; iz < nz_; ++iz) {
          double* ghost = lat.populations_at_halo(0, iy + 1, iz + 1);
          const double* interior = lat.populations_at_halo(1, iy + 1, iz + 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else if (spec.type == DomainBcType::kOutflow) {
            boundary::ApplyOutflowGhost<double, Descriptor>(ghost, interior);
          } else {
            ApplyVelocityGhost(spec, /*ix=*/0, iy, iz, ghost, interior, face);
          }
        }
      }
      break;
    case FaceDir::kXMax:
      for (int iy = 0; iy < ny_; ++iy) {
        for (int iz = 0; iz < nz_; ++iz) {
          double* ghost = lat.populations_at_halo(nx_ + 1, iy + 1, iz + 1);
          const double* interior =
              lat.populations_at_halo(nx_, iy + 1, iz + 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else if (spec.type == DomainBcType::kOutflow) {
            boundary::ApplyOutflowGhost<double, Descriptor>(ghost, interior);
          } else {
            ApplyVelocityGhost(spec, /*ix=*/nx_ - 1, iy, iz, ghost, interior,
                               face);
          }
        }
      }
      break;
    case FaceDir::kYMin:
      for (int ix = 0; ix < nx_; ++ix) {
        for (int iz = 0; iz < nz_; ++iz) {
          double* ghost = lat.populations_at_halo(ix + 1, 0, iz + 1);
          const double* interior = lat.populations_at_halo(ix + 1, 1, iz + 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else if (spec.type == DomainBcType::kOutflow) {
            boundary::ApplyOutflowGhost<double, Descriptor>(ghost, interior);
          } else {
            ApplyVelocityGhost(spec, ix, /*iy=*/0, iz, ghost, interior, face);
          }
        }
      }
      break;
    case FaceDir::kYMax:
      for (int ix = 0; ix < nx_; ++ix) {
        for (int iz = 0; iz < nz_; ++iz) {
          double* ghost = lat.populations_at_halo(ix + 1, ny_ + 1, iz + 1);
          const double* interior = lat.populations_at_halo(ix + 1, ny_, iz + 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else if (spec.type == DomainBcType::kOutflow) {
            boundary::ApplyOutflowGhost<double, Descriptor>(ghost, interior);
          } else {
            ApplyVelocityGhost(spec, ix, /*iy=*/ny_ - 1, iz, ghost, interior,
                               face);
          }
        }
      }
      break;
    case FaceDir::kZMin:
      for (int ix = 0; ix < nx_; ++ix) {
        for (int iy = 0; iy < ny_; ++iy) {
          double* ghost = lat.populations_at_halo(ix + 1, iy + 1, 0);
          const double* interior = lat.populations_at_halo(ix + 1, iy + 1, 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else if (spec.type == DomainBcType::kOutflow) {
            boundary::ApplyOutflowGhost<double, Descriptor>(ghost, interior);
          } else {
            ApplyVelocityGhost(spec, ix, iy, /*iz=*/0, ghost, interior, face);
          }
        }
      }
      break;
    case FaceDir::kZMax:
      for (int ix = 0; ix < nx_; ++ix) {
        for (int iy = 0; iy < ny_; ++iy) {
          double* ghost = lat.populations_at_halo(ix + 1, iy + 1, nz_ + 1);
          const double* interior = lat.populations_at_halo(ix + 1, iy + 1, nz_);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else if (spec.type == DomainBcType::kOutflow) {
            boundary::ApplyOutflowGhost<double, Descriptor>(ghost, interior);
          } else {
            ApplyVelocityGhost(spec, ix, iy, /*iz=*/nz_ - 1, ghost, interior,
                               face);
          }
        }
      }
      break;
  }
}

void ConcreteDomainBoundaryHandler::ApplyVelocityGhost(
    const DomainBcSpec& spec, int ix, int iy, int iz, double* ghost,
    const double* interior, int face) {
  double u[3] = {};
  PrescribedVelocity(spec, ix, iy, iz, current_time_, u);
  boundary::ApplyZouHeVelocityGhost<double, Descriptor>(ghost, interior, face,
                                                         u);
}

void ConcreteDomainBoundaryHandler::apply(CollideRhoStats* rho_stats,
                                          double average_rho,
                                          bool use_const_rho_bgk) {
  if (UsesInterpolatedVelocity() && boundary_lattice_mode_) {
    // One octant has six TreeBoundaryFaces; lattice-mode BC updates all
    // boundary cells in the block, so apply at most once per octant_id.
    std::vector<bool> seen(static_cast<std::size_t>(blocks_.size()), false);
    for (const TreeBoundaryFace& face : faces_) {
      if (seen[static_cast<std::size_t>(face.octant_id)]) {
        continue;
      }
      seen[static_cast<std::size_t>(face.octant_id)] = true;
      DomainBoundaryLattice& lat = blocks_[face.octant_id];
      boundary::CollideDirichletBoundaryCells<double, Descriptor>(
          lat, nx_, ny_, nz_, omega_, specs_, rho_stats, average_rho,
          use_const_rho_bgk);
    }
    return;
  }

  for (const TreeBoundaryFace& face : faces_) {
    DomainBoundaryLattice& lat = blocks_[face.octant_id];
    const DomainBcSpec spec = FindSpec(specs_, face.face_dir);
    ApplyLegacyFaceBc(lat, face.face_dir, spec);
  }
}

bool ConcreteDomainBoundaryHandler::collide_boundary_before_bulk() const {
  return UsesInterpolatedVelocity() && boundary_lattice_mode_;
}

void ConcreteDomainBoundaryHandler::collide_interleaved_with(
    BlockLattice<double, olb::descriptors::D3Q19<>>& lat,
    CollideRhoStats* rho_stats, double average_rho, bool use_const_rho_bgk) {
  if (!UsesInterpolatedVelocity() || !boundary_lattice_mode_) {
    return;
  }
  const double omega = omega_;
  const double avg = average_rho;
  for (int ix = 0; ix < nx_; ++ix) {
    for (int iy = 0; iy < ny_; ++iy) {
      for (int iz = 0; iz < nz_; ++iz) {
        const CellKind kind = lat.cell_kind(ix, iy, iz);
        if (kind == CellKind::kFluid) {
          if (use_const_rho_bgk) {
            lat.collide_const_rho_at(ix, iy, iz, omega, avg, rho_stats);
          } else {
            auto cell = lat.get(ix, iy, iz);
            double rho = 0.0;
            double u[3]{};
            cell.computeRhoU(rho, u);
            olb::lbm<Descriptor>::bgkCollision(cell, rho, u, omega);
          }
        } else if (kind == CellKind::kBoundary) {
          boundary::CollideDirichletBoundaryCellAt<double, Descriptor>(
              lat, ix, iy, iz, nx_, ny_, nz_, omega, specs_, rho_stats);
        }
      }
    }
  }
  lat.collide_overlap_padding_bgk(omega);
}

void ConcreteDomainBoundaryHandler::apply_post_stream() {
  if (!UsesInterpolatedVelocity() || !boundary_lattice_mode_) {
    return;
  }
  std::vector<bool> seen(static_cast<std::size_t>(blocks_.size()), false);
  for (const TreeBoundaryFace& face : faces_) {
    if (seen[static_cast<std::size_t>(face.octant_id)]) {
      continue;
    }
    seen[static_cast<std::size_t>(face.octant_id)] = true;
    DomainBoundaryLattice& lat = blocks_[face.octant_id];
    boundary::ApplyInterpolatedVelocityBoundaryCells<double, Descriptor>(
        lat, nx_, ny_, nz_, omega_, specs_, padding_mode_);
  }
}

}  // namespace octlb
