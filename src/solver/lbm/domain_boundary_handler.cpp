#include "src/solver/lbm/domain_boundary_handler.h"

#include "src/solver/lbm/boundary/bounce_back.h"
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

void ApplyFaceBc(DomainBoundaryLattice& lat, FaceDir dir, const DomainBcSpec& spec,
                 int nx, int ny, int nz) {
  const int face = static_cast<int>(dir);
  switch (dir) {
    case FaceDir::kXMin:
      for (int iy = 0; iy < ny; ++iy) {
        for (int iz = 0; iz < nz; ++iz) {
          double* ghost = lat.populations_at_halo(0, iy + 1, iz + 1);
          const double* interior = lat.populations_at_halo(1, iy + 1, iz + 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else {
            boundary::ApplyZouHeVelocityGhost<double, Descriptor>(
                ghost, interior, face, spec.u_wall.data());
          }
        }
      }
      break;
    case FaceDir::kXMax:
      for (int iy = 0; iy < ny; ++iy) {
        for (int iz = 0; iz < nz; ++iz) {
          double* ghost = lat.populations_at_halo(nx + 1, iy + 1, iz + 1);
          const double* interior =
              lat.populations_at_halo(nx, iy + 1, iz + 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else {
            boundary::ApplyZouHeVelocityGhost<double, Descriptor>(
                ghost, interior, face, spec.u_wall.data());
          }
        }
      }
      break;
    case FaceDir::kYMin:
      for (int ix = 0; ix < nx; ++ix) {
        for (int iz = 0; iz < nz; ++iz) {
          double* ghost = lat.populations_at_halo(ix + 1, 0, iz + 1);
          const double* interior = lat.populations_at_halo(ix + 1, 1, iz + 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else {
            boundary::ApplyZouHeVelocityGhost<double, Descriptor>(
                ghost, interior, face, spec.u_wall.data());
          }
        }
      }
      break;
    case FaceDir::kYMax:
      for (int ix = 0; ix < nx; ++ix) {
        for (int iz = 0; iz < nz; ++iz) {
          double* ghost = lat.populations_at_halo(ix + 1, ny + 1, iz + 1);
          const double* interior = lat.populations_at_halo(ix + 1, ny, iz + 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else {
            boundary::ApplyZouHeVelocityGhost<double, Descriptor>(
                ghost, interior, face, spec.u_wall.data());
          }
        }
      }
      break;
    case FaceDir::kZMin:
      for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
          double* ghost = lat.populations_at_halo(ix + 1, iy + 1, 0);
          const double* interior = lat.populations_at_halo(ix + 1, iy + 1, 1);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else {
            boundary::ApplyZouHeVelocityGhost<double, Descriptor>(
                ghost, interior, face, spec.u_wall.data());
          }
        }
      }
      break;
    case FaceDir::kZMax:
      for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
          double* ghost = lat.populations_at_halo(ix + 1, iy + 1, nz + 1);
          const double* interior = lat.populations_at_halo(ix + 1, iy + 1, nz);
          if (spec.type == DomainBcType::kNoSlip) {
            boundary::ApplyNoSlipGhost<double, Descriptor>(ghost, interior);
          } else {
            boundary::ApplyZouHeVelocityGhost<double, Descriptor>(
                ghost, interior, face, spec.u_wall.data());
          }
        }
      }
      break;
  }
}

}  // namespace

ConcreteDomainBoundaryHandler::ConcreteDomainBoundaryHandler(
    BlockCollection<DomainBoundaryLattice>& blocks,
    const std::vector<TreeBoundaryFace>& faces,
    const std::vector<DomainBcSpec>& specs, int nx, int ny, int nz)
    : blocks_(blocks),
      faces_(faces),
      specs_(specs),
      nx_(nx),
      ny_(ny),
      nz_(nz) {}

void ConcreteDomainBoundaryHandler::apply() {
  for (const TreeBoundaryFace& face : faces_) {
    DomainBoundaryLattice& lat = blocks_[face.octant_id];
    const DomainBcSpec spec = FindSpec(specs_, face.face_dir);
    ApplyFaceBc(lat, face.face_dir, spec, nx_, ny_, nz_);
  }
}

}  // namespace octlb
