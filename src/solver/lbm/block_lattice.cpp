#include "block_lattice.h"

#include "src/solver/lbm/boundary/bouzidi_pull.h"
#include "src/solver/lbm/bouzidi_link_data.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <stdexcept>

namespace octlb {

template <typename T, typename DESCRIPTOR>
int BlockLattice<T, DESCRIPTOR>::face_buffer_count(int nx, int ny, int nz,
                                                   FaceDir dir) {
  switch (dir) {
    case FaceDir::kXMin:
    case FaceDir::kXMax:
      return ny * nz * kQ;
    case FaceDir::kYMin:
    case FaceDir::kYMax:
      return nx * nz * kQ;
    case FaceDir::kZMin:
    case FaceDir::kZMax:
      return nx * ny * kQ;
  }
  return 0;
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::pack_face(FaceDir dir, T* buffer,
                                            int count) const {
  const int expected = face_buffer_count(nx_, ny_, nz_, dir);
  if (count < expected) {
    throw std::runtime_error("BlockLattice::pack_face: buffer too small");
  }
  int out = 0;
  const auto copy_cell = [&](int hx, int hy, int hz) {
    const int base = halo_idx(hx, hy, hz);
    for (int iPop = 0; iPop < kQ; ++iPop) {
      buffer[out++] = populations_[base + iPop];
    }
  };

  switch (dir) {
    case FaceDir::kXMin:
      for (int iy = 0; iy < ny_; ++iy)
        for (int iz = 0; iz < nz_; ++iz) {
          copy_cell(h_, iy + h_, iz + h_);
        }
      break;
    case FaceDir::kXMax:
      for (int iy = 0; iy < ny_; ++iy)
        for (int iz = 0; iz < nz_; ++iz) {
          copy_cell(nx_ + h_ - 1, iy + h_, iz + h_);
        }
      break;
    case FaceDir::kYMin:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iz = 0; iz < nz_; ++iz) {
          copy_cell(ix + h_, h_, iz + h_);
        }
      break;
    case FaceDir::kYMax:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iz = 0; iz < nz_; ++iz) {
          copy_cell(ix + h_, ny_ + h_ - 1, iz + h_);
        }
      break;
    case FaceDir::kZMin:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iy = 0; iy < ny_; ++iy) {
          copy_cell(ix + h_, iy + h_, h_);
        }
      break;
    case FaceDir::kZMax:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iy = 0; iy < ny_; ++iy) {
          copy_cell(ix + h_, iy + h_, nz_ + h_ - 1);
        }
      break;
  }
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::read_ghost_face(FaceDir dir, T* buffer,
                                                  int count) const {
  const int expected = face_buffer_count(nx_, ny_, nz_, dir);
  if (count < expected) {
    throw std::runtime_error("BlockLattice::read_ghost_face: buffer too small");
  }
  int in = 0;
  const auto read_cell = [&](int hx, int hy, int hz) {
    const int base = halo_idx(hx, hy, hz);
    for (int iPop = 0; iPop < kQ; ++iPop) {
      buffer[in++] = populations_[base + iPop];
    }
  };

  switch (dir) {
    case FaceDir::kXMin:
      for (int iy = 0; iy < ny_; ++iy)
        for (int iz = 0; iz < nz_; ++iz) {
          read_cell(h_ - 1, iy + h_, iz + h_);
        }
      break;
    case FaceDir::kXMax:
      for (int iy = 0; iy < ny_; ++iy)
        for (int iz = 0; iz < nz_; ++iz) {
          read_cell(nx_ + h_, iy + h_, iz + h_);
        }
      break;
    case FaceDir::kYMin:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iz = 0; iz < nz_; ++iz) {
          read_cell(ix + h_, h_ - 1, iz + h_);
        }
      break;
    case FaceDir::kYMax:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iz = 0; iz < nz_; ++iz) {
          read_cell(ix + h_, ny_ + h_, iz + h_);
        }
      break;
    case FaceDir::kZMin:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iy = 0; iy < ny_; ++iy) {
          read_cell(ix + h_, iy + h_, h_ - 1);
        }
      break;
    case FaceDir::kZMax:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iy = 0; iy < ny_; ++iy) {
          read_cell(ix + h_, iy + h_, nz_ + h_);
        }
      break;
  }
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::unpack_face(FaceDir dir, const T* buffer,
                                              int count) {
  const int expected = face_buffer_count(nx_, ny_, nz_, dir);
  if (count < expected) {
    throw std::runtime_error("BlockLattice::unpack_face: buffer too small");
  }
  int in = 0;
  const auto paste_cell = [&](int hx, int hy, int hz) {
    const int base = halo_idx(hx, hy, hz);
    for (int iPop = 0; iPop < kQ; ++iPop) {
      populations_[base + iPop] = buffer[in++];
    }
  };

  switch (dir) {
    case FaceDir::kXMin:
      for (int iy = 0; iy < ny_; ++iy)
        for (int iz = 0; iz < nz_; ++iz) {
          paste_cell(h_ - 1, iy + h_, iz + h_);
        }
      break;
    case FaceDir::kXMax:
      for (int iy = 0; iy < ny_; ++iy)
        for (int iz = 0; iz < nz_; ++iz) {
          paste_cell(nx_ + h_, iy + h_, iz + h_);
        }
      break;
    case FaceDir::kYMin:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iz = 0; iz < nz_; ++iz) {
          paste_cell(ix + h_, h_ - 1, iz + h_);
        }
      break;
    case FaceDir::kYMax:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iz = 0; iz < nz_; ++iz) {
          paste_cell(ix + h_, ny_ + h_, iz + h_);
        }
      break;
    case FaceDir::kZMin:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iy = 0; iy < ny_; ++iy) {
          paste_cell(ix + h_, iy + h_, h_ - 1);
        }
      break;
    case FaceDir::kZMax:
      for (int ix = 0; ix < nx_; ++ix)
        for (int iy = 0; iy < ny_; ++iy) {
          paste_cell(ix + h_, iy + h_, nz_ + h_);
        }
      break;
  }
}

template <typename T, typename DESCRIPTOR>
BlockLattice<T, DESCRIPTOR>::BlockLattice(int nx, int ny, int nz, int halo)
    : nx_(nx), ny_(ny), nz_(nz), h_(halo),
      populations_((nx + 2*halo) * (ny + 2*halo) * (nz + 2*halo) * kQ,
                   T{0}),
      cell_kinds_(static_cast<std::size_t>(nx * ny * nz), CellKind::kFluid) {}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::initialize(T rho0, const T* u0) {
  T uSqr = T{0};
  for (int d = 0; d < DESCRIPTOR::d; ++d) uSqr += u0[d] * u0[d];

  for (int ix = 0; ix < nx_; ++ix)
    for (int iy = 0; iy < ny_; ++iy)
      for (int iz = 0; iz < nz_; ++iz) {
        T* f = &populations_[idx(ix, iy, iz)];
        for (int iPop = 0; iPop < kQ; ++iPop)
          f[iPop] = olb::equilibrium<DESCRIPTOR>::secondOrder(
              iPop, rho0, u0, uSqr);
      }
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::set_cell_kind(int ix, int iy, int iz,
                                                  CellKind kind) {
  cell_kinds_[static_cast<std::size_t>((ix * ny_ + iy) * nz_ + iz)] = kind;
}

template <typename T, typename DESCRIPTOR>
CellKind BlockLattice<T, DESCRIPTOR>::cell_kind(int ix, int iy, int iz) const {
  return cell_kinds_[static_cast<std::size_t>((ix * ny_ + iy) * nz_ + iz)];
}

template <typename T, typename DESCRIPTOR>
T* BlockLattice<T, DESCRIPTOR>::populations_at_halo(int hx, int hy, int hz) {
  return &populations_[halo_idx(hx, hy, hz)];
}

template <typename T, typename DESCRIPTOR>
const T* BlockLattice<T, DESCRIPTOR>::populations_at_halo(int hx, int hy,
                                                          int hz) const {
  return &populations_[halo_idx(hx, hy, hz)];
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::collide(T omega) {
  for (int ix = 0; ix < nx_; ++ix)
    for (int iy = 0; iy < ny_; ++iy)
      for (int iz = 0; iz < nz_; ++iz) {
        const CellKind kind = cell_kind(ix, iy, iz);
        if (kind == CellKind::kSolid || kind == CellKind::kBoundary) {
          continue;
        }
        auto cell = get(ix, iy, iz);
        T rho{}, u[DESCRIPTOR::d]{};
        cell.computeRhoU(rho, u);

        T uSqr = T{0};
        for (int d = 0; d < DESCRIPTOR::d; ++d) uSqr += u[d] * u[d];

        for (int iPop = 0; iPop < kQ; ++iPop) {
          const T fEq = olb::equilibrium<DESCRIPTOR>::secondOrder(
              iPop, rho, u, uSqr);
          cell[iPop] += omega * (fEq - cell[iPop]);
        }
      }
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::stream() {
  const int NX = nx_ + 2*h_;
  const int NY = ny_ + 2*h_;
  const int NZ = nz_ + 2*h_;

  std::vector<T> tmp(populations_.size());

  auto kind_at_halo = [&](int hx, int hy, int hz) -> CellKind {
    const int ix = hx - h_;
    const int iy = hy - h_;
    const int iz = hz - h_;
    if (ix < 0 || ix >= nx_ || iy < 0 || iy >= ny_ || iz < 0 || iz >= nz_) {
      return CellKind::kFluid;
    }
    return cell_kind(ix, iy, iz);
  };

  for (int ix = h_; ix < nx_ + h_; ++ix)
    for (int iy = h_; iy < ny_ + h_; ++iy)
      for (int iz = h_; iz < nz_ + h_; ++iz) {
        const int ix_int = ix - h_;
        const int iy_int = iy - h_;
        const int iz_int = iz - h_;
        const CellKind kind = cell_kind(ix_int, iy_int, iz_int);
        if (kind == CellKind::kSolid || kind == CellKind::kBoundary) {
          continue;
        }
        const int dst = (ix * NY * NZ + iy * NZ + iz) * kQ;
        for (int iPop = 0; iPop < kQ; ++iPop) {
          const int sx = ix - olb::descriptors::c<DESCRIPTOR>(iPop, 0);
          const int sy = iy - olb::descriptors::c<DESCRIPTOR>(iPop, 1);
          const int sz = iz - olb::descriptors::c<DESCRIPTOR>(iPop, 2);
          const int src = (sx * NY * NZ + sy * NZ + sz) * kQ + iPop;

          const CellKind src_kind = kind_at_halo(sx, sy, sz);
          if (bouzidi_ != nullptr &&
              (src_kind == CellKind::kSolid || src_kind == CellKind::kBoundary)) {
            const double q_frac =
                bouzidi_->q_frac(octant_id_, ix_int, iy_int, iz_int, iPop);
            if (q_frac > 1.0e-10 && q_frac < 1.0 - 1.0e-10) {
              const int opp = olb::descriptors::opposite<DESCRIPTOR>(iPop);
              const T f_bb_opp = populations_[dst + opp];
              T f_interior_q = populations_[dst + iPop];
              const int nx_h = ix + olb::descriptors::c<DESCRIPTOR>(iPop, 0);
              const int ny_h = iy + olb::descriptors::c<DESCRIPTOR>(iPop, 1);
              const int nz_h = iz + olb::descriptors::c<DESCRIPTOR>(iPop, 2);
              const int nbase =
                  (nx_h * NY * NZ + ny_h * NZ + nz_h) * kQ + iPop;
              if (kind_at_halo(nx_h, ny_h, nz_h) == CellKind::kFluid) {
                f_interior_q = populations_[nbase];
              }
              tmp[dst + iPop] = boundary::BouzidiPostCollisionPull(
                  f_bb_opp, f_interior_q, populations_[dst + iPop], q_frac);
              continue;
            }
            const int opp = olb::descriptors::opposite<DESCRIPTOR>(iPop);
            tmp[dst + iPop] = populations_[dst + opp];
            continue;
          }

          tmp[dst + iPop] = populations_[src];
        }
      }

  // Copy streamed interior back (leave ghost cells unchanged).
  for (int ix = h_; ix < nx_ + h_; ++ix)
    for (int iy = h_; iy < ny_ + h_; ++iy)
      for (int iz = h_; iz < nz_ + h_; ++iz) {
        const int base = (ix * NY * NZ + iy * NZ + iz) * kQ;
        for (int iPop = 0; iPop < kQ; ++iPop)
          populations_[base + iPop] = tmp[base + iPop];
      }
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::fill_periodic_halo() {
  // Copy the outermost real layer into the adjacent ghost layer on each face.
  // After this, stream() can pull from ghost positions and get periodic values.
  const int NX = nx_ + 2*h_;
  const int NY = ny_ + 2*h_;
  const int NZ = nz_ + 2*h_;

  auto cell_base = [&](int ix, int iy, int iz) -> int {
    return (ix * NY * NZ + iy * NZ + iz) * kQ;
  };

  // x faces
  for (int iy = h_; iy < ny_ + h_; ++iy)
    for (int iz = h_; iz < nz_ + h_; ++iz) {
      // ghost at 0 ← real at nx_ (right boundary wraps to left)
      std::memcpy(&populations_[cell_base(0, iy, iz)],
                  &populations_[cell_base(nx_, iy, iz)],
                  kQ * sizeof(T));
      // ghost at nx_+1 ← real at 1 (left boundary wraps to right)
      std::memcpy(&populations_[cell_base(nx_ + 1, iy, iz)],
                  &populations_[cell_base(1, iy, iz)],
                  kQ * sizeof(T));
    }

  // y faces (after x halo is populated)
  for (int ix = 0; ix < NX; ++ix)
    for (int iz = h_; iz < nz_ + h_; ++iz) {
      std::memcpy(&populations_[cell_base(ix, 0, iz)],
                  &populations_[cell_base(ix, ny_, iz)],
                  kQ * sizeof(T));
      std::memcpy(&populations_[cell_base(ix, ny_ + 1, iz)],
                  &populations_[cell_base(ix, 1, iz)],
                  kQ * sizeof(T));
    }

  // z faces (after x and y halos are populated)
  for (int ix = 0; ix < NX; ++ix)
    for (int iy = 0; iy < NY; ++iy) {
      std::memcpy(&populations_[cell_base(ix, iy, 0)],
                  &populations_[cell_base(ix, iy, nz_)],
                  kQ * sizeof(T));
      std::memcpy(&populations_[cell_base(ix, iy, nz_ + 1)],
                  &populations_[cell_base(ix, iy, 1)],
                  kQ * sizeof(T));
    }
}

template <typename T, typename DESCRIPTOR>
CellProxy<T, DESCRIPTOR>
BlockLattice<T, DESCRIPTOR>::get(int ix, int iy, int iz) {
  return CellProxy<T, DESCRIPTOR>(&populations_[idx(ix, iy, iz)]);
}

// ── Explicit instantiations ───────────────────────────────────────────────────
// All TUs that use BlockLattice must match one of these combinations.
template class BlockLattice<double, olb::descriptors::D3Q19<>>;
template class BlockLattice<float,  olb::descriptors::D3Q19<>>;

}  // namespace octlb
