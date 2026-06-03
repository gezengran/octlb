#include "block_lattice.h"

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
                   T{0}) {}

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
void BlockLattice<T, DESCRIPTOR>::collide(T omega) {
  for (int ix = 0; ix < nx_; ++ix)
    for (int iy = 0; iy < ny_; ++iy)
      for (int iz = 0; iz < nz_; ++iz) {
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
  // Pull-scheme: f_new[x][iPop] = f_old[x - c[iPop]][iPop].
  // Reads from [h_-1 .. Nx+h_] (includes ghost halo on both sides).
  const int NX = nx_ + 2*h_;
  const int NY = ny_ + 2*h_;
  const int NZ = nz_ + 2*h_;

  std::vector<T> tmp(populations_.size());

  for (int ix = h_; ix < nx_ + h_; ++ix)
    for (int iy = h_; iy < ny_ + h_; ++iy)
      for (int iz = h_; iz < nz_ + h_; ++iz) {
        const int dst = (ix * NY * NZ + iy * NZ + iz) * kQ;
        for (int iPop = 0; iPop < kQ; ++iPop) {
          const int sx = ix - olb::descriptors::c<DESCRIPTOR>(iPop, 0);
          const int sy = iy - olb::descriptors::c<DESCRIPTOR>(iPop, 1);
          const int sz = iz - olb::descriptors::c<DESCRIPTOR>(iPop, 2);
          const int src = (sx * NY * NZ + sy * NZ + sz) * kQ + iPop;
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
