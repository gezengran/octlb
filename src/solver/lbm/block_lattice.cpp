#include "block_lattice.h"

#include <cassert>
#include <cstring>

namespace octlb {

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
