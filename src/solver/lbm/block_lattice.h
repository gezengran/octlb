#pragma once

// platform.h must come before any OpenLB header that uses the any_platform macro.
// The macro is defined in platform.h but genericVector.h (via scalarVector.h)
// uses it before oalgorithm.h brings it in — so we force the correct order here.
#include "core/platform/platform.h"
#include "dynamics/lbm.h"
#include "descriptor/descriptor.h"

#include <vector>

#include "src/common/types.h"

namespace octlb {

// Proxy for a single lattice cell; satisfies olb::concepts::MinimalCell.
// Points directly into BlockLattice::populations_, no allocation.
template <typename T, typename DESCRIPTOR>
struct CellProxy {
  using value_t      = T;
  using descriptor_t = DESCRIPTOR;

  explicit CellProxy(T* data) : data_(data) {}

  T&       operator[](unsigned iPop)       { return data_[iPop]; }
  const T& operator[](unsigned iPop) const { return data_[iPop]; }

  void computeRhoU(T& rho, T* u) const {
    // OpenLB stores f_i - t_i in cell[iPop]; rho = sum(f_i) = sum(cell) + 1.
    // Momentum J = sum(f_i * c_i) = sum(cell[iPop] * c_i) because
    // sum(t_i * c_i) == 0 for symmetric lattices.
    rho = T{1};
    for (int d = 0; d < DESCRIPTOR::d; ++d) u[d] = T{0};
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      rho += data_[iPop];
      for (int d = 0; d < DESCRIPTOR::d; ++d)
        u[d] += data_[iPop] * olb::descriptors::c<DESCRIPTOR>(iPop, d);
    }
    for (int d = 0; d < DESCRIPTOR::d; ++d) u[d] /= rho;
  }

 private:
  T* data_;
};

// ── BlockLattice ──────────────────────────────────────────────────────────────
//
// One LBM block per p4est octant.  Memory layout:
//   populations_[flat_idx(ix,iy,iz) * Q + iPop]
// where flat_idx includes the ghost halo of width h_ on all sides, so the
// physical domain occupies indices [h_, Nx+h_) × [h_, Ny+h_) × [h_, Nz+h_).
//
// stream() reads from the ghost halo; callers must fill it before streaming:
//   - In production: GhostSchedule<BlockLattice> fills ghost cells via MPI.
//   - In unit tests:  fill_periodic_halo() copies boundary rows into ghost.
template <typename T, typename DESCRIPTOR>
class BlockLattice {
 public:
  using face_value_t = T;

  static constexpr int kQ = DESCRIPTOR::q;

  BlockLattice(int nx, int ny, int nz, int halo = 1);

  static int face_buffer_count(int nx, int ny, int nz, FaceDir dir);

  // Set every interior cell to the Maxwell equilibrium for (rho0, u0).
  void initialize(T rho0, const T* u0);

  // BGK collision on every interior cell.
  void collide(T omega);

  // Pull-scheme streaming.  Reads from the ghost halo — caller must fill it
  // with post-collision values BEFORE calling stream().
  //
  // Standard per-step order:
  //   1. collide()              — interior ← f*
  //   2. fill_periodic_halo()   — ghost ← f* (or GhostSchedule MPI exchange)
  //   3. stream()               — interior ← f*[x - c]
  void stream();

  // Fill ghost halo with periodic boundary values (for unit tests).
  // In production, GhostSchedule<BlockLattice> performs the equivalent MPI
  // exchange after collide() and before stream().
  void fill_periodic_halo();

  void pack_face(FaceDir dir, T* buffer, int count) const;
  void unpack_face(FaceDir dir, const T* buffer, int count);
  void read_ghost_face(FaceDir dir, T* buffer, int count) const;

  // Cell proxy at physical coordinates (0-based, interior only).
  CellProxy<T, DESCRIPTOR> get(int ix, int iy, int iz);

  int nx() const { return nx_; }
  int ny() const { return ny_; }
  int nz() const { return nz_; }

 private:
  int halo_idx(int hx, int hy, int hz) const {
    return (hx * (ny_ + 2 * h_) * (nz_ + 2 * h_) + hy * (nz_ + 2 * h_) + hz) *
           kQ;
  }

  // Index into populations_ (includes halo offset).
  int idx(int ix, int iy, int iz) const {
    return ((ix + h_) * (ny_ + 2*h_) * (nz_ + 2*h_)
          + (iy + h_) * (nz_ + 2*h_)
          + (iz + h_)) * kQ;
  }

  int nx_, ny_, nz_, h_;
  std::vector<T> populations_;
};

}  // namespace octlb
