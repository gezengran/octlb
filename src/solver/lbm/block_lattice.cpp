#include "block_lattice.h"

#include "src/solver/lbm/boundary/bouzidi_pull.h"
#include "src/solver/lbm/bouzidi_link_data.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <utility>
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

namespace {

template <typename T, typename DESCRIPTOR>
void CollideBgkAt(CellProxy<T, DESCRIPTOR>& cell, T omega) {
  T rho{}, u[DESCRIPTOR::d]{};
  cell.computeRhoU(rho, u);
  T u_sqr = T{0};
  for (int d = 0; d < DESCRIPTOR::d; ++d) {
    u_sqr += u[d] * u[d];
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    const T f_eq =
        olb::equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, u_sqr);
    cell[iPop] += omega * (f_eq - cell[iPop]);
  }
}

// OpenLB collision::ConstRhoBGK: ratioRho = 1 + (1 - avg_rho) / rho_cell.
// Rescales each cell's equilibrium toward unit density, suppressing the O(Ma²)
// compressibility errors that regular BGK accumulates at Ma ≈ 0.17.
// Only applied to physical fluid cells — overlap padding cells should not be
// corrected this way (OpenLB leaves them at NoDynamics / material-0).
template <typename T, typename DESCRIPTOR>
void CollideConstRhoBgkAt(CellProxy<T, DESCRIPTOR>& cell, T omega,
                          T average_rho) {
  T rho = T{0};
  T u[DESCRIPTOR::d]{};
  cell.computeRhoU(rho, u);
  CollideConstRhoBgkWithMacroscopic(cell, omega, average_rho, rho, u);
}

// Overlap padding cells are outside the physical domain.  OpenLB uses
// NoDynamics (material 0) for them, so they are never collided.  We keep
// a regular BGK for overlap cells so they stay numerically stable without
// applying the per-cell density rescaling that would disturb the halo.
template <typename T, typename DESCRIPTOR>
void CollideBgkOverlapAt(CellProxy<T, DESCRIPTOR>& cell, T omega) {
  T rho{}, u[DESCRIPTOR::d]{};
  cell.computeRhoU(rho, u);
  T u_sqr = T{0};
  for (int d = 0; d < DESCRIPTOR::d; ++d) u_sqr += u[d] * u[d];
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    const T f_eq =
        olb::equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, u_sqr);
    cell[iPop] += omega * (f_eq - cell[iPop]);
  }
}

}  // namespace

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::collide(T omega) {
  for (int ix = 0; ix < nx_; ++ix) {
    for (int iy = 0; iy < ny_; ++iy) {
      for (int iz = 0; iz < nz_; ++iz) {
        const CellKind kind = cell_kind(ix, iy, iz);
        if (kind == CellKind::kSolid || kind == CellKind::kBoundary) {
          continue;
        }
        auto cell = get(ix, iy, iz);
        CollideBgkAt<T, DESCRIPTOR>(cell, omega);
      }
    }
  }
  if (h_ == 0 ||
      overlap_padding_collide_mode_ == OverlapPaddingCollideMode::kNoDynamics) {
    return;
  }
  for (int hx = 0; hx < nx_ + 2 * h_; ++hx) {
    for (int hy = 0; hy < ny_ + 2 * h_; ++hy) {
      for (int hz = 0; hz < nz_ + 2 * h_; ++hz) {
        const int ix = hx - h_;
        const int iy = hy - h_;
        const int iz = hz - h_;
        const bool outside_core =
            ix < 0 || ix >= nx_ || iy < 0 || iy >= ny_ || iz < 0 || iz >= nz_;
        if (!outside_core) {
          continue;
        }
        if (ix < -h_ || ix >= nx_ + h_ || iy < -h_ || iy >= ny_ + h_ ||
            iz < -h_ || iz >= nz_ + h_) {
          continue;
        }
        if (!OverlapPaddingMaterialNonZero(ix, iy, iz, nx_, ny_, nz_)) {
          continue;
        }
        CellProxy<T, DESCRIPTOR> cell(populations_at_halo(hx, hy, hz));
        CollideBgkAt<T, DESCRIPTOR>(cell, omega);
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
T BlockLattice<T, DESCRIPTOR>::average_fluid_rho() const {
  T sum = T{0};
  int count = 0;
  for (int ix = 0; ix < nx_; ++ix) {
    for (int iy = 0; iy < ny_; ++iy) {
      for (int iz = 0; iz < nz_; ++iz) {
        if (cell_kind(ix, iy, iz) != CellKind::kFluid) {
          continue;
        }
        T rho = T{0};
        T u[DESCRIPTOR::d]{};
        get(ix, iy, iz).computeRhoU(rho, u);
        sum += rho;
        ++count;
      }
    }
  }
  return count > 0 ? sum / static_cast<T>(count) : T{1};
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::collide_const_rho_at(int ix, int iy, int iz,
                                                       T omega, T average_rho,
                                                       CollideRhoStats* rho_stats) {
  if (cell_kind(ix, iy, iz) != CellKind::kFluid) {
    return;
  }
  auto cell = get(ix, iy, iz);
  T rho = T{0};
  T u[DESCRIPTOR::d]{};
  cell.computeRhoU(rho, u);
  if (rho_stats != nullptr) {
    rho_stats->add(static_cast<double>(rho + (T{1} - average_rho)));
  }
  CollideConstRhoBgkAt<T, DESCRIPTOR>(cell, omega, average_rho);
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::collide_overlap_padding_bgk(T omega) {
  if (h_ == 0 ||
      overlap_padding_collide_mode_ == OverlapPaddingCollideMode::kNoDynamics) {
    return;
  }
  for (int hx = 0; hx < nx_ + 2 * h_; ++hx) {
    for (int hy = 0; hy < ny_ + 2 * h_; ++hy) {
      for (int hz = 0; hz < nz_ + 2 * h_; ++hz) {
        const int ix = hx - h_;
        const int iy = hy - h_;
        const int iz = hz - h_;
        const bool outside_core =
            ix < 0 || ix >= nx_ || iy < 0 || iy >= ny_ || iz < 0 || iz >= nz_;
        if (!outside_core) {
          continue;
        }
        if (ix < -h_ || ix >= nx_ + h_ || iy < -h_ || iy >= ny_ + h_ ||
            iz < -h_ || iz >= nz_ + h_) {
          continue;
        }
        if (!OverlapPaddingMaterialNonZero(ix, iy, iz, nx_, ny_, nz_)) {
          continue;
        }
        CellProxy<T, DESCRIPTOR> cell(populations_at_halo(hx, hy, hz));
        CollideBgkOverlapAt<T, DESCRIPTOR>(cell, omega);
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::collide_const_rho(T omega, T average_rho,
                                                    CollideRhoStats* rho_stats) {
  for (int ix = 0; ix < nx_; ++ix) {
    for (int iy = 0; iy < ny_; ++iy) {
      for (int iz = 0; iz < nz_; ++iz) {
        collide_const_rho_at(ix, iy, iz, omega, average_rho, rho_stats);
      }
    }
  }
  collide_overlap_padding_bgk(omega);
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::stream() {
  const int NX = nx_ + 2*h_;
  const int NY = ny_ + 2*h_;
  const int NZ = nz_ + 2*h_;

  stream_tmp_.resize(populations_.size());
  T* tmp = stream_tmp_.data();

  auto kind_at_halo = [&](int hx, int hy, int hz) -> CellKind {
    const int ix = hx - h_;
    const int iy = hy - h_;
    const int iz = hz - h_;
    if (ix < 0 || ix >= nx_ || iy < 0 || iy >= ny_ || iz < 0 || iz >= nz_) {
      return CellKind::kFluid;
    }
    return cell_kind(ix, iy, iz);
  };

  auto in_halo = [&](int hx, int hy, int hz) {
    return hx >= 0 && hx < NX && hy >= 0 && hy < NY && hz >= 0 && hz < NZ;
  };

  const auto in_overlap_padding_lattice = [&](int ix_int, int iy_int,
                                              int iz_int) {
    const bool outside_core =
        ix_int < 0 || ix_int >= nx_ || iy_int < 0 || iy_int >= ny_ ||
        iz_int < 0 || iz_int >= nz_;
    if (!outside_core) {
      return false;
    }
    return ix_int >= -h_ && ix_int < nx_ + h_ && iy_int >= -h_ &&
           iy_int < ny_ + h_ && iz_int >= -h_ && iz_int < nz_ + h_;
  };

  const auto should_stream_halo = [&](int hx, int hy, int hz) {
    const int ix_int = hx - h_;
    const int iy_int = hy - h_;
    const int iz_int = hz - h_;
    if (in_overlap_padding_lattice(ix_int, iy_int, iz_int)) {
      return true;
    }
    if (ix_int >= 0 && ix_int < nx_ && iy_int >= 0 && iy_int < ny_ &&
        iz_int >= 0 && iz_int < nz_) {
      return cell_kind(ix_int, iy_int, iz_int) != CellKind::kSolid;
    }
    return false;
  };

  const int num_cells = NX * NY * NZ;
  // OpenLB ConcreteBlockLattice::stream(): column rotate with wrap on the padded
  // block for every storage slot (incl. solid / NoDynamics). Pull+bounce is kept
  // only for Bouzidi cut links and ymin/ymax padding diagnostic modes.
  const bool use_openlb_rotate =
      bouzidi_ == nullptr &&
      ymin_ymax_out_of_halo_mode_ ==
          YminYmaxPaddingOutOfHaloMode::kOpenLbRotateWrap;

  const auto is_ymin_ymax_stream_only_padding =
      [&](int ix_int, int iy_int, int iz_int) {
        return IsYminYmaxStreamOnlyPadding(ix_int, iy_int, iz_int, nx_, ny_,
                                           nz_);
      };

  const auto pull_openlb_padded_rotate = [&](int dst_linear, int iPop) {
    const int src_linear = OpenLbPaddedBlockStreamSourceLinear<DESCRIPTOR>(
        dst_linear, iPop, NY, NZ, num_cells);
    return populations_[src_linear * kQ + iPop];
  };

  const auto stream_openlb_padded_block_rotate = [&](int hx, int hy, int hz) {
    const int dst_linear = hx * NY * NZ + hy * NZ + hz;
    const int dst = dst_linear * kQ;
    for (int iPop = 0; iPop < kQ; ++iPop) {
      tmp[dst + iPop] = pull_openlb_padded_rotate(dst_linear, iPop);
    }
  };

  const auto stream_one_cell = [&](int hx, int hy, int hz, CellKind kind,
                                   int ix_int, int iy_int, int iz_int) {
    if (kind == CellKind::kSolid) {
      return;
    }
    const int dst = (hx * NY * NZ + hy * NZ + hz) * kQ;
    const int dst_linear = hx * NY * NZ + hy * NZ + hz;
    const bool ymin_ymax_mat0 =
        is_ymin_ymax_stream_only_padding(ix_int, iy_int, iz_int);
    const bool overlap_padding =
        in_overlap_padding_lattice(ix_int, iy_int, iz_int);
    const bool openlb_padding_rotate =
        ymin_ymax_out_of_halo_mode_ ==
            YminYmaxPaddingOutOfHaloMode::kOpenLbRotateWrap &&
        (ymin_ymax_mat0 || overlap_padding);
    for (int iPop = 0; iPop < kQ; ++iPop) {
      if (openlb_padding_rotate) {
        const int cx = olb::descriptors::c<DESCRIPTOR>(iPop, 0);
        const int cy = olb::descriptors::c<DESCRIPTOR>(iPop, 1);
        const int cz = olb::descriptors::c<DESCRIPTOR>(iPop, 2);
        const int sx = hx - cx;
        const int sy = hy - cy;
        const int sz = hz - cz;
        // Mat-0 padding: in-halo neighbors use direct pull (same as rotate when
        // no wrap); exterior links use padded-block torus wrap (seeds f7/f15).
        if (in_halo(sx, sy, sz)) {
          const int src = (sx * NY * NZ + sy * NZ + sz) * kQ + iPop;
          tmp[dst + iPop] = populations_[src];
        } else {
          tmp[dst + iPop] = pull_openlb_padded_rotate(dst_linear, iPop);
        }
        continue;
      }

      const int sx = hx - olb::descriptors::c<DESCRIPTOR>(iPop, 0);
      const int sy = hy - olb::descriptors::c<DESCRIPTOR>(iPop, 1);
      const int sz = hz - olb::descriptors::c<DESCRIPTOR>(iPop, 2);
      if (!in_halo(sx, sy, sz)) {
        if (ymin_ymax_mat0) {
          if (ymin_ymax_out_of_halo_mode_ ==
              YminYmaxPaddingOutOfHaloMode::kKeepSelf) {
            tmp[dst + iPop] = populations_[dst + iPop];
          } else {
            tmp[dst + iPop] = T{0};
          }
        } else if (kind == CellKind::kBoundary) {
          const int opp = olb::descriptors::opposite<DESCRIPTOR>(iPop);
          tmp[dst + iPop] = populations_[dst + opp];
        } else {
          tmp[dst + iPop] = populations_[dst + iPop];
        }
        continue;
      }

      const int sx_int = sx - h_;
      const int sy_int = sy - h_;
      const int sz_int = sz - h_;
      const bool src_in_interior =
          sx_int >= 0 && sx_int < nx_ && sy_int >= 0 && sy_int < ny_ &&
          sz_int >= 0 && sz_int < nz_;
      if (kind == CellKind::kBoundary && !src_in_interior) {
        const int opp = olb::descriptors::opposite<DESCRIPTOR>(iPop);
        tmp[dst + iPop] = populations_[dst + opp];
        continue;
      }

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
          const int nx_h = hx + olb::descriptors::c<DESCRIPTOR>(iPop, 0);
          const int ny_h = hy + olb::descriptors::c<DESCRIPTOR>(iPop, 1);
          const int nz_h = hz + olb::descriptors::c<DESCRIPTOR>(iPop, 2);
          const int nbase = (nx_h * NY * NZ + ny_h * NZ + nz_h) * kQ + iPop;
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
  };

  if (use_openlb_rotate) {
    // OpenLB rotates every iPop column across the full padded block storage.
    for (int hx = 0; hx < NX; ++hx) {
      for (int hy = 0; hy < NY; ++hy) {
        for (int hz = 0; hz < NZ; ++hz) {
          stream_openlb_padded_block_rotate(hx, hy, hz);
        }
      }
    }
    std::memcpy(populations_.data(), tmp,
                populations_.size() * sizeof(T));
    return;
  }

  // Bouzidi / diagnostic stream modes: per-cell pull with mat-0 padding rotate.
  for (int hx = 0; hx < NX; ++hx) {
    for (int hy = 0; hy < NY; ++hy) {
      for (int hz = 0; hz < NZ; ++hz) {
        if (!should_stream_halo(hx, hy, hz)) {
          continue;
        }
        const int ix_int = hx - h_;
        const int iy_int = hy - h_;
        const int iz_int = hz - h_;
        const CellKind kind = in_overlap_padding_lattice(ix_int, iy_int, iz_int)
                                  ? CellKind::kFluid
                                  : cell_kind(ix_int, iy_int, iz_int);
        stream_one_cell(hx, hy, hz, kind, ix_int, iy_int, iz_int);
      }
    }
  }

  for (int hx = 0; hx < NX; ++hx) {
    for (int hy = 0; hy < NY; ++hy) {
      for (int hz = 0; hz < NZ; ++hz) {
        if (!should_stream_halo(hx, hy, hz)) {
          continue;
        }
        const int base = (hx * NY * NZ + hy * NZ + hz) * kQ;
        for (int iPop = 0; iPop < kQ; ++iPop) {
          populations_[base + iPop] = tmp[base + iPop];
        }
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::commit_overlap_padding_stream() {
  if (h_ == 0 || stream_tmp_.empty()) {
    return;
  }
  const int NX = nx_ + 2 * h_;
  const int NY = ny_ + 2 * h_;
  const int NZ = nz_ + 2 * h_;

  const auto in_overlap_padding_lattice = [&](int ix_int, int iy_int,
                                              int iz_int) {
    const bool outside_core =
        ix_int < 0 || ix_int >= nx_ || iy_int < 0 || iy_int >= ny_ ||
        iz_int < 0 || iz_int >= nz_;
    if (!outside_core) {
      return false;
    }
    return ix_int >= -h_ && ix_int < nx_ + h_ && iy_int >= -h_ &&
           iy_int < ny_ + h_ && iz_int >= -h_ && iz_int < nz_ + h_;
  };

  for (int hx = 0; hx < NX; ++hx) {
    for (int hy = 0; hy < NY; ++hy) {
      for (int hz = 0; hz < NZ; ++hz) {
        const int ix_int = hx - h_;
        const int iy_int = hy - h_;
        const int iz_int = hz - h_;
        if (!in_overlap_padding_lattice(ix_int, iy_int, iz_int)) {
          continue;
        }
        const int base = (hx * NY * NZ + hy * NZ + hz) * kQ;
        for (int iPop = 0; iPop < kQ; ++iPop) {
          populations_[base + iPop] = stream_tmp_[base + iPop];
        }
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::stream_overlap_padding_shell() {
  if (h_ == 0) {
    return;
  }
  const int NX = nx_ + 2 * h_;
  const int NY = ny_ + 2 * h_;
  const int NZ = nz_ + 2 * h_;

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

  auto in_halo = [&](int hx, int hy, int hz) {
    return hx >= 0 && hx < NX && hy >= 0 && hy < NY && hz >= 0 && hz < NZ;
  };

  const auto in_overlap_padding_lattice = [&](int ix_int, int iy_int,
                                              int iz_int) {
    const bool outside_core =
        ix_int < 0 || ix_int >= nx_ || iy_int < 0 || iy_int >= ny_ ||
        iz_int < 0 || iz_int >= nz_;
    if (!outside_core) {
      return false;
    }
    return ix_int >= -h_ && ix_int < nx_ + h_ && iy_int >= -h_ &&
           iy_int < ny_ + h_ && iz_int >= -h_ && iz_int < nz_ + h_;
  };

  for (int hx = 0; hx < NX; ++hx) {
    for (int hy = 0; hy < NY; ++hy) {
      for (int hz = 0; hz < NZ; ++hz) {
        const int ix_int = hx - h_;
        const int iy_int = hy - h_;
        const int iz_int = hz - h_;
        if (!in_overlap_padding_lattice(ix_int, iy_int, iz_int)) {
          continue;
        }
        const int dst = (hx * NY * NZ + hy * NZ + hz) * kQ;
        const int dst_linear = hx * NY * NZ + hy * NZ + hz;
        const int num_cells = NX * NY * NZ;
        const bool ymin_ymax_mat0 =
            IsYminYmaxStreamOnlyPadding(ix_int, iy_int, iz_int, nx_, ny_, nz_);
        const bool openlb_padding_rotate =
            ymin_ymax_out_of_halo_mode_ ==
                YminYmaxPaddingOutOfHaloMode::kOpenLbRotateWrap &&
            ymin_ymax_mat0;
        for (int iPop = 0; iPop < kQ; ++iPop) {
          if (openlb_padding_rotate) {
            const int cx = olb::descriptors::c<DESCRIPTOR>(iPop, 0);
            const int cy = olb::descriptors::c<DESCRIPTOR>(iPop, 1);
            const int cz = olb::descriptors::c<DESCRIPTOR>(iPop, 2);
            const int sx = hx - cx;
            const int sy = hy - cy;
            const int sz = hz - cz;
            if (in_halo(sx, sy, sz)) {
              const int src = (sx * NY * NZ + sy * NZ + sz) * kQ + iPop;
              tmp[dst + iPop] = populations_[src];
            } else {
              const int src_linear = OpenLbPaddedBlockStreamSourceLinear<DESCRIPTOR>(
                  dst_linear, iPop, NY, NZ, num_cells);
              tmp[dst + iPop] = populations_[src_linear * kQ + iPop];
            }
            continue;
          }

          const int sx = hx - olb::descriptors::c<DESCRIPTOR>(iPop, 0);
          const int sy = hy - olb::descriptors::c<DESCRIPTOR>(iPop, 1);
          const int sz = hz - olb::descriptors::c<DESCRIPTOR>(iPop, 2);
          if (!in_halo(sx, sy, sz)) {
            tmp[dst + iPop] = populations_[dst + iPop];
            continue;
          }
          const int src = (sx * NY * NZ + sy * NZ + sz) * kQ + iPop;
          (void)kind_at_halo;
          tmp[dst + iPop] = populations_[src];
        }
      }
    }
  }

  for (int hx = 0; hx < NX; ++hx) {
    for (int hy = 0; hy < NY; ++hy) {
      for (int hz = 0; hz < NZ; ++hz) {
        const int ix_int = hx - h_;
        const int iy_int = hy - h_;
        const int iz_int = hz - h_;
        if (!in_overlap_padding_lattice(ix_int, iy_int, iz_int)) {
          continue;
        }
        const int base = (hx * NY * NZ + hy * NZ + hz) * kQ;
        for (int iPop = 0; iPop < kQ; ++iPop) {
          populations_[base + iPop] = tmp[base + iPop];
        }
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void BlockLattice<T, DESCRIPTOR>::fill_overlap_padding_from_core() {
  if (h_ == 0) {
    return;
  }
  const int nx = nx_ + 2 * h_;
  const int ny = ny_ + 2 * h_;
  const int nz = nz_ + 2 * h_;

  const auto mirror_halo = [](int h, int n, int hi) {
    const int core_lo = h;
    const int core_hi = h + n - 1;
    if (hi < core_lo) {
      return 2 * core_lo - hi;
    }
    if (hi > core_hi) {
      return 2 * core_hi - hi;
    }
    return hi;
  };

  for (int hx = 0; hx < nx; ++hx) {
    for (int hy = 0; hy < ny; ++hy) {
      for (int hz = 0; hz < nz; ++hz) {
        const bool in_core = hx >= h_ && hx < h_ + nx_ && hy >= h_ && hy < h_ + ny_ &&
                             hz >= h_ && hz < h_ + nz_;
        if (in_core) {
          continue;
        }
        const int ix_int = hx - h_;
        const int iy_int = hy - h_;
        const int iz_int = hz - h_;
        // OpenLB material 0: stream-only padding (no mirror, no addPoints2CommBC).
        // Includes ymin/ymax slabs and xmin/xmax/zmin/zmax interior-face padding.
        if (!OverlapPaddingMaterialNonZero(ix_int, iy_int, iz_int, nx_, ny_,
                                           nz_)) {
          continue;
        }
        const int cx = mirror_halo(h_, nx_, hx);
        const int cy = mirror_halo(h_, ny_, hy);
        const int cz = mirror_halo(h_, nz_, hz);
        T* dst = populations_at_halo(hx, hy, hz);
        const T* src = populations_at_halo(cx, cy, cz);
        std::memcpy(dst, src, static_cast<std::size_t>(kQ) * sizeof(T));
      }
    }
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

template <typename T, typename DESCRIPTOR>
CellProxy<T, DESCRIPTOR>
BlockLattice<T, DESCRIPTOR>::get(int ix, int iy, int iz) const {
  return CellProxy<T, DESCRIPTOR>(
      const_cast<T*>(&populations_[idx(ix, iy, iz)]));
}

// ── Explicit instantiations ───────────────────────────────────────────────────
// All TUs that use BlockLattice must match one of these combinations.
template class BlockLattice<double, olb::descriptors::D3Q19<>>;
template class BlockLattice<float,  olb::descriptors::D3Q19<>>;

}  // namespace octlb
