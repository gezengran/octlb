#pragma once

// platform.h must come before any OpenLB header that uses the any_platform macro.
// The macro is defined in platform.h but genericVector.h (via scalarVector.h)
// uses it before oalgorithm.h brings it in — so we force the correct order here.
#include "core/platform/platform.h"
#include "dynamics/lbm.h"
#include "descriptor/descriptor.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <optional>
#include <vector>

#include "src/common/types.h"
#include "src/solver/lbm/bc_kind.h"

namespace octlb {

struct CollideRhoStats {
  double sum = 0.0;
  int count = 0;

  void add(double rho) {
    sum += rho;
    ++count;
  }

  double average() const { return count > 0 ? sum / static_cast<double>(count) : 1.0; }
};

// OpenLB overlap geometry: ymin/ymax exterior padding slabs are stream-only
// (material 0: no collide, no mirror fill, no addPoints2CommBC partner copy).
// Full ix x iz extent on each exterior y slab (including iz face corners).
inline bool IsYminYmaxStreamOnlyPadding(int ix, int iy, int iz, int nx, int ny,
                                        int nz) {
  (void)ix;
  (void)iz;
  (void)nx;
  (void)nz;
  return iy < 0 || iy >= ny;
}

inline bool OverlapPaddingMaterialNonZero(int ix, int iy, int iz, int nx, int ny,
                                        int nz) {
  // OpenLB PostStream addPoints2CommBC: only geometry material!=0 padding cells.
  // cavity3d clean() leaves all exterior overlap at material 0 -> no PostStream
  // partner copy. Deeper padding (incl. ymin/ymax slabs) is stream-rotate only.
  return ix >= 0 && ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz;
}

/** Overlap padding collide policy (cavity3d C/D diagnostic). */
enum class OverlapPaddingCollideMode {
  /** Corner/edge padding: plain BGK (current OctLB default). */
  kBgkOnMaterialNonZero,
  /** OpenLB material 0: no collide on any overlap padding cell. */
  kNoDynamics,
};

/** Pull-stream link when upstream leaves halo on ymin/ymax exterior padding. */
enum class YminYmaxPaddingOutOfHaloMode {
  /** OpenLB mat-0: CyclicColumn.rotate wrap within padded block storage. */
  kOpenLbRotateWrap,
  kZero,      // diagnostic: exterior link -> 0
  kKeepSelf,  // diagnostic: keep pre-stream f_i at destination
};

// Euclidean modulo for OpenLB CyclicColumn::rotate pull indices.
inline int OpenLbPaddedBlockFlooredMod(int value, int modulus) {
  int wrapped = value % modulus;
  if (wrapped < 0) {
    wrapped += modulus;
  }
  return wrapped;
}

// OpenLB ConcreteBlockLattice::stream() rotates each iPop column by
// getNeighborDistance(c_i) modulo the **padded** block (core + 2*OVERLAP).
// OVERLAP defaults to 3 in SuperLattice mesh; InterpolatedVelocity FD radius is 2.
template <typename DESCRIPTOR>
inline int OpenLbPaddedBlockStreamSourceLinear(int dst_linear, int iPop, int NY,
                                               int NZ, int num_cells) {
  const int off = olb::descriptors::c<DESCRIPTOR>(iPop, 0) * NY * NZ +
                  olb::descriptors::c<DESCRIPTOR>(iPop, 1) * NZ +
                  olb::descriptors::c<DESCRIPTOR>(iPop, 2);
  return OpenLbPaddedBlockFlooredMod(dst_linear - off, num_cells);
}

// ymin/ymax mat-0 padding alias (audit helpers use this name).
template <typename DESCRIPTOR>
inline int YminYmaxMat0StreamSourceLinear(int dst_linear, int iPop, int NY, int NZ,
                                          int num_cells) {
  return OpenLbPaddedBlockStreamSourceLinear<DESCRIPTOR>(dst_linear, iPop, NY, NZ,
                                                         num_cells);
}

class BouzidiLinkData;

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

// OpenLB collision::ConstRhoBGK with explicit (rho, u) macroscopic fields.
// Boundary InterpolatedVelocity uses FixedVelocityMomentum (prescribed u).
template <typename T, typename DESCRIPTOR>
inline void CollideConstRhoBgkWithMacroscopic(CellProxy<T, DESCRIPTOR>& cell,
                                              T omega, T average_rho, T rho,
                                              const T u[DESCRIPTOR::d]) {
  T u_sqr = T{0};
  for (int d = 0; d < DESCRIPTOR::d; ++d) {
    u_sqr += u[d] * u[d];
  }
  const T delta_rho = T{1} - average_rho;
  const T ratio_rho =
      rho > T{1e-12} ? T{1} + delta_rho / rho : T{1};
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    const T f_eq =
        olb::equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, u_sqr);
    const T t_i = olb::descriptors::t<T, DESCRIPTOR>(iPop);
    cell[iPop] = ratio_rho * (f_eq + t_i) - t_i +
                   (T{1} - omega) * (cell[iPop] - f_eq);
  }
}

inline double ConstRhoStatisticRho(double rho, double average_rho) {
  return rho + (1.0 - average_rho);
}

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

  // ② edge-ghost: an edge is the intersection of two orthogonal faces. The edge
  // line runs along the third axis; buffer length is that axis's cell count.
  static int edge_buffer_count(int nx, int ny, int nz, FaceDir d1, FaceDir d2);

  // Set every interior cell to the Maxwell equilibrium for (rho0, u0).
  void initialize(T rho0, const T* u0);

  // BGK collision on every interior fluid cell.
  void collide(T omega);

  // OpenLB ConstRhoBGKdynamics: global average-rho correction (cavity3d bulk).
  void collide_const_rho(T omega, T average_rho, CollideRhoStats* rho_stats = nullptr);

  // Single fluid cell ConstRhoBGK (OpenLB Dominant collide spatial order).
  void collide_const_rho_at(int ix, int iy, int iz, T omega, T average_rho,
                            CollideRhoStats* rho_stats = nullptr);

  // Plain BGK on overlap padding only (after core collide passes).
  void collide_overlap_padding_bgk(T omega);

  T average_fluid_rho() const;

  // Pull-scheme streaming on core and overlap padding (OpenLB rotate semantics).
  // MPI ghost faces are still filled before stream() via GhostSchedule.
  void stream();

  void set_bouzidi_links(const BouzidiLinkData* links) { bouzidi_ = links; }
  const BouzidiLinkData* bouzidi_links() const { return bouzidi_; }

  void set_octant_id(OctantId id) { octant_id_ = id; }
  OctantId octant_id() const { return octant_id_; }

  // Physical origin (lower corner) + cell width of this block in domain
  // coordinates, set from the owning octant's quadrant_bounds. Used by
  // physical-coordinate inlet fields (e.g. a channel Poiseuille profile) so a
  // spatially varying BC can be evaluated from physical position rather than
  // block-local indices. Defaults keep index-based BCs (origin unset) correct.
  void set_phys_origin(std::array<double, 3> origin, double cell_width) {
    phys_origin_ = origin;
    phys_cell_width_ = cell_width;
  }
  const std::array<double, 3>& phys_origin() const { return phys_origin_; }
  double phys_cell_width() const { return phys_cell_width_; }

  void set_bc_kind(int ix, int iy, int iz, BcKind kind);
  BcKind bc_kind(int ix, int iy, int iz) const;

  // Halo coordinates: hx in [0, nx+2*h_), hy, hz likewise.
  T* populations_at_halo(int hx, int hy, int hz);
  const T* populations_at_halo(int hx, int hy, int hz) const;

  // T11 MEM-timing diagnostic: a snapshot of populations_ taken right before
  // stream() (i.e. post-collide + post-BC + post-ghost, the state the bounce-back
  // actually operates on). The standard momentum-exchange method must read the
  // OUTGOING (pre-stream) f_i on a boundary link, but after stream() that slot is
  // overwritten by the bounced-back value, so the live populations_ no longer hold
  // it. This snapshot lets a post-processor recompute the MEM force at the correct
  // (pre-stream) timing. Allocated lazily; empty until take_post_collide_snapshot().
  void take_post_collide_snapshot();
  bool has_post_collide_snapshot() const { return !post_collide_snapshot_.empty(); }
  const T* post_collide_populations_at_halo(int hx, int hy, int hz) const;

  // Fill ghost halo with periodic boundary values (for unit tests).
  // In production, GhostSchedule<BlockLattice> performs the equivalent MPI
  // exchange after collide() and before stream().
  void fill_periodic_halo();

  void pack_face(FaceDir dir, T* buffer, int count) const;
  void unpack_face(FaceDir dir, const T* buffer, int count);
  void read_ghost_face(FaceDir dir, T* buffer, int count) const;

  // ② edge-ghost line: pack reads the interior edge line (N cells along the
  // third axis at the d1,d2 corner); unpack writes the edge ghost line.
  void pack_edge(FaceDir d1, FaceDir d2, T* buffer, int count) const;
  void unpack_edge(FaceDir d1, FaceDir d2, const T* buffer, int count);

  // Cell proxy at physical coordinates (0-based, interior only).
  CellProxy<T, DESCRIPTOR> get(int ix, int iy, int iz);
  CellProxy<T, DESCRIPTOR> get(int ix, int iy, int iz) const;

  int nx() const { return nx_; }
  int ny() const { return ny_; }
  int nz() const { return nz_; }
  int halo_width() const { return h_; }

  void set_overlap_padding_collide_mode(OverlapPaddingCollideMode mode) {
    overlap_padding_collide_mode_ = mode;
  }
  OverlapPaddingCollideMode overlap_padding_collide_mode() const {
    return overlap_padding_collide_mode_;
  }

  void set_ymin_ymax_padding_out_of_halo_mode(
      YminYmaxPaddingOutOfHaloMode mode) {
    ymin_ymax_out_of_halo_mode_ = mode;
  }
  YminYmaxPaddingOutOfHaloMode ymin_ymax_padding_out_of_halo_mode() const {
    return ymin_ymax_out_of_halo_mode_;
  }

  // Mirror-copy core populations into exterior padding (fallback / tests).
  void fill_overlap_padding_from_core();

  // Pull stream into overlap padding shell (OpenLB rotate on padding only).
  void stream_overlap_padding_shell();

  // Commit overlap-padding stream computed in the last stream() call (PostStream).
  void commit_overlap_padding_stream();

  // OpenLB PostStream: addPoints2CommBC registration + POPULATION communicate.
  template <typename IsBoundaryIndicator, typename IsMaterialNonZero>
  void fill_overlap_padding_bc_post_stream(
      IsBoundaryIndicator&& is_bc_indicator,
      IsMaterialNonZero&& is_material_nonzero);

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
  std::vector<T> stream_tmp_;
  std::vector<T> post_collide_snapshot_;
  std::vector<BcKind> bc_kinds_;
  const BouzidiLinkData* bouzidi_ = nullptr;
  OctantId octant_id_ = 0;
  std::array<double, 3> phys_origin_{0.0, 0.0, 0.0};
  double phys_cell_width_ = 1.0;
  OverlapPaddingCollideMode overlap_padding_collide_mode_ =
      OverlapPaddingCollideMode::kNoDynamics;
  YminYmaxPaddingOutOfHaloMode ymin_ymax_out_of_halo_mode_ =
      YminYmaxPaddingOutOfHaloMode::kOpenLbRotateWrap;
};

template <typename T, typename DESCRIPTOR>
template <typename IsBoundaryIndicator, typename IsMaterialNonZero>
void BlockLattice<T, DESCRIPTOR>::fill_overlap_padding_bc_post_stream(
    IsBoundaryIndicator&& is_bc_indicator,
    IsMaterialNonZero&& is_material_nonzero) {
  const int overlap = h_;
  // InterpolatedVelocity::getNeighborhoodRadius() = 2 (addPoints2CommBC stencil).
  const int bc_overlap = std::min(overlap, 2);
  if (overlap == 0) {
    return;
  }

  const auto in_stencil = [&](int ix, int iy, int iz) {
    return ix >= -overlap && ix < nx_ + overlap && iy >= -overlap &&
           iy < ny_ + overlap && iz >= -overlap && iz < nz_ + overlap;
  };

  const auto in_core = [&](int ix, int iy, int iz) {
    return ix >= 0 && ix < nx_ && iy >= 0 && iy < ny_ && iz >= 0 && iz < nz_;
  };

  const auto copy_pops = [&](int ix_dst, int iy_dst, int iz_dst, int ix_src,
                             int iy_src, int iz_src) {
    if (ix_dst == ix_src && iy_dst == iy_src && iz_dst == iz_src) {
      return;
    }
    T* dst = populations_at_halo(ix_dst + overlap, iy_dst + overlap,
                                 iz_dst + overlap);
    const T* src = nullptr;
    if (in_core(ix_src, iy_src, iz_src)) {
      src = &get(ix_src, iy_src, iz_src)[0];
    } else if (in_stencil(ix_src, iy_src, iz_src)) {
      src = populations_at_halo(ix_src + overlap, iy_src + overlap,
                                iz_src + overlap);
    }
    if (src == nullptr) {
      return;
    }
    std::memcpy(dst, src, static_cast<std::size_t>(kQ) * sizeof(T));
  };

  // OpenLB addPoints2CommBC: only padding within neighborhood radius (2 for
  // InterpolatedVelocity), first boundary neighbor in (dx,dy,dz) loop order.
  for (int ix = -bc_overlap; ix < nx_ + bc_overlap; ++ix) {
    for (int iy = -bc_overlap; iy < ny_ + bc_overlap; ++iy) {
      for (int iz = -bc_overlap; iz < nz_ + bc_overlap; ++iz) {
        if (in_core(ix, iy, iz)) {
          continue;
        }
        if (!is_material_nonzero(ix, iy, iz)) {
          continue;
        }
        bool registered = false;
        int src_ix = 0;
        int src_iy = 0;
        int src_iz = 0;
        for (int dx = -bc_overlap; dx <= bc_overlap && !registered; ++dx) {
          for (int dy = -bc_overlap; dy <= bc_overlap && !registered; ++dy) {
            for (int dz = -bc_overlap; dz <= bc_overlap && !registered; ++dz) {
              const int nx_n = ix + dx;
              const int ny_n = iy + dy;
              const int nz_n = iz + dz;
              if (nx_n < -bc_overlap || nx_n >= nx_ + bc_overlap ||
                  ny_n < -bc_overlap || ny_n >= ny_ + bc_overlap ||
                  nz_n < -bc_overlap || nz_n >= nz_ + bc_overlap) {
                continue;
              }
              if (!is_bc_indicator(nx_n, ny_n, nz_n)) {
                continue;
              }
              registered = true;
              src_ix = nx_n;
              src_iy = ny_n;
              src_iz = nz_n;
            }
          }
        }
        if (!registered) {
          continue;
        }
        copy_pops(ix, iy, iz, src_ix, src_iy, src_iz);
      }
    }
  }
}

}  // namespace octlb
