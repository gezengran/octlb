#ifndef OCTLB_EXAMPLES_CAVITY3D_CAVITY3D_CASE_H_
#define OCTLB_EXAMPLES_CAVITY3D_CAVITY3D_CASE_H_

#include <cstdio>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <numeric>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

#include "core/platform/platform.h"
#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/ghost_schedule.h"
#include "src/solver/io/vtk_writer/amr_vtk_writer.h"
#include "src/solver/io/vtk_writer/vtk_cell_field.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/interpolated_velocity.h"
#include "src/solver/lbm/domain_boundary_handler.h"
#include "src/solver/lbm/level_coupler.h"
#include "src/solver/lbm/time_loop/time_loop.h"
#include "src/solver/lbm/unit_converter/unit_converter.h"
#include "src/solver/lbm/vtk_lbm_fields.h"

namespace octlb {

/** OpenLB cavity3d ValueTracer convergence step (mpirun -n 1 ./cavity3d). */
inline constexpr int kOpenLbCavity3dConvergedSteps = 5269;

/** OpenLB InterpolatedVelocity::getNeighborhoodRadius() (FD stencil width). */
inline constexpr int kCavityInterpolatedVelocityOverlap =
    boundary::kInterpolatedVelocityOverlap;

/** OpenLB SuperLattice block padding (parameters::OVERLAP default = 3).
 *  Stream column-rotate wrap uses the full padded block (core + 2*padding). */
inline constexpr int kOpenLbCavityBlockPadding = 3;


using CavityLattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

inline BoundingBox CavityDomain(const UnitConverter& converter) {
  const scalar L = converter.char_phys_length();
  return {0.0, 0.0, 0.0, L, L, L};
}

/** OpenLB lattice nodes per axis: RESOLUTION=30 -> 31 nodes (0..L). */
inline int CavityLatticeNodesPerAxis(const UnitConverter& converter) {
  return converter.resolution() + 1;
}

/** Fluid material=1 cells per axis inside the boundary shell. */
inline int CavityFluidCellsPerAxis(const UnitConverter& converter) {
  return converter.resolution() - 1;
}

// OpenLB cavity3d setInitialValues: bulk u=0 equilibrium; moving-lid face
// (material 3) gets u_lid in the VELOCITY field only (pops stay u=0 eq until
// collide + PostStream BC). ymin/ymax mat-0 padding (e.g. iy=-2) in-plane
// f7/f15 are seeded only by padded-block stream rotate-wrap over many steps;
// wrong lid pop IC shifts the ymax torus entry and leaves the diagonal ring
// near zero instead of O(1e-3).
inline void InitializeCavityLattice(CavityLattice& lat, int n_lat,
                                    const UnitConverter& converter) {
  (void)converter;
  boundary::MarkDomainBoundaryCellKinds(lat, n_lat, n_lat, n_lat);
  const double u_zero[3] = {0.0, 0.0, 0.0};
  lat.initialize(1.0, u_zero);
}

inline OctreeForest MakePartitionedForest(MPI_Comm comm,
                                         const BoundingBox& domain) {
  OctreeForest forest(comm, domain);
  forest.partition();
  return forest;
}

inline std::vector<DomainBcSpec> CavityBoundarySpecs(
    const UnitConverter& converter) {
  const double u_lid = converter.char_lattice_velocity();
  std::vector<DomainBcSpec> specs;
  specs.push_back({FaceDir::kXMin, DomainBcType::kInterpolatedVelocity,
                   {0.0, 0.0, 0.0}});
  specs.push_back({FaceDir::kXMax, DomainBcType::kInterpolatedVelocity,
                   {0.0, 0.0, 0.0}});
  specs.push_back({FaceDir::kYMin, DomainBcType::kInterpolatedVelocity,
                   {0.0, 0.0, 0.0}});
  specs.push_back({FaceDir::kZMin, DomainBcType::kInterpolatedVelocity,
                   {0.0, 0.0, 0.0}});
  specs.push_back({FaceDir::kZMax, DomainBcType::kInterpolatedVelocity,
                   {0.0, 0.0, 0.0}});
  specs.push_back({FaceDir::kYMax, DomainBcType::kInterpolatedVelocity,
                   {u_lid, 0.0, 0.0}});
  return specs;
}

/** Ghia et al. 1982 centerline u/u_lid at Re=100 (OpenLB cavity2d table). */
inline constexpr std::array<double, 17> kGhiaRe100_YOverH = {
    1.0,    0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172,
    0.5,    0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547,
    0.0};

/** OpenLB cavity3d defaults, converged iT=5269 (mpirun -n 1 ./cavity3d).
 *  Verified against fresh GhiaDump output in
 *  tests/data/openlb_cavity3d_centerline_dump.txt (2026-06-08). */
inline constexpr std::array<double, 17> kOpenLbCavity3d_UOverLid = {
    1.0,     0.629597, 0.506129, 0.429369, 0.371218, 0.110531, 0.0733846,
    0.0481395, 0.00778279, -0.0115753, -0.152547, -0.213294, -0.167816,
    -0.124812, -0.112694, -0.0995179, 0.0};

inline constexpr std::array<double, 17> kGhiaRe100_UOverLid = {
    1.0,     0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332,
    -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775,
    -0.04192, -0.03717, 0.0};

inline double RelativeL2(const std::vector<double>& simulated,
                         const std::vector<double>& reference) {
  double sum_sq_diff = 0.0;
  double sum_sq_ref = 0.0;
  for (std::size_t i = 0; i < simulated.size(); ++i) {
    const double diff = simulated[i] - reference[i];
    sum_sq_diff += diff * diff;
    sum_sq_ref += reference[i] * reference[i];
  }
  return std::sqrt(sum_sq_diff / sum_sq_ref);
}

struct Cavity3dCase {
  UnitConverter converter;
  BoundingBox domain;
  OctreeForest forest;
  FacePairList face_pairs;
  int n_lat;
  int n_fluid;
  BlockCollection<CavityLattice> blocks;
  GhostSchedule<CavityLattice> ghosts;
  LevelCoupler coupler;
  ConcreteDomainBoundaryHandler domain_bc;
  TimeLoop loop;

  Cavity3dCase(MPI_Comm comm, UnitConverter converter_in,
               OverlapPaddingMode padding_mode = OverlapPaddingMode::kHybrid,
               ConstRhoStatsScope const_rho_stats_scope =
                   ConstRhoStatsScope::kFluidAndBoundary,
               OverlapPaddingCollideMode padding_collide_mode =
                   OverlapPaddingCollideMode::kNoDynamics,
               bool use_const_rho_bgk = true,
               YminYmaxPaddingOutOfHaloMode ymin_ymax_out_of_halo_mode =
                   YminYmaxPaddingOutOfHaloMode::kOpenLbRotateWrap,
               int halo_overlap = kOpenLbCavityBlockPadding)
      : converter(std::move(converter_in)),
        domain(CavityDomain(converter)),
        forest(MakePartitionedForest(comm, domain)),
        face_pairs(forest),
        n_lat(CavityLatticeNodesPerAxis(converter)),
        n_fluid(CavityFluidCellsPerAxis(converter)),
        blocks(forest.local_num_octants(),
               [n = n_lat, conv = converter,
                pcm = padding_collide_mode,
                oo_mode = ymin_ymax_out_of_halo_mode,
                halo = halo_overlap](OctantId) {
                 CavityLattice lat(n, n, n, halo);
                 lat.set_overlap_padding_collide_mode(pcm);
                 lat.set_ymin_ymax_padding_out_of_halo_mode(oo_mode);
                 InitializeCavityLattice(lat, n, conv);
                 return lat;
               }),
        ghosts(comm, face_pairs, blocks, n_lat, n_lat, n_lat),
        coupler(comm, face_pairs, forest, blocks, n_lat, n_lat, n_lat,
                converter.omega()),
        domain_bc(blocks, face_pairs.tree_boundary_faces(),
                  CavityBoundarySpecs(converter), n_lat, n_lat, n_lat,
                  converter.omega(), /*boundary_lattice_mode=*/true,
                  padding_mode),
        loop(forest, blocks, ghosts, coupler, domain_bc, converter.omega(),
             use_const_rho_bgk, const_rho_stats_scope) {}

  void advance_steps(int num_steps) {
    for (int step = 0; step < num_steps; ++step) {
      loop.advance_one();
    }
  }

  CavityLattice& lattice() { return blocks[0]; }
  const CavityLattice& lattice() const { return blocks[0]; }

  double lattice_ux_at(int ix, int iy, int iz) const {
    if (lattice().cell_kind(ix, iy, iz) == CellKind::kBoundary) {
      double u_wall[3] = {};
      boundary::detail::PrescribedBoundaryU(ix, iy, iz, n_lat, n_lat, n_lat,
                                            CavityBoundarySpecs(converter),
                                            u_wall);
      return u_wall[0];
    }
    double rho = 0.0;
    double u[3] = {};
    lattice().get(ix, iy, iz).computeRhoU(rho, u);
    return u[0];
  }

  // Velocity from populations only (OpenLB bulk computeRhoU on all cells).
  double lattice_ux_raw_at(int ix, int iy, int iz) const {
    double rho = 0.0;
    double u[3] = {};
    lattice().get(ix, iy, iz).computeRhoU(rho, u);
    return u[0];
  }

  static const char* CellKindLabel(CellKind kind) {
    switch (kind) {
      case CellKind::kFluid:
        return "fluid";
      case CellKind::kBoundary:
        return "boundary";
      case CellKind::kSolid:
        return "solid";
    }
    return "unknown";
  }

  // Trilinear interpolation aligned with OpenLB AnalyticalFfromSuperF3D:
  // floor(phys/dx) base node, weights from phys - origin - latticeR*dx.
  double sample_ux(scalar x, scalar y, scalar z) const {
    const scalar dx = converter.phys_delta_x();

    auto floor_lattice = [this, dx](scalar coord) {
      int idx = static_cast<int>(std::floor(coord / dx));
      if (idx < 0) {
        idx = 0;
      }
      if (idx > n_lat - 2) {
        idx = n_lat - 2;
      }
      return idx;
    };

    const int ix0 = floor_lattice(x);
    const int iy0 = floor_lattice(y);
    const int iz0 = floor_lattice(z);
    const int ix1 = ix0 + 1;
    const int iy1 = iy0 + 1;
    const int iz1 = iz0 + 1;

    const double wx = (x - ix0 * dx) / dx;
    const double wy = (y - iy0 * dx) / dx;
    const double wz = (z - iz0 * dx) / dx;

    auto ux_at = [this](int ix, int iy, int iz) {
      return lattice_ux_at(ix, iy, iz);
    };

    const double c000 = ux_at(ix0, iy0, iz0);
    const double c100 = ux_at(ix1, iy0, iz0);
    const double c010 = ux_at(ix0, iy1, iz0);
    const double c110 = ux_at(ix1, iy1, iz0);
    const double c001 = ux_at(ix0, iy0, iz1);
    const double c101 = ux_at(ix1, iy0, iz1);
    const double c011 = ux_at(ix0, iy1, iz1);
    const double c111 = ux_at(ix1, iy1, iz1);

    const double c00 = c000 * (1.0 - wx) + c100 * wx;
    const double c10 = c010 * (1.0 - wx) + c110 * wx;
    const double c01 = c001 * (1.0 - wx) + c101 * wx;
    const double c11 = c011 * (1.0 - wx) + c111 * wx;
    const double c0 = c00 * (1.0 - wy) + c10 * wy;
    const double c1 = c01 * (1.0 - wy) + c11 * wy;
    return c0 * (1.0 - wz) + c1 * wz;
  }

  int centerline_ix() const {
    const scalar dx = converter.phys_delta_x();
    const scalar x = 0.5 * converter.char_phys_length();
    return static_cast<int>(std::floor(x / dx));
  }

  int centerline_iz() const {
    const scalar dx = converter.phys_delta_x();
    const scalar z = 0.5 * converter.char_phys_length();
    return static_cast<int>(std::floor(z / dx));
  }

  int floor_lattice_index(scalar coord) const {
    const scalar dx = converter.phys_delta_x();
    int idx = static_cast<int>(std::floor(coord / dx));
    if (idx < 0) {
      idx = 0;
    }
    if (idx > n_lat - 2) {
      idx = n_lat - 2;
    }
    return idx;
  }

  int nearest_lattice_index(scalar coord) const {
    const scalar dx = converter.phys_delta_x();
    int idx = static_cast<int>(std::lround(coord / dx));
    if (idx < 0) {
      idx = 0;
    }
    if (idx >= n_lat) {
      idx = n_lat - 1;
    }
    return idx;
  }

  std::vector<double> sample_ghia_centerline() const {
    const scalar L = converter.char_phys_length();
    const scalar x = 0.5 * L;
    const scalar z = 0.5 * L;
    const double u_lid = converter.char_lattice_velocity();

    std::vector<double> profile;
    profile.reserve(kGhiaRe100_YOverH.size());
    for (double y_over_h : kGhiaRe100_YOverH) {
      const double ux = sample_ux(x, y_over_h * L, z);
      profile.push_back(ux / u_lid);
    }
    return profile;
  }

  bool has_non_finite_velocity() const {
    for (int ix = 1; ix < n_lat - 1; ++ix) {
      for (int iy = 1; iy < n_lat - 1; ++iy) {
        for (int iz = 1; iz < n_lat - 1; ++iz) {
          double rho = 0.0;
          double u[3] = {};
          lattice().get(ix, iy, iz).computeRhoU(rho, u);
          for (int d = 0; d < 3; ++d) {
            if (!std::isfinite(u[d])) {
              return true;
            }
          }
        }
      }
    }
    return false;
  }

  bool any_ux_above(double threshold) const {
    for (int ix = 1; ix < n_lat - 1; ++ix) {
      for (int iy = 1; iy < n_lat - 1; ++iy) {
        for (int iz = 1; iz < n_lat - 1; ++iz) {
          double rho = 0.0;
          double u[3] = {};
          lattice().get(ix, iy, iz).computeRhoU(rho, u);
          if (u[0] > threshold) {
            return true;
          }
        }
      }
    }
    return false;
  }

  void write_vtk_timestep(MPI_Comm comm, int iT, const std::string& output_dir,
                          const std::string& base_name = "cavity3d") const {
    AmrVtkWriter writer(comm, forest, n_lat, n_lat, n_lat, output_dir, base_name);
    VtkVelocityField<double, olb::descriptors::D3Q19<>> velocity(lattice());
    VtkPhysVelocityField<double, olb::descriptors::D3Q19<>> phys_velocity(
        lattice(), converter);
    const VtkCellFieldView velocity_field = VtkCellFieldView::From(velocity);
    const VtkCellFieldView phys_field = VtkCellFieldView::From(phys_velocity);
    const std::array<VtkCellFieldView, 2> fields = {velocity_field, phys_field};
    const std::vector<std::string> paths = writer.WriteTimestep(
        iT, blocks,
        std::span<const VtkCellFieldView>(fields.data(), fields.size()));
    writer.WriteVtmAndPvd(iT, paths);
  }

  void write_centerline_csv(const std::string& path) const {
    const std::vector<double> simulated = sample_ghia_centerline();
    std::ofstream out(path);
    out << "y_over_H,u_over_lid,ghia\n";
    for (std::size_t i = 0; i < kGhiaRe100_YOverH.size(); ++i) {
      out << kGhiaRe100_YOverH[i] << ',' << simulated[i] << ','
          << kGhiaRe100_UOverLid[i] << '\n';
    }
  }
};

struct CavityCenterlinePointRow {
  double y_over_h = 0.0;
  int ix = 0;
  int iy_floor = 0;
  int iy_nearest = 0;
  int iz = 0;
  const char* kind_floor = "unknown";
  const char* kind_nearest = "unknown";
  double trilinear = 0.0;
  double ux_at_floor = 0.0;
  double raw_at_floor = 0.0;
  double ux_at_nearest = 0.0;
  double raw_at_nearest = 0.0;
  double openlb = 0.0;
};

struct CavityCenterlinePointComparison {
  int ix_center = 0;
  int iz_center = 0;
  std::vector<CavityCenterlinePointRow> rows;
  double l2_trilinear = 0.0;
  double l2_ux_at_floor = 0.0;
  double l2_raw_at_floor = 0.0;
  double l2_ux_at_nearest = 0.0;
  double l2_raw_at_nearest = 0.0;
  double l2_trilinear_minus_ux_floor = 0.0;
};

inline CavityCenterlinePointComparison AnalyzeCenterlinePoints(
    const Cavity3dCase& cavity) {
  CavityCenterlinePointComparison cmp;
  cmp.ix_center = cavity.centerline_ix();
  cmp.iz_center = cavity.centerline_iz();
  const double u_lid = cavity.converter.char_lattice_velocity();
  const scalar L = cavity.converter.char_phys_length();
  const scalar x = 0.5 * L;
  const scalar z = 0.5 * L;

  std::vector<double> trilinear;
  std::vector<double> ux_at_floor;
  std::vector<double> raw_at_floor;
  std::vector<double> ux_at_nearest;
  std::vector<double> raw_at_nearest;
  std::vector<double> openlb;

  cmp.rows.reserve(kGhiaRe100_YOverH.size());
  trilinear.reserve(kGhiaRe100_YOverH.size());
  ux_at_floor.reserve(kGhiaRe100_YOverH.size());
  raw_at_floor.reserve(kGhiaRe100_YOverH.size());
  ux_at_nearest.reserve(kGhiaRe100_YOverH.size());
  raw_at_nearest.reserve(kGhiaRe100_YOverH.size());
  openlb.reserve(kGhiaRe100_YOverH.size());

  for (std::size_t i = 0; i < kGhiaRe100_YOverH.size(); ++i) {
    const double y_over_h = kGhiaRe100_YOverH[i];
    const scalar y = y_over_h * L;
    const int iy_floor = cavity.floor_lattice_index(y);
    const int iy_nearest = cavity.nearest_lattice_index(y);

    CavityCenterlinePointRow row;
    row.y_over_h = y_over_h;
    row.ix = cmp.ix_center;
    row.iy_floor = iy_floor;
    row.iy_nearest = iy_nearest;
    row.iz = cmp.iz_center;
    row.kind_floor =
        Cavity3dCase::CellKindLabel(cavity.lattice().cell_kind(
            row.ix, iy_floor, row.iz));
    row.kind_nearest =
        Cavity3dCase::CellKindLabel(cavity.lattice().cell_kind(
            row.ix, iy_nearest, row.iz));

    row.trilinear = cavity.sample_ux(x, y, z) / u_lid;
    row.ux_at_floor =
        cavity.lattice_ux_at(row.ix, iy_floor, row.iz) / u_lid;
    row.raw_at_floor =
        cavity.lattice_ux_raw_at(row.ix, iy_floor, row.iz) / u_lid;
    row.ux_at_nearest =
        cavity.lattice_ux_at(row.ix, iy_nearest, row.iz) / u_lid;
    row.raw_at_nearest =
        cavity.lattice_ux_raw_at(row.ix, iy_nearest, row.iz) / u_lid;
    row.openlb = kOpenLbCavity3d_UOverLid[i];

    trilinear.push_back(row.trilinear);
    ux_at_floor.push_back(row.ux_at_floor);
    raw_at_floor.push_back(row.raw_at_floor);
    ux_at_nearest.push_back(row.ux_at_nearest);
    raw_at_nearest.push_back(row.raw_at_nearest);
    openlb.push_back(row.openlb);
    cmp.rows.push_back(row);
  }

  cmp.l2_trilinear = RelativeL2(trilinear, openlb);
  cmp.l2_ux_at_floor = RelativeL2(ux_at_floor, openlb);
  cmp.l2_raw_at_floor = RelativeL2(raw_at_floor, openlb);
  cmp.l2_ux_at_nearest = RelativeL2(ux_at_nearest, openlb);
  cmp.l2_raw_at_nearest = RelativeL2(raw_at_nearest, openlb);
  cmp.l2_trilinear_minus_ux_floor = RelativeL2(trilinear, ux_at_floor);
  return cmp;
}

inline void PrintCenterlinePointComparison(
    const UnitConverter& converter, int steps,
    const CavityCenterlinePointComparison& cmp,
    std::ostream& out = std::cout) {
  out << "OctLB cavity3d centerline point comparison (ix=" << cmp.ix_center
      << " iz=" << cmp.iz_center << ")\n";
  out << "  N=" << converter.resolution()
      << " tau=" << converter.lattice_relaxation_time()
      << " omega=" << converter.omega() << " iT=" << steps << "\n";
  out << "  L2 vs OpenLB: trilinear=" << cmp.l2_trilinear
      << " ux@floor=" << cmp.l2_ux_at_floor
      << " raw@floor=" << cmp.l2_raw_at_floor
      << " ux@nearest=" << cmp.l2_ux_at_nearest
      << " raw@nearest=" << cmp.l2_raw_at_nearest << "\n";
  out << "  L2(trilinear - ux@floor)=" << cmp.l2_trilinear_minus_ux_floor
      << "  (sampling-only artifact)\n";
  out << "  columns: y/H iyF kindF iyN kindN | tri uxF rawF uxN rawN | "
         "openlb | err_tri err_uxF err_rawF\n";
  for (const CavityCenterlinePointRow& row : cmp.rows) {
    out << "  y/H=" << row.y_over_h << " iyF=" << row.iy_floor << ' '
        << row.kind_floor << " iyN=" << row.iy_nearest << ' '
        << row.kind_nearest << " | tri=" << row.trilinear
        << " uxF=" << row.ux_at_floor << " rawF=" << row.raw_at_floor
        << " uxN=" << row.ux_at_nearest << " rawN=" << row.raw_at_nearest
        << " | olb=" << row.openlb << " | d_tri=" << (row.trilinear - row.openlb)
        << " d_uxF=" << (row.ux_at_floor - row.openlb)
        << " d_rawF=" << (row.raw_at_floor - row.openlb) << '\n';
  }
}

// Piecewise-linear OpenLB centerline reference between Ghia sample heights.
// kGhiaRe100_YOverH is descending (1.0 -> 0.0).
inline double InterpolateOpenLbCenterline(double y_over_h) {
  if (y_over_h >= kGhiaRe100_YOverH.front()) {
    return kOpenLbCavity3d_UOverLid.front();
  }
  if (y_over_h <= kGhiaRe100_YOverH.back()) {
    return kOpenLbCavity3d_UOverLid.back();
  }
  for (std::size_t i = 0; i + 1 < kGhiaRe100_YOverH.size(); ++i) {
    const double y_hi = kGhiaRe100_YOverH[i];
    const double y_lo = kGhiaRe100_YOverH[i + 1];
    if (y_over_h <= y_hi && y_over_h >= y_lo) {
      const double u_hi = kOpenLbCavity3d_UOverLid[i];
      const double u_lo = kOpenLbCavity3d_UOverLid[i + 1];
      const double t = (y_over_h - y_hi) / (y_lo - y_hi);
      return u_hi + t * (u_lo - u_hi);
    }
  }
  return 0.0;
}

inline bool IsGhiaCenterlineHeight(double y_over_h) {
  for (double y : kGhiaRe100_YOverH) {
    if (std::abs(y - y_over_h) < 1e-6) {
      return true;
    }
  }
  return false;
}

struct CenterlineColumnRow {
  int iy = 0;
  double y_over_h = 0.0;
  const char* kind = "unknown";
  double ux_over_lid = 0.0;
  double rho = 0.0;
  double openlb_interp = 0.0;
  double err = 0.0;
  bool ghia_sample = false;
};

struct CenterlineColumn {
  int ix = 0;
  int iz = 0;
  std::vector<CenterlineColumnRow> rows;
  double l2_node_vs_openlb_interp = 0.0;
  double avg_rho_fluid = 0.0;
  double avg_rho_all = 0.0;
};

inline CenterlineColumn AnalyzeCenterlineColumn(const Cavity3dCase& cavity) {
  CenterlineColumn col;
  col.ix = cavity.centerline_ix();
  col.iz = cavity.centerline_iz();
  const double u_lid = cavity.converter.char_lattice_velocity();
  const scalar L = cavity.converter.char_phys_length();
  const scalar dx = cavity.converter.phys_delta_x();

  double sum_rho_fluid = 0.0;
  int count_fluid = 0;
  double sum_rho_all = 0.0;
  int count_all = 0;

  std::vector<double> sim;
  std::vector<double> ref;
  col.rows.reserve(static_cast<std::size_t>(cavity.n_lat));

  for (int iy = 0; iy < cavity.n_lat; ++iy) {
    const double y_over_h = (static_cast<double>(iy) * dx) / L;
    CenterlineColumnRow row;
    row.iy = iy;
    row.y_over_h = y_over_h;
    row.kind = Cavity3dCase::CellKindLabel(
        cavity.lattice().cell_kind(col.ix, iy, col.iz));

    double u[3] = {};
    row.rho = 0.0;
    cavity.lattice().get(col.ix, iy, col.iz).computeRhoU(row.rho, u);
    row.ux_over_lid = cavity.lattice_ux_at(col.ix, iy, col.iz) / u_lid;
    row.openlb_interp = InterpolateOpenLbCenterline(y_over_h);
    row.err = row.ux_over_lid - row.openlb_interp;
    row.ghia_sample = IsGhiaCenterlineHeight(y_over_h);

    sum_rho_all += row.rho;
    ++count_all;
    if (cavity.lattice().cell_kind(col.ix, iy, col.iz) == CellKind::kFluid) {
      sum_rho_fluid += row.rho;
      ++count_fluid;
    }

    if (cavity.lattice().cell_kind(col.ix, iy, col.iz) == CellKind::kFluid) {
      sim.push_back(row.ux_over_lid);
      ref.push_back(row.openlb_interp);
    }
    col.rows.push_back(row);
  }

  col.avg_rho_fluid =
      count_fluid > 0 ? sum_rho_fluid / static_cast<double>(count_fluid) : 1.0;
  col.avg_rho_all =
      count_all > 0 ? sum_rho_all / static_cast<double>(count_all) : 1.0;
  col.l2_node_vs_openlb_interp = RelativeL2(sim, ref);
  return col;
}

inline double BandFluidL2VsOpenLbInterp(const CenterlineColumn& col,
                                        int iy_min, int iy_max) {
  std::vector<double> sim;
  std::vector<double> ref;
  for (const CenterlineColumnRow& row : col.rows) {
    if (row.iy < iy_min || row.iy > iy_max) {
      continue;
    }
    if (std::strcmp(row.kind, "fluid") != 0) {
      continue;
    }
    sim.push_back(row.ux_over_lid);
    ref.push_back(row.openlb_interp);
  }
  return sim.empty() ? 0.0 : RelativeL2(sim, ref);
}

struct CavityCenterlineDiagnostics {
  double relative_l2 = 0.0;
  std::vector<double> simulated;
  std::vector<double> point_errors;
};

inline CavityCenterlineDiagnostics AnalyzeGhiaCenterline(
    const std::vector<double>& simulated) {
  CavityCenterlineDiagnostics diag;
  diag.simulated = simulated;
  diag.point_errors.resize(kGhiaRe100_UOverLid.size());
  for (std::size_t i = 0; i < simulated.size(); ++i) {
    diag.point_errors[i] = simulated[i] - kGhiaRe100_UOverLid[i];
  }
  diag.relative_l2 = RelativeL2(simulated, std::vector<double>(
      kGhiaRe100_UOverLid.begin(), kGhiaRe100_UOverLid.end()));
  return diag;
}

inline double RelativeL2VsOpenLb(const std::vector<double>& simulated) {
  return RelativeL2(simulated, std::vector<double>(
      kOpenLbCavity3d_UOverLid.begin(), kOpenLbCavity3d_UOverLid.end()));
}

inline void PrintGhiaDiagnostics(const UnitConverter& converter, int steps,
                                 const CavityCenterlineDiagnostics& diag,
                                 std::ostream& out = std::cout) {
  out << "OctLB cavity3d diagnostics\n";
  out << "  N=" << converter.resolution()
      << " tau=" << converter.lattice_relaxation_time()
      << " omega=" << converter.omega() << "\n";
  out << "  Re=" << converter.reynolds()
      << " u_lid_lattice=" << converter.char_lattice_velocity()
      << " iT=" << steps << "\n";
  out << "  relative_L2_ghia=" << diag.relative_l2
      << " relative_L2_openlb=" << RelativeL2VsOpenLb(diag.simulated) << "\n";
  for (std::size_t i = 0; i < kGhiaRe100_YOverH.size(); ++i) {
    out << "  y/H=" << kGhiaRe100_YOverH[i] << " sim=" << diag.simulated[i]
        << " ghia=" << kGhiaRe100_UOverLid[i]
        << " err=" << diag.point_errors[i] << "\n";
  }
}

inline bool CavityWriteVtkEnabled() {
  const char* flag = std::getenv("OCTLB_CAVITY_WRITE_VTK");
  return flag != nullptr && flag[0] == '1' && flag[1] == '\0';
}

}  // namespace octlb

#endif  // OCTLB_EXAMPLES_CAVITY3D_CAVITY3D_CASE_H_
