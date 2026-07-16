#ifndef OCTLB_EXAMPLES_CYLINDER3D_CYLINDER3D_CASE_H_
#define OCTLB_EXAMPLES_CYLINDER3D_CYLINDER3D_CASE_H_

#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/geometry/geometry_config.h"
#include "src/mesh/geometry/geometry_engine.h"
#include "src/mesh/geometry/geometry_types.h"
#include "src/mesh/geometry/material_field.h"
#include "src/mesh/io/stl_reader/stl_reader.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/ghost_schedule.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/inlet_velocity_field.h"
#include "src/solver/lbm/bouzidi_link_data.h"
#include "src/solver/lbm/bc_installer.h"
#include "src/solver/lbm/domain_boundary_handler.h"
#include "src/solver/lbm/drag_postprocessor.h"
#include "src/solver/lbm/lattice_material_init.h"
#include "src/solver/lbm/level_coupler.h"
#include "src/solver/lbm/time_loop/time_loop.h"

namespace octlb {

// W2 (T11) uniform cylinder3d: STL cylinder in a uniform (no-AMR) octree grid,
// flow in +x. Sanity oracle (runs / mass conserved / drag finite & positive /
// no block-edge artifact). AMR + strict Cd<1% gate come in W3/W4.
//
// BC mode is legacy (boundary_lattice_mode=false): inlet = Zou-He velocity
// (kMovingLid u=(U,0,0)), outlet = zero-gradient outflow (kOutflow), walls =
// no-slip (kNoSlip), cylinder surface = Bouzidi (kBoundary, driven by the
// MaterialField). Legacy mode is required because boundary_lattice_mode=true
// (InterpolatedVelocity) would collide *all* kBoundary cells -- including the
// cylinder surface -- with a tree-face Dirichlet velocity. The overlap-padding
// machinery (cavity3d's InterpolatedVelocity calibration) is not used in W2.

using CylinderLattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

// Domain encloses the cylinder fixture (axis z, R=0.5, centred (0,0),
// z in [0,1]): inlet at x=-1, cylinder at x=0, outlet at x=2 (1.5 wake).
inline BoundingBox CylinderDomain() {
  return {-1.0, -1.0, 0.0, 2.0, 1.0, 1.0};
}

// Uniform grid: refine every octant to max_level -> 8^max_level same-level
// octants, balance, partition. W2 uses max_level=1 (8 octants) so multi-block
// ghost exchange (and the ② edge-ghost path) is exercised.
inline OctreeForest MakeUniformForest(MPI_Comm comm, const BoundingBox& domain,
                                      int max_level = 1) {
  OctreeForest forest(comm, domain);
  forest.refine([](OctantId) { return true; }, max_level);
  forest.balance();
  forest.partition();
  return forest;
}

inline GeometryAssembly CylinderAssembly(const std::string& stl_path) {
  GeometryAssembly assembly;
  GeometryPart part;
  part.name = "cylinder";
  part.role = GeometryPartRole::kExternalObstacle;
  part.priority = 0;
  part.soup = read_stl_file(stl_path);
  assembly.parts.push_back(std::move(part));
  return assembly;
}

inline std::vector<DomainBcSpec> CylinderBcSpecs(double u_inlet) {
  std::vector<DomainBcSpec> specs;
  DomainBcSpec inlet;
  inlet.face = FaceDir::kXMin;
  inlet.type = DomainBcType::kInterpolatedVelocity;  // velocity Dirichlet inlet
  inlet.inlet_field = std::make_shared<boundary::UniformInletProfile>(
      std::array<double, 3>{u_inlet, 0.0, 0.0});
  specs.push_back(inlet);
  DomainBcSpec outlet;
  outlet.face = FaceDir::kXMax;
  outlet.type = DomainBcType::kInterpolatedPressure;  // pressure outlet p=0
  outlet.rho_target = 1.0;
  specs.push_back(outlet);
  for (FaceDir wall : {FaceDir::kYMin, FaceDir::kYMax, FaceDir::kZMin,
                       FaceDir::kZMax}) {
    DomainBcSpec s;
    s.face = wall;
    s.type = DomainBcType::kInterpolatedVelocity;  // FD no-slip (u=0), cavity3d-proven
    specs.push_back(s);
  }
  return specs;
}

struct Cylinder3dCase {
  BoundingBox domain;
  int n;             // lattice cells per octant per axis (cell_width)
  double omega;
  double u_inlet;
  double rho0;
  OctreeForest forest;
  GeometryAssembly assembly;
  MaterialField material;
  BouzidiLinkData bouzidi;
  FacePairList face_pairs;
  BlockCollection<CylinderLattice> blocks;
  GhostSchedule<CylinderLattice> ghosts;
  LevelCoupler coupler;
  ConcreteDomainBoundaryHandler domain_bc;
  TimeLoop loop;

  Cylinder3dCase(MPI_Comm comm, const std::string& stl_path, int n_cells,
                 double omega_in, double u_inlet_in, double rho0_in,
                 int max_level = 1)
      : domain(CylinderDomain()),
        n(n_cells),
        omega(omega_in),
        u_inlet(u_inlet_in),
        rho0(rho0_in),
        forest(MakeUniformForest(comm, domain, max_level)),
        assembly(CylinderAssembly(stl_path)),
        material(GeometryEngine{}.build(
            forest, assembly,
            GeometryConfig{.max_level = 0,
                           .resolve_surface_times = 0,
                           .cell_width = n_cells})),
        bouzidi(BouzidiLinkData::Build(forest, material, assembly)),
        face_pairs(forest),
        blocks(forest.local_num_octants(),
               [this](OctantId id) {
                 CylinderLattice lat(n, n, n, /*halo=*/1);
                 const double u0[3] = {0.0, 0.0, 0.0};
                 initialize_from_material(id, lat, material, rho0, u0, omega);
                 lat.set_bouzidi_links(&bouzidi);
                 return lat;
               }),
        ghosts(comm, face_pairs, blocks, n, n, n),
        coupler(comm, face_pairs, forest, blocks, n, n, n, omega),
        domain_bc(blocks, face_pairs.tree_boundary_faces(),
                  CylinderBcSpecs(u_inlet), n, n, n, omega),
        loop(forest, blocks, ghosts, coupler, domain_bc, omega,
             /*use_const_rho_bgk=*/false) {
      // Per-cell dispatch: stamp the domain-outer face cells (inlet/outlet/
      // walls) per spec. The cylinder surface keeps kBouzidi from
      // initialize_from_material; the fluid interior keeps kBulk.
      bc::StampTreeBoundaryCells(blocks, face_pairs.tree_boundary_faces(),
                                 CylinderBcSpecs(u_inlet), n, n, n);
    }

  void advance_steps(int num_steps) {
    for (int step = 0; step < num_steps; ++step) {
      loop.advance_one();
    }
  }

  label num_blocks() const { return forest.local_num_octants(); }

  // Momentum-exchange drag aggregated over all local blocks and reduced across
  // ranks. Returns force on the fluid (per W1 MomentumExchangeDrag convention).
  std::array<double, 3> aggregate_force_on_fluid(MPI_Comm comm) const {
    std::array<double, 3> local{0.0, 0.0, 0.0};
    for (OctantId id = 0; id < forest.local_num_octants(); ++id) {
      const MomentumExchangeDrag drag(blocks[id]);
      const std::array<double, 3> f = drag.force_on_fluid();
      local[0] += f[0];
      local[1] += f[1];
      local[2] += f[2];
    }
    std::array<double, 3> global{0.0, 0.0, 0.0};
    MPI_Allreduce(local.data(), global.data(), 3, MPI_DOUBLE, MPI_SUM, comm);
    return global;
  }

  // Cd = 2 * (-force_on_fluid . flow_dir) / (rho * u^2 * area).
  double drag_coefficient(MPI_Comm comm, double area,
                          std::array<double, 3> flow_dir) const {
    const std::array<double, 3> f = aggregate_force_on_fluid(comm);
    const double fdot =
        f[0] * flow_dir[0] + f[1] * flow_dir[1] + f[2] * flow_dir[2];
    return 2.0 * (-fdot) / (rho0 * u_inlet * u_inlet * area);
  }

  // Total fluid mass, summed over local fluid cells and reduced across ranks.
  double total_mass(MPI_Comm comm) const {
    double local = 0.0;
    for (OctantId id = 0; id < forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      for (int iz = 0; iz < n; ++iz) {
        for (int iy = 0; iy < n; ++iy) {
          for (int ix = 0; ix < n; ++ix) {
            if (lat.bc_kind(ix, iy, iz) != BcKind::kBulk) {
              continue;
            }
            double rho = 0.0;
            double u[3] = {};
            lat.get(ix, iy, iz).computeRhoU(rho, u);
            local += rho;
          }
        }
      }
    }
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, comm);
    return global;
  }

  bool has_non_finite_velocity() const {
    for (OctantId id = 0; id < forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      for (int iz = 0; iz < n; ++iz) {
        for (int iy = 0; iy < n; ++iy) {
          for (int ix = 0; ix < n; ++ix) {
            if (lat.bc_kind(ix, iy, iz) == BcKind::kSolid) {
              continue;
            }
            double rho = 0.0;
            double u[3] = {};
            lat.get(ix, iy, iz).computeRhoU(rho, u);
            for (int d = 0; d < 3; ++d) {
              if (!std::isfinite(u[d])) {
                return true;
              }
            }
          }
        }
      }
    }
    return false;
  }

  // Map a physical point to the local octant + cell that contains it and read
  // (rho, u_x). Used by the W2 ② probe: compare the multi-block field (with
  // inter-block ghost) against a single-block reference (no inter-block ghost)
  // at the same physical cell centres -- same dx, so a mismatch localizes the
  // edge-ghost (②) corruption.
  std::pair<double, double> sample_rho_ux(scalar x, scalar y, scalar z) const {
    const auto d = sample_diag(x, y, z);
    return {d.rho, d.ux};
  }

  struct SampleDiag {
    double rho = 0.0;
    double ux = 0.0;
    BcKind kind = BcKind::kSolid;
    int edge_count = -1;  // # of local coords on a block boundary (0=interior,
                          // 1=face, 2=edge, 3=corner); -1 = not found locally
  };

  // Same as sample_rho_ux but also returns the cell kind and how many of the
  // cell's local coords sit on a block boundary -- used to categorize which
  // cells diverge in the ② probe.
  SampleDiag sample_diag(scalar x, scalar y, scalar z) const {
    for (OctantId id = 0; id < forest.local_num_octants(); ++id) {
      const BoundingBox b = forest.quadrant_bounds(id);
      if (x < b.x_min || x > b.x_max || y < b.y_min || y > b.y_max ||
          z < b.z_min || z > b.z_max) {
        continue;
      }
      const double ddx = (b.x_max - b.x_min) / static_cast<double>(n);
      const double ddy = (b.y_max - b.y_min) / static_cast<double>(n);
      const double ddz = (b.z_max - b.z_min) / static_cast<double>(n);
      int ix = static_cast<int>((x - b.x_min) / ddx);
      int iy = static_cast<int>((y - b.y_min) / ddy);
      int iz = static_cast<int>((z - b.z_min) / ddz);
      if (ix < 0) ix = 0;
      if (ix > n - 1) ix = n - 1;
      if (iy < 0) iy = 0;
      if (iy > n - 1) iy = n - 1;
      if (iz < 0) iz = 0;
      if (iz > n - 1) iz = n - 1;
      const CylinderLattice& lat = blocks[id];
      double rho = 0.0;
      double u[3] = {};
      lat.get(ix, iy, iz).computeRhoU(rho, u);
      const int ec = (ix == 0 || ix == n - 1) + (iy == 0 || iy == n - 1) +
                     (iz == 0 || iz == n - 1);
      return {rho, u[0], lat.bc_kind(ix, iy, iz), ec};
    }
    return {0.0, 0.0, BcKind::kSolid, -1};
  }
};

}  // namespace octlb

#endif  // OCTLB_EXAMPLES_CYLINDER3D_CYLINDER3D_CASE_H_