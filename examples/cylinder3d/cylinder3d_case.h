#ifndef OCTLB_EXAMPLES_CYLINDER3D_CYLINDER3D_CASE_H_
#define OCTLB_EXAMPLES_CYLINDER3D_CYLINDER3D_CASE_H_

#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
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
#include "src/mesh/geometry/triangle_soup.h"
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
#include "src/solver/lbm/vtk_lbm_fields.h"
#include "src/solver/io/vtk_writer/structured_grid_writer.h"
#include "src/solver/io/vtk_writer/vtk_block_meta.h"
#include "src/solver/io/vtk_writer/vtk_cell_field.h"
#include "src/solver/io/vtk_writer/parallel_vtk_index.h"
#include <fstream>
#include <iomanip>

namespace octlb {

// W3 (T11) Schäfer-Turek cylinder3d, option A (cubic domain + carved channel):
// the computation domain is a CUBE [0,2.5]^3 (so the single-tree octree yields
// isotropic cubic cells), and the 2.5 x 0.41 x 0.41 channel is carved into it
// via a kInternalChannel part (outer surface = the cube, inner cavity = the
// channel duct): the channel lumen is fluid, the cube-minus-channel wall is
// solid. The cylinder (D=0.1, axis z) sits in the channel as a kExternalObstacle
// (Bouzidi surface). Inlet is the channel cross-section patch on the x=0 cube
// face (velocity Dirichlet); outlet the patch on x=2.5 (pressure Dirichlet);
// the channel walls (y,z = 1.045,1.455) are internal Bouzidi; the four lateral
// cube faces are deep solid (skipped by the geometry-aware stamp).
//
// W3 uses a uniform inlet (UniformInletProfile); the Schäfer-Turek Poiseuille
// inlet requires a physical-coordinate inlet field (the inlet face is a patch
// on the larger cube face, not the full channel cross-section), landed in W4.

using CylinderLattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

// Channel centered in the cube so all four channel walls are internal (Bouzidi)
// and the four lateral cube faces are deep solid. Cylinder keeps its
// Schäfer-Turek position relative to the walls (0.2 off the y=wall, 0.205 off
// the z=wall): cylinder centre (0.45, 1.245, 1.25), D=0.1, z in [1.045,1.455].
inline constexpr double kCubeSide = 2.5;
inline constexpr double kChannelHalf = 0.205;       // channel half-width 0.41/2
inline constexpr double kChannelCenterYZ = 1.25;    // = cube center
inline constexpr double kChannelLo = kChannelCenterYZ - kChannelHalf;  // 1.045
inline constexpr double kChannelHi = kChannelCenterYZ + kChannelHalf;  // 1.455

inline BoundingBox CylinderDomain() {
  return {0.0, 0.0, 0.0, kCubeSide, kCubeSide, kCubeSide};
}

// Build an axis-aligned box triangle soup (12 triangles) with per-axis bounds.
// inward_normals = true gives normals pointing into the box (used for the
// kInternalChannel inner cavity so is_inside marks the lumen).
inline TriangleSoup BoxSoup(scalar xmin, scalar ymin, scalar zmin, scalar xmax,
                            scalar ymax, scalar zmax, bool inward_normals) {
  TriangleSoup soup;
  const auto p000 = std::array<scalar, 3>{xmin, ymin, zmin};
  const auto p100 = std::array<scalar, 3>{xmax, ymin, zmin};
  const auto p010 = std::array<scalar, 3>{xmin, ymax, zmin};
  const auto p110 = std::array<scalar, 3>{xmax, ymax, zmin};
  const auto p001 = std::array<scalar, 3>{xmin, ymin, zmax};
  const auto p101 = std::array<scalar, 3>{xmax, ymin, zmax};
  const auto p011 = std::array<scalar, 3>{xmin, ymax, zmax};
  const auto p111 = std::array<scalar, 3>{xmax, ymax, zmax};
  auto add = [&](const std::array<scalar, 3>& a,
                 const std::array<scalar, 3>& b,
                 const std::array<scalar, 3>& c) {
    Triangle tri;
    tri.v0 = a;
    tri.v1 = b;
    tri.v2 = c;
    const scalar ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
    const scalar vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2];
    tri.normal[0] = uy * vz - uz * vy;
    tri.normal[1] = uz * vx - ux * vz;
    tri.normal[2] = ux * vy - uy * vx;
    const scalar len = std::sqrt(tri.normal[0] * tri.normal[0] +
                                 tri.normal[1] * tri.normal[1] +
                                 tri.normal[2] * tri.normal[2]);
    if (len > 0.0) {
      tri.normal[0] /= len;
      tri.normal[1] /= len;
      tri.normal[2] /= len;
    }
    if (inward_normals) {
      // Reverse winding -> outward normal flips to inward.
      tri.v1 = c;
      tri.v2 = b;
      tri.normal[0] = -tri.normal[0];
      tri.normal[1] = -tri.normal[1];
      tri.normal[2] = -tri.normal[2];
    }
    soup.add_triangle(tri);
  };
  add(p000, p100, p010);
  add(p100, p110, p010);
  add(p100, p110, p101);
  add(p100, p101, p110);
  add(p110, p011, p101);
  add(p110, p011, p111);
  add(p010, p110, p011);
  add(p010, p011, p001);
  add(p001, p011, p101);
  add(p001, p101, p100);
  add(p000, p010, p001);
  add(p000, p001, p100);
  soup.update_bounding_box();
  return soup;
}

inline GeometryAssembly CylinderAssembly(const std::string& stl_path,
                                          bool include_cylinder = true) {
  GeometryAssembly assembly;
  // Channel duct: outer surface = the cube (in_wall = cube-minus-channel =
  // solid), inner cavity = the channel box (anisotropic: 2.5 x 0.41 x 0.41).
  GeometryPart channel;
  channel.name = "channel";
  channel.role = GeometryPartRole::kInternalChannel;
  channel.priority = 0;
  channel.soup = BoxSoup(0.0, 0.0, 0.0, kCubeSide, kCubeSide, kCubeSide, false);
  channel.inner_cavity_soup =
      BoxSoup(0.0, kChannelLo, kChannelLo, kCubeSide, kChannelHi, kChannelHi,
              true);
  assembly.parts.push_back(std::move(channel));
  // Cylinder obstacle (higher priority overlays the channel fluid). Skipped for
  // the ② clean re-probe (straight duct, no surface voxelization noise).
  if (include_cylinder) {
    GeometryPart cylinder;
    cylinder.name = "cylinder";
    cylinder.role = GeometryPartRole::kExternalObstacle;
    cylinder.priority = 1;
    cylinder.soup = read_stl_file(stl_path);
    assembly.parts.push_back(std::move(cylinder));
  }
  return assembly;
}

// Schäfer-Turek face roles. The four lateral cube faces (y,z = 0, kCubeSide) are
// deep solid (no spec -- the geometry-aware stamp skips their kSolid cells).
// Inlet = channel patch on x=0 (velocity Dirichlet); outlet = patch on x=2.5
// (pressure Dirichlet, p=0). Walls of the channel are internal (Bouzidi via
// StampFromMaterial), not domain faces.
inline std::vector<DomainBcSpec> CylinderBcSpecs(double u_inlet,
                                                 double ramp_end_t) {
  std::vector<DomainBcSpec> specs;
  DomainBcSpec inlet;
  inlet.face = FaceDir::kXMin;
  inlet.type = DomainBcType::kInterpolatedVelocity;  // velocity Dirichlet inlet
  inlet.inlet_field = std::make_shared<boundary::UniformInletProfile>(
      std::array<double, 3>{u_inlet, 0.0, 0.0}, ramp_end_t);
  specs.push_back(inlet);
  DomainBcSpec outlet;
  outlet.face = FaceDir::kXMax;
  outlet.type = DomainBcType::kInterpolatedPressure;  // pressure outlet p=0
  outlet.rho_target = 1.0;
  specs.push_back(outlet);
  return specs;
}

enum class Cylinder3dMode { kUniform, kAmr };

struct Cylinder3dConfig {
  Cylinder3dMode mode = Cylinder3dMode::kUniform;
  std::string stl_path;
  int n = 8;                   // lattice cells per octant per axis
  int max_level = 2;           // uniform refine level (kUniform) / AMR finest (kAmr)
  double omega = 1.0 / 0.53;   // ~1.887, tau=0.53 (Schäfer-Turek Re=20)
  double u_inlet = 0.02;       // latticeU mean (Poiseuille peak 2.25x in W4)
  double rho0 = 1.0;
  double ramp_end_t = 0.0;     // 0 disables the inlet ramp
  double charU = 0.2;          // physical char velocity (for drag Cd)
  // ② W3 clean re-probe: include_cylinder=false builds a straight channel duct
  // (no obstacle) so block edge/corner artifacts can be isolated from the
  // cylinder-surface voxelization noise. enable_edge_exchange=false builds the
  // no-edge-exchange baseline (Stage B cross-rank edge ghosts stay unfilled).
  bool include_cylinder = true;
  bool enable_edge_exchange = true;
};

// Forest + material built together: GeometryEngine.build mutates the forest
// (refine/balance) AND returns the material, so they must be produced in one
// place. partition() runs after the geometry refines the forest.
struct ForestAndMaterial {
  OctreeForest forest;
  MaterialField material;
};

inline ForestAndMaterial MakeUniformForestAndMaterial(
    MPI_Comm comm, const BoundingBox& domain, const GeometryAssembly& assembly,
    const Cylinder3dConfig& cfg) {
  OctreeForest forest(comm, domain);
  forest.refine([](OctantId) { return true; }, cfg.max_level);
  forest.balance();
  forest.partition();
  GeometryConfig gcfg;
  gcfg.max_level = 0;
  gcfg.resolve_surface_times = 0;
  gcfg.cell_width = cfg.n;
  MaterialField material = GeometryEngine{}.build(forest, assembly, gcfg);
  return {std::move(forest), std::move(material)};
}

inline ForestAndMaterial MakeAmrForestAndMaterial(
    MPI_Comm comm, const BoundingBox& domain, const GeometryAssembly& assembly,
    const Cylinder3dConfig& cfg) {
  OctreeForest forest(comm, domain);
  GeometryConfig gcfg;
  gcfg.max_level = cfg.max_level;
  gcfg.resolve_surface_times = 3;
  gcfg.bound_width = 0.005;
  gcfg.wake_length = 2.0;  // cylinder (~x=0.5) to outlet x=2.5
  gcfg.wake_direction = 0;
  gcfg.cell_width = cfg.n;
  // GeometryEngine.build refines (resolve_bounding + resolve_surface) and
  // partitions (resolve_surface partitions on its first iteration), then
  // voxelizes the partitioned forest -- so the returned material matches the
  // partitioned local octant ids. Do NOT re-partition after build (it would
  // invalidate the material layout). But build's later resolve_surface
  // iterations refine WITHOUT re-partitioning, so the face-ghost is stale; rebuild
  // it (same partition) so FacePairList sees remote octants at their current
  // level (stale ghost -> asymmetric same-level/coarse-fine enumeration across
  // ranks -> MPI deadlock in GhostSchedule).
  MaterialField material = GeometryEngine{}.build(forest, assembly, gcfg);
  forest.RebuildGhost();
  return {std::move(forest), std::move(material)};
}

struct Cylinder3dCase {
  Cylinder3dConfig cfg;
  BoundingBox domain;
  GeometryAssembly assembly;
  ForestAndMaterial fam;
  BouzidiLinkData bouzidi;
  FacePairList face_pairs;
  BlockCollection<CylinderLattice> blocks;
  GhostSchedule<CylinderLattice> ghosts;
  LevelCoupler coupler;
  ConcreteDomainBoundaryHandler domain_bc;
  TimeLoop loop;

  Cylinder3dCase(MPI_Comm comm, Cylinder3dConfig config)
      : cfg(std::move(config)),
        domain(CylinderDomain()),
        assembly(CylinderAssembly(cfg.stl_path, cfg.include_cylinder)),
        fam(cfg.mode == Cylinder3dMode::kUniform
                ? MakeUniformForestAndMaterial(comm, domain, assembly, cfg)
                : MakeAmrForestAndMaterial(comm, domain, assembly, cfg)),
        bouzidi(BouzidiLinkData::Build(fam.forest, fam.material, assembly)),
        face_pairs(fam.forest),
        blocks(fam.forest.local_num_octants(),
               [this](OctantId id) {
                 CylinderLattice lat(cfg.n, cfg.n, cfg.n, /*halo=*/1);
                 const double u0[3] = {0.0, 0.0, 0.0};
                 initialize_from_material(id, lat, fam.material, cfg.rho0, u0,
                                          cfg.omega);
                 lat.set_bouzidi_links(&bouzidi);
                 return lat;
               }),
        ghosts(comm, face_pairs, blocks, cfg.n, cfg.n, cfg.n,
               cfg.enable_edge_exchange),
        coupler(comm, face_pairs, fam.forest, blocks, cfg.n, cfg.n, cfg.n,
                cfg.omega),
        domain_bc(blocks, face_pairs.tree_boundary_faces(),
                  CylinderBcSpecs(cfg.u_inlet, cfg.ramp_end_t), cfg.n, cfg.n,
                  cfg.n, cfg.omega),
        loop(fam.forest, blocks, ghosts, coupler, domain_bc, cfg.omega,
             /*use_const_rho_bgk=*/false) {
    // Per-cell dispatch: stamp geometry (lumen kBulk, cube-minus-channel kSolid,
    // channel walls + cylinder kBouzidi) then the domain-outer face cells per
    // spec, geometry-aware (skipping the solid exterior on the inlet/outlet).
    bc::StampFromMaterial(blocks, fam.material);
    bc::StampTreeBoundaryCells(blocks, face_pairs.tree_boundary_faces(),
                               CylinderBcSpecs(cfg.u_inlet, cfg.ramp_end_t),
                               cfg.n, cfg.n, cfg.n);
    // Diagnostic (T11 W3 AMR deadlock hunt): mark construct-done so the CI log
    // distinguishes a construct-phase hang (collective Allreduce) from a
    // step-loop hang (LevelCoupler/GhostSchedule exchange). Remove once green.
    int dbg_rank = 0;
    MPI_Comm_rank(comm, &dbg_rank);
    if (dbg_rank == 0) {
      std::fprintf(stderr,
                   "[cyl] constructed mode=%d max_level=%d n=%d u_inlet=%g "
                   "blocks=%lld faces=%zu cf=%zu edges=%zu\n",
                   static_cast<int>(cfg.mode), cfg.max_level, cfg.n, cfg.u_inlet,
                   static_cast<long long>(fam.forest.local_num_octants()),
                   face_pairs.same_level_faces().size(),
                   face_pairs.coarse_fine_faces().size(),
                   face_pairs.cross_rank_edges().size());
      std::fflush(stderr);
    }
  }

  void advance_steps(int num_steps) {
    for (int step = 0; step < num_steps; ++step) {
      loop.advance_one();
    }
  }

  label num_blocks() const { return fam.forest.local_num_octants(); }

  std::array<double, 3> aggregate_force_on_fluid(MPI_Comm comm) const {
    std::array<double, 3> local{0.0, 0.0, 0.0};
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
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

  // Cd = 2 * (-force_on_fluid . flow_dir) / (rho * u^2 * A). A is the cylinder
  // frontal area (D * z_span); its exact value matters for W4 accuracy.
  double drag_coefficient(MPI_Comm comm, double area,
                          std::array<double, 3> flow_dir) const {
    const std::array<double, 3> f = aggregate_force_on_fluid(comm);
    const double fdot =
        f[0] * flow_dir[0] + f[1] * flow_dir[1] + f[2] * flow_dir[2];
    return 2.0 * (-fdot) / (cfg.rho0 * cfg.u_inlet * cfg.u_inlet * area);
  }

  double total_mass(MPI_Comm comm) const {
    double local = 0.0;
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      for (int iz = 0; iz < cfg.n; ++iz) {
        for (int iy = 0; iy < cfg.n; ++iy) {
          for (int ix = 0; ix < cfg.n; ++ix) {
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
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      for (int iz = 0; iz < cfg.n; ++iz) {
        for (int iy = 0; iy < cfg.n; ++iy) {
          for (int ix = 0; ix < cfg.n; ++ix) {
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

  // ② W3 clean re-probe: aggregate the worst block edge/corner cell of every
  // local block. Corner cells are those touching a block boundary in all three
  // axes (ix,iy,iz in {0, n-1}); edge cells touch a boundary in exactly two
  // axes. These are the cells whose stream pulls from edge/corner ghosts -- the
  // site of the ② artifact when cross-rank edge ghosts are not exchanged.
  // Returns local (per-rank) maxima; the test Allreduces across ranks.
  struct EdgeCornerStats {
    double max_u = 0.0;
    double max_rho_dev = 0.0;
  };
  EdgeCornerStats local_edge_corner_stats() const {
    EdgeCornerStats s;
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      for (int iz = 0; iz < cfg.n; ++iz) {
        for (int iy = 0; iy < cfg.n; ++iy) {
          for (int ix = 0; ix < cfg.n; ++ix) {
            const int bx = (ix == 0 || ix == cfg.n - 1) ? 1 : 0;
            const int by = (iy == 0 || iy == cfg.n - 1) ? 1 : 0;
            const int bz = (iz == 0 || iz == cfg.n - 1) ? 1 : 0;
            if (bx + by + bz < 2) continue;  // edge (2) or corner (3) only
            // Only interior fluid cells at a block edge: kBulk cells whose
            // stream pulls from a neighbour block's edge ghost (the ② site).
            // Domain BC cells (inlet/outlet Dirichlet) and Bouzidi walls
            // legitimately deviate from rho0 -- they are not ② artifacts.
            if (lat.bc_kind(ix, iy, iz) != BcKind::kBulk) continue;
            double rho = 0.0;
            double u[3] = {};
            lat.get(ix, iy, iz).computeRhoU(rho, u);
            const double umag =
                std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
            if (std::isfinite(umag)) s.max_u = std::max(s.max_u, umag);
            const double rd = std::abs(rho - cfg.rho0);
            if (std::isfinite(rd)) s.max_rho_dev = std::max(s.max_rho_dev, rd);
          }
        }
      }
    }
    return s;
  }

  // W4 VTK output: one StructuredGrid .vts per local octant with per-block
  // velocity + pressure fields, plus a .vtm/.pvd index for ParaView. Each block
  // gets its OWN field view bound to that block's lattice (the AmrVtkWriter
  // shared-field loop would alias every block to block 0; we loop blocks here
  // and bind per block so multi-block AMR VTK is correct).
  void write_vtk_timestep(MPI_Comm comm, int iT, const std::string& output_dir,
                          const std::string& base_name = "cylinder3d") const {
    std::vector<std::string> local_paths;
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      VtkVelocityField<double, olb::descriptors::D3Q19<>> velocity(lat);
      VtkPressureField<double, olb::descriptors::D3Q19<>> pressure(lat);
      const VtkCellFieldView vel_field = VtkCellFieldView::From(velocity);
      const VtkCellFieldView pres_field = VtkCellFieldView::From(pressure);
      const std::array<VtkCellFieldView, 2> fields = {vel_field, pres_field};
      VtkBlockMeta meta;
      meta.id = id;
      meta.nx = cfg.n;
      meta.ny = cfg.n;
      meta.nz = cfg.n;
      meta.bounds = fam.forest.quadrant_bounds(id);
      const std::string path =
          output_dir + "/" + base_name + "_block" + std::to_string(id) + "_T" +
          std::to_string(iT) + ".vts";
      WriteStructuredGridVts(
          path, meta,
          std::span<const VtkCellFieldView>(fields.data(), fields.size()));
      local_paths.push_back(path);
    }
    const std::vector<std::string> rel_paths =
        CollectVtsRelativePaths(comm, local_paths, output_dir);
    WriteVtmIndex(comm, iT, output_dir, base_name, rel_paths);
    if (GetMpiRank(comm) == 0) {
      AppendPvdTimestep(iT, output_dir + "/" + base_name + ".pvd",
                        base_name + "_T" + std::to_string(iT) + ".vtm");
    }
  }
};

}  // namespace octlb

#endif  // OCTLB_EXAMPLES_CYLINDER3D_CYLINDER3D_CASE_H_