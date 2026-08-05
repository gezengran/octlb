#ifndef OCTLB_EXAMPLES_CYLINDER3D_CYLINDER3D_CASE_H_
#define OCTLB_EXAMPLES_CYLINDER3D_CYLINDER3D_CASE_H_

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
  // Watertight axis-aligned box: 6 faces x 2 triangles, each face wound CCW
  // viewed from OUTSIDE so normals point outward (the inward_normals branch
  // below reverses winding + flips the normal). The earlier hand-written
  // triangulation had a duplicate triangle on the +x face, missed the second
  // triangle on the +z face, and added a spurious interior diagonal triangle
  // (p110,p011,p101) -- the cube was non-watertight, so CgalSurfaceMesh::
  // is_inside's ray parity returned false for interior points and the
  // kInternalChannel voxelization fell back to kFluid for the cube-minus-duct
  // "earth" (158k spurious fluid cells at max_level=4), corrupting the flow.
  add(p000, p001, p011);  // x=xmin face
  add(p000, p011, p010);
  add(p100, p110, p111);  // x=xmax face
  add(p100, p111, p101);
  add(p000, p100, p101);  // y=ymin face
  add(p000, p101, p001);
  add(p010, p011, p111);  // y=ymax face
  add(p010, p111, p110);
  add(p000, p010, p110);  // z=zmin face
  add(p000, p110, p100);
  add(p001, p101, p111);  // z=zmax face
  add(p001, p111, p011);
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

// Bounding box of the cylinder obstacle (the kExternalObstacle part's STL
// soup), expanded by `pad` on every side. Used to restrict the momentum-
// exchange drag sum to the cylinder surface only: in this ducted cylinder3d
// both the curved cylinder surface AND the flat duct walls are stamped kBouzidi,
// so an unfiltered MEM sum would include the duct-wall shear (a long strip of
// Bouzidi cells along the whole duct) and inflate Cd by ~100x. The duct walls
// lie outside the cylinder's x/y projection, so a bbox test on the obstacle
// neighbour cleanly isolates the cylinder. `pad` should cover the surface cell
// ring (~2-3 finest cells) so no curved-surface cell just outside the STL bbox
// is dropped. Returns a huge box if there is no obstacle (uniform-duct probe),
// so the filter is a no-op there.
inline BoundingBox CylinderObstacleBbox(const GeometryAssembly& assembly,
                                        double pad) {
  for (const GeometryPart& part : assembly.parts) {
    if (part.role != GeometryPartRole::kExternalObstacle) continue;
    BoundingBox b = part.soup.bounding_box();
    b.x_min -= pad;
    b.y_min -= pad;
    b.z_min -= pad;
    b.x_max += pad;
    b.y_max += pad;
    b.z_max += pad;
    return b;
  }
  BoundingBox all;
  all.x_min = all.y_min = all.z_min = -1e30;
  all.x_max = all.y_max = all.z_max = 1e30;
  return all;
}

// Schäfer-Turek face roles. The four lateral cube faces (y,z = 0, kCubeSide) are
// deep solid (no spec -- the geometry-aware stamp skips their kSolid cells).
// Inlet = channel patch on x=0 (velocity Dirichlet); outlet = patch on x=2.5
// (pressure Dirichlet, p=0). Walls of the channel are internal (Bouzidi via
// StampFromMaterial), not domain faces.
inline std::vector<DomainBcSpec> CylinderBcSpecs(double u_inlet,
                                                 double ramp_end_t,
                                                 bool poiseuille_inlet) {
  std::vector<DomainBcSpec> specs;
  DomainBcSpec inlet;
  inlet.face = FaceDir::kXMin;
  inlet.type = DomainBcType::kInterpolatedVelocity;  // velocity Dirichlet inlet
  if (poiseuille_inlet) {
    // Schäfer-Turek Poiseuille duct profile on the channel cross-section patch
    // (y,z in [kChannelLo, kChannelHi]). The rectangular-duct product of two
    // parabolas has mean/peak = (2/3)^2 = 4/9. Schaefer-Turek defines BOTH Re and
    // Cd with the peak inlet velocity U_max, so cfg.u_inlet IS U_max (the profile
    // peak); the duct mean is u_inlet * 4/9. Re_lat = U_max * (D/dx) / nu_lat =
    // 0.02 * 10.24 / 0.01 = 20.5 (Schaefer-Turek Re=20). Evaluated from physical
    // cell position (BlockLattice::phys_origin) so the patch on the larger cube
    // face is correct regardless of the AMR octant layout.
    inlet.inlet_field =
        std::make_shared<boundary::ChannelPoiseuilleInletProfile>(
            FaceDir::kXMin, kChannelLo, kChannelHi, kChannelLo, kChannelHi,
            u_inlet, std::array<double, 3>{1.0, 0.0, 0.0}, ramp_end_t);
  } else {
    inlet.inlet_field = std::make_shared<boundary::UniformInletProfile>(
        std::array<double, 3>{u_inlet, 0.0, 0.0}, ramp_end_t);
  }
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
  double u_inlet = 0.02;       // Schäfer-Turek U_max (Poiseuille peak); Re=U_max*D/dx/nu=20
  double rho0 = 1.0;
  double ramp_end_t = 0.0;     // 0 disables the inlet ramp
  double charU = 0.2;          // physical char velocity (for drag Cd)
  // ② W3 clean re-probe: include_cylinder=false builds a straight channel duct
  // (no obstacle) so block edge/corner artifacts can be isolated from the
  // cylinder-surface voxelization noise. enable_edge_exchange=false builds the
  // no-edge-exchange baseline (Stage B cross-rank edge ghosts stay unfilled).
  bool include_cylinder = true;
  bool enable_edge_exchange = true;
  // W4: Schäfer-Turek Poiseuille inlet (channel duct profile) instead of the
  // uniform inlet used for W2/W3 sanity. Enabled for the alignment (magnitude
  // + 1% Cd) gates and the example; left false for the sanity cases.
  bool poiseuille_inlet = false;
  // T11 alignment: pre-refine the channel duct (inlet + development + around
  // cylinder + wake) to max_level BEFORE geometry build, so the material
  // voxelization sees the fine channel. Without this the geometry engine only
  // refines the *cylinder* bbox (resolve_bounding targets geom_bbox, not the
  // channel) -> channel bulk + inlet sit at max_level-2 (~10 cells across the
  // 0.41 duct) and the Poiseuille inlet cannot resolve -> Cd far off the
  // OpenLB reference. Opt-in: the sanity cases keep the cheap coarse mesh.
  bool refine_channel_to_finest = false;
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
  // T11 alignment: pre-refine the channel duct to the finest level so the
  // Poiseuille inlet resolves and the bulk channel matches the OpenLB uniform-dx
  // reference. The geometry engine otherwise refines only the *cylinder* bbox
  // (resolve_bounding) + the cylinder surface ring (resolve_surface), leaving
  // the channel bulk + inlet at max_level-2 (~10 cells across the 0.41 duct).
  // Done BEFORE build so the material voxelization matches the refined forest.
  if (cfg.refine_channel_to_finest) {
    auto overlap = [](const BoundingBox& a, const BoundingBox& b) {
      return !(a.x_max < b.x_min || a.x_min > b.x_max ||
               a.y_max < b.y_min || a.y_min > b.y_max ||
               a.z_max < b.z_min || a.z_min > b.z_max);
    };
    const BoundingBox channel_box{0.0, kChannelLo, kChannelLo, kCubeSide,
                                  kChannelHi, kChannelHi};
    // Recursive refine (p8est_refine_ext is recursive): a bounds-based
    // criterion cascades the WHOLE channel to max_level in one call. An
    // id-indexed mark would NOT work -- during recursive refine the first
    // child inherits the parent's local id, so only a 1-octant-wide chain
    // cascades (the original coarse-mesh root cause).
    forest.refine(
        [&](const BoundingBox& box, int level) {
          if (level >= cfg.max_level) return false;
          return overlap(box, channel_box);
        },
        cfg.max_level);
    forest.balance();
    // Distribute the refined channel across ranks BEFORE build. Without this
    // the pre-refine (run only on the root-owner rank) leaves every octant on
    // one rank, and GeometryEngine.build never repartitions (resolve_surface
    // breaks immediately once the channel is already at max_level) -> the sim
    // runs near-serial. build does not re-partition, so this partition sticks.
    forest.partition();
    if (std::getenv("OCTLB_CYL_DEBUG") != nullptr) {
      int local_n = static_cast<int>(forest.local_num_octants());
      int global_n = 0;
      MPI_Allreduce(&local_n, &global_n, 1, MPI_INT, MPI_SUM, comm);
      std::vector<int> hist(cfg.max_level + 1, 0);
      for (label id = 0; id < forest.local_num_octants(); ++id) {
        ++hist[forest.quadrant_level(id)];
      }
      std::vector<int> ghist(cfg.max_level + 1, 0);
      MPI_Allreduce(hist.data(), ghist.data(), cfg.max_level + 1, MPI_INT,
                    MPI_SUM, comm);
      int r = 0;
      MPI_Comm_rank(comm, &r);
      if (r == 0) {
        std::fprintf(stderr, "[cyl] after channel pre-refine: global=%d ",
                     global_n);
        for (int lv = 0; lv <= cfg.max_level; ++lv)
          std::fprintf(stderr, "L%d=%d ", lv, ghist[lv]);
        std::fprintf(stderr, "\n");
        std::fflush(stderr);
      }
    }
  }
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
  int my_rank_ = 0;
  // Bounding box of the cylinder obstacle (STL bbox + pad), used to restrict
  // the MEM drag sum to the cylinder surface and exclude the duct-wall Bouzidi.
  BoundingBox cylinder_obstacle_bbox_{};

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
                 // Physical origin + cell width so a physical-coordinate inlet
                 // field (channel Poiseuille) can be evaluated from domain
                 // position. The octant is cubic, so cell width is uniform.
                 const BoundingBox ob = fam.forest.quadrant_bounds(id);
                 const double cw =
                     (ob.x_max - ob.x_min) / static_cast<double>(cfg.n);
                 lat.set_phys_origin({ob.x_min, ob.y_min, ob.z_min}, cw);
                 lat.set_bouzidi_links(&bouzidi);
                 return lat;
               }),
        ghosts(comm, face_pairs, blocks, cfg.n, cfg.n, cfg.n,
               cfg.enable_edge_exchange),
        coupler(comm, face_pairs, fam.forest, blocks, cfg.n, cfg.n, cfg.n,
                cfg.omega),
        domain_bc(blocks, face_pairs.tree_boundary_faces(),
                  CylinderBcSpecs(cfg.u_inlet, cfg.ramp_end_t,
                                  cfg.poiseuille_inlet),
                  cfg.n, cfg.n, cfg.n, cfg.omega),
        loop(fam.forest, blocks, ghosts, coupler, domain_bc, cfg.omega,
             /*use_const_rho_bgk=*/false) {
    MPI_Comm_rank(comm, &my_rank_);
    // Cylinder obstacle bbox + a ~3-finest-cell pad, so the MEM drag filter
    // keeps the whole curved-surface ring while still excluding the duct walls
    // (which sit well outside the cylinder's x/y projection).
    const double finest_dx = kCubeSide / std::pow(2.0, cfg.max_level) /
                             static_cast<double>(cfg.n);
    cylinder_obstacle_bbox_ = CylinderObstacleBbox(assembly, 3.0 * finest_dx);
    // Per-cell dispatch: stamp geometry (lumen kBulk, cube-minus-channel kSolid,
    // channel walls + cylinder kBouzidi) then the domain-outer face cells per
    // spec, geometry-aware (skipping the solid exterior on the inlet/outlet).
    bc::StampFromMaterial(blocks, fam.material);
    bc::StampTreeBoundaryCells(blocks, face_pairs.tree_boundary_faces(),
                               CylinderBcSpecs(cfg.u_inlet, cfg.ramp_end_t,
                                               cfg.poiseuille_inlet),
                               cfg.n, cfg.n, cfg.n);
    if (std::getenv("OCTLB_CYL_DEBUG") != nullptr) {
      int dbg_rank = 0;
      MPI_Comm_rank(comm, &dbg_rank);
      std::fprintf(stderr,
                   "[cyl] rank=%d mode=%d max_level=%d u_inlet=%g "
                   "blocks=%lld faces=%zu cf=%zu edges=%zu\n",
                   dbg_rank, static_cast<int>(cfg.mode), cfg.max_level,
                   cfg.u_inlet,
                   static_cast<long long>(fam.forest.local_num_octants()),
                   face_pairs.same_level_faces().size(),
                   face_pairs.coarse_fine_faces().size(),
                   face_pairs.cross_rank_edges().size());
      std::fflush(stderr);
    }
  }

  void advance_steps(int num_steps) {
    // T11 MEM timing: snapshot every block's post-collide state right before
    // each stream so a subsequent drag_coefficient() reads the correct outgoing
    // f_i (pre-stream snapshot) + bounced f_bar (post-stream live). Without this
    // the MEM falls back to the post-stream live array and overestimates Cd
    // ~6.6x. Always-on: the memcpy per step is negligible vs collide/stream,
    // and it makes drag_coefficient correct by default (no caller footgun).
    loop.set_snapshot_post_collide(true);
    for (int step = 0; step < num_steps; ++step) {
      loop.advance_one();
    }
  }

  label num_blocks() const { return fam.forest.local_num_octants(); }

  std::array<double, 3> aggregate_force_on_fluid(MPI_Comm comm,
                                                   bool legacy = false) const {
    std::array<double, 3> local{0.0, 0.0, 0.0};
    // Restrict the MEM sum to the cylinder surface: a kBouzidi/kSolid neighbour
    // only contributes when its physical centre lies in the cylinder obstacle
    // bbox. The duct walls (also kBouzidi) and the cube-minus-duct earth
    // (kSolid) are outside this bbox, so their (large) shear is excluded and
    // only the cylinder drag is measured -- the Schäfer-Turek Cd target.
    const BoundingBox bbox = cylinder_obstacle_bbox_;
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      const std::array<double, 3> org = lat.phys_origin();
      const double cw = lat.phys_cell_width();
      auto in_cylinder = [org, cw, bbox](int i, int j, int k) {
        const double x = org[0] + (i + 0.5) * cw;
        const double y = org[1] + (j + 0.5) * cw;
        const double z = org[2] + (k + 0.5) * cw;
        return x >= bbox.x_min && x <= bbox.x_max && y >= bbox.y_min &&
               y <= bbox.y_max && z >= bbox.z_min && z <= bbox.z_max;
      };
      const MomentumExchangeDrag drag(lat);
      const std::array<double, 3> f = legacy
                                         ? drag.force_on_fluid_if_legacy(in_cylinder)
                                         : drag.force_on_fluid_if(in_cylinder);
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
  // legacy (T11 timing diagnostic): use the legacy 2*live MEM (wrong timing) to
  // quantify the correction; the default uses the correct (snapshot+live) MEM.
  double drag_coefficient(MPI_Comm comm, double area,
                          std::array<double, 3> flow_dir,
                          bool legacy = false) const {
    const std::array<double, 3> f = aggregate_force_on_fluid(comm, legacy);
    const double fdot =
        f[0] * flow_dir[0] + f[1] * flow_dir[1] + f[2] * flow_dir[2];
    return 2.0 * (-fdot) / (cfg.rho0 * cfg.u_inlet * cfg.u_inlet * area);
  }

  // T11 drag-isolation diagnostic: count kBouzidi/kSolid boundary cells inside
  // vs outside the cylinder obstacle bbox, and report the MEM force both
  // filtered (cylinder-only, used for Cd) and unfiltered (all boundaries, the
  // pre-fix value). A large unfiltered/filtered ratio means the duct walls
  // still dominate the sum; a near-1 ratio means the bbox already isolates the
  // cylinder. Also reports the per-rank max |u| cell so a frozen field is
  // obvious. Printed once per call; caller gates with OCTLB_CYL_DRAG_DEBUG.
  void drag_breakdown(MPI_Comm comm, std::ostream& os) const {
    const BoundingBox bbox = cylinder_obstacle_bbox_;
    long long nb_in = 0, nb_out = 0, ns_in = 0, ns_out = 0;
    std::array<double, 3> f_all{0.0, 0.0, 0.0};
    std::array<double, 3> f_cyl{0.0, 0.0, 0.0};
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      const std::array<double, 3> org = lat.phys_origin();
      const double cw = lat.phys_cell_width();
      auto phys = [&](int i, int j, int k, double& x, double& y, double& z) {
        x = org[0] + (i + 0.5) * cw;
        y = org[1] + (j + 0.5) * cw;
        z = org[2] + (k + 0.5) * cw;
      };
      auto in_bbox = [&](double x, double y, double z) {
        return x >= bbox.x_min && x <= bbox.x_max && y >= bbox.y_min &&
               y <= bbox.y_max && z >= bbox.z_min && z <= bbox.z_max;
      };
      for (int k = 0; k < cfg.n; ++k) {
        for (int j = 0; j < cfg.n; ++j) {
          for (int i = 0; i < cfg.n; ++i) {
            const BcKind bk = lat.bc_kind(i, j, k);
            if (bk == BcKind::kBulk) continue;
            double x, y, z;
            phys(i, j, k, x, y, z);
            const bool in = in_bbox(x, y, z);
            if (bk == BcKind::kBouzidi) {
              if (in) ++nb_in; else ++nb_out;
            } else if (bk == BcKind::kSolid) {
              if (in) ++ns_in; else ++ns_out;
            }
          }
        }
      }
      // Unfiltered (all kSolid/kBouzidi) and filtered (cylinder bbox) MEM force.
      MomentumExchangeDrag drag(lat);
      const std::array<double, 3> fa = drag.force_on_fluid();
      const std::array<double, 3> fc =
          drag.force_on_fluid_if([&](int i, int j, int k) {
            double x, y, z;
            phys(i, j, k, x, y, z);
            return in_bbox(x, y, z);
          });
      f_all[0] += fa[0]; f_all[1] += fa[1]; f_all[2] += fa[2];
      f_cyl[0] += fc[0]; f_cyl[1] += fc[1]; f_cyl[2] += fc[2];
    }
    long long g[4] = {0, 0, 0, 0};
    MPI_Allreduce(&nb_in, &g[0], 1, MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&nb_out, &g[1], 1, MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&ns_in, &g[2], 1, MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&ns_out, &g[3], 1, MPI_LONG_LONG, MPI_SUM, comm);
    double ga[3], gc[3];
    MPI_Allreduce(f_all.data(), ga, 3, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(f_cyl.data(), gc, 3, MPI_DOUBLE, MPI_SUM, comm);
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      const double area = (0.1 / (kCubeSide / std::pow(2.0, cfg.max_level) /
                                  static_cast<double>(cfg.n))) *
                          ((kChannelHi - kChannelLo) /
                           (kCubeSide / std::pow(2.0, cfg.max_level) /
                            static_cast<double>(cfg.n)));
      const double denom = 0.5 * cfg.rho0 * cfg.u_inlet * cfg.u_inlet * area;
      os << "[drag] Bouzidi in_bbox=" << g[0] << " out_bbox=" << g[1]
         << " | kSolid in_bbox=" << g[2] << " out_bbox=" << g[3]
         << "\n[drag] F_all=(" << ga[0] << "," << ga[1] << "," << ga[2]
         << ") Cd_all=" << (-ga[0] / denom)
         << "\n[drag] F_cyl=(" << gc[0] << "," << gc[1] << "," << gc[2]
         << ") Cd_cyl=" << (-gc[0] / denom)
         << " area=" << area << " bbox_x[" << bbox.x_min << "," << bbox.x_max
         << "] y[" << bbox.y_min << "," << bbox.y_max << "] z[" << bbox.z_min
         << "," << bbox.z_max << "]\n" << std::flush;
    }
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

  // T11 diagnostic: mean streamwise (ux) + max |u| over bulk fluid cells, and
  // the net x-force on the fluid (sign sanity). Used to compare -n 1 vs -n 4.
  struct FlowStats { double mean_ux = 0.0; double mean_rho = 0.0; double max_umag = 0.0; long long n = 0; };
  FlowStats flow_stats(MPI_Comm comm) const {
    double local_sum_u = 0.0;
    double local_sum_rho = 0.0;
    double local_max = 0.0;
    long long local_n = 0;
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      for (int iz = 0; iz < cfg.n; ++iz) {
        for (int iy = 0; iy < cfg.n; ++iy) {
          for (int ix = 0; ix < cfg.n; ++ix) {
            if (lat.bc_kind(ix, iy, iz) != BcKind::kBulk) continue;
            double rho = 0.0;
            double u[3] = {};
            lat.get(ix, iy, iz).computeRhoU(rho, u);
            local_sum_u += u[0];
            local_sum_rho += rho;
            const double m = std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
            if (std::isfinite(m)) local_max = std::max(local_max, m);
            ++local_n;
          }
        }
      }
    }
    FlowStats s;
    double gx = 0.0, grho = 0.0;
    long long gn = 0;
    MPI_Allreduce(&local_sum_u, &gx, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&local_sum_rho, &grho, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&local_n, &gn, 1, MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&local_max, &s.max_umag, 1, MPI_DOUBLE, MPI_MAX, comm);
    s.mean_ux = gn > 0 ? gx / static_cast<double>(gn) : 0.0;
    s.mean_rho = gn > 0 ? grho / static_cast<double>(gn) : 0.0;
    s.n = gn;
    return s;
  }

  // T11 mass-conservation diagnostic: cross-sectional x-flux = sum(rho*u_x) over
  // the inlet-adjacent bulk slice (the first interior x-plane of duct fluid).
  // For an inlet-driven incompressible-ish flow this should match the prescribed
  // inlet mass flux rho0*u_inlet*inlet_area up to O(Mach^2) compressibility. A much
  // smaller flux signals the inlet is under-injecting or mass is leaking.
  double inlet_adjacent_mass_flux(MPI_Comm comm) const {
    double local = 0.0;
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      const std::array<double, 3> org = lat.phys_origin();
      const double cw = lat.phys_cell_width();
      for (int k = 0; k < cfg.n; ++k) {
        for (int j = 0; j < cfg.n; ++j) {
          for (int i = 0; i < cfg.n; ++i) {
            if (lat.bc_kind(i, j, k) != BcKind::kBulk) continue;
            const double x = org[0] + (i + 0.5) * cw;
            // first interior plane: just downstream of the x=0 inlet face
            if (x > 2.0 * cw) continue;
            double rho = 0.0;
            double u[3] = {};
            lat.get(i, j, k).computeRhoU(rho, u);
            local += rho * u[0];
          }
        }
      }
    }
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, comm);
    return global;
  }

  // T11 diagnostic: locate the global max-|u| cell. Each rank finds its local
  // max over NON-solid cells (incl. BC/CF cells), then we Allreduce-locmax and
  // the owning rank prints (rank, local block id, ijk, bc_kind, rho, u, and the
  // block's world-space bounding box so we can tell if it sits on a partition /
  // coarse-fine boundary). Used to pin down the stuck high-velocity cell that
  // drives the cross-rank Cd sign flip.
  void report_max_umag_cell(MPI_Comm comm, std::ostream& os) const {
    double local_max = -1.0;
    int local_rank = my_rank_;
    OctantId local_block = 0;
    int lix = 0, liy = 0, liz = 0;
    int local_bc = -1;
    double local_rho = 0.0, local_ux = 0.0, local_uy = 0.0, local_uz = 0.0;
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const CylinderLattice& lat = blocks[id];
      for (int iz = 0; iz < cfg.n; ++iz) {
        for (int iy = 0; iy < cfg.n; ++iy) {
          for (int ix = 0; ix < cfg.n; ++ix) {
            if (lat.bc_kind(ix, iy, iz) == BcKind::kSolid) continue;
            double rho = 0.0;
            double u[3] = {};
            lat.get(ix, iy, iz).computeRhoU(rho, u);
            const double m = std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
            if (!std::isfinite(m)) continue;
            if (m > local_max) {
              local_max = m; local_block = id; lix = ix; liy = iy; liz = iz;
              local_bc = static_cast<int>(lat.bc_kind(ix, iy, iz));
              local_rho = rho; local_ux = u[0]; local_uy = u[1]; local_uz = u[2];
            }
          }
        }
      }
    }
    struct { double val; int rank; } in, out;
    in.val = local_max; in.rank = local_rank;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    if (out.rank == local_rank) {
      const BoundingBox bb = fam.forest.quadrant_bounds(local_block);
      os << "[maxu] rank=" << local_rank << " block=" << local_block
         << " ijk=(" << lix << "," << liy << "," << liz << ")"
         << " bc=" << local_bc << " rho=" << local_rho
         << " u=(" << local_ux << "," << local_uy << "," << local_uz << ")"
         << " umag=" << local_max
         << " box=[" << bb.x_min << "," << bb.y_min << "," << bb.z_min
         << "]-[" << bb.x_max << "," << bb.y_max << "," << bb.z_max << "]"
         << " lvl=" << fam.forest.quadrant_level(local_block) << std::endl;
    }
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

  // T11 audit: detect spurious FLUID classified OUTSIDE the channel duct. The
  // channel duct is y,z in [kChannelLo, kChannelHi]; the surrounding "earth"
  // (cube minus duct) must be kSolid. A kBulk cell whose physical center lies
  // outside the duct is a material-classification defect -- spurious fluid that
  // evolves unconstrained and corrupts the drag. Reports per-rank + global
  // counts and the max |u| among the spurious cells. Guarded by env so it does
  // not run in the tests.
  void audit_spurious_fluid(MPI_Comm comm, std::ostream& os) const {
    long long local_bulk = 0;       // kBulk cells inside duct (legitimate)
    long long local_spurious = 0;   // kBulk cells outside duct (defect)
    long long local_solid_out = 0;  // kSolid cells outside duct (correct)
    double local_spur_maxu = 0.0;
    double worst_y = 0.0, worst_z = 0.0;
    for (OctantId id = 0; id < fam.forest.local_num_octants(); ++id) {
      const BoundingBox oct = fam.forest.quadrant_bounds(id);
      const double dx = (oct.x_max - oct.x_min) / cfg.n;
      const double dy = (oct.y_max - oct.y_min) / cfg.n;
      const double dz = (oct.z_max - oct.z_min) / cfg.n;
      const CylinderLattice& lat = blocks[id];
      for (int k = 0; k < cfg.n; ++k)
        for (int j = 0; j < cfg.n; ++j)
          for (int i = 0; i < cfg.n; ++i) {
            const double cy = oct.y_min + (j + 0.5) * dy;
            const double cz = oct.z_min + (k + 0.5) * dz;
            const bool in_duct =
                cy >= kChannelLo && cy <= kChannelHi &&
                cz >= kChannelLo && cz <= kChannelHi;
            const BcKind bk = lat.bc_kind(i, j, k);
            if (bk == BcKind::kBulk) {
              if (in_duct) {
                ++local_bulk;
              } else {
                ++local_spurious;
                double rho = 0.0, u[3] = {};
                lat.get(i, j, k).computeRhoU(rho, u);
                const double m =
                    std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
                if (std::isfinite(m) && m > local_spur_maxu) {
                  local_spur_maxu = m;
                  worst_y = cy;
                  worst_z = cz;
                }
              }
            } else if (bk == BcKind::kSolid && !in_duct) {
              ++local_solid_out;
            }
          }
    }
    long long g_bulk = 0, g_spur = 0, g_solid_out = 0;
    double g_spur_maxu = 0.0;
    MPI_Allreduce(&local_bulk, &g_bulk, 1, MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&local_spurious, &g_spur, 1, MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&local_solid_out, &g_solid_out, 1, MPI_LONG_LONG, MPI_SUM,
                  comm);
    MPI_Allreduce(&local_spur_maxu, &g_spur_maxu, 1, MPI_DOUBLE, MPI_MAX, comm);
    if (my_rank_ == 0) {
      os << "[audit] bulk_in_duct=" << g_bulk
         << " solid_out_duct=" << g_solid_out
         << " SPURIOUS_bulk_out_duct=" << g_spur
         << " spurious_maxu=" << g_spur_maxu;
      if (g_spur > 0) os << " (worst y=" << worst_y << " z=" << worst_z << ")";
      os << "\n";
    }
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