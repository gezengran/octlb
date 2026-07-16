#include <gtest/gtest.h>
#include <mpi.h>

#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/geometry/geometry_engine.h"
#include "src/mesh/geometry/geometry_types.h"
#include "src/mesh/geometry/material_field.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/bc_installer.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/domain_boundary_handler.h"
#include "tests/unit/mesh/geometry_fixtures.h"

namespace octlb {
namespace {

using Lattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

BoundingBox UnitCubeDomain() { return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}; }

// A centered solid cube stands in for the Schäfer-Turek cylinder: the geometry
// only needs a kBoundary surface / kSolid interior so BcInstaller can be
// observed splitting it into kBouzidi / kSolid while the six domain faces take
// their per-spec BcKind.
GeometryAssembly CenteredObstacleAssembly() {
  GeometryAssembly assembly;
  GeometryPart part;
  part.name = "obstacle";
  part.role = GeometryPartRole::kExternalObstacle;
  part.priority = 0;
  part.soup = test_geom::SolidCube(0.3, 0.3, 0.3, 0.7, 0.7, 0.7);
  assembly.parts.push_back(std::move(part));
  return assembly;
}

// Schäfer-Turek face roles: velocity inlet, pressure outlet, no-slip walls.
std::vector<DomainBcSpec> SchaeferTurekSpecs(double u_inlet) {
  std::vector<DomainBcSpec> specs;
  DomainBcSpec inlet;
  inlet.face = FaceDir::kXMin;
  inlet.type = DomainBcType::kInterpolatedVelocity;
  inlet.u_wall = {u_inlet, 0.0, 0.0};
  specs.push_back(inlet);
  DomainBcSpec outlet;
  outlet.face = FaceDir::kXMax;
  outlet.type = DomainBcType::kInterpolatedPressure;
  outlet.rho_target = 1.0;
  specs.push_back(outlet);
  for (FaceDir wall : {FaceDir::kYMin, FaceDir::kYMax, FaceDir::kZMin,
                       FaceDir::kZMax}) {
    DomainBcSpec s;
    s.face = wall;
    s.type = DomainBcType::kInterpolatedVelocity;  // FD no-slip (u=0)
    specs.push_back(s);
  }
  return specs;
}

// DomainBcType -> BcKind mapping covers every spec type (the Schäfer-Turek
// face config above only exercises kInterpolatedVelocity/kInterpolatedPressure).
TEST(BcInstaller, BcKindFromSpecType_MapsAllTypes) {
  EXPECT_EQ(bc::BcKindFromSpecType(DomainBcType::kNoSlip), BcKind::kBounceBack);
  EXPECT_EQ(bc::BcKindFromSpecType(DomainBcType::kMovingLid),
            BcKind::kMovingBounceBack);
  EXPECT_EQ(bc::BcKindFromSpecType(DomainBcType::kInterpolatedVelocity),
            BcKind::kVelocityDirichlet);
  EXPECT_EQ(bc::BcKindFromSpecType(DomainBcType::kOutflow), BcKind::kOutflow);
  EXPECT_EQ(bc::BcKindFromSpecType(DomainBcType::kInterpolatedPressure),
            BcKind::kPressureDirichlet);
}

// R3: BcInstaller resolves a Schäfer-Turek-style setup into the right per-cell
// BcKind distribution -- inlet kVelocityDirichlet, outlet kPressureDirichlet,
// walls kBounceBack, obstacle surface kBouzidi, obstacle interior kSolid, and
// fluid interior kBulk -- with no global mode flag.
TEST(BcInstaller, SchaeferTurek_StampingCorrect) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();
  const FacePairList pairs(forest);

  GeometryConfig cfg = test_geom::DefaultConfig();
  cfg.cell_width = 8;
  cfg.max_level = 0;
  cfg.resolve_surface_times = 0;
  const MaterialField material =
      GeometryEngine{}.build(forest, CenteredObstacleAssembly(), cfg);
  const int nx = material.nx();
  const int ny = material.ny();
  const int nz = material.nz();

  BlockCollection<Lattice> blocks(material.num_octants(), [&](OctantId) {
    return Lattice(nx, ny, nz, 1);
  });
  for (label oid = 0; oid < material.num_octants(); ++oid) {
    blocks[static_cast<OctantId>(oid)].set_octant_id(
        static_cast<OctantId>(oid));
  }

  const double u_inlet = 0.05;
  const auto specs = SchaeferTurekSpecs(u_inlet);
  bc::StampFromMaterial(blocks, material);
  bc::StampTreeBoundaryCells(blocks, pairs.tree_boundary_faces(), specs, nx, ny,
                              nz);

  Lattice& lat = blocks[0];

  // Domain face cells (obstacle is centered, so the outer slabs are fluid):
  // inlet, outlet, and the four walls take their per-spec BcKind.
  EXPECT_EQ(lat.bc_kind(0, ny / 2, nz / 2), BcKind::kVelocityDirichlet)
      << "inlet (kXMin) face cell";
  EXPECT_EQ(lat.bc_kind(nx - 1, ny / 2, nz / 2), BcKind::kPressureDirichlet)
      << "outlet (kXMax) face cell";
  EXPECT_EQ(lat.bc_kind(nx / 2, 0, nz / 2), BcKind::kVelocityDirichlet)
      << "kYMin wall cell (FD no-slip)";
  EXPECT_EQ(lat.bc_kind(nx / 2, ny - 1, nz / 2), BcKind::kVelocityDirichlet)
      << "kYMax wall cell (FD no-slip)";
  EXPECT_EQ(lat.bc_kind(nx / 2, ny / 2, 0), BcKind::kVelocityDirichlet)
      << "kZMin wall cell (FD no-slip)";
  EXPECT_EQ(lat.bc_kind(nx / 2, ny / 2, nz - 1), BcKind::kVelocityDirichlet)
      << "kZMax wall cell (FD no-slip)";

  // Obstacle present: surface -> kBouzidi, interior -> kSolid; fluid bulk ->
  // kBulk. Count across the block so the assertion does not depend on the exact
  // voxelized surface cells.
  int n_bouzidi = 0;
  int n_solid = 0;
  int n_bulk = 0;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        const BcKind kind = lat.bc_kind(i, j, k);
        if (kind == BcKind::kBouzidi) ++n_bouzidi;
        else if (kind == BcKind::kSolid) ++n_solid;
        else if (kind == BcKind::kBulk) ++n_bulk;
      }
    }
  }
  EXPECT_GT(n_bouzidi, 0) << "obstacle surface must be stamped kBouzidi";
  EXPECT_GT(n_solid, 0) << "obstacle interior must be stamped kSolid";
  EXPECT_GT(n_bulk, 0) << "fluid interior must remain kBulk";
}

// R4: cavity3d migration. The single-block closed cavity has all six domain
// faces as tree boundaries; StampTreeBoundaryCells with six kInterpolatedVelocity
// specs stamps the whole outer shell kVelocityDirichlet and leaves the interior
// kBulk -- equivalent to the legacy geometric MarkDomainBoundaryBcKinds, so
// test_cavity3d_serial #12 (L2<2%) does not regress.
TEST(BcInstaller, Cavity_LegacyEquivalent) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();
  const FacePairList pairs(forest);

  constexpr int kN = 8;
  BlockCollection<Lattice> blocks(1, [](OctantId) {
    return Lattice(kN, kN, kN, 1);
  });

  std::vector<DomainBcSpec> specs;
  for (FaceDir face : {FaceDir::kXMin, FaceDir::kXMax, FaceDir::kYMin,
                       FaceDir::kYMax, FaceDir::kZMin, FaceDir::kZMax}) {
    DomainBcSpec s;
    s.face = face;
    s.type = DomainBcType::kInterpolatedVelocity;
    specs.push_back(s);
  }
  bc::StampTreeBoundaryCells(blocks, pairs.tree_boundary_faces(), specs, kN,
                              kN, kN);

  Lattice& lat = blocks[0];
  for (int i = 0; i < kN; ++i) {
    for (int j = 0; j < kN; ++j) {
      for (int k = 0; k < kN; ++k) {
        const bool on_shell = i == 0 || i == kN - 1 || j == 0 ||
                              j == kN - 1 || k == 0 || k == kN - 1;
        if (on_shell) {
          EXPECT_EQ(lat.bc_kind(i, j, k), BcKind::kVelocityDirichlet)
              << "cavity shell cell (" << i << ',' << j << ',' << k << ')';
        } else {
          EXPECT_EQ(lat.bc_kind(i, j, k), BcKind::kBulk)
              << "cavity interior cell (" << i << ',' << j << ',' << k << ')';
        }
      }
    }
  }
}

}  // namespace
}  // namespace octlb