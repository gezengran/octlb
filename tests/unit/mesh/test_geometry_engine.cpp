#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/geometry/cgal_surface_mesh.h"
#include "src/mesh/geometry/geometry_build_error.h"
#include "src/mesh/geometry/geometry_engine.h"
#include "src/mesh/topology/face_pair_list.h"
#include "tests/unit/mesh/geometry_fixtures.h"

namespace octlb {
namespace {

BoundingBox UnitCubeDomain() { return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}; }

OctantId FindFinestOctantContaining(const OctreeForest& forest, scalar x,
                                    scalar y, scalar z) {
  OctantId best = 0;
  int best_level = -1;
  for (label id = 0; id < forest.local_num_octants(); ++id) {
    const BoundingBox b = forest.quadrant_bounds(id);
    if (x >= b.x_min && x <= b.x_max && y >= b.y_min && y <= b.y_max &&
        z >= b.z_min && z <= b.z_max) {
      const int lvl = forest.quadrant_level(id);
      if (lvl > best_level) {
        best_level = lvl;
        best = id;
      }
    }
  }
  return best;
}

MaterialKind MaterialAtPoint(const OctreeForest& forest,
                             const MaterialField& field, scalar x, scalar y,
                             scalar z) {
  const OctantId id = FindFinestOctantContaining(forest, x, y, z);
  const BoundingBox b = forest.quadrant_bounds(id);
  const int n = field.nx();
  const scalar fx = (x - b.x_min) / (b.x_max - b.x_min) * static_cast<scalar>(n);
  const scalar fy = (y - b.y_min) / (b.y_max - b.y_min) * static_cast<scalar>(n);
  const scalar fz = (z - b.z_min) / (b.z_max - b.z_min) * static_cast<scalar>(n);
  const int i =
      std::clamp(static_cast<int>(fx), 0, n - 1);
  const int j =
      std::clamp(static_cast<int>(fy), 0, n - 1);
  const int k =
      std::clamp(static_cast<int>(fz), 0, n - 1);
  return field.at(id, i, j, k);
}

bool FieldContainsKind(const MaterialField& field, MaterialKind want) {
  for (label id = 0; id < field.num_octants(); ++id) {
    for (int k = 0; k < field.nz(); ++k) {
      for (int j = 0; j < field.ny(); ++j) {
        for (int i = 0; i < field.nx(); ++i) {
          if (field.at(id, i, j, k) == want) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

GeometryAssembly SolidObstacleAssembly() {
  GeometryAssembly assembly;
  GeometryPart part;
  part.name = "cube";
  part.role = GeometryPartRole::kExternalObstacle;
  part.priority = 0;
  part.soup = test_geom::SolidCube(0.4, 0.4, 0.4, 0.6, 0.6, 0.6);
  assembly.parts.push_back(part);
  return assembly;
}

GeometryAssembly ChannelAssembly() {
  GeometryAssembly assembly;
  GeometryPart part;
  part.name = "duct";
  part.role = GeometryPartRole::kInternalChannel;
  part.priority = 0;
  test_geom::FillHollowChannelPart(&part, 0.2, 0.8, 0.3, 0.7);
  assembly.parts.push_back(part);
  return assembly;
}

GeometryAssembly WindTunnelAssembly(int model_priority) {
  GeometryAssembly assembly;
  GeometryPart channel;
  channel.name = "tunnel";
  channel.role = GeometryPartRole::kInternalChannel;
  channel.priority = 0;
  test_geom::FillHollowChannelPart(&channel, 0.15, 0.85, 0.25, 0.75);
  GeometryPart model;
  model.name = "sphere";
  model.role = GeometryPartRole::kExternalObstacle;
  model.priority = model_priority;
  model.soup = test_geom::CenteredSphere(0.5, 0.5, 0.5, 0.12, 8);
  assembly.parts.push_back(channel);
  assembly.parts.push_back(model);
  return assembly;
}

TEST(GeometryEngineTest, CubeSoup_CenterIsInside) {
  const CgalSurfaceMesh mesh =
      CgalSurfaceMesh::from_soup(test_geom::SolidCube(0.4, 0.4, 0.4, 0.6, 0.6, 0.6));
  EXPECT_TRUE(mesh.is_inside(0.5, 0.5, 0.5));
  EXPECT_FALSE(mesh.is_inside(0.1, 0.1, 0.1));
}

TEST(GeometryEngineTest, S1_ExternalObstacle_SolidFluidBoundary) {
  MPI_Comm comm = MPI_COMM_WORLD;
  OctreeForest forest(comm, UnitCubeDomain());
  GeometryConfig cfg = test_geom::DefaultConfig();
  GeometryEngine engine;
  const MaterialField field =
      engine.build(forest, SolidObstacleAssembly(), cfg);

  EXPECT_TRUE(FieldContainsKind(field, MaterialKind::kSolid));
  EXPECT_TRUE(FieldContainsKind(field, MaterialKind::kBoundary));
  EXPECT_EQ(MaterialAtPoint(forest, field, 0.05, 0.05, 0.05),
            MaterialKind::kFluid);
  const MaterialKind center = MaterialAtPoint(forest, field, 0.5, 0.5, 0.5);
  EXPECT_TRUE(center == MaterialKind::kSolid ||
              center == MaterialKind::kBoundary);
}

TEST(GeometryEngineTest, S3_InternalChannel_LumenWallBoundary) {
  MPI_Comm comm = MPI_COMM_WORLD;
  OctreeForest forest(comm, UnitCubeDomain());
  GeometryConfig cfg = test_geom::DefaultConfig();
  GeometryEngine engine;
  const MaterialField field = engine.build(forest, ChannelAssembly(), cfg);

  EXPECT_EQ(MaterialAtPoint(forest, field, 0.5, 0.5, 0.5), MaterialKind::kFluid);
  EXPECT_TRUE(FieldContainsKind(field, MaterialKind::kSolid));
  EXPECT_TRUE(FieldContainsKind(field, MaterialKind::kBoundary));
}

TEST(GeometryEngineTest, S5_WindTunnel_SphereSolidLumenFluid) {
  MPI_Comm comm = MPI_COMM_WORLD;
  OctreeForest forest(comm, UnitCubeDomain());
  GeometryConfig cfg = test_geom::DefaultConfig();
  GeometryEngine engine;
  const MaterialField field =
      engine.build(forest, WindTunnelAssembly(/*model_priority=*/10), cfg);

  EXPECT_EQ(MaterialAtPoint(forest, field, 0.5, 0.5, 0.5), MaterialKind::kSolid);
  EXPECT_TRUE(FieldContainsKind(field, MaterialKind::kFluid));
}

TEST(GeometryEngineTest, S6_PriorityOverlay_SphereDominatesObstacle) {
  MPI_Comm comm = MPI_COMM_WORLD;
  GeometryConfig cfg = test_geom::DefaultConfig();
  GeometryEngine engine;

  OctreeForest forest_high(comm, UnitCubeDomain());
  const MaterialField high =
      engine.build(forest_high, WindTunnelAssembly(10), cfg);


  const MaterialKind kh =
      MaterialAtPoint(forest_high, high, 0.5, 0.5, 0.5);
  GeometryAssembly low_pri = WindTunnelAssembly(0);
  low_pri.parts[0].priority = 10;  // channel wins over model at equal geometry
  low_pri.parts[1].priority = 0;
  OctreeForest forest_low2(MPI_COMM_WORLD, UnitCubeDomain());
  const MaterialField low =
      engine.build(forest_low2, low_pri, cfg);
  EXPECT_EQ(kh, MaterialKind::kSolid);
  EXPECT_EQ(MaterialAtPoint(forest_low2, low, 0.5, 0.5, 0.5),
            MaterialKind::kFluid);
}

TEST(GeometryEngineTest, S7_IllegalInternalChannel_SolidBoxFails) {
  MPI_Comm comm = MPI_COMM_WORLD;
  OctreeForest forest(comm, UnitCubeDomain());
  GeometryAssembly assembly;
  GeometryPart part;
  part.name = "solid_wall";
  part.role = GeometryPartRole::kInternalChannel;
  part.soup = test_geom::SolidCube(0.3, 0.3, 0.3, 0.7, 0.7, 0.7);
  assembly.parts.push_back(part);

  GeometryEngine engine;
  EXPECT_THROW(engine.build(forest, assembly, test_geom::DefaultConfig()),
               GeometryBuildError);
}

// W3-d (T11) option A: a kInternalChannel whose OUTER surface is the cubic
// domain itself and whose INNER cavity is the channel duct carves the channel
// as fluid and the cube-minus-channel region as solid (in_wall) -- so an
// anisotropic channel fits inside a cubic domain with isotropic cells. The
// cube-minus-channel point must be kSolid (not kFluid, which the "else" branch
// would give if the outer were smaller than the domain).
TEST(GeometryEngineTest, InternalChannel_CubicOuter_CarvesChannelSolidExterior) {
  MPI_Comm comm = MPI_COMM_WORLD;
  OctreeForest forest(comm, UnitCubeDomain());
  GeometryAssembly assembly;
  GeometryPart part;
  part.name = "channel";
  part.role = GeometryPartRole::kInternalChannel;
  part.priority = 0;
  // outer = the whole domain cube [0,1]^3; inner = the channel [0.3,0.7]^3.
  test_geom::FillHollowChannelPart(&part, /*outer_min=*/0.0, /*outer_max=*/1.0,
                                   /*inner_min=*/0.3, /*inner_max=*/0.7);
  assembly.parts.push_back(part);

  GeometryConfig cfg = test_geom::DefaultConfig();
  GeometryEngine engine;
  const MaterialField field = engine.build(forest, assembly, cfg);

  // Lumen (inside the inner cavity) is fluid.
  EXPECT_EQ(MaterialAtPoint(forest, field, 0.5, 0.5, 0.5), MaterialKind::kFluid);
  // Cube-minus-channel (between inner and outer = the wall) is solid.
  EXPECT_EQ(MaterialAtPoint(forest, field, 0.1, 0.5, 0.5), MaterialKind::kSolid)
      << "cube-minus-channel region must be solid (in_wall)";
  EXPECT_TRUE(FieldContainsKind(field, MaterialKind::kBoundary));
}

TEST(GeometryEngineTest, Build_ThenFacePairList_NonEmpty) {
  MPI_Comm comm = MPI_COMM_WORLD;
  OctreeForest forest(comm, UnitCubeDomain());
  GeometryEngine engine;
  (void)engine.build(forest, WindTunnelAssembly(10), test_geom::DefaultConfig());
  const FacePairList pairs(forest);
  EXPECT_GT(pairs.same_level_faces().size(), 0u);
  EXPECT_GT(pairs.coarse_fine_faces().size(), 0u);
}

}  // namespace
}  // namespace octlb
