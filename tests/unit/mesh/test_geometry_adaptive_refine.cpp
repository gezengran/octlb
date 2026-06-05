#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/geometry/geometry_adaptive_refinement.h"
#include "src/mesh/geometry/geometry_config.h"
#include "src/mesh/geometry/geometry_assembly.h"
#include "src/mesh/geometry/geometry_engine.h"
#include "tests/unit/mesh/geometry_fixtures.h"

namespace octlb {
namespace {

BoundingBox UnitCubeDomain() { return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}; }

int MaxLevelNear(const OctreeForest& forest, const BoundingBox& region) {
  int max_lvl = 0;
  for (label id = 0; id < forest.local_num_octants(); ++id) {
    const BoundingBox b = forest.quadrant_bounds(id);
    if (!(b.x_max < region.x_min || b.x_min > region.x_max ||
          b.y_max < region.y_min || b.y_min > region.y_max ||
          b.z_max < region.z_min || b.z_min > region.z_max)) {
      max_lvl = std::max(max_lvl, forest.quadrant_level(id));
    }
  }
  return max_lvl;
}

int GlobalMaxLevel(const OctreeForest& forest, MPI_Comm comm) {
  int local = 0;
  for (label id = 0; id < forest.local_num_octants(); ++id) {
    local = std::max(local, forest.quadrant_level(id));
  }
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, comm);
  return global;
}

TEST(GeometryAdaptiveRefineTest, SurfaceRefinement_IncreasesLevelNearGeometry) {
  MPI_Comm comm = MPI_COMM_WORLD;
  OctreeForest forest(comm, UnitCubeDomain());
  TriangleSoup soup =
      test_geom::SolidCube(0.45, 0.45, 0.45, 0.55, 0.55, 0.55);
  GeometryConfig cfg;
  cfg.max_level = 5;
  cfg.resolve_surface_times = 4;
  cfg.bound_width = 0.0;

  resolve_surface(forest, soup, cfg);

  const BoundingBox near_geom{0.44, 0.44, 0.44, 0.56, 0.56, 0.56};
  const BoundingBox far_corner{0.01, 0.01, 0.01, 0.12, 0.12, 0.12};
  const int near_level = MaxLevelNear(forest, near_geom);
  const int far_level = MaxLevelNear(forest, far_corner);
  EXPECT_GT(near_level, far_level);
  EXPECT_GE(GlobalMaxLevel(forest, comm), 2);
}

TEST(GeometryAdaptiveRefineTest, ResolveBounding_RefinesInsideExtendedBBox) {
  MPI_Comm comm = MPI_COMM_WORLD;
  OctreeForest forest(comm, UnitCubeDomain());
  TriangleSoup soup =
      test_geom::SolidCube(0.4, 0.4, 0.4, 0.6, 0.6, 0.6);
  GeometryConfig cfg;
  cfg.max_level = 4;
  cfg.bound_width = 0.05;
  const BoundingBox geom_bbox = extended_geometry_bbox(soup, cfg);
  resolve_bounding(forest, geom_bbox, cfg);
  forest.balance();

  const BoundingBox inside{0.35, 0.35, 0.35, 0.65, 0.65, 0.65};
  const BoundingBox corner{0.0, 0.0, 0.0, 0.12, 0.12, 0.12};
  EXPECT_GT(MaxLevelNear(forest, inside), MaxLevelNear(forest, corner));
}

TEST(GeometryAdaptiveRefineTest, TwoRank_RefineCollective) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    GTEST_SKIP() << "Requires at least 2 MPI ranks";
  }

  MPI_Comm comm = MPI_COMM_WORLD;
  OctreeForest forest(comm, UnitCubeDomain());
  GeometryAssembly assembly;
  GeometryPart part;
  part.soup = test_geom::SolidCube(0.3, 0.3, 0.3, 0.7, 0.7, 0.7);
  part.role = GeometryPartRole::kExternalObstacle;
  assembly.parts.push_back(part);

  GeometryConfig cfg = test_geom::DefaultConfig();
  GeometryEngine engine;
  const MaterialField field = engine.build(forest, assembly, cfg);
  EXPECT_GT(forest.local_num_octants(), 0);
  EXPECT_EQ(field.num_octants(), forest.local_num_octants());
}

}  // namespace
}  // namespace octlb
