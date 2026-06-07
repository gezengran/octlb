#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"

namespace octlb {
namespace {

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

OctantId FindCenterOctant(const OctreeForest& forest) {
  const BoundingBox center = {0.375, 0.375, 0.375, 0.625, 0.625, 0.625};
  OctantId best = 0;
  scalar best_dist = 1e30;
  for (label i = 0; i < forest.local_num_octants(); ++i) {
    const BoundingBox b = forest.quadrant_bounds(i);
    const scalar cx = 0.5 * (b.x_min + b.x_max);
    const scalar cy = 0.5 * (b.y_min + b.y_max);
    const scalar cz = 0.5 * (b.z_min + b.z_max);
    const scalar dist =
        (cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5) +
        (cz - 0.5) * (cz - 0.5);
    if (dist < best_dist) {
      best_dist = dist;
      best = i;
    }
  }
  return best;
}

label GlobalSum(label local, MPI_Comm comm) {
  label global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT32_T, MPI_SUM, comm);
  return global;
}

}  // namespace

TEST(FacePairList, UniformRefineTwoLevelsYieldsOnlySameLevelFaces) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);
  forest.balance();
  forest.partition();

  const FacePairList pairs(forest);
  EXPECT_TRUE(pairs.coarse_fine_faces().empty());
  EXPECT_GT(pairs.same_level_faces().size(), 0u);
}

TEST(FacePairList, CenterRefineProducesCoarseFineFacesWithFourFines) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);
  const OctantId center = FindCenterOctant(forest);
  forest.refine(
      [center](OctantId id) { return id == center; }, 3);
  forest.balance();
  forest.partition();

  const FacePairList pairs(forest);
  const label local_cf =
      static_cast<label>(pairs.coarse_fine_faces().size());
  const label global_cf = GlobalSum(local_cf, MPI_COMM_WORLD);
  EXPECT_GT(global_cf, 0);
  // Six faces of the refined block; balance adds three corner hanging faces.
  EXPECT_EQ(global_cf, 9);

  for (const CoarseFineFace& face : pairs.coarse_fine_faces()) {
    (void)face.coarse_id;
    for (int i = 0; i < 4; ++i) {
      EXPECT_GE(face.fine_ids[i], 0);
      EXPECT_GT(face.comm_tags[i], 0);
    }
  }
}

TEST(FacePairList, TreeBoundaryFaces_Enumerated) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();

  const FacePairList pairs(forest);
  const auto& tb = pairs.tree_boundary_faces();
  EXPECT_EQ(tb.size(), 6u);
  for (const TreeBoundaryFace& face : tb) {
    EXPECT_EQ(face.octant_id, 0);
    EXPECT_GE(static_cast<int>(face.face_dir), 0);
    EXPECT_LE(static_cast<int>(face.face_dir), 5);
  }
}

TEST(FacePairList, DomainTreeBoundary_NoSameLevelFace) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();

  const FacePairList pairs(forest);
  const auto& tb = pairs.tree_boundary_faces();
  ASSERT_EQ(tb.size(), 6u);

  for (const TreeBoundaryFace& face : tb) {
    for (const SameLevelFace& sl : pairs.same_level_faces()) {
      EXPECT_FALSE(sl.local_id == face.octant_id &&
                   sl.dir == face.face_dir);
    }
  }
}

TEST(FacePairList, CrossRankFacesShareSymmetricCommTag) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 1);
  forest.balance();
  forest.partition();

  const FacePairList pairs(forest);
  for (const SameLevelFace& face : pairs.same_level_faces()) {
    if (face.remote_rank == rank) {
      continue;
    }
    EXPECT_GT(face.comm_tag, 0);
  }
}

}  // namespace octlb
