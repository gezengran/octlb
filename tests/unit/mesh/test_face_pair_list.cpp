#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <vector>

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

// Stronger symmetry guard (requires exactly 2 ranks): the two ranks sharing a
// partition boundary must assign the SAME multiset of comm_tags for their
// cross-rank same-level faces (per-rank-pair sorted-index tags are symmetric
// and unique). Run via the test_face_pair_list_symmetry_mpi2 ctest entry.
TEST(FacePairList, CrossRankTagsAgreeBetweenRanks) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "requires exactly 2 ranks";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);  // 64 octants -> boundary
  forest.balance();
  forest.partition();

  const FacePairList pairs(forest);
  std::vector<int> tags;
  for (const SameLevelFace& face : pairs.same_level_faces()) {
    if (face.remote_rank != rank) {  // cross-rank
      tags.push_back(face.comm_tag);
    }
  }
  std::sort(tags.begin(), tags.end());
  ASSERT_FALSE(tags.empty()) << "no cross-rank faces; partition too coarse";

  // Per-pair uniqueness: no duplicate comm_tag on this rank.
  for (std::size_t i = 1; i < tags.size(); ++i) {
    ASSERT_NE(tags[i], tags[i - 1]) << "duplicate cross-rank comm_tag on rank "
                                    << rank;
  }

  // Symmetry: the peer rank assigns the same multiset of tags for the same
  // shared faces. Exchange counts, then the sorted tag vectors.
  const int peer = 1 - rank;
  int my_count = static_cast<int>(tags.size());
  int peer_count = 0;
  MPI_Sendrecv(&my_count, 1, MPI_INT, peer, 0, &peer_count, 1, MPI_INT, peer, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  ASSERT_EQ(my_count, peer_count)
      << "cross-rank face count differs between ranks";

  std::vector<int> other(static_cast<std::size_t>(peer_count), 0);
  MPI_Sendrecv(tags.data(), my_count, MPI_INT, peer, 1, other.data(),
               peer_count, MPI_INT, peer, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::sort(other.begin(), other.end());
  EXPECT_EQ(tags, other)
      << "cross-rank comm_tag multiset differs between ranks";
}

// ② W3 Stage B (T11): cross-rank edge-diagonal pairs are enumerated by the
// p8est edge callback and given symmetric comm_tags (per-rank-pair canonical
// edge key), so GhostSchedule Stage B can Isend/Irecv edge halos without tag
// mismatch. Requires exactly 2 ranks; uniform refine to a single level keeps
// the forest edge-balanced (= face-balanced, no 2:1) so the edge callback fires.
TEST(FacePairList, CrossRankEdges_EnumeratedSymmetric) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "requires exactly 2 ranks";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);  // 64 octants -> boundary
  forest.balance();
  forest.partition();

  const FacePairList pairs(forest);

  // Cross-rank edges exist (diagonal blocks across the partition boundary).
  std::vector<int> tags;
  for (const CrossRankEdge& e : pairs.cross_rank_edges()) {
    EXPECT_NE(e.remote_rank, rank) << "cross_rank_edges must be cross-rank only";
    EXPECT_GE(static_cast<int>(e.d1), 0);
    EXPECT_LE(static_cast<int>(e.d1), 5);
    EXPECT_GE(static_cast<int>(e.d2), 0);
    EXPECT_LE(static_cast<int>(e.d2), 5);
    EXPECT_NE(static_cast<int>(e.d1) / 2, static_cast<int>(e.d2) / 2)
        << "edge faces must be on different axes";
    tags.push_back(e.comm_tag);
  }
  ASSERT_FALSE(tags.empty())
      << "no cross-rank edges; partition too coarse or edge callback did not fire";

  // Per-peer uniqueness: no duplicate comm_tag on this rank.
  std::sort(tags.begin(), tags.end());
  for (std::size_t i = 1; i < tags.size(); ++i) {
    ASSERT_NE(tags[i], tags[i - 1])
        << "duplicate cross-rank edge comm_tag on rank " << rank;
  }

  // Symmetry: the peer rank assigns the same multiset of edge tags for the
  // shared edges (canonical edge key is symmetric).
  const int peer = 1 - rank;
  int my_count = static_cast<int>(tags.size());
  int peer_count = 0;
  MPI_Sendrecv(&my_count, 1, MPI_INT, peer, 0, &peer_count, 1, MPI_INT, peer, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  ASSERT_EQ(my_count, peer_count)
      << "cross-rank edge count differs between ranks";

  std::vector<int> other(static_cast<std::size_t>(peer_count), 0);
  MPI_Sendrecv(tags.data(), my_count, MPI_INT, peer, 1, other.data(),
               peer_count, MPI_INT, peer, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::sort(other.begin(), other.end());
  EXPECT_EQ(tags, other) << "cross-rank edge comm_tag multiset differs between ranks";
}

}  // namespace octlb
