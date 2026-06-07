#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"

namespace octlb {
namespace {

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

label GlobalSum(label local, MPI_Comm comm) {
  label global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT32_T, MPI_SUM, comm);
  return global;
}

scalar BoxVolume(const BoundingBox& b) {
  return (b.x_max - b.x_min) * (b.y_max - b.y_min) * (b.z_max - b.z_min);
}

bool BoxesOverlap(const BoundingBox& a, const BoundingBox& b) {
  constexpr scalar kTol = 1e-12;
  const bool sep_x = a.x_max <= b.x_min + kTol || b.x_max <= a.x_min + kTol;
  const bool sep_y = a.y_max <= b.y_min + kTol || b.y_max <= a.y_min + kTol;
  const bool sep_z = a.z_max <= b.z_min + kTol || b.z_max <= a.z_min + kTol;
  return !(sep_x || sep_y || sep_z);
}

int FaceLevelGap(const BoundingBox& a, int level_a, const BoundingBox& b,
                 int level_b) {
  constexpr scalar kTol = 1e-10;
  const bool share_x =
      std::abs(a.x_max - b.x_min) < kTol || std::abs(b.x_max - a.x_min) < kTol;
  const bool share_y =
      std::abs(a.y_max - b.y_min) < kTol || std::abs(b.y_max - a.y_min) < kTol;
  const bool share_z =
      std::abs(a.z_max - b.z_min) < kTol || std::abs(b.z_max - a.z_min) < kTol;
  const bool overlap_y =
      a.y_min < b.y_max - kTol && b.y_min < a.y_max - kTol;
  const bool overlap_z =
      a.z_min < b.z_max - kTol && b.z_min < a.z_max - kTol;
  const bool overlap_x =
      a.x_min < b.x_max - kTol && b.x_min < a.x_max - kTol;
  const bool overlap_z2 =
      a.z_min < b.z_max - kTol && b.z_min < a.z_max - kTol;
  const bool overlap_y2 =
      a.y_min < b.y_max - kTol && b.y_min < a.y_max - kTol;

  const bool face_x = share_x && overlap_y && overlap_z;
  const bool face_y = share_y && overlap_x && overlap_z2;
  const bool face_z = share_z && overlap_x && overlap_y2;
  if (!face_x && !face_y && !face_z) {
    return 0;
  }
  return std::abs(level_a - level_b);
}

struct GlobalOctant {
  BoundingBox bounds;
  int level = 0;
};

void GatherGlobalOctants(const OctreeForest& forest, MPI_Comm comm,
                         std::vector<GlobalOctant>* out) {
  const label local_n = forest.local_num_octants();
  std::vector<double> flat(static_cast<std::size_t>(local_n) * 7);
  for (label i = 0; i < local_n; ++i) {
    const BoundingBox b = forest.quadrant_bounds(i);
    const std::size_t o = static_cast<std::size_t>(i) * 7;
    flat[o + 0] = b.x_min;
    flat[o + 1] = b.y_min;
    flat[o + 2] = b.z_min;
    flat[o + 3] = b.x_max;
    flat[o + 4] = b.y_max;
    flat[o + 5] = b.z_max;
    flat[o + 6] = static_cast<double>(forest.quadrant_level(i));
  }

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  std::vector<int> counts(size, 0);
  MPI_Allgather(&local_n, 1, MPI_INT32_T, counts.data(), 1, MPI_INT32_T, comm);

  std::vector<int> displs(size, 0);
  int total = 0;
  for (int r = 0; r < size; ++r) {
    displs[r] = total;
    total += counts[r];
  }

  std::vector<double> gathered(static_cast<std::size_t>(total) * 7);
  std::vector<int> recv_counts(size);
  std::vector<int> recv_displs(size);
  for (int r = 0; r < size; ++r) {
    recv_counts[r] = counts[r] * 7;
    recv_displs[r] = displs[r] * 7;
  }

  MPI_Allgatherv(flat.data(), local_n * 7, MPI_DOUBLE, gathered.data(),
                 recv_counts.data(), recv_displs.data(), MPI_DOUBLE, comm);

  out->clear();
  out->reserve(static_cast<std::size_t>(total));
  for (int i = 0; i < total; ++i) {
    const std::size_t o = static_cast<std::size_t>(i) * 7;
    GlobalOctant g;
    g.bounds = {gathered[o + 0], gathered[o + 1], gathered[o + 2],
                gathered[o + 3], gathered[o + 4], gathered[o + 5]};
    g.level = static_cast<int>(gathered[o + 6]);
    out->push_back(g);
  }
}

}  // namespace

TEST(OctreeForest, Brick2x1x1_HasTwoRootOctantsSplitAlongX) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain(), 2, 1, 1);
  const label global_n =
      GlobalSum(forest.local_num_octants(), MPI_COMM_WORLD);
  EXPECT_EQ(global_n, 2);

  std::vector<GlobalOctant> octants;
  GatherGlobalOctants(forest, MPI_COMM_WORLD, &octants);
  ASSERT_EQ(octants.size(), 2u);

  scalar x_min = octants[0].bounds.x_min;
  scalar x_max = octants[0].bounds.x_max;
  for (const GlobalOctant& o : octants) {
    EXPECT_NEAR(o.bounds.y_min, 0.0, 1e-12);
    EXPECT_NEAR(o.bounds.y_max, 1.0, 1e-12);
    EXPECT_NEAR(o.bounds.z_min, 0.0, 1e-12);
    EXPECT_NEAR(o.bounds.z_max, 1.0, 1e-12);
    EXPECT_EQ(o.level, 0);
    EXPECT_NEAR(o.bounds.x_max - o.bounds.x_min, 0.5, 1e-12);
    x_min = std::min(x_min, o.bounds.x_min);
    x_max = std::max(x_max, o.bounds.x_max);
  }
  EXPECT_NEAR(x_min, 0.0, 1e-12);
  EXPECT_NEAR(x_max, 1.0, 1e-12);
}

TEST(OctreeForest, UniformRefineTwoLevelsYields64OctantsGlobally) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);
  const label total = GlobalSum(forest.local_num_octants(), MPI_COMM_WORLD);
  EXPECT_EQ(total, 64);
}

TEST(OctreeForest, QuadrantLevelIsTwoAfterUniformRefine) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);
  for (label i = 0; i < forest.local_num_octants(); ++i) {
    EXPECT_EQ(forest.quadrant_level(i), 2);
  }
}

TEST(OctreeForest, QuadrantBoundsTileUnitCubeWithoutOverlap) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);

  std::vector<GlobalOctant> octants;
  GatherGlobalOctants(forest, MPI_COMM_WORLD, &octants);
  ASSERT_EQ(octants.size(), 64u);

  const BoundingBox domain = UnitCubeDomain();
  scalar volume_sum = 0;
  for (const auto& o : octants) {
    volume_sum += BoxVolume(o.bounds);
    EXPECT_GE(o.bounds.x_min, domain.x_min - 1e-12);
    EXPECT_GE(o.bounds.y_min, domain.y_min - 1e-12);
    EXPECT_GE(o.bounds.z_min, domain.z_min - 1e-12);
    EXPECT_LE(o.bounds.x_max, domain.x_max + 1e-12);
    EXPECT_LE(o.bounds.y_max, domain.y_max + 1e-12);
    EXPECT_LE(o.bounds.z_max, domain.z_max + 1e-12);
  }
  EXPECT_NEAR(volume_sum, BoxVolume(domain), 1e-10);

  for (std::size_t i = 0; i < octants.size(); ++i) {
    for (std::size_t j = i + 1; j < octants.size(); ++j) {
      EXPECT_FALSE(BoxesOverlap(octants[i].bounds, octants[j].bounds))
          << "octants " << i << " and " << j << " overlap";
    }
  }
}

TEST(OctreeForest, BalanceEnforcesTwoToOneLevelJumpAcrossFaces) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 1);
  forest.refine([](OctantId id) { return (id % 2) == 0; }, 2);
  forest.balance();

  std::vector<GlobalOctant> octants;
  GatherGlobalOctants(forest, MPI_COMM_WORLD, &octants);
  for (std::size_t i = 0; i < octants.size(); ++i) {
    for (std::size_t j = i + 1; j < octants.size(); ++j) {
      const int gap = FaceLevelGap(octants[i].bounds, octants[i].level,
                                   octants[j].bounds, octants[j].level);
      if (gap > 0) {
        EXPECT_LE(gap, 1) << "face neighbors " << i << " and " << j;
      }
    }
  }
}

TEST(OctreeForest, PartitionBalancesLocalOctantCounts) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);
  forest.partition();

  const label local = forest.local_num_octants();
  label min_local = local;
  label max_local = local;
  MPI_Allreduce(&local, &min_local, 1, MPI_INT32_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local, &max_local, 1, MPI_INT32_T, MPI_MAX, MPI_COMM_WORLD);
  EXPECT_LE(max_local - min_local, 1);
}

}  // namespace octlb
