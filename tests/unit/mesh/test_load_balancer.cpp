#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/load_balance/weighted_load_balancer.h"

namespace octlb {
namespace {

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

int64_t LocalWeightSum(const OctreeForest& forest) {
  int64_t sum = 0;
  for (label i = 0; i < forest.local_num_octants(); ++i) {
    sum += static_cast<int64_t>(1) << forest.quadrant_level(i);
  }
  return sum;
}

int64_t GlobalWeightSum(const OctreeForest& forest, MPI_Comm comm) {
  const int64_t local = LocalWeightSum(forest);
  int64_t global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT64_T, MPI_SUM, comm);
  return global;
}

void BuildMixedLevelForest(OctreeForest* forest) {
  forest->refine([](OctantId) { return true; }, 2);
  const OctantId center = [&]() {
    OctantId best = 0;
    scalar best_dist = 1e30;
    for (label i = 0; i < forest->local_num_octants(); ++i) {
      const BoundingBox b = forest->quadrant_bounds(i);
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
  }();
  forest->refine([center](OctantId id) { return id == center; }, 3);
  forest->balance();
}

}  // namespace

TEST(LoadBalancer, WeightedPartitionBalancesTotalWeight) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  BuildMixedLevelForest(&forest);

  const int64_t target = GlobalWeightSum(forest, MPI_COMM_WORLD);
  forest.partition(make_level_weight_fn(forest));

  const int64_t local_w = LocalWeightSum(forest);
  int64_t min_w = local_w;
  int64_t max_w = local_w;
  MPI_Allreduce(&local_w, &min_w, 1, MPI_INT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_w, &max_w, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);

  const double imbalance =
      static_cast<double>(max_w - min_w) / static_cast<double>(target);
  EXPECT_LT(imbalance, 0.05);
}

double WeightImbalanceFraction(const OctreeForest& forest, MPI_Comm comm) {
  const int64_t local = LocalWeightSum(forest);
  const int64_t global = GlobalWeightSum(forest, comm);
  int64_t min_w = local;
  int64_t max_w = local;
  MPI_Allreduce(&local, &min_w, 1, MPI_INT64_T, MPI_MIN, comm);
  MPI_Allreduce(&local, &max_w, 1, MPI_INT64_T, MPI_MAX, comm);
  if (global == 0) {
    return 0.0;
  }
  return static_cast<double>(max_w - min_w) / static_cast<double>(global);
}

TEST(LoadBalancer, WeightedPartitionImprovesWeightBalanceVersusUniform) {
  OctreeForest uniform(MPI_COMM_WORLD, UnitCubeDomain());
  BuildMixedLevelForest(&uniform);
  uniform.partition();

  OctreeForest weighted(MPI_COMM_WORLD, UnitCubeDomain());
  BuildMixedLevelForest(&weighted);
  weighted.partition(make_level_weight_fn(weighted));

  const double uniform_imbalance =
      WeightImbalanceFraction(uniform, MPI_COMM_WORLD);
  const double weighted_imbalance =
      WeightImbalanceFraction(weighted, MPI_COMM_WORLD);

  EXPECT_GT(uniform_imbalance, weighted_imbalance);
  EXPECT_LT(weighted_imbalance, 0.05);
}

}  // namespace octlb
