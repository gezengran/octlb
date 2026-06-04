#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/ghost_schedule.h"
#include "src/solver/lbm/level_coupler.h"
#include "src/solver/lbm/time_loop/time_loop.h"

namespace octlb {
namespace {

constexpr int kN = 4;
constexpr double kOmega = 1.0;

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

OctantId FindCenterOctant(const OctreeForest& forest) {
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

OctreeForest MakeThreeLevelForest() {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 1);
  OctantId center = FindCenterOctant(forest);
  forest.refine([center](OctantId id) { return id == center; }, 2);
  forest.balance();
  forest.partition();
  return forest;
}

int GlobalMaxLevel(const OctreeForest& forest) {
  int local_max = 0;
  for (label i = 0; i < forest.local_num_octants(); ++i) {
    local_max = std::max(local_max, forest.quadrant_level(i));
  }
  int global_max = 0;
  MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return global_max;
}

}  // namespace

TEST(TimeLoopLevels, ThreeLevels_StepCountRatio_1_2_4) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest = MakeThreeLevelForest();
  ASSERT_EQ(GlobalMaxLevel(forest), 2);

  const FacePairList pairs(forest);
  BlockCollection<TimeLoopLattice> blocks(
      forest.local_num_octants(), [](OctantId) {
        TimeLoopLattice lat(kN, kN, kN, 1);
        const double u0[3] = {0.0, 0.0, 0.0};
        lat.initialize(1.0, u0);
        return lat;
      });

  GhostSchedule<TimeLoopLattice> ghosts(MPI_COMM_WORLD, pairs, blocks, kN, kN,
                                        kN);
  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);
  TimeLoop loop(forest, blocks, ghosts, coupler, kOmega);

  loop.advance_one();

  EXPECT_EQ(loop.counters().collide[0], 1);
  EXPECT_EQ(loop.counters().collide[1], 2);
  EXPECT_EQ(loop.counters().collide[2], 4);
  EXPECT_EQ(loop.counters().stream[0], 1);
  EXPECT_EQ(loop.counters().stream[1], 2);
  EXPECT_EQ(loop.counters().stream[2], 4);
}

TEST(TimeLoopLevels, ThreeLevels_CouplerCallOrder) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest = MakeThreeLevelForest();
  const FacePairList pairs(forest);
  BlockCollection<TimeLoopLattice> blocks(
      forest.local_num_octants(), [](OctantId) {
        TimeLoopLattice lat(kN, kN, kN, 1);
        const double u0[3] = {0.0, 0.0, 0.0};
        lat.initialize(1.0, u0);
        return lat;
      });

  GhostSchedule<TimeLoopLattice> ghosts(MPI_COMM_WORLD, pairs, blocks, kN, kN,
                                        kN);
  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);
  TimeLoop loop(forest, blocks, ghosts, coupler, kOmega);
  loop.advance_one();

  const auto& calls = loop.counters().coupler_calls;
  ASSERT_FALSE(calls.empty());
  EXPECT_EQ(calls.front().phase, CouplerPhase::kHalfTime);
  EXPECT_EQ(calls.front().coarse_level, 0);

  std::size_t full0 = calls.size();
  std::size_t restrict0 = calls.size();
  for (std::size_t i = 1; i < calls.size(); ++i) {
    if (full0 == calls.size() && calls[i].phase == CouplerPhase::kFullTime &&
        calls[i].coarse_level == 0) {
      full0 = i;
    }
    if (restrict0 == calls.size() &&
        calls[i].phase == CouplerPhase::kRestrict &&
        calls[i].coarse_level == 0) {
      restrict0 = i;
    }
  }
  EXPECT_LT(full0, calls.size());
  EXPECT_LT(restrict0, calls.size());
  EXPECT_LT(0u, full0);
  EXPECT_LT(full0, restrict0);
}

TEST(TimeLoopLevels, AdvanceOne_WithBlockLattice_NoThrow) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest = MakeThreeLevelForest();
  const FacePairList pairs(forest);
  BlockCollection<TimeLoopLattice> blocks(
      forest.local_num_octants(), [](OctantId) {
        TimeLoopLattice lat(kN, kN, kN, 1);
        const double u0[3] = {0.0, 0.0, 0.0};
        lat.initialize(1.0, u0);
        return lat;
      });

  GhostSchedule<TimeLoopLattice> ghosts(MPI_COMM_WORLD, pairs, blocks, kN, kN,
                                        kN);
  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);
  TimeLoop loop(forest, blocks, ghosts, coupler, kOmega);

  for (int step = 0; step < 3; ++step) {
    loop.advance_one();
  }

  for (label id = 0; id < blocks.size(); ++id) {
    auto& lat = blocks[static_cast<OctantId>(id)];
    for (int ix = 0; ix < kN; ++ix) {
      for (int iy = 0; iy < kN; ++iy) {
        for (int iz = 0; iz < kN; ++iz) {
          auto cell = lat.get(ix, iy, iz);
          double rho = 0.0;
          double u[3] = {};
          cell.computeRhoU(rho, u);
          EXPECT_TRUE(std::isfinite(rho));
          EXPECT_TRUE(std::isfinite(u[0]));
        }
      }
    }
  }
}

}  // namespace octlb
