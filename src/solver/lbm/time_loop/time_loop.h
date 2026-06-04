#ifndef OCTLB_SRC_SOLVER_LBM_TIME_LOOP_TIME_LOOP_H_
#define OCTLB_SRC_SOLVER_LBM_TIME_LOOP_TIME_LOOP_H_

#include <string>
#include <vector>

#include "src/mesh/forest/octree_forest.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/ghost_schedule.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/level_coupler.h"

namespace octlb {

enum class CouplerPhase {
  kHalfTime,
  kFullTime,
  kRestrict,
};

struct CouplerCall {
  CouplerPhase phase;
  int coarse_level;
};

struct TimeLoopCounters {
  std::vector<int> collide;
  std::vector<int> stream;
  std::vector<CouplerCall> coupler_calls;
};

using TimeLoopLattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

/** Recursive subcycling time integrator for static AMR (T06). */
class TimeLoop {
 public:
  TimeLoop(const OctreeForest& forest, BlockCollection<TimeLoopLattice>& blocks,
           GhostSchedule<TimeLoopLattice>& ghosts, LevelCoupler& coupler,
           double omega);

  void advance_one();

  int max_level() const { return max_level_; }

  const TimeLoopCounters& counters() const { return counters_; }
  void reset_counters();

 private:
  void advance(int level);
  void collide_level(int level);
  void stream_level(int level);

  BlockCollection<TimeLoopLattice>& blocks_;
  GhostSchedule<TimeLoopLattice>& ghosts_;
  LevelCoupler& coupler_;
  double omega_;
  int max_level_;

  std::vector<std::vector<OctantId>> octants_by_level_;
  TimeLoopCounters counters_;
};

inline TimeLoop::TimeLoop(const OctreeForest& forest,
                          BlockCollection<TimeLoopLattice>& blocks,
                          GhostSchedule<TimeLoopLattice>& ghosts,
                          LevelCoupler& coupler, double omega)
    : blocks_(blocks),
      ghosts_(ghosts),
      coupler_(coupler),
      omega_(omega),
      max_level_(0) {
  for (label i = 0; i < forest.local_num_octants(); ++i) {
    const int lvl = forest.quadrant_level(static_cast<OctantId>(i));
    max_level_ = std::max(max_level_, lvl);
  }

  octants_by_level_.assign(static_cast<std::size_t>(max_level_ + 1), {});
  for (label i = 0; i < forest.local_num_octants(); ++i) {
    const int lvl = forest.quadrant_level(static_cast<OctantId>(i));
    octants_by_level_[static_cast<std::size_t>(lvl)].push_back(
        static_cast<OctantId>(i));
  }

  counters_.collide.assign(static_cast<std::size_t>(max_level_ + 1), 0);
  counters_.stream.assign(static_cast<std::size_t>(max_level_ + 1), 0);
}

inline void TimeLoop::reset_counters() {
  std::fill(counters_.collide.begin(), counters_.collide.end(), 0);
  std::fill(counters_.stream.begin(), counters_.stream.end(), 0);
  counters_.coupler_calls.clear();
}

inline void TimeLoop::collide_level(int level) {
  for (OctantId id : octants_by_level_[static_cast<std::size_t>(level)]) {
    blocks_[id].collide(static_cast<double>(omega_));
  }
  ++counters_.collide[static_cast<std::size_t>(level)];
}

inline void TimeLoop::stream_level(int level) {
  for (OctantId id : octants_by_level_[static_cast<std::size_t>(level)]) {
    blocks_[id].stream();
  }
  ++counters_.stream[static_cast<std::size_t>(level)];
}

inline void TimeLoop::advance(int level) {
  collide_level(level);
  // Global same-level halo exchange (v1: all faces each call).
  ghosts_.exchange();
  stream_level(level);

  if (level >= max_level_) {
    return;
  }

  coupler_.apply_half_time(level);
  counters_.coupler_calls.push_back({CouplerPhase::kHalfTime, level});

  advance(level + 1);

  coupler_.apply_full_time(level);
  counters_.coupler_calls.push_back({CouplerPhase::kFullTime, level});

  advance(level + 1);

  coupler_.restrict(level + 1);
  counters_.coupler_calls.push_back({CouplerPhase::kRestrict, level});
}

inline void TimeLoop::advance_one() {
  advance(0);
}

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_TIME_LOOP_TIME_LOOP_H_
