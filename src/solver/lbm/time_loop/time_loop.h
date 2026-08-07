#ifndef OCTLB_SRC_SOLVER_LBM_TIME_LOOP_TIME_LOOP_H_
#define OCTLB_SRC_SOLVER_LBM_TIME_LOOP_TIME_LOOP_H_

#include <functional>
#include <string>
#include <vector>

#include <mpi.h>

#include "src/mesh/forest/octree_forest.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/ghost_schedule.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/level_coupler.h"
#include "src/solver/lbm/domain_boundary_handler.h"

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

/** ConstRhoBGK statistics aggregation scope (cavity3d A/B). */
enum class ConstRhoStatsScope {
  /** OpenLB-like: fluid corrected rho + boundary raw rho in one average. */
  kFluidAndBoundary,
  /** Option B: only kFluid ConstRho cells contribute to average_rho. */
  kFluidOnly,
};

/** Recursive subcycling time integrator for static AMR (T06). */
class TimeLoop {
 public:
  TimeLoop(const OctreeForest& forest, BlockCollection<TimeLoopLattice>& blocks,
           GhostSchedule<TimeLoopLattice>& ghosts, LevelCoupler& coupler,
           DomainBoundaryHandler& domain_bc, double omega,
           bool use_const_rho_bgk = false,
           ConstRhoStatsScope const_rho_stats_scope =
               ConstRhoStatsScope::kFluidAndBoundary);

  void advance_one();

  // Single-level diagnostic hooks (cavity3d flat loop).
  enum class FlatPhase {
    kAfterCollide,
    kAfterBoundaryCollide,
    kAfterGhostExchange,
    kAfterStream,
    kAfterPostStream,
  };

  void advance_one_flat_with_hooks(
      const std::function<void(FlatPhase)>& hook = {});

  // Same phase order as advance_one() / advance(0) (interleaved BC when enabled).
  void advance_one_with_hooks(
      const std::function<void(FlatPhase)>& hook = {});

  int max_level() const { return max_level_; }

  const TimeLoopCounters& counters() const { return counters_; }
  void reset_counters();

  // T11 MEM-timing diagnostic: when enabled, snapshot every block's populations_
  // right before stream_level (the post-collide + post-BC + post-ghost pre-stream
  // state) so a drag post-processor can recompute the MEM force at the correct
  // (outgoing-f_i) timing instead of the post-stream (bounced-back) timing. Off
  // by default (extra memcpy per block per step).
  void set_snapshot_post_collide(bool on) { snapshot_post_collide_ = on; }
  bool snapshot_post_collide() const { return snapshot_post_collide_; }

  /** ConstRhoBGK average rho from the previous collide (OpenLB statistics lag). */
  double average_rho() const { return average_rho_; }

 private:
  void advance(int level);
  void collide_level(int level, CollideRhoStats* rho_stats = nullptr,
                     bool openlb_spatial_collide_order = false);
  void stream_level(int level);

  BlockCollection<TimeLoopLattice>& blocks_;
  GhostSchedule<TimeLoopLattice>& ghosts_;
  LevelCoupler& coupler_;
  DomainBoundaryHandler& domain_bc_;
  double omega_;
  bool use_const_rho_bgk_;
  ConstRhoStatsScope const_rho_stats_scope_;
  double average_rho_;
  int max_level_;
  // 0-based coarse-step counter threaded to the boundary handler as time (for
  // inlet velocity ramp-up). OpenLB iT semantics: the first step uses t=0.
  int step_ = 0;

  CollideRhoStats* boundary_rho_stats(CollideRhoStats* fluid_stats) const {
    if (fluid_stats == nullptr) {
      return nullptr;
    }
    return const_rho_stats_scope_ == ConstRhoStatsScope::kFluidAndBoundary
               ? fluid_stats
               : nullptr;
  }

  std::vector<std::vector<OctantId>> octants_by_level_;
  TimeLoopCounters counters_;
  bool snapshot_post_collide_ = false;
};

inline TimeLoop::TimeLoop(const OctreeForest& forest,
                          BlockCollection<TimeLoopLattice>& blocks,
                          GhostSchedule<TimeLoopLattice>& ghosts,
                          LevelCoupler& coupler,
                          DomainBoundaryHandler& domain_bc, double omega,
                          bool use_const_rho_bgk,
                          ConstRhoStatsScope const_rho_stats_scope)
    : blocks_(blocks),
      ghosts_(ghosts),
      coupler_(coupler),
      domain_bc_(domain_bc),
      omega_(omega),
      use_const_rho_bgk_(use_const_rho_bgk),
      const_rho_stats_scope_(const_rho_stats_scope),
      average_rho_(1.0),
      max_level_(0) {
  for (label i = 0; i < forest.local_num_octants(); ++i) {
    const int lvl = forest.quadrant_level(static_cast<OctantId>(i));
    max_level_ = std::max(max_level_, lvl);
  }
  // max_level_ must be GLOBAL (the max level across ALL ranks), not per-rank.
  // advance(level) recurses to max_level_ and calls coupler.apply_half_time /
  // apply_full_time at each level. If ranks have different per-rank max levels
  // (AMR partition: some ranks own the refined octants, others only coarse
  // level-0), a rank with the deeper max would call apply_half_time(L) and
  // post a cross-rank Irecv at a level L the shallower peer never reaches (its
  // advance returns at its own lower max_level_), so the peer never posts the
  // matching Isend -> Waitall hangs. Allreduce so every rank recurses to the
  // same depth; ranks with no coarse-fine at level L post nothing (empty
  // Waitall), keeping the coupler exchange symmetric per level.
  MPI_Allreduce(MPI_IN_PLACE, &max_level_, 1, MPI_INT, MPI_MAX, forest.comm());

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
  step_ = 0;
}

inline void TimeLoop::collide_level(int level, CollideRhoStats* rho_stats,
                                   bool openlb_spatial_collide_order) {
  if (openlb_spatial_collide_order) {
    for (OctantId id : octants_by_level_[static_cast<std::size_t>(level)]) {
      domain_bc_.collide_interleaved_with(
          blocks_[id], rho_stats, average_rho_, use_const_rho_bgk_);
    }
  } else {
    for (OctantId id : octants_by_level_[static_cast<std::size_t>(level)]) {
      if (use_const_rho_bgk_) {
        blocks_[id].collide_const_rho(static_cast<double>(omega_),
                                     static_cast<double>(average_rho_),
                                     rho_stats);
      } else {
        blocks_[id].collide(static_cast<double>(omega_));
      }
    }
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
  CollideRhoStats rho_stats;
  CollideRhoStats* stats_ptr =
      (use_const_rho_bgk_ && level >= max_level_) ? &rho_stats : nullptr;
  const bool openlb_order = domain_bc_.collide_boundary_before_bulk();
  collide_level(level, stats_ptr, openlb_order);
  if (!openlb_order) {
    domain_bc_.apply(boundary_rho_stats(stats_ptr), average_rho_,
                     use_const_rho_bgk_);
  }
  ghosts_.exchange();
  if (snapshot_post_collide_) {
    for (OctantId id : octants_by_level_[static_cast<std::size_t>(level)]) {
      blocks_[id].take_post_collide_snapshot();
    }
  }
  stream_level(level);
  domain_bc_.apply_post_stream();
  // OpenLB collectStatistics() runs after PostStream; stats are from collide only.
  if (stats_ptr != nullptr) {
    average_rho_ = rho_stats.average();
  }

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
  domain_bc_.set_time(static_cast<double>(step_));
  advance(0);
  ++step_;
}

inline void TimeLoop::advance_one_flat_with_hooks(
    const std::function<void(FlatPhase)>& hook) {
  CollideRhoStats rho_stats;
  CollideRhoStats* stats_ptr =
      (use_const_rho_bgk_ && max_level_ >= 0) ? &rho_stats : nullptr;
  collide_level(0, stats_ptr, /*openlb_spatial_collide_order=*/false);
  if (hook) {
    hook(FlatPhase::kAfterCollide);
  }
  domain_bc_.apply(boundary_rho_stats(stats_ptr), average_rho_,
                   use_const_rho_bgk_);
  if (hook) {
    hook(FlatPhase::kAfterBoundaryCollide);
  }
  ghosts_.exchange();
  if (hook) {
    hook(FlatPhase::kAfterGhostExchange);
  }
  stream_level(0);
  if (hook) {
    hook(FlatPhase::kAfterStream);
  }
  domain_bc_.apply_post_stream();
  if (stats_ptr != nullptr) {
    average_rho_ = rho_stats.average();
  }
  if (hook) {
    hook(FlatPhase::kAfterPostStream);
  }
}

inline void TimeLoop::advance_one_with_hooks(
    const std::function<void(FlatPhase)>& hook) {
  CollideRhoStats rho_stats;
  CollideRhoStats* stats_ptr =
      (use_const_rho_bgk_ && max_level_ >= 0) ? &rho_stats : nullptr;
  const bool openlb_order = domain_bc_.collide_boundary_before_bulk();
  collide_level(0, stats_ptr, openlb_order);
  if (hook) {
    hook(FlatPhase::kAfterCollide);
  }
  if (!openlb_order) {
    domain_bc_.apply(boundary_rho_stats(stats_ptr), average_rho_,
                     use_const_rho_bgk_);
    if (hook) {
      hook(FlatPhase::kAfterBoundaryCollide);
    }
  }
  ghosts_.exchange();
  if (hook) {
    hook(FlatPhase::kAfterGhostExchange);
  }
  stream_level(0);
  if (hook) {
    hook(FlatPhase::kAfterStream);
  }
  domain_bc_.apply_post_stream();
  if (stats_ptr != nullptr) {
    average_rho_ = rho_stats.average();
  }
  if (hook) {
    hook(FlatPhase::kAfterPostStream);
  }
}

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_TIME_LOOP_TIME_LOOP_H_
