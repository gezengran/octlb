#ifndef OCTLB_SRC_SOLVER_LBM_LEVEL_COUPLER_H_
#define OCTLB_SRC_SOLVER_LBM_LEVEL_COUPLER_H_

#include <cstddef>
#include <vector>

#include <mpi.h>

#include "block_lattice.h"
#include "src/common/types.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/face_iterator.h"

namespace octlb {

struct CouplingPoint {
  OctantId coarse_id;
  int ci, cj, ck;
  OctantId fine_id;
  int fi, fj, fk;
  FaceDir normal;
  int coarse_level;
  int fine_level;
  int remote_rank;
  int comm_tag;
  int fine_slot;
  int coarse_remote_rank;
};

using LatticeD3Q19 = BlockLattice<double, olb::descriptors::D3Q19<>>;

/** Lagrava coarse-fine coupling for static AMR (T06). */
class LevelCoupler {
 public:
  LevelCoupler(MPI_Comm comm, const FacePairList& faces,
               const OctreeForest& forest, BlockCollection<LatticeD3Q19>& blocks,
               int nx, int ny, int nz, double omega);

  void apply_half_time(int coarse_level);
  void apply_full_time(int coarse_level);
  void restrict(int fine_level);

  const std::vector<CouplingPoint>& coupling_plan() const { return plan_; }

 private:
  struct MacroState {
    double rho{1.0};
    double u[3]{};
    double f_neq[LatticeD3Q19::kQ]{};
  };

  static constexpr int kMacroDoubles = 1 + 3 + LatticeD3Q19::kQ;

  void PackMacro(const MacroState& m, double* buf) const;
  void UnpackMacro(const double* buf, MacroState* m) const;
  void ReadMacro(const LatticeD3Q19& lat, int ix, int iy, int iz,
                 MacroState* out) const;
  void WriteProlongation(LatticeD3Q19& fine, int fi, int fj, int fk,
                         const MacroState& macro, double scaling) const;
  void WriteRestriction(LatticeD3Q19& coarse, int ci, int cj, int ck,
                        const MacroState& fine_macro,
                        const MacroState& neighbor_fneq_sum,
                        int neighbor_count, double scaling) const;

  void ExchangeCoarseMacros(int coarse_level, bool include_prev);
  void ApplyProlongation(std::size_t begin, std::size_t end, bool half_time);

  std::size_t LevelRangeBegin(int coarse_level) const;
  std::size_t LevelRangeEnd(int coarse_level) const;

  MPI_Comm comm_;
  BlockCollection<LatticeD3Q19>& blocks_;
  int nx_, ny_, nz_;
  double omega_;
  double tau_;
  double prolong_scale_;
  double restrict_scale_;
  int my_rank_;
  label num_local_octants_;

  std::vector<CouplingPoint> plan_;
  std::vector<MacroState> prev_macro_;
  std::vector<MacroState> curr_macro_;
  std::vector<MacroState> recv_macro_;
  std::vector<MacroState> recv_prev_macro_;

  std::vector<std::size_t> level_begin_;
  std::vector<std::size_t> level_end_;

  struct MpiBatch {
    int peer_rank;
    int comm_tag;
    bool coarse_sends;
    std::vector<std::size_t> plan_indices;
  };
  std::vector<MpiBatch> mpi_batches_;
  std::vector<double> send_buf_;
  std::vector<double> recv_buf_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_LEVEL_COUPLER_H_
