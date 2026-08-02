#include "level_coupler.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <unordered_set>

#include "core/util.h"
#include "dynamics/lbm.h"
#include "src/common/bounding_box.h"

namespace octlb {
namespace {

struct AxisInfo {
  int normal_axis;
  int tan0;
  int tan1;
  int coarse_normal_idx;
  int fine_normal_idx;
};

AxisInfo MakeAxisInfo(FaceDir normal, int nx, int ny, int nz) {
  const int face = static_cast<int>(normal);
  const int normal_axis = face / 2;
  const bool is_max = (face % 2) == 1;
  const int n[3] = {nx, ny, nz};

  AxisInfo info{};
  info.normal_axis = normal_axis;
  info.tan0 = (normal_axis + 1) % 3;
  info.tan1 = (normal_axis + 2) % 3;
  info.coarse_normal_idx = is_max ? n[normal_axis] - 1 : 0;
  info.fine_normal_idx = is_max ? 0 : n[normal_axis] - 1;
  return info;
}

int QuadrantFromCenter(scalar center, scalar ref) {
  return center < ref ? 0 : 1;
}

struct CouplingKey {
  OctantId coarse_id;
  int ci, cj, ck;
  OctantId fine_id;
  int fi, fj, fk;

  bool operator==(const CouplingKey& o) const {
    return coarse_id == o.coarse_id && ci == o.ci && cj == o.cj && ck == o.ck &&
           fine_id == o.fine_id && fi == o.fi && fj == o.fj && fk == o.fk;
  }
};

struct CouplingKeyHash {
  std::size_t operator()(const CouplingKey& k) const {
    std::size_t h = static_cast<std::size_t>(k.coarse_id);
    h = h * 131u + static_cast<std::size_t>(k.ci);
    h = h * 131u + static_cast<std::size_t>(k.cj);
    h = h * 131u + static_cast<std::size_t>(k.ck);
    h = h * 131u + static_cast<std::size_t>(k.fine_id);
    h = h * 131u + static_cast<std::size_t>(k.fi);
    h = h * 131u + static_cast<std::size_t>(k.fj);
    h = h * 131u + static_cast<std::size_t>(k.fk);
    return h;
  }
};

void AppendFaceCouplingPoints(const CoarseFineFace& face, int nx, int ny,
                              int nz, std::vector<CouplingPoint>* out) {
  if ((nx % 2) != 0 || (ny % 2) != 0 || (nz % 2) != 0) {
    throw std::runtime_error("LevelCoupler: block dimensions must be even");
  }

  const BoundingBox coarse_bounds = face.coarse_bounds;
  const scalar coarse_center[3] = {
      0.5 * (coarse_bounds.x_min + coarse_bounds.x_max),
      0.5 * (coarse_bounds.y_min + coarse_bounds.y_max),
      0.5 * (coarse_bounds.z_min + coarse_bounds.z_max)};

  const AxisInfo axes = MakeAxisInfo(face.normal, nx, ny, nz);
  const int n[3] = {nx, ny, nz};
  const int half[3] = {nx / 2, ny / 2, nz / 2};

  scalar fine_centers[4][3];
  for (int i = 0; i < 4; ++i) {
    // Defect 5: use the geometry snapshotted by FacePairList at iterate time,
    // not forest.quadrant_bounds(fine_ids[i]) -- the remote fine side's quadid
    // is a transient ghost-array index, invalid after FacePairList construction.
    const BoundingBox fb = face.fine_bounds[i];
    fine_centers[i][0] = 0.5 * (fb.x_min + fb.x_max);
    fine_centers[i][1] = 0.5 * (fb.y_min + fb.y_max);
    fine_centers[i][2] = 0.5 * (fb.z_min + fb.z_max);
  }

  for (int slot = 0; slot < 4; ++slot) {
    const int qu0 = QuadrantFromCenter(fine_centers[slot][axes.tan0],
                                       coarse_center[axes.tan0]);
    const int qu1 = QuadrantFromCenter(fine_centers[slot][axes.tan1],
                                       coarse_center[axes.tan1]);

    const int t0_begin = qu0 * half[axes.tan0];
    const int t0_end = t0_begin + half[axes.tan0];
    const int t1_begin = qu1 * half[axes.tan1];
    const int t1_end = t1_begin + half[axes.tan1];

    for (int t0 = t0_begin; t0 < t0_end; ++t0) {
      for (int t1 = t1_begin; t1 < t1_end; ++t1) {
        int coarse_idx[3] = {0, 0, 0};
        int fine_idx[3] = {0, 0, 0};
        coarse_idx[axes.normal_axis] = axes.coarse_normal_idx;
        fine_idx[axes.normal_axis] = axes.fine_normal_idx;
        coarse_idx[axes.tan0] = t0;
        coarse_idx[axes.tan1] = t1;
        fine_idx[axes.tan0] = t0 - t0_begin;
        fine_idx[axes.tan1] = t1 - t1_begin;

        CouplingPoint pt{};
        pt.coarse_id = face.coarse_id;
        pt.ci = coarse_idx[0];
        pt.cj = coarse_idx[1];
        pt.ck = coarse_idx[2];
        pt.fine_id = face.fine_ids[slot];
        pt.fi = fine_idx[0];
        pt.fj = fine_idx[1];
        pt.fk = fine_idx[2];
        pt.normal = face.normal;
        pt.coarse_level = face.coarse_level;
        pt.fine_level = face.fine_level;
        pt.fine_slot = slot;
        pt.comm_tag = face.comm_tags[slot];
        pt.remote_rank = face.remote_ranks[slot];
        pt.coarse_remote_rank = face.coarse_remote_rank;
        out->push_back(pt);
      }
    }
  }
}

std::vector<CouplingPoint> BuildCouplingPlan(const FacePairList& faces,
                                            int nx, int ny, int nz) {
  std::vector<CouplingPoint> plan;
  FaceIterator it(faces);
  for (const CoarseFineFace& face : it) {
    AppendFaceCouplingPoints(face, nx, ny, nz, &plan);
  }
  return plan;
}

// Ownership must NOT use "id < num_local": for a cross-rank coarse-fine face the
// remote side's id is a p4est ghost-array index, which routinely falls in
// [0, num_local) and would look "owned". Use the ranks FacePairList recorded.
bool OwnsCoarse(const CouplingPoint& pt, int my_rank) {
  return pt.coarse_remote_rank == my_rank;
}
bool OwnsFine(const CouplingPoint& pt, int my_rank) {
  return pt.remote_rank == my_rank;
}

}  // namespace

void LevelCoupler::PackMacro(const MacroState& m, double* buf) const {
  buf[0] = m.rho;
  buf[1] = m.u[0];
  buf[2] = m.u[1];
  buf[3] = m.u[2];
  for (int i = 0; i < LatticeD3Q19::kQ; ++i) {
    buf[4 + i] = m.f_neq[i];
  }
}

void LevelCoupler::UnpackMacro(const double* buf, MacroState* m) const {
  m->rho = buf[0];
  m->u[0] = buf[1];
  m->u[1] = buf[2];
  m->u[2] = buf[3];
  for (int i = 0; i < LatticeD3Q19::kQ; ++i) {
    m->f_neq[i] = buf[4 + i];
  }
}

void LevelCoupler::ReadMacro(const LatticeD3Q19& lat, int ix, int iy, int iz,
                             MacroState* out) const {
  auto cell = const_cast<LatticeD3Q19&>(lat).get(ix, iy, iz);
  olb::lbm<olb::descriptors::D3Q19<>>::computeRhoU(cell, out->rho, out->u);
  olb::lbm<olb::descriptors::D3Q19<>>::computeFneq(cell, out->f_neq, out->rho,
                                                   out->u);
}

void LevelCoupler::WriteProlongation(LatticeD3Q19& fine, int fi, int fj, int fk,
                                     const MacroState& macro,
                                     double scaling) const {
  auto cell = fine.get(fi, fj, fk);
  const double uSqr = olb::util::normSqr<double, 3>(macro.u);
  for (int iPop = 0; iPop < LatticeD3Q19::kQ; ++iPop) {
    cell[iPop] = olb::equilibrium<olb::descriptors::D3Q19<>>::secondOrder(
                     iPop, macro.rho, macro.u, uSqr) +
                 scaling * macro.f_neq[iPop];
  }
}

void LevelCoupler::WriteRestriction(LatticeD3Q19& coarse, int ci, int cj, int ck,
                                    const MacroState& fine_macro,
                                    const MacroState& neighbor_fneq_sum,
                                    int neighbor_count, double scaling) const {
  MacroState avg = fine_macro;
  for (int iPop = 0; iPop < LatticeD3Q19::kQ; ++iPop) {
    avg.f_neq[iPop] =
        (fine_macro.f_neq[iPop] + neighbor_fneq_sum.f_neq[iPop]) /
        static_cast<double>(neighbor_count + 1);
  }
  WriteProlongation(coarse, ci, cj, ck, avg, scaling);
}

LevelCoupler::LevelCoupler(MPI_Comm comm, const FacePairList& faces,
                           const OctreeForest& forest,
                           BlockCollection<LatticeD3Q19>& blocks, int nx,
                           int ny, int nz, double omega)
    : comm_(comm),
      blocks_(blocks),
      nx_(nx),
      ny_(ny),
      nz_(nz),
      omega_(omega),
      tau_(1.0 / omega),
      prolong_scale_((tau_ - 0.25) / tau_),
      restrict_scale_(tau_ / (tau_ - 0.25)),
      num_local_octants_(forest.local_num_octants()) {
  MPI_Comm_rank(comm_, &my_rank_);

  plan_ = BuildCouplingPlan(faces, nx, ny, nz);
  prev_macro_.resize(plan_.size());
  curr_macro_.resize(plan_.size());
  recv_macro_.resize(plan_.size());
  recv_prev_macro_.resize(plan_.size());

  int max_level = 0;
  for (const CouplingPoint& pt : plan_) {
    max_level = std::max(max_level, pt.coarse_level);
  }
  level_begin_.assign(static_cast<std::size_t>(max_level + 1), plan_.size());
  level_end_.assign(static_cast<std::size_t>(max_level + 1), 0);
  for (std::size_t i = 0; i < plan_.size(); ++i) {
    const int lvl = plan_[i].coarse_level;
    level_begin_[static_cast<std::size_t>(lvl)] =
        std::min(level_begin_[static_cast<std::size_t>(lvl)], i);
    level_end_[static_cast<std::size_t>(lvl)] =
        std::max(level_end_[static_cast<std::size_t>(lvl)], i + 1);
  }

  auto find_batch = [&](int peer, int tag, bool sends) -> MpiBatch* {
    for (MpiBatch& batch : mpi_batches_) {
      if (batch.peer_rank == peer && batch.comm_tag == tag &&
          batch.coarse_sends == sends) {
        return &batch;
      }
    }
    return nullptr;
  };

  for (std::size_t i = 0; i < plan_.size(); ++i) {
    const CouplingPoint& pt = plan_[i];
    const bool owns_coarse = OwnsCoarse(pt, my_rank_);
    const bool owns_fine = OwnsFine(pt, my_rank_);

    if (owns_coarse) {
      ReadMacro(blocks_[pt.coarse_id], pt.ci, pt.cj, pt.ck, &prev_macro_[i]);
      curr_macro_[i] = prev_macro_[i];
    }

    if (owns_coarse && !owns_fine) {
      MpiBatch* batch = find_batch(pt.remote_rank, pt.comm_tag, true);
      if (batch == nullptr) {
        mpi_batches_.push_back(
            MpiBatch{pt.remote_rank, pt.comm_tag, true, {}});
        batch = &mpi_batches_.back();
      }
      batch->plan_indices.push_back(i);
    } else if (owns_fine && !owns_coarse) {
      MpiBatch* batch =
          find_batch(pt.coarse_remote_rank, pt.comm_tag, false);
      if (batch == nullptr) {
        mpi_batches_.push_back(
            MpiBatch{pt.coarse_remote_rank, pt.comm_tag, false, {}});
        batch = &mpi_batches_.back();
      }
      batch->plan_indices.push_back(i);
    }
  }

  std::size_t max_payload = 0;
  for (const MpiBatch& batch : mpi_batches_) {
    max_payload = std::max(max_payload, batch.plan_indices.size());
  }
  const std::size_t buf_stride =
      max_payload * static_cast<std::size_t>(kMacroDoubles) * 2;
  send_buf_.resize(mpi_batches_.size() * buf_stride);
  recv_buf_.resize(mpi_batches_.size() * buf_stride);
}

std::size_t LevelCoupler::LevelRangeBegin(int coarse_level) const {
  if (coarse_level < 0 ||
      static_cast<std::size_t>(coarse_level) >= level_begin_.size()) {
    return plan_.size();
  }
  return level_begin_[static_cast<std::size_t>(coarse_level)];
}

std::size_t LevelCoupler::LevelRangeEnd(int coarse_level) const {
  if (coarse_level < 0 ||
      static_cast<std::size_t>(coarse_level) >= level_end_.size()) {
    return plan_.size();
  }
  return level_end_[static_cast<std::size_t>(coarse_level)];
}

void LevelCoupler::ExchangeCoarseMacros(int coarse_level, bool include_prev) {
  const std::size_t begin = LevelRangeBegin(coarse_level);
  const std::size_t end = LevelRangeEnd(coarse_level);

  for (std::size_t i = begin; i < end; ++i) {
    const CouplingPoint& pt = plan_[i];
    if (OwnsCoarse(pt, my_rank_)) {
      ReadMacro(blocks_[pt.coarse_id], pt.ci, pt.cj, pt.ck, &curr_macro_[i]);
    }
  }

  std::size_t max_payload = 0;
  for (const MpiBatch& batch : mpi_batches_) {
    max_payload = std::max(max_payload, batch.plan_indices.size());
  }
  const std::size_t stride =
      max_payload * static_cast<std::size_t>(kMacroDoubles) * 2;
  const int macro_block =
      static_cast<int>(kMacroDoubles * (include_prev ? 2 : 1));

  std::vector<MPI_Request> requests;
  requests.reserve(mpi_batches_.size() * 2);

  const bool dbg = std::getenv("OCTLB_COUPLE_DEBUG") != nullptr;
  const bool verbose = std::getenv("OCTLB_MPI_VERBOSE") != nullptr;
  if (dbg) {
    std::fprintf(stderr,
                 "[couple r%d] ExchangeCoarseMacros L=%d prev=%d batches=%zu\n",
                 my_rank_, coarse_level, include_prev ? 1 : 0,
                 mpi_batches_.size());
    if (verbose) {
      for (std::size_t b = 0; b < mpi_batches_.size(); ++b) {
        const MpiBatch& batch = mpi_batches_[b];
        int c = 0;
        for (std::size_t pi : batch.plan_indices) {
          if (plan_[pi].coarse_level == coarse_level) ++c;
        }
        std::fprintf(stderr,
                     "[couple r%d]   batch[%zu] peer=%d tag=%d sends=%d "
                     "plan_idx=%zu count@L=%d\n",
                     my_rank_, b, batch.peer_rank, batch.comm_tag,
                     batch.coarse_sends ? 1 : 0, batch.plan_indices.size(), c);
      }
    }
    std::fflush(stderr);
  }

  for (std::size_t b = 0; b < mpi_batches_.size(); ++b) {
    const MpiBatch& batch = mpi_batches_[b];
    double* base = &send_buf_[b * stride];

    int count = 0;
    for (std::size_t plan_index : batch.plan_indices) {
      if (plan_[plan_index].coarse_level != coarse_level) {
        continue;
      }
      PackMacro(curr_macro_[plan_index], base + count * kMacroDoubles);
      if (include_prev) {
        PackMacro(prev_macro_[plan_index],
                  base + count * kMacroDoubles + kMacroDoubles);
      }
      ++count;
    }
    if (count == 0) {
      continue;
    }

    if (batch.coarse_sends) {
      MPI_Request req{};
      MPI_Isend(base, count * macro_block, MPI_DOUBLE, batch.peer_rank,
                batch.comm_tag, comm_, &req);
      requests.push_back(req);
    } else {
      MPI_Request req{};
      MPI_Irecv(&recv_buf_[b * stride], count * macro_block, MPI_DOUBLE,
                batch.peer_rank, batch.comm_tag, comm_, &req);
      requests.push_back(req);
    }
  }

  if (dbg) {
    std::fprintf(stderr, "[couple r%d] Waitall %zu requests\n", my_rank_,
                 requests.size());
    std::fflush(stderr);
  }

  if (!requests.empty()) {
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(),
                MPI_STATUSES_IGNORE);
  }

  for (std::size_t b = 0; b < mpi_batches_.size(); ++b) {
    const MpiBatch& batch = mpi_batches_[b];
    double* base = batch.coarse_sends ? &send_buf_[b * stride]
                                      : &recv_buf_[b * stride];

    int count = 0;
    for (std::size_t plan_index : batch.plan_indices) {
      if (plan_[plan_index].coarse_level != coarse_level) {
        continue;
      }
      if (!batch.coarse_sends) {
        UnpackMacro(base + count * kMacroDoubles, &recv_macro_[plan_index]);
        if (include_prev) {
          UnpackMacro(base + count * kMacroDoubles + kMacroDoubles,
                      &recv_prev_macro_[plan_index]);
        }
      }
      ++count;
    }
  }
}

void LevelCoupler::ApplyProlongation(std::size_t begin, std::size_t end,
                                     bool half_time) {
  for (std::size_t i = begin; i < end; ++i) {
    const CouplingPoint& pt = plan_[i];
    if (!OwnsFine(pt, my_rank_)) {
      continue;
    }

    MacroState macro{};
    const bool owns_coarse = OwnsCoarse(pt, my_rank_);

    if (owns_coarse) {
      ReadMacro(blocks_[pt.coarse_id], pt.ci, pt.cj, pt.ck, &macro);
    } else {
      macro = recv_macro_[i];
    }

    if (half_time) {
      const MacroState& prev =
          owns_coarse ? prev_macro_[i] : recv_prev_macro_[i];
      macro.rho = 0.5 * (prev.rho + macro.rho);
      for (int d = 0; d < 3; ++d) {
        macro.u[d] = 0.5 * (prev.u[d] + macro.u[d]);
      }
      for (int p = 0; p < LatticeD3Q19::kQ; ++p) {
        macro.f_neq[p] = 0.5 * (prev.f_neq[p] + macro.f_neq[p]);
      }
    }

    WriteProlongation(blocks_[pt.fine_id], pt.fi, pt.fj, pt.fk, macro,
                      prolong_scale_);

    if (!half_time) {
      if (owns_coarse) {
        ReadMacro(blocks_[pt.coarse_id], pt.ci, pt.cj, pt.ck, &prev_macro_[i]);
      } else {
        prev_macro_[i] = macro;
      }
    }
  }
}

void LevelCoupler::apply_half_time(int coarse_level) {
  ExchangeCoarseMacros(coarse_level, true);
  ApplyProlongation(LevelRangeBegin(coarse_level), LevelRangeEnd(coarse_level),
                    true);
}

void LevelCoupler::apply_full_time(int coarse_level) {
  ExchangeCoarseMacros(coarse_level, false);
  ApplyProlongation(LevelRangeBegin(coarse_level), LevelRangeEnd(coarse_level),
                    false);
}

void LevelCoupler::restrict(int fine_level) {
  for (std::size_t i = 0; i < plan_.size(); ++i) {
    const CouplingPoint& pt = plan_[i];
    if (pt.fine_level != fine_level) {
      continue;
    }
    // Restriction only when both sides are local (cross-rank restrict is a
    // separate exchange path; ghost ids must not be treated as local).
    if (!OwnsCoarse(pt, my_rank_) || !OwnsFine(pt, my_rank_)) {
      continue;
    }

    MacroState fine_macro{};
    ReadMacro(blocks_[pt.fine_id], pt.fi, pt.fj, pt.fk, &fine_macro);

    MacroState neighbor_sum{};
    int neighbor_count = 0;
    using D = olb::descriptors::D3Q19<>;
    for (int jPop = 1; jPop < D::q; ++jPop) {
      const int nix = pt.fi + olb::descriptors::c<D>(jPop, 0);
      const int niy = pt.fj + olb::descriptors::c<D>(jPop, 1);
      const int niz = pt.fk + olb::descriptors::c<D>(jPop, 2);
      if (nix < 0 || nix >= nx_ || niy < 0 || niy >= ny_ || niz < 0 ||
          niz >= nz_) {
        continue;
      }
      MacroState nmacro{};
      ReadMacro(blocks_[pt.fine_id], nix, niy, niz, &nmacro);
      for (int p = 0; p < LatticeD3Q19::kQ; ++p) {
        neighbor_sum.f_neq[p] += nmacro.f_neq[p];
      }
      ++neighbor_count;
    }

    WriteRestriction(blocks_[pt.coarse_id], pt.ci, pt.cj, pt.ck, fine_macro,
                     neighbor_sum, neighbor_count, restrict_scale_);
  }
}

}  // namespace octlb
