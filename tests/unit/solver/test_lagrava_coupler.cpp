#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>
#include <set>
#include <unordered_set>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/level_coupler.h"
#include "dynamics/lbm.h"

namespace octlb {
namespace {

constexpr int kN = 4;
constexpr double kOmega = 1.0;

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

OctreeForest MakeCenterRefinedForest(int extra_levels) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);
  OctantId center = FindCenterOctant(forest);
  for (int lvl = 0; lvl < extra_levels; ++lvl) {
    forest.refine([center](OctantId id) { return id == center; }, 3 + lvl);
    center = FindCenterOctant(forest);
  }
  forest.balance();
  forest.partition();
  return forest;
}

BlockCollection<LatticeD3Q19> MakeUniformBlocks(const OctreeForest& forest,
                                                double rho0 = 1.0) {
  return BlockCollection<LatticeD3Q19>(
      forest.local_num_octants(), [rho0](OctantId) {
        LatticeD3Q19 lat(kN, kN, kN, 1);
        const double u0[3] = {0.0, 0.0, 0.0};
        lat.initialize(rho0, u0);
        return lat;
      });
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

int NormalAxis(FaceDir normal) { return static_cast<int>(normal) / 2; }

bool IsInterfaceLayer(FaceDir normal, int ci, int cj, int ck) {
  const int axis = NormalAxis(normal);
  const bool is_max = (static_cast<int>(normal) % 2) == 1;
  const int idx[3] = {ci, cj, ck};
  return idx[axis] == (is_max ? kN - 1 : 0);
}

double TotalMass(BlockCollection<LatticeD3Q19>& blocks) {
  double sum = 0.0;
  for (label id = 0; id < blocks.size(); ++id) {
    auto& lat = blocks[static_cast<OctantId>(id)];
    for (int ix = 0; ix < kN; ++ix) {
      for (int iy = 0; iy < kN; ++iy) {
        for (int iz = 0; iz < kN; ++iz) {
          auto cell = lat.get(ix, iy, iz);
          double rho = 1.0;
          double u[3] = {};
          cell.computeRhoU(rho, u);
          sum += rho;
        }
      }
    }
  }
  return sum;
}

}  // namespace

TEST(LagravaCoupler, CouplingPlan_MatchesTopology) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest = MakeCenterRefinedForest(1);
  const FacePairList pairs(forest);
  auto blocks = MakeUniformBlocks(forest);
  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);

  const auto& plan = coupler.coupling_plan();
  ASSERT_GT(plan.size(), 0u);

  std::unordered_set<CouplingKey, CouplingKeyHash> seen;
  for (const CouplingPoint& pt : plan) {
    CouplingKey key{pt.coarse_id, pt.ci, pt.cj, pt.ck,
                    pt.fine_id,    pt.fi, pt.fj, pt.fk};
    EXPECT_TRUE(seen.insert(key).second) << "duplicate coupling point";
  }

  for (const CoarseFineFace& face : pairs.coarse_fine_faces()) {
    for (int slot = 0; slot < 4; ++slot) {
      const OctantId expected_fine = face.fine_ids[slot];
      bool found = false;
      for (const CouplingPoint& pt : plan) {
        if (pt.coarse_id == face.coarse_id && pt.fine_id == expected_fine) {
          found = true;
          EXPECT_EQ(pt.fine_slot, slot);
        }
      }
      EXPECT_TRUE(found) << "missing plan entries for fine slot " << slot;
    }
  }
}

TEST(LagravaCoupler, CouplingPlan_InterfaceIndicesInRange) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest = MakeCenterRefinedForest(1);
  const FacePairList pairs(forest);
  auto blocks = MakeUniformBlocks(forest);
  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);

  for (const CouplingPoint& pt : coupler.coupling_plan()) {
    EXPECT_GE(pt.ci, 0);
    EXPECT_GE(pt.cj, 0);
    EXPECT_GE(pt.ck, 0);
    EXPECT_LT(pt.ci, kN);
    EXPECT_LT(pt.cj, kN);
    EXPECT_LT(pt.ck, kN);
    EXPECT_GE(pt.fi, 0);
    EXPECT_GE(pt.fj, 0);
    EXPECT_GE(pt.fk, 0);
    EXPECT_LT(pt.fi, kN);
    EXPECT_LT(pt.fj, kN);
    EXPECT_LT(pt.fk, kN);
    EXPECT_TRUE(IsInterfaceLayer(pt.normal, pt.ci, pt.cj, pt.ck));
  }
}

TEST(LagravaCoupler, OneRank_CoarseFine_InterfaceContinuous) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest = MakeCenterRefinedForest(1);
  const FacePairList pairs(forest);
  auto blocks = MakeUniformBlocks(forest, 1.0);
  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);

  const int coarse_level = 1;
  coupler.apply_full_time(coarse_level);

  for (const CouplingPoint& pt : coupler.coupling_plan()) {
    if (pt.coarse_level != coarse_level) {
      continue;
    }
    auto coarse = blocks[pt.coarse_id].get(pt.ci, pt.cj, pt.ck);
    auto fine = blocks[pt.fine_id].get(pt.fi, pt.fj, pt.fk);
    double rho_c = 0.0, rho_f = 0.0;
    double u_c[3] = {}, u_f[3] = {};
    coarse.computeRhoU(rho_c, u_c);
    fine.computeRhoU(rho_f, u_f);
    EXPECT_NEAR(rho_c, rho_f, 1e-6);
    EXPECT_NEAR(u_c[0], u_f[0], 1e-6);
    EXPECT_NEAR(u_c[1], u_f[1], 1e-6);
    EXPECT_NEAR(u_c[2], u_f[2], 1e-6);
  }
}

TEST(LagravaCoupler, OneRank_ProlongRestrict_MassConserved) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest = MakeCenterRefinedForest(1);
  const FacePairList pairs(forest);
  auto blocks = MakeUniformBlocks(forest, 1.0);
  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);

  const double mass_before = TotalMass(blocks);
  const int coarse_level = 1;
  const int fine_level = 2;
  coupler.apply_full_time(coarse_level);
  coupler.restrict(fine_level);
  const double mass_after = TotalMass(blocks);

  EXPECT_NEAR(mass_before, mass_after, 1e-10);
}

TEST(LagravaCoupler, TwoRank_CoarseFine_MacroExchange) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "two-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 1);
  OctantId center = FindCenterOctant(forest);
  forest.refine([center](OctantId id) { return id == center; }, 2);
  forest.balance();
  forest.partition();

  const FacePairList pairs(forest);
  auto blocks = MakeUniformBlocks(forest, 1.0);

  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);

  bool has_cross = false;
  double expected_rho = 1.0;
  int rho_source_rank = -1;
  for (const CouplingPoint& pt : coupler.coupling_plan()) {
    const bool owns_coarse =
        pt.coarse_id >= 0 && pt.coarse_id < forest.local_num_octants();
    const bool owns_fine =
        pt.fine_id >= 0 && pt.fine_id < forest.local_num_octants();
    if (owns_coarse != owns_fine) {
      has_cross = true;
      if (owns_coarse) {
        auto cell = blocks[pt.coarse_id].get(pt.ci, pt.cj, pt.ck);
        expected_rho = 1.2;
        rho_source_rank = rank;
        const double u0[3] = {0.0, 0.0, 0.0};
        const double uSqr = 0.0;
        for (int iPop = 0; iPop < LatticeD3Q19::kQ; ++iPop) {
          cell[iPop] = olb::equilibrium<olb::descriptors::D3Q19<>>::secondOrder(
              iPop, expected_rho, u0, uSqr);
        }
      }
    }
  }

  int bcast_root = rho_source_rank;
  MPI_Allreduce(MPI_IN_PLACE, &bcast_root, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  int local_has = has_cross ? 1 : 0;
  int global_has = 0;
  MPI_Allreduce(&local_has, &global_has, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (global_has == 0) {
    GTEST_SKIP() << "partition did not split coarse/fine across ranks";
  }
  if (bcast_root >= 0) {
    MPI_Bcast(&expected_rho, 1, MPI_DOUBLE, bcast_root, MPI_COMM_WORLD);
  }

  const int coarse_level = 0;
  coupler.apply_full_time(coarse_level);

  for (const CouplingPoint& pt : coupler.coupling_plan()) {
    if (pt.coarse_level != coarse_level) {
      continue;
    }
    const bool owns_coarse =
        pt.coarse_id >= 0 && pt.coarse_id < forest.local_num_octants();
    const bool owns_fine =
        pt.fine_id >= 0 && pt.fine_id < forest.local_num_octants();
    if (!owns_fine || owns_coarse) {
      continue;
    }
    auto fine = blocks[pt.fine_id].get(pt.fi, pt.fj, pt.fk);
    double rho_f = 0.0;
    double u_f[3] = {};
    fine.computeRhoU(rho_f, u_f);
    EXPECT_NEAR(rho_f, expected_rho, 1e-5);
  }
}

TEST(LagravaCoupler, Rebuild_RecreateCoupler) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest = MakeCenterRefinedForest(1);
  FacePairList pairs(forest);
  auto blocks = MakeUniformBlocks(forest);
  LevelCoupler coupler(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);
  const std::size_t before = coupler.coupling_plan().size();

  pairs.rebuild(forest);
  LevelCoupler rebuilt(MPI_COMM_WORLD, pairs, forest, blocks, kN, kN, kN,
                       kOmega);
  EXPECT_EQ(rebuilt.coupling_plan().size(), before);
}

TEST(LagravaCoupler, CoarseFine_CommTagsUnique) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    GTEST_SKIP() << "two-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 1);
  OctantId center = FindCenterOctant(forest);
  forest.refine([center](OctantId id) { return id == center; }, 2);
  forest.balance();
  forest.partition();

  const FacePairList pairs(forest);
  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // Only slots this rank actually exchanges over MPI need a valid, unique tag
  // (LevelCoupler reads comm_tags[slot] solely for these). A coarse-owning rank
  // uses every slot with peer = remote_ranks[i] (cross-rank when fine_i is
  // remote); a fine-owning rank uses only its owned slot (fine_i local) with
  // peer = coarse_remote_rank. Non-owned slots are not exchanged here, so their
  // comm_tag is unused and left 0 by FacePairList.
  std::set<std::pair<int, int>> tags;
  for (const CoarseFineFace& face : pairs.coarse_fine_faces()) {
    const bool owns_coarse = face.coarse_remote_rank == my_rank;
    for (int i = 0; i < 4; ++i) {
      int peer = -1;
      if (owns_coarse) {
        if (face.remote_ranks[i] == my_rank) {
          continue;  // same-rank coarse-fine, no MPI
        }
        peer = face.remote_ranks[i];
      } else {
        if (face.remote_ranks[i] != my_rank) {
          continue;  // fine_i not owned by this rank -> slot unused here
        }
        peer = face.coarse_remote_rank;  // != my_rank since coarse is a ghost
      }
      EXPECT_GT(face.comm_tags[i], 0);
      const std::pair<int, int> key{peer, face.comm_tags[i]};
      EXPECT_TRUE(tags.insert(key).second)
          << "duplicate coarse-fine comm tag (peer,tag) on rank " << my_rank;
    }
  }
}

}  // namespace octlb
