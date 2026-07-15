#include <gtest/gtest.h>
#include <mpi.h>

#include <set>
#include <unordered_map>
#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/face_iterator.h"
#include "src/solver/field/face_packable.h"
#include "src/solver/field/ghost_schedule.h"
#include "src/solver/lbm/block_lattice.h"

namespace octlb {
namespace {

constexpr int kNx = 4;
constexpr int kNy = 4;
constexpr int kNz = 4;

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

struct DummyBlock {
  using face_value_t = int;

  int id = 0;

  static int face_buffer_count(int nx, int ny, int nz, FaceDir dir) {
    return FaceBufferCount(nx, ny, nz, dir);
  }

  void pack_face(FaceDir dir, int* buf, int count) const {
    for (int i = 0; i < count; ++i) {
      buf[i] = Signature(dir, i);
    }
  }

  void unpack_face(FaceDir dir, const int* buf, int count) {
    ghost_[static_cast<int>(dir)].assign(buf, buf + count);
  }

  // ② edge-ghost (test double). No test inspects edge ghosts on DummyBlock, so
  // unpack discards; pack emits a signature so the exchange path is exercised.
  static int edge_buffer_count(int nx, int ny, int nz, FaceDir d1, FaceDir d2) {
    return EdgeBufferCount(nx, ny, nz, d1, d2);
  }
  void pack_edge(FaceDir d1, FaceDir d2, int* buf, int count) const {
    for (int i = 0; i < count; ++i) {
      buf[i] = id * 100000 + static_cast<int>(d1) * 1000 +
               static_cast<int>(d2) * 10 + i;
    }
  }
  void unpack_edge(FaceDir /*d1*/, FaceDir /*d2*/, const int* /*buf*/,
                   int /*count*/) {}

  const std::vector<int>& ghost(FaceDir dir) const {
    return ghost_[static_cast<int>(dir)];
  }

  int Signature(FaceDir dir, int index = 0) const {
    return id * 1000 + static_cast<int>(dir) * 100 + index;
  }

 private:
  std::array<std::vector<int>, 6> ghost_{};
};

OctreeForest MakeUniformForest(int max_level) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, max_level);
  forest.balance();
  forest.partition();
  return forest;
}

const SameLevelFace* FindFace(const FacePairList& pairs, OctantId local_id,
                              FaceDir dir) {
  for (const SameLevelFace& face : pairs.same_level_faces()) {
    if (face.local_id == local_id && face.dir == dir) {
      return &face;
    }
  }
  return nullptr;
}

bool AllRanksTrue(bool local_ok) {
  int flag = local_ok ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&flag, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return global == 1;
}

void ExpectGhostMatchesNeighborPack(const DummyBlock& receiver, FaceDir recv_dir,
                                  const DummyBlock& sender, FaceDir send_dir) {
  const std::vector<int>& ghost = receiver.ghost(recv_dir);
  ASSERT_EQ(ghost.size(),
            static_cast<std::size_t>(
                DummyBlock::face_buffer_count(kNx, kNy, kNz, recv_dir)));
  for (int i = 0; i < static_cast<int>(ghost.size()); ++i) {
    EXPECT_EQ(ghost[static_cast<std::size_t>(i)], sender.Signature(send_dir, i));
  }
}

}  // namespace

TEST(GhostSchedule, EmptySchedule_NoOp) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();

  const FacePairList pairs(forest);
  label local_faces = static_cast<label>(pairs.same_level_faces().size());
  label global_faces = 0;
  MPI_Allreduce(&local_faces, &global_faces, 1, MPI_INT32_T, MPI_SUM,
                MPI_COMM_WORLD);
  EXPECT_EQ(global_faces, 0);

  if (forest.local_num_octants() == 0) {
    return;
  }

  BlockCollection<DummyBlock> blocks(
      forest.local_num_octants(), [](OctantId id) {
        DummyBlock block;
        block.id = static_cast<int>(id);
        return block;
      });

  const int before = blocks[0].Signature(FaceDir::kXMin, 0);
  GhostSchedule schedule(MPI_COMM_WORLD, pairs, blocks, kNx, kNy, kNz);
  schedule.exchange();
  EXPECT_EQ(blocks[0].Signature(FaceDir::kXMin, 0), before);
  EXPECT_TRUE(blocks[0].ghost(FaceDir::kXMin).empty());
}

TEST(GhostSchedule, OneRank_TwoAdjacentBlocks) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  OctreeForest forest = MakeUniformForest(1);
  const FacePairList pairs(forest);
  const SameLevelFace* ab = FindFace(pairs, 0, FaceDir::kXMax);
  if (ab == nullptr || ab->remote_rank != rank) {
    ab = FindFace(pairs, 1, FaceDir::kXMin);
  }
  if (ab == nullptr || ab->remote_rank != rank) {
    return;
  }
  const OctantId a_id = ab->local_id;
  const OctantId b_id = ab->remote_id;

  BlockCollection<DummyBlock> blocks(
      forest.local_num_octants(), [](OctantId id) {
        DummyBlock block;
        block.id = static_cast<int>(id);
        return block;
      });

  GhostSchedule schedule(MPI_COMM_WORLD, pairs, blocks, kNx, kNy, kNz);
  schedule.exchange();

  ExpectGhostMatchesNeighborPack(blocks[a_id], ab->dir, blocks[b_id],
                                 OppositeFace(ab->dir));
  ExpectGhostMatchesNeighborPack(blocks[b_id], OppositeFace(ab->dir), blocks[a_id],
                                 ab->dir);
}

TEST(GhostSchedule, TwoRank_OneBlockEach) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ASSERT_EQ(size, 2);

  OctreeForest forest = MakeUniformForest(1);
  const FacePairList pairs(forest);

  const SameLevelFace* cross = nullptr;
  for (const SameLevelFace& face : pairs.same_level_faces()) {
    if (face.remote_rank != rank) {
      cross = &face;
      break;
    }
  }
  ASSERT_NE(cross, nullptr);

  BlockCollection<DummyBlock> blocks(
      forest.local_num_octants(), [rank](OctantId id) {
        DummyBlock block;
        block.id = rank * 100 + static_cast<int>(id);
        return block;
      });

  GhostSchedule schedule(MPI_COMM_WORLD, pairs, blocks, kNx, kNy, kNz);
  schedule.exchange();

  const int count =
      DummyBlock::face_buffer_count(kNx, kNy, kNz, cross->dir);
  const std::vector<int> ghost_received = blocks[cross->local_id].ghost(cross->dir);

  std::vector<int> outbound(count);
  blocks[cross->local_id].pack_face(cross->dir, outbound.data(), count);

  std::vector<int> peer_outbound(count);
  if (rank == 0) {
    MPI_Recv(peer_outbound.data(), count, MPI_INT, 1, 9100, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (int i = 0; i < count; ++i) {
      EXPECT_EQ(ghost_received[static_cast<std::size_t>(i)],
                peer_outbound[static_cast<std::size_t>(i)]);
    }
  } else {
    MPI_Send(outbound.data(), count, MPI_INT, 0, 9100, MPI_COMM_WORLD);
  }
}

TEST(GhostSchedule, MixedLocalAndRemoteFaces) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ASSERT_GE(size, 2);

  OctreeForest forest = MakeUniformForest(2);
  const FacePairList pairs(forest);

  bool has_local = false;
  bool has_remote = false;
  for (const SameLevelFace& face : pairs.same_level_faces()) {
    if (face.remote_rank == rank) {
      has_local = true;
    } else {
      has_remote = true;
    }
  }
  if (!AllRanksTrue(has_local && has_remote)) {
    return;
  }

  BlockCollection<DummyBlock> blocks(
      forest.local_num_octants(), [rank](OctantId id) {
        DummyBlock block;
        block.id = rank * 1000 + static_cast<int>(id);
        return block;
      });

  GhostSchedule schedule(MPI_COMM_WORLD, pairs, blocks, kNx, kNy, kNz);
  schedule.exchange();

  for (const SameLevelFace& face : pairs.same_level_faces()) {
    DummyBlock& local = blocks[face.local_id];
    if (face.remote_rank == rank) {
      ExpectGhostMatchesNeighborPack(local, face.dir, blocks[face.remote_id],
                                     OppositeFace(face.dir));
    } else {
      ASSERT_FALSE(local.ghost(face.dir).empty());
    }
  }
}

TEST(GhostSchedule, CornerGhostConsistent) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  OctreeForest forest = MakeUniformForest(1);
  ASSERT_GE(forest.local_num_octants(), 2);

  const FacePairList pairs(forest);
  BlockCollection<DummyBlock> blocks(
      forest.local_num_octants(), [](OctantId id) {
        DummyBlock block;
        block.id = static_cast<int>(id);
        return block;
      });

  GhostSchedule schedule(MPI_COMM_WORLD, pairs, blocks, kNx, kNy, kNz);
  schedule.exchange();

  OctantId corner_id = -1;
  const SameLevelFace* x_face = nullptr;
  const SameLevelFace* y_face = nullptr;
  for (label id = 0; id < forest.local_num_octants(); ++id) {
    const SameLevelFace* xf = FindFace(pairs, static_cast<OctantId>(id),
                                      FaceDir::kXMax);
    const SameLevelFace* yf = FindFace(pairs, static_cast<OctantId>(id),
                                      FaceDir::kYMax);
    if (xf != nullptr && yf != nullptr && xf->remote_rank == rank &&
        yf->remote_rank == rank) {
      corner_id = static_cast<OctantId>(id);
      x_face = xf;
      y_face = yf;
      break;
    }
  }
  if (corner_id < 0) {
    return;
  }

  ExpectGhostMatchesNeighborPack(blocks[corner_id], FaceDir::kXMax,
                                 blocks[x_face->remote_id],
                                 OppositeFace(x_face->dir));
  ExpectGhostMatchesNeighborPack(blocks[corner_id], FaceDir::kYMax,
                                 blocks[y_face->remote_id],
                                 OppositeFace(y_face->dir));
}

TEST(GhostSchedule, CommTagUnique) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ASSERT_GE(size, 2);

  OctreeForest forest = MakeUniformForest(2);
  const FacePairList pairs(forest);

  std::set<std::pair<int, int>> keys;
  for (const SameLevelFace& face : pairs.same_level_faces()) {
    if (face.remote_rank == rank) {
      continue;
    }
    const std::pair<int, int> key{face.remote_rank, face.comm_tag};
    EXPECT_TRUE(keys.insert(key).second)
        << "duplicate (remote_rank, comm_tag) " << key.first << ","
        << key.second;
  }
}

TEST(GhostSchedule, OppositeDirPairing) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ASSERT_EQ(size, 2);

  OctreeForest forest = MakeUniformForest(1);
  const FacePairList pairs(forest);

  struct FaceReport {
    int comm_tag;
    int dir;
    int local_id;
  };

  FaceReport local_report{-1, -1, -1};
  for (const SameLevelFace& face : pairs.same_level_faces()) {
    if (face.remote_rank != rank) {
      local_report = {face.comm_tag, static_cast<int>(face.dir),
                      static_cast<int>(face.local_id)};
      break;
    }
  }
  if (!AllRanksTrue(local_report.comm_tag >= 0)) {
    return;
  }

  FaceReport peer_report{-1, -1, -1};
  MPI_Sendrecv(&local_report, sizeof(FaceReport), MPI_BYTE, 1 - rank, 9001,
               &peer_report, sizeof(FaceReport), MPI_BYTE, 1 - rank, 9001,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  EXPECT_EQ(local_report.comm_tag, peer_report.comm_tag)
      << "symmetric FacePairList comm_tag must match";
  EXPECT_EQ(static_cast<FaceDir>(local_report.dir),
            OppositeFace(static_cast<FaceDir>(peer_report.dir)));
}

TEST(GhostSchedule, FaceIterator_MatchesCoarseFineCount) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.refine([](OctantId) { return true; }, 2);
  OctantId center = 0;
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
      center = static_cast<OctantId>(i);
    }
  }
  forest.refine([center](OctantId id) { return id == center; }, 3);
  forest.balance();
  forest.partition();

  const FacePairList pairs(forest);
  FaceIterator it(pairs);
  std::size_t count = 0;
  for (const CoarseFineFace& face : it) {
    ++count;
    for (int i = 0; i < 4; ++i) {
      EXPECT_GE(face.fine_ids[i], 0);
    }
    (void)face.coarse_id;
  }
  EXPECT_EQ(count, pairs.coarse_fine_faces().size());
}

TEST(GhostSchedule, TwoRank_BlockLattice_FaceValuesMatch) {
  using Lattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ASSERT_EQ(size, 2);

  OctreeForest forest = MakeUniformForest(1);
  const FacePairList pairs(forest);

  const SameLevelFace* cross = nullptr;
  for (const SameLevelFace& face : pairs.same_level_faces()) {
    if (face.remote_rank != rank) {
      cross = &face;
      break;
    }
  }
  ASSERT_NE(cross, nullptr);

  BlockCollection<Lattice> blocks(forest.local_num_octants(), [](OctantId id) {
    Lattice lat(kNx, kNy, kNz, 1);
    const double rho0 = 1.0 + 0.01 * static_cast<double>(id);
    const double u0[3] = {0.0, 0.0, 0.0};
    lat.initialize(rho0, u0);
    return lat;
  });

  GhostSchedule<Lattice> schedule(MPI_COMM_WORLD, pairs, blocks, kNx, kNy,
                                  kNz);
  schedule.exchange();

  const int n = Lattice::face_buffer_count(kNx, kNy, kNz, cross->dir);
  std::vector<double> ghost_values(n);
  blocks[cross->local_id].read_ghost_face(cross->dir, ghost_values.data(), n);

  std::vector<double> outbound(n);
  blocks[cross->local_id].pack_face(cross->dir, outbound.data(), n);

  std::vector<double> peer_outbound(n);
  if (rank == 0) {
    MPI_Recv(peer_outbound.data(), n, MPI_DOUBLE, 1, 9200, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (int i = 0; i < n; ++i) {
      EXPECT_NEAR(ghost_values[static_cast<std::size_t>(i)],
                  peer_outbound[static_cast<std::size_t>(i)], 1e-12);
    }
  } else {
    MPI_Send(outbound.data(), n, MPI_DOUBLE, 0, 9200, MPI_COMM_WORLD);
  }
}

}  // namespace octlb
