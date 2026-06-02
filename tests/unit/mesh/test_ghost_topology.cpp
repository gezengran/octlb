#include <gtest/gtest.h>
#include <mpi.h>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/forest/octree_forest_access.h"
#include "src/mesh/topology/face_pair_list.h"

namespace octlb {
namespace {

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

int GhostOwnerRank(const p8est_ghost_t* ghost, p4est_locidx_t ghost_index) {
  for (int r = 0; r < ghost->mpisize; ++r) {
    if (ghost_index >= ghost->proc_offsets[r] &&
        ghost_index < ghost->proc_offsets[r + 1]) {
      return r;
    }
  }
  return -1;
}

void BuildGhostTopologyForest(OctreeForest* forest) {
  forest->refine([](OctantId) { return true; }, 1);
  forest->refine([](OctantId id) { return (id % 3) == 0; }, 2);
  forest->refine([](OctantId id) { return (id % 5) == 0; }, 3);
  forest->balance();
  forest->partition();
}

}  // namespace

TEST(GhostTopology, FacePairRemoteRanksMatchGhostLayer) {
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  BuildGhostTopologyForest(&forest);

  const FacePairList pairs(forest);
  ASSERT_GT(pairs.same_level_faces().size(), 0u);
  ASSERT_GT(pairs.coarse_fine_faces().size(), 0u);

  p8est_ghost_t* ghost = MeshForestAccess::Ghost(forest);
  ASSERT_NE(ghost, nullptr);

  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  for (const SameLevelFace& face : pairs.same_level_faces()) {
    if (face.remote_rank == my_rank) {
      continue;
    }
    const int owner = GhostOwnerRank(ghost, face.remote_id);
    EXPECT_EQ(owner, face.remote_rank)
        << "same-level ghost index " << face.remote_id;
  }

  for (const CoarseFineFace& face : pairs.coarse_fine_faces()) {
    for (int i = 0; i < 4; ++i) {
      if (face.remote_ranks[i] == my_rank) {
        continue;
      }
      const int owner = GhostOwnerRank(ghost, face.fine_ids[i]);
      EXPECT_EQ(owner, face.remote_ranks[i])
          << "coarse-fine ghost index " << face.fine_ids[i];
    }
  }
}

}  // namespace octlb
