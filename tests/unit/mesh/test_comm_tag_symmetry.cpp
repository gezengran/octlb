#include <gtest/gtest.h>

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "src/common/comm_tag_assigner.h"

namespace octlb {
namespace {

// MPI_TAG_UB on MPICH (0x3FFFFFFF); used as a realistic upper bound in tests
// where face counts stay well below it, so per-peer indices stay unique and
// in range.
constexpr int kTagUb = 0x3FFFFFFF;

// Two endpoints of a cross-rank face both observe the same set of faces for a
// given peer, but p4est_iterate may visit them in different orders on each
// rank. The assigned tags must be order-independent so both ranks agree.
TEST(CommTagAssigner, TagsAreOrderIndependentAcrossInsertionOrder) {
  const std::vector<uint64_t> order_a = {17ULL, 5ULL, 99ULL, 42ULL, 8ULL};
  const std::vector<uint64_t> order_b = {99ULL, 8ULL, 17ULL, 42ULL, 5ULL};

  CommTagAssigner assigner_a(kTagUb);
  CommTagAssigner assigner_b(kTagUb);
  for (uint64_t key : order_a) {
    assigner_a.add_face(key, /*peer_rank=*/1);
  }
  for (uint64_t key : order_b) {
    assigner_b.add_face(key, /*peer_rank=*/1);
  }
  assigner_a.assign();
  assigner_b.assign();

  for (uint64_t key : order_a) {
    EXPECT_EQ(assigner_a.tag_for(key, /*peer_rank=*/1),
              assigner_b.tag_for(key, /*peer_rank=*/1))
        << "tag for face " << key << " differs across insertion order";
  }
}

// Within a single rank-pair, MPI matches messages by (src, dst, tag); two faces
// sharing the same pair must therefore carry distinct tags.
TEST(CommTagAssigner, TagsAreUniqueWithinAPeer) {
  CommTagAssigner assigner(kTagUb);
  const std::vector<uint64_t> faces = {3ULL, 1ULL, 4ULL, 1ULL, 5ULL, 9ULL};
  for (uint64_t key : faces) {
    assigner.add_face(key, /*peer_rank=*/2);
  }
  assigner.assign();

  std::unordered_map<int, uint64_t> tag_owner;
  for (uint64_t key : faces) {
    const int tag = assigner.tag_for(key, /*peer_rank=*/2);
    EXPECT_GT(tag, 0);
    const auto it = tag_owner.find(tag);
    EXPECT_TRUE(it == tag_owner.end() || it->second == key)
        << "tag " << tag << " reused for distinct faces in same peer";
    tag_owner[tag] = key;
  }
}

// Every assigned tag must be a legal MPI tag in [1, MPI_TAG_UB].
TEST(CommTagAssigner, TagsStayWithinTagUb) {
  const int small_ub = 64;
  CommTagAssigner assigner(small_ub);
  for (uint64_t key = 1; key <= 100; ++key) {
    assigner.add_face(key * 7919ULL, /*peer_rank=*/3);  // spread keys
  }
  assigner.assign();

  for (uint64_t key = 1; key <= 100; ++key) {
    const int tag = assigner.tag_for(key * 7919ULL, /*peer_rank=*/3);
    EXPECT_GE(tag, 1);
    EXPECT_LE(tag, small_ub);
  }
}

// Regression guard for the original bug: the old global per-rank tag_owner map
// let faces of one rank-pair rehash tags of another pair. Tags of a given pair
// must be independent of faces belonging to a different pair.
TEST(CommTagAssigner, PeerTagsAreIndependentOfOtherPeers) {
  const std::vector<uint64_t> peer_b_faces = {11ULL, 27ULL, 33ULL};
  const std::vector<uint64_t> peer_c_faces = {2ULL, 4ULL, 8ULL, 16ULL,
                                              32ULL, 64ULL, 128ULL};

  CommTagAssigner only_b(kTagUb);
  for (uint64_t key : peer_b_faces) {
    only_b.add_face(key, /*peer_rank=*/1);
  }
  only_b.assign();

  CommTagAssigner b_and_c(kTagUb);
  for (uint64_t key : peer_b_faces) {
    b_and_c.add_face(key, /*peer_rank=*/1);
  }
  for (uint64_t key : peer_c_faces) {
    b_and_c.add_face(key, /*peer_rank=*/5);
  }
  b_and_c.assign();

  for (uint64_t key : peer_b_faces) {
    EXPECT_EQ(only_b.tag_for(key, /*peer_rank=*/1),
              b_and_c.tag_for(key, /*peer_rank=*/1))
        << "peer-B tag changed by presence of peer-C faces (old global-tag bug)";
  }
}

// Regression: the same uint64 key under two peers (hash collision, or a
// same-rank "self" bucket colliding with a cross-rank face) must keep
// independent per-peer tags. A flat key→tag map overwrites and desynchronizes
// GhostSchedule face Waitall across the rank-pair (cylinder3d AMR hang).
TEST(CommTagAssigner, SameKeyUnderTwoPeersKeepsIndependentTags) {
  constexpr uint64_t kShared = 4286180941653651009ULL;
  CommTagAssigner assigner(kTagUb);
  // Peer 3: shared key sorts to index 2 → tag 3 (with two smaller keys).
  assigner.add_face(412347995506250705ULL, /*peer_rank=*/3);
  assigner.add_face(3281948207910144494ULL, /*peer_rank=*/3);
  assigner.add_face(kShared, /*peer_rank=*/3);
  // Peer 2 (self): shared key alone → tag 1. Flat map would overwrite peer-3.
  assigner.add_face(kShared, /*peer_rank=*/2);
  assigner.assign();

  EXPECT_EQ(assigner.tag_for(kShared, /*peer_rank=*/3), 3);
  EXPECT_EQ(assigner.tag_for(kShared, /*peer_rank=*/2), 1);
}

}  // namespace
}  // namespace octlb
