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
    EXPECT_EQ(assigner_a.tag_for(key), assigner_b.tag_for(key))
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

  std::unordered_map<uint64_t, int> seen;
  for (uint64_t key : faces) {
    const int tag = assigner.tag_for(key);
    EXPECT_GT(tag, 0);
    EXPECT_TRUE(seen.find(tag) == seen.end() || seen[tag] == static_cast<int>(key))
        << "tag " << tag << " reused for distinct faces in same peer";
    seen[tag] = static_cast<int>(key);
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
    const int tag = assigner.tag_for(key * 7919ULL);
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
    EXPECT_EQ(only_b.tag_for(key), b_and_c.tag_for(key))
        << "peer-B tag changed by presence of peer-C faces (old global-tag bug)";
  }
}

}  // namespace
}  // namespace octlb