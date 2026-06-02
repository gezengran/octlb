#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "src/common/types.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/block_iterator.h"

namespace octlb {
namespace {

constexpr label kNumOctants = 5;

}  // namespace

TEST(BlockCollection, IntFactoryStoresIdTimesTen) {
  BlockCollection<int> collection(
      kNumOctants, [](OctantId id) { return static_cast<int>(id) * 10; });

  for (OctantId id = 0; id < kNumOctants; ++id) {
    EXPECT_EQ(collection[id], static_cast<int>(id) * 10);
  }
}

TEST(BlockCollection, StringFactoryStoresStringifiedId) {
  BlockCollection<std::string> collection(
      kNumOctants, [](OctantId id) { return std::to_string(id); });

  for (OctantId id = 0; id < kNumOctants; ++id) {
    EXPECT_EQ(collection[id], std::to_string(id));
  }
}

TEST(BlockCollection, SizeMatchesConstructionCount) {
  BlockCollection<int> collection(
      kNumOctants, [](OctantId id) { return static_cast<int>(id); });
  EXPECT_EQ(collection.size(), kNumOctants);
}

TEST(BlockIterator, YieldsSequentialOctantIds) {
  std::vector<OctantId> seen;
  for (OctantId id : BlockIterator(kNumOctants)) {
    seen.push_back(id);
  }

  ASSERT_EQ(seen.size(), static_cast<std::size_t>(kNumOctants));
  for (label i = 0; i < kNumOctants; ++i) {
    EXPECT_EQ(seen[static_cast<std::size_t>(i)], i);
  }
}

TEST(BlockIterator, EmptyRangeProducesNoIds) {
  std::vector<OctantId> seen;
  for (OctantId id : BlockIterator(0)) {
    seen.push_back(id);
  }
  EXPECT_TRUE(seen.empty());
}

}  // namespace octlb
