#include "src/mesh/load_balance/weighted_load_balancer.h"

#include "src/mesh/forest/octree_forest.h"

namespace octlb {

std::function<int(OctantId)> make_level_weight_fn(const OctreeForest& forest) {
  return [&forest](OctantId id) {
    const int level = forest.quadrant_level(id);
    return 1 << level;
  };
}

}  // namespace octlb
