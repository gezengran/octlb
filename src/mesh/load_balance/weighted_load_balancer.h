#ifndef OCTLB_SRC_MESH_LOAD_BALANCE_WEIGHTED_LOAD_BALANCER_H_
#define OCTLB_SRC_MESH_LOAD_BALANCE_WEIGHTED_LOAD_BALANCER_H_

#include <functional>

#include "src/common/types.h"

namespace octlb {

class OctreeForest;

/** Returns weight(octant) = 2^level for \c OctreeForest::partition. */
std::function<int(OctantId)> make_level_weight_fn(const OctreeForest& forest);

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_LOAD_BALANCE_WEIGHTED_LOAD_BALANCER_H_
