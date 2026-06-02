#ifndef OCTLB_SRC_SOLVER_FIELD_BLOCK_COLLECTION_H_
#define OCTLB_SRC_SOLVER_FIELD_BLOCK_COLLECTION_H_

#include <functional>
#include <vector>

#include "src/common/types.h"

namespace octlb {

/** Stores one T per local octant, keyed by OctantId (0…size-1).
 *
 * T need not be default-constructible; factory(id) is called once
 * per octant at construction time.
 */
template <typename T>
class BlockCollection {
 public:
  BlockCollection(label num_octants, std::function<T(OctantId)> factory) {
    blocks_.reserve(static_cast<std::size_t>(num_octants));
    for (label i = 0; i < num_octants; ++i) {
      blocks_.push_back(factory(static_cast<OctantId>(i)));
    }
  }

  T& operator[](OctantId id) {
    return blocks_[static_cast<std::size_t>(id)];
  }

  const T& operator[](OctantId id) const {
    return blocks_[static_cast<std::size_t>(id)];
  }

  label size() const { return static_cast<label>(blocks_.size()); }

 private:
  std::vector<T> blocks_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_FIELD_BLOCK_COLLECTION_H_
