#ifndef OCTLB_SRC_SOLVER_FIELD_BLOCK_ITERATOR_H_
#define OCTLB_SRC_SOLVER_FIELD_BLOCK_ITERATOR_H_

#include "src/common/types.h"

namespace octlb {

/** Iterates OctantId values 0…num_octants-1.
 *
 * Level-agnostic: does not know which octants belong to which level.
 * TimeLoop is responsible for filtering by level using OctreeForest.
 */
class BlockIterator {
 public:
  class Iterator {
   public:
    explicit Iterator(OctantId current) : current_(current) {}

    OctantId operator*() const { return current_; }

    Iterator& operator++() {
      ++current_;
      return *this;
    }

    bool operator!=(const Iterator& other) const {
      return current_ != other.current_;
    }

   private:
    OctantId current_;
  };

  explicit BlockIterator(label num_octants) : num_octants_(num_octants) {}

  Iterator begin() const { return Iterator(0); }

  Iterator end() const { return Iterator(num_octants_); }

 private:
  label num_octants_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_FIELD_BLOCK_ITERATOR_H_
