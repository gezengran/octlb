#ifndef OCTLB_SRC_SOLVER_FIELD_FACE_ITERATOR_H_
#define OCTLB_SRC_SOLVER_FIELD_FACE_ITERATOR_H_

#include <vector>

#include "src/mesh/topology/face_pair_list.h"

namespace octlb {

/** Iterates CoarseFineFace entries from a FacePairList (T06 LevelCoupler input). */
class FaceIterator {
 public:
  class Iterator {
   public:
    Iterator(const CoarseFineFace* ptr, const CoarseFineFace* end)
        : ptr_(ptr), end_(end) {}

    const CoarseFineFace& operator*() const { return *ptr_; }

    Iterator& operator++() {
      ++ptr_;
      return *this;
    }

    bool operator!=(const Iterator& other) const { return ptr_ != other.ptr_; }

   private:
    const CoarseFineFace* ptr_;
    const CoarseFineFace* end_;
  };

  explicit FaceIterator(const FacePairList& faces)
      : faces_(faces.coarse_fine_faces()) {}

  Iterator begin() const {
    return Iterator(faces_.data(), faces_.data() + faces_.size());
  }

  Iterator end() const {
    return Iterator(faces_.data() + faces_.size(),
                    faces_.data() + faces_.size());
  }

  std::size_t size() const { return faces_.size(); }

 private:
  const std::vector<CoarseFineFace>& faces_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_FIELD_FACE_ITERATOR_H_
