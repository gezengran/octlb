#ifndef OCTLB_SRC_SOLVER_FIELD_FACE_PACKABLE_H_
#define OCTLB_SRC_SOLVER_FIELD_FACE_PACKABLE_H_

#include <algorithm>
#include <concepts>

#include "src/common/types.h"

namespace octlb {

inline FaceDir OppositeFace(FaceDir dir) {
  return static_cast<FaceDir>(static_cast<int>(dir) ^ 1);
}

inline int FaceBufferCount(int nx, int ny, int nz, FaceDir dir) {
  switch (dir) {
    case FaceDir::kXMin:
    case FaceDir::kXMax:
      return ny * nz;
    case FaceDir::kYMin:
    case FaceDir::kYMax:
      return nx * nz;
    case FaceDir::kZMin:
    case FaceDir::kZMax:
      return nx * ny;
  }
  return 0;
}

inline int MaxFaceBufferCount(int nx, int ny, int nz,
                              int (*count_fn)(int, int, int, FaceDir)) {
  int max_count = 0;
  for (int d = 0; d < 6; ++d) {
    max_count = std::max(
        max_count, count_fn(nx, ny, nz, static_cast<FaceDir>(d)));
  }
  return max_count;
}

template <typename T>
concept FacePackable = requires(const T& ct, T& t, FaceDir dir,
                                typename T::face_value_t* buf, int n) {
  typename T::face_value_t;
  { ct.pack_face(dir, buf, n) } -> std::same_as<void>;
  { t.unpack_face(dir, buf, n) } -> std::same_as<void>;
  { T::face_buffer_count(0, 0, 0, FaceDir::kXMin) } -> std::same_as<int>;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_FIELD_FACE_PACKABLE_H_
