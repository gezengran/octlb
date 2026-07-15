#ifndef OCTLB_SRC_SOLVER_FIELD_FACE_PACKABLE_H_
#define OCTLB_SRC_SOLVER_FIELD_FACE_PACKABLE_H_

#include <algorithm>
#include <concepts>

#include "src/common/types.h"

namespace octlb {

inline FaceDir OppositeFace(FaceDir dir) {
  return static_cast<FaceDir>(static_cast<int>(dir) ^ 1);
}

// Axis index of a face direction: 0=x, 1=y, 2=z.
inline int FaceAxis(FaceDir dir) {
  return static_cast<int>(dir) / 2;  // kXMin/kXMax->0, kY..->1, kZ..->2
}

// Whether the face is the + (max) side of its axis.
inline bool IsMaxFace(FaceDir dir) {
  return (static_cast<int>(dir) & 1) == 1;  // kXMax, kYMax, kZMax
}

// The two face directions d1, d2 must be orthogonal (different axes). The edge
// runs along the remaining (third) axis; 0+1+2 == 3.
inline int EdgeAxis(FaceDir d1, FaceDir d2) {
  return 3 - FaceAxis(d1) - FaceAxis(d2);
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

// Edge ghost line length (cells along the third axis) for an nx*ny*nz block.
inline int EdgeBufferCount(int nx, int ny, int nz, FaceDir d1, FaceDir d2) {
  const int n[3] = {nx, ny, nz};
  return n[EdgeAxis(d1, d2)];
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

inline int MaxEdgeBufferCount(
    int nx, int ny, int nz,
    int (*count_fn)(int, int, int, FaceDir, FaceDir)) {
  int max_count = 0;
  for (int d1 = 0; d1 < 6; ++d1) {
    for (int d2 = 0; d2 < 6; ++d2) {
      const FaceDir f1 = static_cast<FaceDir>(d1);
      const FaceDir f2 = static_cast<FaceDir>(d2);
      if (FaceAxis(f1) == FaceAxis(f2)) continue;  // not orthogonal -> not an edge
      max_count = std::max(max_count, count_fn(nx, ny, nz, f1, f2));
    }
  }
  return max_count;
}

template <typename T>
concept FacePackable = requires(const T& ct, T& t, FaceDir dir,
                                typename T::face_value_t* buf, int n,
                                FaceDir e1, FaceDir e2) {
  typename T::face_value_t;
  { ct.pack_face(dir, buf, n) } -> std::same_as<void>;
  { t.unpack_face(dir, buf, n) } -> std::same_as<void>;
  { T::face_buffer_count(0, 0, 0, FaceDir::kXMin) } -> std::same_as<int>;
  // Edge ghost line (② fix): an edge is the intersection of two orthogonal
  // faces; pack reads the interior edge line, unpack writes the edge ghost.
  { ct.pack_edge(e1, e2, buf, n) } -> std::same_as<void>;
  { t.unpack_edge(e1, e2, buf, n) } -> std::same_as<void>;
  { T::edge_buffer_count(0, 0, 0, FaceDir::kXMin, FaceDir::kYMin) } ->
      std::same_as<int>;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_FIELD_FACE_PACKABLE_H_
