#ifndef OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_BUILD_ERROR_H_
#define OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_BUILD_ERROR_H_

#include <stdexcept>
#include <string>

namespace octlb {

class GeometryBuildError : public std::runtime_error {
 public:
  explicit GeometryBuildError(const std::string& message)
      : std::runtime_error(message) {}
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_BUILD_ERROR_H_
