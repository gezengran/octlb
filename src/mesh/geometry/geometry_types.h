#ifndef OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_TYPES_H_
#define OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_TYPES_H_

#include <cstdint>
#include <string>
#include <vector>

#include "src/mesh/geometry/triangle_soup.h"

namespace octlb {

enum class GeometryPartRole : std::uint8_t {
  kInternalChannel,
  kExternalObstacle,
};

enum class MaterialKind : std::uint8_t {
  kFluid,
  kSolid,
  kBoundary,
};

struct GeometryPart {
  TriangleSoup soup;
  /** Inner cavity surface for \c kInternalChannel (thin-wall duct). */
  TriangleSoup inner_cavity_soup;
  GeometryPartRole role = GeometryPartRole::kExternalObstacle;
  int priority = 0;
  std::string name;
};

struct GeometryAssembly {
  std::vector<GeometryPart> parts;
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_TYPES_H_
