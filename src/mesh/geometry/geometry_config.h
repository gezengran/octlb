#ifndef OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_CONFIG_H_
#define OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_CONFIG_H_

namespace octlb {

struct GeometryConfig {
  int max_level = 6;
  int resolve_surface_times = 3;
  scalar bound_width = 0.0;
  scalar wake_length = 0.0;
  int wake_direction = 0;  // 0=x, 1=y, 2=z
  int cell_width = 4;
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_GEOMETRY_GEOMETRY_CONFIG_H_
