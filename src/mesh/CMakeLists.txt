add_library(octlb_mesh_stl STATIC
  geometry/triangle_soup.cpp
  io/stl_reader/stl_reader.cpp
)

target_include_directories(octlb_mesh_stl
  PUBLIC
    ${CMAKE_SOURCE_DIR}
)

target_compile_features(octlb_mesh_stl PUBLIC cxx_std_17)

target_sources(octlb_mesh PRIVATE
  forest/octree_forest.cpp
  topology/face_pair_list.cpp
  load_balance/weighted_load_balancer.cpp
  geometry/material_field.cpp
  geometry/cgal_surface_mesh.cpp
  geometry/part_voxelizer.cpp
  geometry/geometry_adaptive_refinement.cpp
  geometry/geometry_engine.cpp
  geometry/voxelization/fast_intersect.cpp
)

target_link_libraries(octlb_mesh PUBLIC octlb_mesh_stl CGAL::CGAL)
