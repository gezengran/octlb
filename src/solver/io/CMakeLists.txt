add_library(octlb_io STATIC
  vtk_writer/vtk_binary_codec.cpp
  vtk_writer/structured_grid_writer.cpp
  vtk_writer/amr_vtk_writer.cpp
  vtk_writer/parallel_vtk_index.cpp
)

target_include_directories(octlb_io
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../..
)

target_compile_features(octlb_io PUBLIC cxx_std_20)

target_link_libraries(octlb_io
  PUBLIC
    octlb_field
    octlb_mesh
    MPI::MPI_CXX
)
