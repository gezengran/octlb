add_library(octlb_lbm STATIC
  block_lattice.cpp
  level_coupler.cpp
)

target_include_directories(octlb_lbm
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../..
)

target_compile_features(octlb_lbm PUBLIC cxx_std_20)

target_link_libraries(octlb_lbm
  PUBLIC
    octlb_mesh
    octlb_field
    MPI::MPI_CXX
)
