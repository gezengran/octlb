add_library(octlb_mesh STATIC)

target_include_directories(octlb_mesh
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_compile_features(octlb_mesh PUBLIC cxx_std_17)

target_link_libraries(octlb_mesh
  PUBLIC
    P4est::p4est
    MPI::MPI_CXX
)

add_subdirectory(mesh)
