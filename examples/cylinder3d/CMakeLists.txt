add_executable(cylinder3d cylinder3d.cpp)

target_compile_features(cylinder3d PRIVATE cxx_std_20)

target_include_directories(cylinder3d
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

# Default STL fixture for --stl (the example accepts --stl PATH to override).
target_compile_definitions(cylinder3d
  PRIVATE
    OCTLB_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)

target_link_libraries(cylinder3d
  PRIVATE
    octlb_lbm
    octlb_io
    octlb_field_schedule
    MPI::MPI_CXX
)