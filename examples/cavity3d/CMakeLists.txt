add_executable(cavity3d cavity3d.cpp)

target_compile_features(cavity3d PRIVATE cxx_std_20)

target_include_directories(cavity3d
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(cavity3d
  PRIVATE
    octlb_lbm
    octlb_io
    octlb_field_schedule
    MPI::MPI_CXX
)
