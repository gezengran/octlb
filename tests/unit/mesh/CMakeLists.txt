add_executable(test_octree_forest
  test_octree_forest.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_octree_forest PRIVATE cxx_std_17)

target_link_libraries(test_octree_forest
  PRIVATE
    octlb_mesh
    GTest::gtest
    MPI::MPI_CXX
)

include(GoogleTest)

set(_octlb_mpi_test_cmd
  ${MPIEXEC_EXECUTABLE}
  ${MPIEXEC_NUMPROC_FLAG} 4
  --oversubscribe
  $<TARGET_FILE:test_octree_forest>
)

add_test(NAME test_octree_forest COMMAND ${_octlb_mpi_test_cmd})
