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

function(octlb_add_mpi_mesh_test target_name)
  add_executable(${target_name}
    ${target_name}.cpp
    ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
  )
  target_compile_features(${target_name} PRIVATE cxx_std_17)
  target_link_libraries(${target_name}
    PRIVATE
      octlb_mesh
      GTest::gtest
      MPI::MPI_CXX
  )
  add_test(NAME ${target_name}
    COMMAND
      ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG} 4
      --oversubscribe
      $<TARGET_FILE:${target_name}>
  )
endfunction()

octlb_add_mpi_mesh_test(test_face_pair_list)
octlb_add_mpi_mesh_test(test_load_balancer)
octlb_add_mpi_mesh_test(test_ghost_topology)
