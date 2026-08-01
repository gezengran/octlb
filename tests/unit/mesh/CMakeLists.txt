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

# 2-rank symmetry guard for cross-rank comm_tag agreement.
add_test(NAME test_face_pair_list_symmetry_mpi2
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2
    --oversubscribe
    $<TARGET_FILE:test_face_pair_list>
    --gtest_filter=FacePairList.CrossRankTagsAgreeBetweenRanks
)

# ② W3 Stage B: 2-rank cross-rank edge enumeration + symmetric comm_tag.
add_test(NAME test_face_pair_list_edges_mpi2
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2
    --oversubscribe
    $<TARGET_FILE:test_face_pair_list>
    --gtest_filter=FacePairList.CrossRankEdges_EnumeratedSymmetric
)

add_executable(test_stl_reader test_stl_reader.cpp)
target_compile_features(test_stl_reader PRIVATE cxx_std_17)
target_compile_definitions(test_stl_reader PRIVATE
  OCTLB_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data")
target_link_libraries(test_stl_reader PRIVATE octlb_mesh_stl GTest::gtest_main)
add_test(NAME test_stl_reader COMMAND test_stl_reader)

# CommTagAssigner is a pure header-only component in src/common/ (no MPI/p4est),
# so this test needs only the project-root include path + gtest_main.
add_executable(test_comm_tag_symmetry test_comm_tag_symmetry.cpp)
target_compile_features(test_comm_tag_symmetry PRIVATE cxx_std_17)
target_include_directories(test_comm_tag_symmetry PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(test_comm_tag_symmetry PRIVATE GTest::gtest_main)
add_test(NAME test_comm_tag_symmetry COMMAND test_comm_tag_symmetry)

function(octlb_add_mpi_mesh_test_1rank target_name)
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
      ${MPIEXEC_NUMPROC_FLAG} 1
      --oversubscribe
      $<TARGET_FILE:${target_name}>
  )
endfunction()

octlb_add_mpi_mesh_test_1rank(test_geometry_adaptive_refine)
octlb_add_mpi_mesh_test_1rank(test_geometry_engine)

add_test(NAME test_geometry_adaptive_refine_two_rank
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2
    --oversubscribe
    $<TARGET_FILE:test_geometry_adaptive_refine>
    --gtest_filter=GeometryAdaptiveRefineTest.TwoRank_RefineCollective)
