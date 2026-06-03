add_executable(test_collision_bgk
  test_collision_bgk.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_collision_bgk PRIVATE cxx_std_20)

target_link_libraries(test_collision_bgk
  PRIVATE
    octlb_lbm
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_collision_bgk
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_collision_bgk>
)

add_executable(test_block_collection
  test_block_collection.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_block_collection PRIVATE cxx_std_17)

target_link_libraries(test_block_collection
  PRIVATE
    octlb_field
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_block_collection
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_block_collection>
)

add_executable(test_ghost_schedule
  test_ghost_schedule.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_ghost_schedule PRIVATE cxx_std_20)

target_link_libraries(test_ghost_schedule
  PRIVATE
    octlb_field_schedule
    octlb_lbm
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_ghost_schedule
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2
    --oversubscribe
    $<TARGET_FILE:test_ghost_schedule>
)
