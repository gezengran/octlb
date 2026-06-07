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

add_executable(test_lagrava_coupler
  test_lagrava_coupler.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_lagrava_coupler PRIVATE cxx_std_20)

target_link_libraries(test_lagrava_coupler
  PRIVATE
    octlb_lbm
    octlb_field_schedule
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_lagrava_coupler
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_lagrava_coupler>
)

add_test(NAME test_lagrava_coupler_mpi2
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2
    --oversubscribe
    $<TARGET_FILE:test_lagrava_coupler>
    --gtest_filter=LagravaCoupler.TwoRank*:LagravaCoupler.CoarseFine_CommTagsUnique*
)

add_executable(test_time_loop_levels
  test_time_loop_levels.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_time_loop_levels PRIVATE cxx_std_20)

target_link_libraries(test_time_loop_levels
  PRIVATE
    octlb_lbm
    octlb_field_schedule
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_time_loop_levels
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_time_loop_levels>
)

add_executable(test_vtk_writer
  test_vtk_writer.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_vtk_writer PRIVATE cxx_std_20)

target_include_directories(test_vtk_writer
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(test_vtk_writer
  PRIVATE
    octlb_io
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_vtk_writer
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_vtk_writer>
)

add_test(NAME test_vtk_writer_two_rank
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2
    --oversubscribe
    $<TARGET_FILE:test_vtk_writer>
    --gtest_filter=VtkWriter.TwoRank_VtmListsAllBlocks
)
