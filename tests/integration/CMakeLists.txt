add_executable(test_cavity3d_serial
  test_cavity3d_serial.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_cavity3d_serial PRIVATE cxx_std_20)

target_include_directories(test_cavity3d_serial
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_compile_definitions(test_cavity3d_serial
  PRIVATE
    OCTLB_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
    OCTLB_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

target_link_libraries(test_cavity3d_serial
  PRIVATE
    octlb_lbm
    octlb_io
    octlb_field_schedule
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_cavity3d_serial
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_cavity3d_serial>
)

add_executable(test_cylinder3d_uniform
  test_cylinder3d_uniform.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_cylinder3d_uniform PRIVATE cxx_std_20)

target_include_directories(test_cylinder3d_uniform
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_compile_definitions(test_cylinder3d_uniform
  PRIVATE
    OCTLB_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)

target_link_libraries(test_cylinder3d_uniform
  PRIVATE
    octlb_lbm
    octlb_field_schedule
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_cylinder3d_uniform
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2
    --oversubscribe
    $<TARGET_FILE:test_cylinder3d_uniform>
)

add_executable(test_cylinder3d_amr
  test_cylinder3d_amr.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_cylinder3d_amr PRIVATE cxx_std_20)

target_include_directories(test_cylinder3d_amr
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_compile_definitions(test_cylinder3d_amr
  PRIVATE
    OCTLB_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)

target_link_libraries(test_cylinder3d_amr
  PRIVATE
    octlb_lbm
    octlb_field_schedule
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_cylinder3d_amr
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 4
    --oversubscribe
    $<TARGET_FILE:test_cylinder3d_amr>
)

# Diagnostic (T11 W3 AMR deadlock hunt): fail fast (600s, not the 40-min job
# timeout) and emit per-exchange debug so the CI log shows exactly which
# LevelCoupler batch / GhostSchedule Waitall hangs. Remove once green.
set_tests_properties(test_cylinder3d_amr PROPERTIES
  TIMEOUT 600
  ENVIRONMENT "OCTLB_COUPLE_DEBUG=1;OCTLB_GHOST_DEBUG=1"
)
