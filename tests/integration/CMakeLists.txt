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
