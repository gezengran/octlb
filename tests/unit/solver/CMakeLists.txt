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
