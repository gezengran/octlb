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

add_executable(test_lattice_cell_kind
  test_lattice_cell_kind.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_lattice_cell_kind PRIVATE cxx_std_20)

target_link_libraries(test_lattice_cell_kind
  PRIVATE
    octlb_lbm
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_lattice_cell_kind
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_lattice_cell_kind>
)

add_executable(test_domain_boundary
  test_domain_boundary.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_domain_boundary PRIVATE cxx_std_20)

target_link_libraries(test_domain_boundary
  PRIVATE
    octlb_lbm
    octlb_field
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_domain_boundary
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_domain_boundary>
)

add_executable(test_interpolated_velocity_compute_rho
  test_interpolated_velocity_compute_rho.cpp
)

target_compile_features(test_interpolated_velocity_compute_rho PRIVATE cxx_std_20)

target_link_libraries(test_interpolated_velocity_compute_rho
  PRIVATE
    octlb_lbm
    GTest::gtest_main
)

add_test(NAME test_interpolated_velocity_compute_rho
  COMMAND test_interpolated_velocity_compute_rho)

add_executable(test_drag_postprocessor
  test_drag_postprocessor.cpp
)

target_compile_features(test_drag_postprocessor PRIVATE cxx_std_20)

target_link_libraries(test_drag_postprocessor
  PRIVATE
    octlb_lbm
    GTest::gtest_main
)

add_test(NAME test_drag_postprocessor COMMAND test_drag_postprocessor)

add_executable(test_outlet_bc
  test_outlet_bc.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_outlet_bc PRIVATE cxx_std_20)

target_link_libraries(test_outlet_bc
  PRIVATE
    octlb_lbm
    octlb_field
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_outlet_bc
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_outlet_bc>
)

add_executable(test_bouzidi_link
  test_bouzidi_link.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_bouzidi_link PRIVATE cxx_std_20)

target_include_directories(test_bouzidi_link
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(test_bouzidi_link
  PRIVATE
    octlb_lbm
    octlb_mesh
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_bouzidi_link
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_bouzidi_link>
)

add_executable(test_inlet_velocity_field test_inlet_velocity_field.cpp)

target_compile_features(test_inlet_velocity_field PRIVATE cxx_std_20)

target_link_libraries(test_inlet_velocity_field
  PRIVATE
    octlb_lbm
    GTest::gtest_main
)

add_test(NAME test_inlet_velocity_field COMMAND test_inlet_velocity_field)

add_executable(test_inlet_bc_integration
  test_inlet_bc_integration.cpp
  ${CMAKE_SOURCE_DIR}/tests/mpi_main.cpp
)

target_compile_features(test_inlet_bc_integration PRIVATE cxx_std_20)

target_link_libraries(test_inlet_bc_integration
  PRIVATE
    octlb_lbm
    octlb_field_schedule
    GTest::gtest
    MPI::MPI_CXX
)

add_test(NAME test_inlet_bc_integration
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 1
    --oversubscribe
    $<TARGET_FILE:test_inlet_bc_integration>
)

add_executable(test_unit_converter test_unit_converter.cpp)

target_compile_features(test_unit_converter PRIVATE cxx_std_20)

target_link_libraries(test_unit_converter
  PRIVATE
    octlb_lbm
    GTest::gtest_main
)

add_test(NAME test_unit_converter COMMAND test_unit_converter)
