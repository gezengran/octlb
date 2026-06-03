add_library(octlb_lbm STATIC
  block_lattice.cpp
)

target_include_directories(octlb_lbm
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../..
)

target_compile_features(octlb_lbm PUBLIC cxx_std_20)

# No PLATFORM_CPU_SISD needed: OctLB calls equilibrium formulas directly,
# bypassing OpenLB's platform dispatch framework.
# No OLB_VERSION, no tinyxml2, no MPI singleton from OpenLB.

target_link_libraries(octlb_lbm
  PUBLIC MPI::MPI_CXX
)
