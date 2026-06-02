add_library(octlb_field INTERFACE)

# Project root (same as octlb_mesh): includes like "src/common/types.h".
target_include_directories(octlb_field
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../..
)

target_compile_features(octlb_field INTERFACE cxx_std_17)

add_subdirectory(field)
add_subdirectory(lbm)
