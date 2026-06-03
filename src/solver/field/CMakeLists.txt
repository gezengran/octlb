# octlb_field is defined in parent solver/CMakeLists.txt.

add_library(octlb_field_schedule INTERFACE)
target_link_libraries(octlb_field_schedule INTERFACE
  octlb_field
  octlb_mesh
  MPI::MPI_CXX
)
