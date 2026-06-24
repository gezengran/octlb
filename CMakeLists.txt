cmake_minimum_required(VERSION 3.15)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project(OctLB
  VERSION 0.1.0
  LANGUAGES CXX
  DESCRIPTION "Octree-based parallel LBM-AMR framework"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(MPI REQUIRED)
find_package(GTest REQUIRED)
find_package(P4est REQUIRED)
find_package(CGAL REQUIRED)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(examples)
