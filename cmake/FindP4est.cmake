#[=======================================================================[.rst:
FindP4est
---------

Locates p4est (3D API: p8est) and exposes imported target ``P4est::p4est``.

Cache / environment
^^^^^^^^^^^^^^^^^^^

``P4EST_ROOT``
  Installation prefix (``include/``, ``lib/``).

Result variables
^^^^^^^^^^^^^^^^

``P4est_FOUND``
``P4EST_INCLUDE_DIRS``
``P4EST_LIBRARIES``  (p4est + sc)
#]=======================================================================]

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

set(_p4est_hints)
if(P4EST_ROOT)
  list(APPEND _p4est_hints
    "${P4EST_ROOT}/include"
    "${P4EST_ROOT}/lib"
    "${P4EST_ROOT}")
endif()

find_path(P4EST_INCLUDE_DIR
  NAMES p4est.h
  HINTS ${_p4est_hints}
  PATH_SUFFIXES include)

find_library(P4EST_LIBRARY
  NAMES p4est
  HINTS ${_p4est_hints}
  PATH_SUFFIXES lib lib64)

find_library(SC_LIBRARY
  NAMES sc
  HINTS ${_p4est_hints}
  PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(P4est
  REQUIRED_VARS P4EST_INCLUDE_DIR P4EST_LIBRARY SC_LIBRARY)

if(P4est_FOUND)
  set(P4EST_INCLUDE_DIRS "${P4EST_INCLUDE_DIR}")
  set(P4EST_LIBRARIES "${P4EST_LIBRARY}" "${SC_LIBRARY}")

  if(NOT TARGET P4est::p4est)
    add_library(P4est::p4est INTERFACE IMPORTED)
    target_include_directories(P4est::p4est SYSTEM INTERFACE
      "${P4EST_INCLUDE_DIR}")
    target_link_libraries(P4est::p4est INTERFACE
      "${P4EST_LIBRARY}"
      "${SC_LIBRARY}"
      MPI::MPI_CXX)
  endif()
endif()

mark_as_advanced(P4EST_INCLUDE_DIR P4EST_LIBRARY SC_LIBRARY)
