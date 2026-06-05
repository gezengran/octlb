#ifndef OCTLB_SRC_MESH_IO_STL_READER_STL_READER_H_
#define OCTLB_SRC_MESH_IO_STL_READER_STL_READER_H_

#include <string>

#include "src/mesh/geometry/triangle_soup.h"

namespace octlb {

/** Loads ASCII or binary STL into a \c TriangleSoup (OpenLB-style parsing). */
TriangleSoup read_stl_file(const std::string& path);

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_IO_STL_READER_STL_READER_H_
