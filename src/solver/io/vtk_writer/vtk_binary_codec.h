#ifndef OCTLB_SRC_SOLVER_IO_VTK_WRITER_VTK_BINARY_CODEC_H_
#define OCTLB_SRC_SOLVER_IO_VTK_WRITER_VTK_BINARY_CODEC_H_

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

namespace octlb {
namespace vtk_binary {

/** VTK XML binary array: base64( uint32 byte_length + raw little-endian data ). */
void WriteFloat64Array(std::ostream& out, const double* data, std::size_t count);

std::vector<double> DecodeFloat64Array(const std::string& base64_payload);

}  // namespace vtk_binary
}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_IO_VTK_WRITER_VTK_BINARY_CODEC_H_
