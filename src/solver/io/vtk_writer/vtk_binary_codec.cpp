#include "src/solver/io/vtk_writer/vtk_binary_codec.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <stdexcept>

namespace octlb {
namespace vtk_binary {
namespace {

constexpr char kEncodeTable[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

std::string EncodeBase64(const unsigned char* data, std::size_t len) {
  std::string out;
  out.reserve(((len + 2) / 3) * 4);
  for (std::size_t i = 0; i < len; i += 3) {
    const unsigned int b0 = data[i];
    const unsigned int b1 = (i + 1 < len) ? data[i + 1] : 0;
    const unsigned int b2 = (i + 2 < len) ? data[i + 2] : 0;
    const unsigned int triple = (b0 << 16) | (b1 << 8) | b2;
    out.push_back(kEncodeTable[(triple >> 18) & 0x3F]);
    out.push_back(kEncodeTable[(triple >> 12) & 0x3F]);
    out.push_back((i + 1 < len) ? kEncodeTable[(triple >> 6) & 0x3F] : '=');
    out.push_back((i + 2 < len) ? kEncodeTable[triple & 0x3F] : '=');
  }
  return out;
}

int DecodeChar(unsigned char c) {
  if (c >= 'A' && c <= 'Z') return c - 'A';
  if (c >= 'a' && c <= 'z') return c - 'a' + 26;
  if (c >= '0' && c <= '9') return c - '0' + 52;
  if (c == '+') return 62;
  if (c == '/') return 63;
  return -1;
}

std::vector<unsigned char> DecodeBase64(const std::string& in) {
  std::vector<unsigned char> out;
  out.reserve(in.size() * 3 / 4);
  unsigned int accum = 0;
  int bits = 0;
  for (unsigned char c : in) {
    if (c == '=' || c == '\n' || c == '\r' || c == ' ' || c == '\t') {
      continue;
    }
    const int v = DecodeChar(c);
    if (v < 0) {
      continue;
    }
    accum = (accum << 6) | static_cast<unsigned int>(v);
    bits += 6;
    if (bits >= 8) {
      bits -= 8;
      out.push_back(static_cast<unsigned char>((accum >> bits) & 0xFF));
    }
  }
  return out;
}

}  // namespace

void WriteFloat64Array(std::ostream& out, const double* data, std::size_t count) {
  const std::size_t payload_bytes = count * sizeof(double);
  std::vector<unsigned char> raw(sizeof(std::uint32_t) + payload_bytes);
  const auto byte_len = static_cast<std::uint32_t>(payload_bytes);
  std::memcpy(raw.data(), &byte_len, sizeof(byte_len));
  std::memcpy(raw.data() + sizeof(std::uint32_t), data, payload_bytes);
  out << EncodeBase64(raw.data(), raw.size());
}

std::vector<double> DecodeFloat64Array(const std::string& base64_payload) {
  const std::vector<unsigned char> raw = DecodeBase64(base64_payload);
  if (raw.size() < sizeof(std::uint32_t)) {
    throw std::runtime_error("VTK binary array too short");
  }
  std::uint32_t byte_len = 0;
  std::memcpy(&byte_len, raw.data(), sizeof(byte_len));
  if (raw.size() < sizeof(std::uint32_t) + byte_len) {
    throw std::runtime_error("VTK binary array truncated");
  }
  if (byte_len % sizeof(double) != 0) {
    throw std::runtime_error("VTK binary array size not multiple of 8");
  }
  const std::size_t count = byte_len / sizeof(double);
  std::vector<double> values(count);
  std::memcpy(values.data(), raw.data() + sizeof(std::uint32_t), byte_len);
  return values;
}

}  // namespace vtk_binary
}  // namespace octlb
