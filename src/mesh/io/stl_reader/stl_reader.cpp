#include "src/mesh/io/stl_reader/stl_reader.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace octlb {
namespace {

void ComputeNormal(Triangle* tri) {
  const scalar ax = tri->v1[0] - tri->v0[0];
  const scalar ay = tri->v1[1] - tri->v0[1];
  const scalar az = tri->v1[2] - tri->v0[2];
  const scalar bx = tri->v2[0] - tri->v0[0];
  const scalar by = tri->v2[1] - tri->v0[1];
  const scalar bz = tri->v2[2] - tri->v0[2];
  tri->normal[0] = ay * bz - az * by;
  tri->normal[1] = az * bx - ax * bz;
  tri->normal[2] = ax * by - ay * bx;
  const scalar len = std::sqrt(tri->normal[0] * tri->normal[0] +
                               tri->normal[1] * tri->normal[1] +
                               tri->normal[2] * tri->normal[2]);
  if (len > 0.0) {
    tri->normal[0] /= len;
    tri->normal[1] /= len;
    tri->normal[2] /= len;
  }
}

bool IsBinaryStl(const std::string& path) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    return false;
  }
  char header[80] = {};
  file.read(header, 80);
  if (file.gcount() != 80) {
    return false;
  }
  std::string header_str(header, 80);
  std::transform(header_str.begin(), header_str.end(), header_str.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return header_str.find("solid") == std::string::npos;
}

TriangleSoup ReadAsciiStl(std::istream& in) {
  TriangleSoup soup;
  Triangle tri;
  bool in_facet = false;
  bool have_normal = false;
  int vertex_count = 0;

  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    std::string token;
    if (!(iss >> token)) {
      continue;
    }
    if (token == "facet") {
      in_facet = true;
      have_normal = false;
      vertex_count = 0;
      std::string normal_token;
      if (iss >> normal_token && normal_token == "normal") {
        iss >> tri.normal[0] >> tri.normal[1] >> tri.normal[2];
        have_normal = true;
      }
    } else if (token == "vertex" && in_facet) {
      std::array<scalar, 3>* v = nullptr;
      if (vertex_count == 0) {
        v = &tri.v0;
      } else if (vertex_count == 1) {
        v = &tri.v1;
      } else if (vertex_count == 2) {
        v = &tri.v2;
      }
      if (v != nullptr) {
        iss >> (*v)[0] >> (*v)[1] >> (*v)[2];
        ++vertex_count;
      }
    } else if (token == "endfacet" && in_facet) {
      if (vertex_count != 3) {
        throw std::runtime_error("STL ASCII facet does not have 3 vertices");
      }
      if (!have_normal) {
        ComputeNormal(&tri);
      }
      soup.add_triangle(tri);
      in_facet = false;
    } else if (token == "endsolid") {
      break;
    }
  }
  return soup;
}

TriangleSoup ReadBinaryStl(std::istream& in) {
  TriangleSoup soup;
  char header[80];
  in.read(header, 80);
  if (!in) {
    throw std::runtime_error("STL binary header read failed");
  }
  std::int32_t n_facets = 0;
  in.read(reinterpret_cast<char*>(&n_facets), sizeof(n_facets));
  if (!in) {
    throw std::runtime_error("STL binary facet count read failed");
  }
  for (std::int32_t f = 0; f < n_facets; ++f) {
    float data[12];
    in.read(reinterpret_cast<char*>(data), sizeof(data));
    if (!in) {
      throw std::runtime_error("STL binary facet data read failed");
    }
    std::uint16_t attr = 0;
    in.read(reinterpret_cast<char*>(&attr), sizeof(attr));
    Triangle tri;
    tri.normal[0] = data[0];
    tri.normal[1] = data[1];
    tri.normal[2] = data[2];
    tri.v0[0] = data[3];
    tri.v0[1] = data[4];
    tri.v0[2] = data[5];
    tri.v1[0] = data[6];
    tri.v1[1] = data[7];
    tri.v1[2] = data[8];
    tri.v2[0] = data[9];
    tri.v2[1] = data[10];
    tri.v2[2] = data[11];
    ComputeNormal(&tri);
    soup.add_triangle(tri);
  }
  return soup;
}

}  // namespace

TriangleSoup read_stl_file(const std::string& path) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    throw std::runtime_error("Cannot open STL file: " + path);
  }
  if (IsBinaryStl(path)) {
    return ReadBinaryStl(file);
  }
  file.clear();
  file.seekg(0);
  return ReadAsciiStl(file);
}

}  // namespace octlb
