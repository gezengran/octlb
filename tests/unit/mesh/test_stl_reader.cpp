#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <string>

#include "src/mesh/io/stl_reader/stl_reader.h"

namespace octlb {
namespace {

std::string MeshDataPath(const char* name) {
  return (std::filesystem::path(OCTLB_TEST_DATA_DIR) / "mesh" / name)
      .string();
}

void ExpectSingleZUpTriangle(const TriangleSoup& soup) {
  ASSERT_EQ(soup.triangles().size(), 1u);
  const Triangle& tri = soup.triangles().front();
  EXPECT_NEAR(tri.normal[0], 0.0, 1e-6);
  EXPECT_NEAR(tri.normal[1], 0.0, 1e-6);
  EXPECT_NEAR(tri.normal[2], 1.0, 1e-6);

  const BoundingBox& bb = soup.bounding_box();
  EXPECT_NEAR(bb.x_min, 0.0, 1e-6);
  EXPECT_NEAR(bb.y_min, 0.0, 1e-6);
  EXPECT_NEAR(bb.z_min, 0.0, 1e-6);
  EXPECT_NEAR(bb.x_max, 1.0, 1e-6);
  EXPECT_NEAR(bb.y_max, 1.0, 1e-6);
  EXPECT_NEAR(bb.z_max, 0.0, 1e-6);
}

TEST(StlReaderTest, AsciiTriangle_ParsesCountNormalAndBBox) {
  const TriangleSoup soup =
      read_stl_file(MeshDataPath("triangle_ascii.stl"));
  ExpectSingleZUpTriangle(soup);
}

TEST(StlReaderTest, BinaryTriangle_ParsesCountNormalAndBBox) {
  const TriangleSoup soup =
      read_stl_file(MeshDataPath("triangle_binary.stl"));
  ExpectSingleZUpTriangle(soup);
}

}  // namespace
}  // namespace octlb
