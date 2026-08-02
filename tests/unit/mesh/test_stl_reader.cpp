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

// W1 (T11) cylinder fixture: a closed cylinder (axis along z, R=0.5, H=1.0,
// 16 segments) used by the cylinder3d case. Verify it parses with the expected
// triangle count and bounding box.
TEST(StlReaderTest, Cylinder_ParsesCountAndBBox) {
  const TriangleSoup soup = read_stl_file(MeshDataPath("cylinder.stl"));
  EXPECT_EQ(soup.triangles().size(), 64u);
  const BoundingBox& bb = soup.bounding_box();
  EXPECT_NEAR(bb.x_min, -0.5, 1e-4);
  EXPECT_NEAR(bb.x_max, 0.5, 1e-4);
  EXPECT_NEAR(bb.y_min, -0.5, 1e-4);
  EXPECT_NEAR(bb.y_max, 0.5, 1e-4);
  EXPECT_NEAR(bb.z_min, 0.0, 1e-4);
  EXPECT_NEAR(bb.z_max, 1.0, 1e-4);
}

// W3 (T11) Schäfer-Turek cylinder fixture: D=0.1 (R=0.05), axis along z,
// centred (0.45, 1.245, 1.25) -- the channel is centred in the cubic domain
// [0,2.5]^3 at y,z in [1.045,1.455], so the cylinder keeps its Schäfer-Turek
// position relative to the channel walls (0.2 off y=wall, 0.205 off z=wall).
// 16 segments -> 64 triangles (same convention as the placeholder cylinder).
TEST(StlReaderTest, CylinderSt_ParsesCountAndBBox) {
  const TriangleSoup soup = read_stl_file(MeshDataPath("cylinder_st.stl"));
  EXPECT_EQ(soup.triangles().size(), 64u);
  const BoundingBox& bb = soup.bounding_box();
  EXPECT_NEAR(bb.x_min, 0.40, 1e-4);
  EXPECT_NEAR(bb.x_max, 0.50, 1e-4);
  EXPECT_NEAR(bb.y_min, 1.195, 1e-4);
  EXPECT_NEAR(bb.y_max, 1.295, 1e-4);
  EXPECT_NEAR(bb.z_min, 1.045, 1e-4);
  EXPECT_NEAR(bb.z_max, 1.455, 1e-4);
}

}  // namespace
}  // namespace octlb
