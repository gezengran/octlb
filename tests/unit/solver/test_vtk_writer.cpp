#include <gtest/gtest.h>
#include <mpi.h>
#include <unistd.h>

#include <filesystem>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/io/vtk_writer/amr_vtk_writer.h"
#include "src/solver/io/vtk_writer/structured_grid_writer.h"
#include "src/solver/io/vtk_writer/vtk_cell_field.h"
#include "tests/unit/solver/vtk_writer_fixtures.h"

namespace octlb {
namespace {

namespace fs = std::filesystem;
using vtk_test::ApproxEq;
using vtk_test::DecodeCellArrayNamed;
using vtk_test::DecodePointsArray;
using vtk_test::DummyScalarField;
using vtk_test::DummyVectorField;
using vtk_test::InteriorMarkedField;
using vtk_test::MpiSharedTempDir;
using vtk_test::ReadFile;

std::string UniqueTempDir(const std::string& tag) {
  const fs::path dir =
      fs::temp_directory_path() /
      ("octlb_vtk_" + tag + "_" + std::to_string(::getpid()));
  fs::create_directories(dir);
  return dir.string();
}

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

VtkBlockMeta MakeMeta(OctantId id, int n, const BoundingBox& bounds) {
  VtkBlockMeta meta;
  meta.id = id;
  meta.nx = n;
  meta.ny = n;
  meta.nz = n;
  meta.bounds = bounds;
  return meta;
}

bool AllRanksTrue(bool local_ok, MPI_Comm comm) {
  int flag = local_ok ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&flag, &global, 1, MPI_INT, MPI_MIN, comm);
  return global == 1;
}

/** 2×1×1 brick: exactly two root octants split along x on one rank. */
OctreeForest MakeTwoOctantBrickForest(MPI_Comm comm) {
  return OctreeForest(comm, UnitCubeDomain(), 2, 1, 1);
}

void ExpectCornerBoundsMatchFile(const std::string& path,
                                 const BoundingBox& expected) {
  const std::vector<double> pts = DecodePointsArray(ReadFile(path));
  ASSERT_FALSE(pts.empty());
  EXPECT_TRUE(ApproxEq(pts[0], expected.x_min));
  EXPECT_TRUE(ApproxEq(pts[1], expected.y_min));
  EXPECT_TRUE(ApproxEq(pts[2], expected.z_min));
  const std::size_t last = pts.size() - 3;
  EXPECT_TRUE(ApproxEq(pts[last + 0], expected.x_max));
  EXPECT_TRUE(ApproxEq(pts[last + 1], expected.y_max));
  EXPECT_TRUE(ApproxEq(pts[last + 2], expected.z_max));
}

}  // namespace

// ── W1: single-block .vts kernel (WriteStructuredGridVts) ─────────────────

TEST(VtkWriter, SingleBlock_ConstantScalar) {
  constexpr int kN = 2;
  const std::string dir = UniqueTempDir("w1_scalar");
  const std::string path = dir + "/block.vts";
  const VtkBlockMeta meta = MakeMeta(0, kN, UnitCubeDomain());
  DummyScalarField field(3.5);
  const VtkCellFieldView view = VtkCellFieldView::From(field);
  WriteStructuredGridVts(path, meta, std::span<const VtkCellFieldView>(&view, 1));

  const std::string xml = ReadFile(path);
  EXPECT_NE(xml.find("StructuredGrid"), std::string::npos);
  EXPECT_NE(xml.find("NumberOfPoints=\"27\""), std::string::npos);
  EXPECT_NE(xml.find("NumberOfCells=\"8\""), std::string::npos);

  const std::vector<double> cells = DecodeCellArrayNamed(xml, "dummy_scalar");
  ASSERT_EQ(cells.size(), static_cast<std::size_t>(kN * kN * kN));
  for (double v : cells) {
    EXPECT_DOUBLE_EQ(v, 3.5);
  }
}

TEST(VtkWriter, SingleBlock_CornerBoundsMatch) {
  constexpr int kN = 4;
  const BoundingBox bounds{0.2, 0.3, 0.4, 0.8, 0.9, 1.0};
  const std::string dir = UniqueTempDir("w1_bounds");
  const std::string path = dir + "/block.vts";
  DummyScalarField field;
  const VtkCellFieldView view = VtkCellFieldView::From(field);
  WriteStructuredGridVts(path, MakeMeta(0, kN, bounds),
                         std::span<const VtkCellFieldView>(&view, 1));

  ExpectCornerBoundsMatchFile(path, bounds);
}

TEST(VtkWriter, SingleBlock_MultiField) {
  constexpr int kN = 2;
  const std::string path = UniqueTempDir("w1_multi") + "/block.vts";
  DummyScalarField scalar(1.0);
  DummyVectorField vector;
  const VtkCellFieldView fields[2] = {VtkCellFieldView::From(scalar),
                                      VtkCellFieldView::From(vector)};
  WriteStructuredGridVts(path, MakeMeta(0, kN, UnitCubeDomain()),
                         std::span<const VtkCellFieldView>(fields, 2));

  const std::string xml = ReadFile(path);
  EXPECT_NE(xml.find("Name=\"dummy_scalar\""), std::string::npos);
  EXPECT_NE(xml.find("Name=\"dummy_vector\""), std::string::npos);
  EXPECT_NE(xml.find("NumberOfComponents=\"3\""), std::string::npos);

  const std::vector<double> vec = DecodeCellArrayNamed(xml, "dummy_vector");
  ASSERT_EQ(vec.size(), static_cast<std::size_t>(kN * kN * kN * 3));
  EXPECT_DOUBLE_EQ(vec[0], 0.0);
  EXPECT_DOUBLE_EQ(vec[1], 0.0);
  EXPECT_DOUBLE_EQ(vec[2], 0.0);
}

// ── W2: AmrVtkWriter + BlockCollection (PRD #11) ──────────────────────────

TEST(VtkWriter, TwoOctants_DistinctBounds) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank != 0) {
    return;
  }

  constexpr int kN = 4;
  OctreeForest forest = MakeTwoOctantBrickForest(comm);
  ASSERT_EQ(forest.local_num_octants(), 2);

  const BoundingBox b0 = forest.quadrant_bounds(0);
  const BoundingBox b1 = forest.quadrant_bounds(1);
  EXPECT_TRUE(ApproxEq(b0.x_max, 0.5));
  EXPECT_TRUE(ApproxEq(b1.x_min, 0.5));
  EXPECT_NE(b0.x_min, b1.x_min);

  const std::string dir = UniqueTempDir("w2_two");
  AmrVtkWriter writer(comm, forest, kN, kN, kN, dir, "amr");
  DummyScalarField field;
  const VtkCellFieldView view = VtkCellFieldView::From(field);
  BlockCollection<int> blocks(2, [](OctantId id) { return static_cast<int>(id); });
  const auto paths =
      writer.WriteTimestep(0, blocks, std::span<const VtkCellFieldView>(&view, 1));
  ASSERT_EQ(paths.size(), 2u);

  ExpectCornerBoundsMatchFile(paths[0], b0);
  ExpectCornerBoundsMatchFile(paths[1], b1);
}

TEST(VtkWriter, CellData_CountMatchesN) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank != 0) {
    return;
  }

  constexpr int kN = 8;
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  const std::string dir = UniqueTempDir("w2_count");
  AmrVtkWriter writer(comm, forest, kN, kN, kN, dir, "count");
  DummyScalarField field;
  const VtkCellFieldView view = VtkCellFieldView::From(field);
  BlockCollection<int> blocks(forest.local_num_octants(),
                              [](OctantId) { return 0; });
  const auto paths =
      writer.WriteTimestep(7, blocks, std::span<const VtkCellFieldView>(&view, 1));
  ASSERT_FALSE(paths.empty());

  const std::vector<double> cells =
      DecodeCellArrayNamed(ReadFile(paths[0]), "dummy_scalar");
  EXPECT_EQ(cells.size(), static_cast<std::size_t>(kN * kN * kN));
}

TEST(VtkWriter, InteriorOnly_NoGhost) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank != 0) {
    return;
  }

  constexpr int kN = 4;
  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  const std::string dir = UniqueTempDir("w2_interior");
  AmrVtkWriter writer(comm, forest, kN, kN, kN, dir, "interior");
  InteriorMarkedField field(kN, kN, kN);
  const VtkCellFieldView view = VtkCellFieldView::From(field);
  BlockCollection<int> blocks(forest.local_num_octants(),
                              [](OctantId) { return 0; });
  const auto paths =
      writer.WriteTimestep(0, blocks, std::span<const VtkCellFieldView>(&view, 1));
  ASSERT_EQ(paths.size(), 1u);

  const std::vector<double> cells =
      DecodeCellArrayNamed(ReadFile(paths[0]), "interior_marker");
  ASSERT_EQ(cells.size(), static_cast<std::size_t>(kN * kN * kN));
  for (double v : cells) {
    EXPECT_NE(v, field.HaloMarker());
  }
  EXPECT_DOUBLE_EQ(cells[0], 0.0);
}

// ── W3: parallel .vtm / .pvd (2 rank) ─────────────────────────────────────

TEST(VtkWriter, TwoRank_VtmListsAllBlocks) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  if (size < 2) {
    GTEST_SKIP() << "Requires at least 2 MPI ranks";
  }

  OctreeForest forest = MakeTwoOctantBrickForest(comm);
  forest.partition();
  ASSERT_EQ(forest.local_num_octants(), 1);

  constexpr int kN = 2;
  const std::string dir = MpiSharedTempDir("w3_mpi");
  AmrVtkWriter writer(comm, forest, kN, kN, kN, dir, "parallel");
  DummyScalarField field;
  const VtkCellFieldView view = VtkCellFieldView::From(field);
  BlockCollection<int> blocks(1, [](OctantId id) { return static_cast<int>(id); });
  const auto paths =
      writer.WriteTimestep(0, blocks, std::span<const VtkCellFieldView>(&view, 1));
  ASSERT_EQ(paths.size(), 1u);
  writer.WriteVtmAndPvd(0, paths);

  bool ok = true;
  if (rank == 0) {
    const std::string vtm_xml = ReadFile(dir + "/parallel_T00000.vtm");
    const std::string pvd_xml = ReadFile(dir + "/parallel.pvd");

    const std::regex block_re("<Block index=\"\\d+\"");
    const auto begin =
        std::sregex_iterator(vtm_xml.begin(), vtm_xml.end(), block_re);
    const auto end = std::sregex_iterator();
    ok = std::distance(begin, end) == 2;

    const std::regex file_re("file=\"([^\"]+)\"");
    std::set<std::string> vts_refs;
    for (auto it = std::sregex_iterator(vtm_xml.begin(), vtm_xml.end(), file_re);
         it != std::sregex_iterator(); ++it) {
      vts_refs.insert((*it)[1].str());
    }
    ok = ok && vts_refs.size() == 2u;
    for (const std::string& rel : vts_refs) {
      ok = ok && fs::exists(fs::path(dir) / rel);
    }

    ok = ok && pvd_xml.find("<DataSet timestep=\"0\"") != std::string::npos;
    ok = ok && pvd_xml.find("parallel_T00000.vtm") != std::string::npos;
  }

  EXPECT_TRUE(AllRanksTrue(ok, comm));
}

}  // namespace octlb
