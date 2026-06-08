#include <gtest/gtest.h>
#include <mpi.h>

#include <array>

#include "src/common/bounding_box.h"
#include "src/mesh/forest/octree_forest.h"
#include "src/mesh/geometry/geometry_engine.h"
#include "src/mesh/geometry/material_field.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/block_lattice.h"
#include "src/solver/lbm/boundary/bouzidi_pull.h"
#include "src/solver/lbm/bouzidi_link_data.h"
#include "src/solver/lbm/lattice_material_init.h"
#include "tests/unit/mesh/geometry_fixtures.h"

namespace octlb {
namespace {

using Lattice = BlockLattice<double, olb::descriptors::D3Q19<>>;
using Descriptor = olb::descriptors::D3Q19<>;

constexpr int kN = 8;

BoundingBox UnitCubeDomain() {
  return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

GeometryAssembly SolidObstacleAssembly() {
  GeometryAssembly assembly;
  GeometryPart part;
  part.name = "cube";
  part.role = GeometryPartRole::kExternalObstacle;
  part.priority = 0;
  part.soup = test_geom::SolidCube(0.4, 0.4, 0.4, 0.6, 0.6, 0.6);
  assembly.parts.push_back(part);
  return assembly;
}

MaterialField MakeFluidSolidField(int nx, int ny, int nz) {
  MaterialField field(1, nx, ny, nz);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        const MaterialKind kind =
            (i >= nx / 2) ? MaterialKind::kSolid : MaterialKind::kFluid;
        field.set(0, i, j, k, kind);
      }
    }
  }
  return field;
}

GeometryAssembly WallAtMidXAssembly() {
  GeometryAssembly assembly;
  GeometryPart part;
  part.name = "wall";
  part.role = GeometryPartRole::kExternalObstacle;
  part.priority = 0;
  Triangle tri;
  tri.v0 = {0.5, 0.0, 0.0};
  tri.v1 = {0.5, 1.0, 0.0};
  tri.v2 = {0.5, 0.0, 1.0};
  tri.normal = {1.0, 0.0, 0.0};
  part.soup.add_triangle(tri);
  assembly.parts.push_back(part);
  return assembly;
}

CellKind MaterialToCellKind(MaterialKind mk) {
  switch (mk) {
    case MaterialKind::kFluid:
      return CellKind::kFluid;
    case MaterialKind::kSolid:
      return CellKind::kSolid;
    case MaterialKind::kBoundary:
      return CellKind::kBoundary;
  }
  return CellKind::kFluid;
}

}  // namespace

TEST(BouzidiLink, LinkData_CutLinkQFrac) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  forest.partition();
  const MaterialField material = MakeFluidSolidField(kN, kN, kN);
  const BouzidiLinkData links =
      BouzidiLinkData::Build(forest, material, WallAtMidXAssembly());

  bool found = false;
  for (int k = 0; k < kN; ++k) {
    for (int j = 0; j < kN; ++j) {
      const int i = kN / 2 - 1;
      for (int iPop = 1; iPop < Descriptor::q; ++iPop) {
        const double q = links.q_frac(0, i, j, k, iPop);
        if (q > 0.0 && q < 1.0) {
          found = true;
          EXPECT_GT(q, 0.0);
          EXPECT_LT(q, 1.0);
        }
      }
    }
  }
  EXPECT_TRUE(found);
}

TEST(BouzidiLink, Stream_BouzidiReplacesSolidPull) {
  Lattice lat(kN, kN, kN, 1);
  const double u0[3] = {0.0, 0.0, 0.0};
  lat.initialize(1.0, u0);

  for (int k = 0; k < kN; ++k) {
    for (int j = 0; j < kN; ++j) {
      for (int i = kN / 2; i < kN; ++i) {
        lat.set_cell_kind(i, j, k, CellKind::kSolid);
      }
    }
  }

  BouzidiLinkData links(1, kN, kN, kN, Descriptor::q);
  const int i = kN / 2 - 1;
  const int iPop = 1;  // pulls from +x solid neighbor
  constexpr double kQ = 0.35;
  links.set_q_frac(0, i, 0, 0, iPop, kQ);

  lat.set_bouzidi_links(&links);
  lat.collide(1.0);

  auto cell = lat.get(i, 0, 0);
  cell[3] += 0.05;
  const double f_bb = cell[olb::descriptors::opposite<Descriptor>(iPop)];
  const double f_same = cell[iPop];
  const double expected =
      boundary::BouzidiPostCollisionPull(f_bb, f_same, f_same, kQ);

  lat.stream();
  cell = lat.get(i, 0, 0);
  EXPECT_NEAR(cell[iPop], expected, 1e-12);

  lat.set_bouzidi_links(nullptr);
  lat.initialize(1.0, u0);
  for (int k = 0; k < kN; ++k) {
    for (int j = 0; j < kN; ++j) {
      for (int ii = kN / 2; ii < kN; ++ii) {
        lat.set_cell_kind(ii, j, k, CellKind::kSolid);
      }
    }
  }
  lat.collide(1.0);
  cell = lat.get(i, 0, 0);
  cell[3] += 0.05;
  auto solid = lat.get(kN / 2, 0, 0);
  solid[iPop] = 0.25;
  lat.stream();
  cell = lat.get(i, 0, 0);
  EXPECT_NE(cell[iPop], expected);
}

TEST(BouzidiLink, MaterialField_MapsCellKind) {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) {
    GTEST_SKIP() << "single-rank test";
  }

  OctreeForest forest(MPI_COMM_WORLD, UnitCubeDomain());
  GeometryConfig cfg = test_geom::DefaultConfig();
  GeometryEngine engine;
  const MaterialField field = engine.build(forest, SolidObstacleAssembly(), cfg);

  BlockCollection<Lattice> blocks(field.num_octants(), [&](OctantId) {
    return Lattice(field.nx(), field.ny(), field.nz(), 1);
  });
  const double u0[3] = {0.0, 0.0, 0.0};
  initialize_from_material(blocks, field, 1.0, u0, 1.0);

  for (label oid = 0; oid < field.num_octants(); ++oid) {
    const OctantId id = static_cast<OctantId>(oid);
    Lattice& lat = blocks[id];
    for (int k = 0; k < field.nz(); ++k) {
      for (int j = 0; j < field.ny(); ++j) {
        for (int i = 0; i < field.nx(); ++i) {
          EXPECT_EQ(lat.cell_kind(i, j, k),
                    MaterialToCellKind(field.at(id, i, j, k)));
        }
      }
    }
  }
}

}  // namespace octlb
