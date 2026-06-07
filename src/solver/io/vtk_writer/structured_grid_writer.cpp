#include "src/solver/io/vtk_writer/structured_grid_writer.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "src/solver/io/vtk_writer/vtk_binary_codec.h"

namespace octlb {
namespace {

double Lerp(scalar a, scalar b, int i, int n) {
  if (n <= 0) {
    return a;
  }
  return a + (b - a) * static_cast<double>(i) / static_cast<double>(n);
}

int CellCount(int nx, int ny, int nz) { return nx * ny * nz; }

int PointCount(int nx, int ny, int nz) { return (nx + 1) * (ny + 1) * (nz + 1); }

void BuildPoints(const VtkBlockMeta& meta, std::vector<double>* xyz) {
  const int nx = meta.nx;
  const int ny = meta.ny;
  const int nz = meta.nz;
  const int npts = PointCount(nx, ny, nz);
  xyz->resize(static_cast<std::size_t>(npts) * 3);
  std::size_t idx = 0;
  for (int k = 0; k <= nz; ++k) {
    for (int j = 0; j <= ny; ++j) {
      for (int i = 0; i <= nx; ++i) {
        (*xyz)[idx++] = Lerp(meta.bounds.x_min, meta.bounds.x_max, i, nx);
        (*xyz)[idx++] = Lerp(meta.bounds.y_min, meta.bounds.y_max, j, ny);
        (*xyz)[idx++] = Lerp(meta.bounds.z_min, meta.bounds.z_max, k, nz);
      }
    }
  }
}

void BuildCellField(const VtkBlockMeta& meta, const VtkCellFieldView& field,
                    std::vector<double>* values) {
  const int nx = meta.nx;
  const int ny = meta.ny;
  const int nz = meta.nz;
  const int ncomp = field.components;
  const int ncells = CellCount(nx, ny, nz);
  values->resize(static_cast<std::size_t>(ncells * ncomp));
  std::array<double, 3> sample{};
  std::size_t out = 0;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        field.sample(i, j, k, sample.data());
        for (int c = 0; c < ncomp; ++c) {
          (*values)[out++] = sample[static_cast<std::size_t>(c)];
        }
      }
    }
  }
}

}  // namespace

void WriteStructuredGridVts(const std::string& filename, const VtkBlockMeta& meta,
                            std::span<const VtkCellFieldView> fields) {
  if (meta.nx <= 0 || meta.ny <= 0 || meta.nz <= 0) {
    throw std::invalid_argument("WriteStructuredGridVts: nx, ny, nz must be positive");
  }

  std::vector<double> points;
  BuildPoints(meta, &points);

  std::ofstream out(filename, std::ios::trunc);
  if (!out) {
    throw std::runtime_error("WriteStructuredGridVts: failed to open " + filename);
  }

  const int nx = meta.nx;
  const int ny = meta.ny;
  const int nz = meta.nz;
  const int npts = PointCount(nx, ny, nz);
  const int ncells = CellCount(nx, ny, nz);

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" "
         "byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
  out << "  <StructuredGrid WholeExtent=\"0 " << nx << " 0 " << ny << " 0 " << nz
      << "\">\n";
  out << "    <Piece Extent=\"0 " << nx << " 0 " << ny << " 0 " << nz
      << "\" NumberOfPoints=\"" << npts << "\" NumberOfCells=\"" << ncells
      << "\">\n";
  out << "      <PointData>\n";
  out << "      </PointData>\n";
  out << "      <CellData>\n";
  for (const VtkCellFieldView& field : fields) {
    std::vector<double> cell_values;
    BuildCellField(meta, field, &cell_values);
    out << "        <DataArray type=\"Float64\" Name=\"" << field.name
        << "\" NumberOfComponents=\"" << field.components
        << "\" format=\"binary\" encoding=\"base64\">\n          ";
    vtk_binary::WriteFloat64Array(out, cell_values.data(), cell_values.size());
    out << "\n        </DataArray>\n";
  }
  out << "      </CellData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         "format=\"binary\" encoding=\"base64\">\n          ";
  vtk_binary::WriteFloat64Array(out, points.data(), points.size());
  out << "\n        </DataArray>\n";
  out << "      </Points>\n";
  out << "    </Piece>\n";
  out << "  </StructuredGrid>\n";
  out << "</VTKFile>\n";
}

}  // namespace octlb
