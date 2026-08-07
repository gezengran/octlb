#include "src/mesh/geometry/geometry_engine.h"

#include <algorithm>

#include "src/mesh/geometry/cgal_surface_mesh.h"
#include "src/mesh/geometry/geometry_adaptive_refinement.h"
#include "src/mesh/geometry/geometry_build_error.h"
#include "src/mesh/geometry/part_voxelizer.h"

namespace octlb {

MaterialField GeometryEngine::build(OctreeForest& forest,
                                    const GeometryAssembly& assembly,
                                    const GeometryConfig& config) const {
  if (assembly.parts.empty()) {
    throw GeometryBuildError("GeometryAssembly has no parts");
  }

  for (const GeometryPart& part : assembly.parts) {
    validate_geometry_part(part);
  }

  const TriangleSoup surface_soup = merge_assembly_soup(assembly);
  const BoundingBox geom_bbox = extended_geometry_bbox(surface_soup, config);
  resolve_bounding(forest, geom_bbox, config);
  resolve_surface(forest, surface_soup, config);

  std::vector<GeometryPart> ordered = assembly.parts;
  std::sort(ordered.begin(), ordered.end(),
            [](const GeometryPart& a, const GeometryPart& b) {
              return a.priority < b.priority;
            });

  MaterialField merged(forest.local_num_octants(), config.cell_width,
                       config.cell_width, config.cell_width);
  for (const GeometryPart& part : ordered) {
    MaterialField part_field = voxelize_part(forest, part, config);
    merge_material_field(&merged, part_field, part.role);
  }
  return merged;
}

}  // namespace octlb
