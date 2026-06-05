#include "src/mesh/geometry/geometry_adaptive_refinement.h"

#include <algorithm>
#include <functional>
#include <unordered_set>
#include <vector>

#include <mpi.h>

#include "src/mesh/geometry/voxelization/fast_intersect.h"
#include "src/mesh/geometry/cgal_surface_mesh.h"
#include "src/mesh/load_balance/weighted_load_balancer.h"

namespace octlb {
namespace {

bool BoxesOverlap(const BoundingBox& a, const BoundingBox& b) {
  return !(a.x_max < b.x_min || a.x_min > b.x_max || a.y_max < b.y_min ||
           a.y_min > b.y_max || a.z_max < b.z_min || a.z_min > b.z_max);
}

void RefineCollective(OctreeForest& forest,
                      const std::function<bool(OctantId)>& criterion,
                      int max_level, MPI_Comm comm) {
  const label n_before = forest.local_num_octants();
  forest.refine(criterion, max_level);
  const label n_after = forest.local_num_octants();
  int local_changed = (n_after != n_before) ? 1 : 0;
  int global_changed = 0;
  MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_MAX, comm);
  if (global_changed > 0) {
    forest.balance();
  }
}

int GlobalFinestLevel(const OctreeForest& forest, MPI_Comm comm) {
  int local_max = 0;
  for (label id = 0; id < forest.local_num_octants(); ++id) {
    local_max = std::max(local_max, forest.quadrant_level(id));
  }
  int global_max = 0;
  MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, comm);
  return global_max;
}

}  // namespace

BoundingBox extended_geometry_bbox(const TriangleSoup& soup,
                                   const GeometryConfig& config) {
  BoundingBox bb = soup.bounding_box();
  bb.x_min -= config.bound_width;
  bb.y_min -= config.bound_width;
  bb.z_min -= config.bound_width;
  bb.x_max += config.bound_width;
  bb.y_max += config.bound_width;
  bb.z_max += config.bound_width;
  if (config.wake_length > 0.0) {
    switch (config.wake_direction) {
      case 0:
        bb.x_max += config.wake_length;
        break;
      case 1:
        bb.y_max += config.wake_length;
        break;
      case 2:
        bb.z_max += config.wake_length;
        break;
      default:
        break;
    }
  }
  return bb;
}

void resolve_bounding(OctreeForest& forest, const BoundingBox& geom_bbox,
                      const GeometryConfig& config) {
  const int target_level = std::max(0, config.max_level - 2);
  MPI_Comm comm = MPI_COMM_WORLD;
  for (int level = 0; level <= target_level; ++level) {
    std::unordered_set<OctantId> mark;
    for (label id = 0; id < forest.local_num_octants(); ++id) {
      if (forest.quadrant_level(id) == level &&
          BoxesOverlap(forest.quadrant_bounds(id), geom_bbox)) {
        mark.insert(id);
      }
    }
    RefineCollective(
        forest, [&](OctantId id) { return mark.count(id) > 0; },
        config.max_level, comm);
  }
}

void resolve_surface(OctreeForest& forest, const TriangleSoup& surface_soup,
                     const GeometryConfig& config) {
  const CgalTriangleSurface surface = CgalTriangleSurface::from_soup(surface_soup);
  MPI_Comm comm = MPI_COMM_WORLD;

  for (int iteration = 0; iteration < config.resolve_surface_times; ++iteration) {
    const int finest = GlobalFinestLevel(forest, comm);
    std::vector<BoundingBox> leaf_boxes;
    std::vector<OctantId> leaf_ids;
    for (label id = 0; id < forest.local_num_octants(); ++id) {
      const int lvl = forest.quadrant_level(id);
      if (lvl == finest && lvl < config.max_level) {
        leaf_boxes.push_back(forest.quadrant_bounds(id));
        leaf_ids.push_back(id);
      }
    }

    int local_has = leaf_boxes.empty() ? 0 : 1;
    int global_has = 0;
    MPI_Allreduce(&local_has, &global_has, 1, MPI_INT, MPI_MAX, comm);
    if (global_has == 0) {
      break;
    }

    const std::vector<bool> crossed =
        fast_surface_intersect_octants(leaf_boxes, surface, config.bound_width);

    bool any_marked = false;
    for (const bool flag : crossed) {
      if (flag) {
        any_marked = true;
        break;
      }
    }
    int local_refine = any_marked ? 1 : 0;
    int global_refine = 0;
    MPI_Allreduce(&local_refine, &global_refine, 1, MPI_INT, MPI_MAX, comm);
    if (global_refine == 0 || finest >= config.max_level) {
      break;
    }

    std::unordered_set<OctantId> mark;
    for (std::size_t i = 0; i < leaf_ids.size(); ++i) {
      if (crossed[i]) {
        mark.insert(leaf_ids[i]);
      }
    }
    RefineCollective(
        forest, [&](OctantId id) { return mark.count(id) > 0; },
        config.max_level, comm);

    if (iteration == 0) {
      forest.partition(make_level_weight_fn(forest));
    }
  }
}

}  // namespace octlb
