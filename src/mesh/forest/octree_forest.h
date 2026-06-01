#ifndef OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_H_
#define OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_H_

#include <functional>
#include <memory>

#include <mpi.h>

#include "src/common/bounding_box.h"
#include "src/common/types.h"

namespace octlb {

/** Parallel octree forest (p4est backend); topology only, no physics. */
class OctreeForest {
 public:
  OctreeForest(MPI_Comm comm, BoundingBox domain);
  ~OctreeForest();

  OctreeForest(const OctreeForest&) = delete;
  OctreeForest& operator=(const OctreeForest&) = delete;
  OctreeForest(OctreeForest&&) noexcept;
  OctreeForest& operator=(OctreeForest&&) noexcept;

  void refine(std::function<bool(OctantId)> criterion, int max_level);
  void balance();
  void partition();

  label local_num_octants() const;
  BoundingBox quadrant_bounds(OctantId id) const;
  int quadrant_level(OctantId id) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_H_
