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
  friend struct MeshForestAccess;
 public:
  OctreeForest(MPI_Comm comm, BoundingBox domain);
  ~OctreeForest();

  OctreeForest(const OctreeForest&) = delete;
  OctreeForest& operator=(const OctreeForest&) = delete;
  OctreeForest(OctreeForest&&) noexcept;
  OctreeForest& operator=(OctreeForest&&) noexcept;

  void refine(std::function<bool(OctantId)> criterion, int max_level);
  void balance();
  /** Uniform partition when \a weight_fn is null; else \c p8est_partition_ext. */
  void partition(std::function<int(OctantId)> weight_fn = nullptr);

  label local_num_octants() const;
  BoundingBox quadrant_bounds(OctantId id) const;
  int quadrant_level(OctantId id) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;

  static void RebuildGhostLayer(Impl* impl);
};

}  // namespace octlb

#endif  // OCTLB_SRC_MESH_FOREST_OCTREE_FOREST_H_
