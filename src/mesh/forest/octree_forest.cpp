#include "src/mesh/forest/octree_forest.h"

#include <mutex>
#include <stdexcept>
#include <utility>

#include <p8est.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <sc.h>

#include "src/mesh/forest/octree_forest_access.h"

namespace octlb {
namespace {

void EnsureP4estInitialized(MPI_Comm comm) {
  static std::once_flag once;
  std::call_once(once, [comm]() {
    sc_init(comm, 0, 1, nullptr, SC_LP_ERROR);
    p4est_init(nullptr, SC_LP_ERROR);
  });
}

OctantId LocalOctantIndex(const p8est_t* forest, p4est_topidx_t which_tree,
                          const p8est_quadrant_t* quadrant) {
  p8est_tree_t* tree =
      p8est_tree_array_index(forest->trees, which_tree);
  p8est_quadrant_t* first =
      p8est_quadrant_array_index(&tree->quadrants, 0);
  return static_cast<OctantId>(quadrant - first);
}

BoundingBox ScaleVertexBox(const BoundingBox& domain, const double v_lo[3],
                           const double v_hi[3]) {
  const scalar sx = domain.x_max - domain.x_min;
  const scalar sy = domain.y_max - domain.y_min;
  const scalar sz = domain.z_max - domain.z_min;
  BoundingBox out;
  out.x_min = domain.x_min + static_cast<scalar>(v_lo[0]) * sx;
  out.y_min = domain.y_min + static_cast<scalar>(v_lo[1]) * sy;
  out.z_min = domain.z_min + static_cast<scalar>(v_lo[2]) * sz;
  out.x_max = domain.x_min + static_cast<scalar>(v_hi[0]) * sx;
  out.y_max = domain.y_min + static_cast<scalar>(v_hi[1]) * sy;
  out.z_max = domain.z_min + static_cast<scalar>(v_hi[2]) * sz;
  return out;
}

struct RefineBridge {
  std::function<bool(OctantId)> criterion;
  int max_level = 0;
};

struct PartitionBridge {
  std::function<int(OctantId)> weight_fn;
};

int PartitionWeightCallback(p8est_t* forest, p4est_topidx_t which_tree,
                            p8est_quadrant_t* quadrant) {
  auto* bridge = static_cast<PartitionBridge*>(forest->user_pointer);
  const OctantId id = LocalOctantIndex(forest, which_tree, quadrant);
  return bridge->weight_fn(id);
}

int RefineCallback(p8est_t* forest, p4est_topidx_t which_tree,
                   p8est_quadrant_t* quadrant) {
  auto* bridge = static_cast<RefineBridge*>(forest->user_pointer);
  if (quadrant->level >= static_cast<int8_t>(bridge->max_level)) {
    return 0;
  }
  const OctantId id = LocalOctantIndex(forest, which_tree, quadrant);
  return bridge->criterion(id) ? 1 : 0;
}

}  // namespace

struct OctreeForest::Impl {
  MPI_Comm comm = MPI_COMM_NULL;
  BoundingBox domain{};
  p8est_connectivity_t* connectivity = nullptr;
  p8est_t* forest = nullptr;
  p8est_ghost_t* ghost = nullptr;
  RefineBridge refine_bridge{};
  PartitionBridge partition_bridge{};

  ~Impl() {
    if (ghost != nullptr) {
      p8est_ghost_destroy(ghost);
      ghost = nullptr;
    }
    if (forest != nullptr) {
      p8est_destroy(forest);
      forest = nullptr;
    }
    if (connectivity != nullptr) {
      p8est_connectivity_destroy(connectivity);
      connectivity = nullptr;
    }
  }
};

void OctreeForest::RebuildGhostLayer(Impl* impl) {
  if (impl->ghost != nullptr) {
    p8est_ghost_destroy(impl->ghost);
    impl->ghost = nullptr;
  }
  impl->ghost = p8est_ghost_new(impl->forest, P8EST_CONNECT_FACE);
  if (impl->ghost == nullptr) {
    throw std::runtime_error("p8est_ghost_new failed");
  }
}

OctreeForest::OctreeForest(MPI_Comm comm, BoundingBox domain)
    : impl_(std::make_unique<Impl>()) {
  EnsureP4estInitialized(comm);
  impl_->comm = comm;
  impl_->domain = domain;
  impl_->connectivity = p8est_connectivity_new_unitcube();
  if (impl_->connectivity == nullptr) {
    throw std::runtime_error("p8est_connectivity_new_unitcube failed");
  }
  impl_->forest = p8est_new(comm, impl_->connectivity, 0, nullptr, nullptr);
  if (impl_->forest == nullptr) {
    throw std::runtime_error("p8est_new failed");
  }
}

OctreeForest::~OctreeForest() = default;

OctreeForest::OctreeForest(OctreeForest&&) noexcept = default;
OctreeForest& OctreeForest::operator=(OctreeForest&&) noexcept = default;

void OctreeForest::refine(std::function<bool(OctantId)> criterion,
                          int max_level) {
  if (!criterion) {
    throw std::invalid_argument("refine criterion must not be null");
  }
  impl_->refine_bridge.criterion = std::move(criterion);
  impl_->refine_bridge.max_level = max_level;
  impl_->forest->user_pointer = &impl_->refine_bridge;
  p8est_refine_ext(impl_->forest, 1, max_level, RefineCallback, nullptr,
                   nullptr);
}

void OctreeForest::balance() {
  p8est_balance(impl_->forest, P8EST_CONNECT_FACE, nullptr);
}

void OctreeForest::partition(std::function<int(OctantId)> weight_fn) {
  if (weight_fn) {
    impl_->partition_bridge.weight_fn = std::move(weight_fn);
    impl_->forest->user_pointer = &impl_->partition_bridge;
    p8est_partition_ext(impl_->forest, 0, PartitionWeightCallback);
  } else {
    p8est_partition(impl_->forest, 0, nullptr);
  }
  OctreeForest::RebuildGhostLayer(impl_.get());
}

label OctreeForest::local_num_octants() const {
  return impl_->forest->local_num_quadrants;
}

BoundingBox OctreeForest::quadrant_bounds(OctantId id) const {
  if (id < 0 || id >= local_num_octants()) {
    throw std::out_of_range("OctantId out of range");
  }
  const p4est_topidx_t treeid = impl_->forest->first_local_tree;
  p8est_tree_t* tree = p8est_tree_array_index(impl_->forest->trees, treeid);
  p8est_quadrant_t* q = p8est_quadrant_array_index(&tree->quadrants, id);
  const p4est_qcoord_t len = P8EST_QUADRANT_LEN(q->level);
  double v_lo[3];
  double v_hi[3];
  p8est_qcoord_to_vertex(impl_->connectivity, treeid, q->x, q->y, q->z, v_lo);
  p8est_qcoord_to_vertex(impl_->connectivity, treeid,
                         static_cast<p4est_qcoord_t>(q->x + len),
                         static_cast<p4est_qcoord_t>(q->y + len),
                         static_cast<p4est_qcoord_t>(q->z + len), v_hi);
  return ScaleVertexBox(impl_->domain, v_lo, v_hi);
}

int OctreeForest::quadrant_level(OctantId id) const {
  if (id < 0 || id >= local_num_octants()) {
    throw std::out_of_range("OctantId out of range");
  }
  p8est_tree_t* tree = p8est_tree_array_index(
      impl_->forest->trees, impl_->forest->first_local_tree);
  p8est_quadrant_t* q = p8est_quadrant_array_index(&tree->quadrants, id);
  return static_cast<int>(q->level);
}

p8est_t* MeshForestAccess::Forest(const OctreeForest& forest) {
  return forest.impl_->forest;
}

p8est_ghost_t* MeshForestAccess::Ghost(const OctreeForest& forest) {
  return forest.impl_->ghost;
}

}  // namespace octlb
