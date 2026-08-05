#include "src/mesh/forest/octree_forest.h"

#include <algorithm>
#include <limits>
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

struct LocalQuadrantRef {
  p4est_topidx_t treeid = 0;
  p4est_locidx_t quad_idx = 0;
};

LocalQuadrantRef ResolveLocalQuadrant(const p8est_t* forest, OctantId id) {
  p4est_locidx_t cursor = 0;
  for (p4est_topidx_t tt = forest->first_local_tree; tt <= forest->last_local_tree;
       ++tt) {
    p8est_tree_t* tree = p8est_tree_array_index(forest->trees, tt);
    const p4est_locidx_t n =
        static_cast<p4est_locidx_t>(tree->quadrants.elem_count);
    if (id < cursor + n) {
      return {tt, id - cursor};
    }
    cursor += n;
  }
  throw std::out_of_range("OctantId out of range");
}

BoundingBox ConnectivityVertexExtents(const p8est_connectivity_t* conn) {
  BoundingBox ext;
  ext.x_min = ext.y_min = ext.z_min = std::numeric_limits<scalar>::infinity();
  ext.x_max = ext.y_max = ext.z_max = -std::numeric_limits<scalar>::infinity();
  for (p4est_topidx_t i = 0; i < conn->num_vertices; ++i) {
    const double* v = conn->vertices + 3 * static_cast<std::size_t>(i);
    ext.x_min = std::min(ext.x_min, static_cast<scalar>(v[0]));
    ext.y_min = std::min(ext.y_min, static_cast<scalar>(v[1]));
    ext.z_min = std::min(ext.z_min, static_cast<scalar>(v[2]));
    ext.x_max = std::max(ext.x_max, static_cast<scalar>(v[0]));
    ext.y_max = std::max(ext.y_max, static_cast<scalar>(v[1]));
    ext.z_max = std::max(ext.z_max, static_cast<scalar>(v[2]));
  }
  return ext;
}

scalar NormalizeVertexCoord(double v, scalar ext_min, scalar ext_max) {
  const scalar span = ext_max - ext_min;
  if (span <= std::numeric_limits<scalar>::epsilon()) {
    return scalar{0};
  }
  return (static_cast<scalar>(v) - ext_min) / span;
}

BoundingBox ScaleVertexBox(const BoundingBox& domain, const BoundingBox& vtx_ext,
                           const double v_lo[3], const double v_hi[3]) {
  const scalar sx = domain.x_max - domain.x_min;
  const scalar sy = domain.y_max - domain.y_min;
  const scalar sz = domain.z_max - domain.z_min;
  BoundingBox out;
  out.x_min =
      domain.x_min +
      NormalizeVertexCoord(v_lo[0], vtx_ext.x_min, vtx_ext.x_max) * sx;
  out.y_min =
      domain.y_min +
      NormalizeVertexCoord(v_lo[1], vtx_ext.y_min, vtx_ext.y_max) * sy;
  out.z_min =
      domain.z_min +
      NormalizeVertexCoord(v_lo[2], vtx_ext.z_min, vtx_ext.z_max) * sz;
  out.x_max =
      domain.x_min +
      NormalizeVertexCoord(v_hi[0], vtx_ext.x_min, vtx_ext.x_max) * sx;
  out.y_max =
      domain.y_min +
      NormalizeVertexCoord(v_hi[1], vtx_ext.y_min, vtx_ext.y_max) * sy;
  out.z_max =
      domain.z_min +
      NormalizeVertexCoord(v_hi[2], vtx_ext.z_min, vtx_ext.z_max) * sz;
  return out;
}

struct RefineBridge {
  std::function<bool(OctantId)> criterion;
  int max_level = 0;
};

struct RefineBoundsBridge {
  std::function<bool(const BoundingBox&, int)> criterion;
  int max_level = 0;
  BoundingBox domain{};
  BoundingBox vtx_ext{};
};

int RefineBoundsCallback(p8est_t* forest, p4est_topidx_t which_tree,
                         p8est_quadrant_t* quadrant) {
  auto* bridge = static_cast<RefineBoundsBridge*>(forest->user_pointer);
  if (quadrant->level >= static_cast<int8_t>(bridge->max_level)) {
    return 0;
  }
  const p4est_qcoord_t len = P8EST_QUADRANT_LEN(quadrant->level);
  double v_lo[3];
  double v_hi[3];
  p8est_qcoord_to_vertex(forest->connectivity, which_tree, quadrant->x,
                         quadrant->y, quadrant->z, v_lo);
  p8est_qcoord_to_vertex(forest->connectivity, which_tree,
                         static_cast<p4est_qcoord_t>(quadrant->x + len),
                         static_cast<p4est_qcoord_t>(quadrant->y + len),
                         static_cast<p4est_qcoord_t>(quadrant->z + len), v_hi);
  const BoundingBox box =
      ScaleVertexBox(bridge->domain, bridge->vtx_ext, v_lo, v_hi);
  return bridge->criterion(box, static_cast<int>(quadrant->level)) ? 1 : 0;
}

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
  RefineBoundsBridge refine_bounds_bridge{};
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

void OctreeForest::RebuildGhost() { RebuildGhostLayer(impl_.get()); }

MPI_Comm OctreeForest::comm() const { return impl_->comm; }

OctreeForest::OctreeForest(MPI_Comm comm, BoundingBox domain, int bricks_x,
                           int bricks_y, int bricks_z)
    : impl_(std::make_unique<Impl>()) {
  if (bricks_x < 1 || bricks_y < 1 || bricks_z < 1) {
    throw std::invalid_argument("brick dimensions must be >= 1");
  }
  EnsureP4estInitialized(comm);
  impl_->comm = comm;
  impl_->domain = domain;
  if (bricks_x == 1 && bricks_y == 1 && bricks_z == 1) {
    impl_->connectivity = p8est_connectivity_new_unitcube();
    if (impl_->connectivity == nullptr) {
      throw std::runtime_error("p8est_connectivity_new_unitcube failed");
    }
  } else {
    impl_->connectivity =
        p8est_connectivity_new_brick(bricks_x, bricks_y, bricks_z, 0, 0, 0);
    if (impl_->connectivity == nullptr) {
      throw std::runtime_error("p8est_connectivity_new_brick failed");
    }
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

void OctreeForest::refine(std::function<bool(const BoundingBox&, int)> criterion,
                          int max_level) {
  if (!criterion) {
    throw std::invalid_argument("refine criterion must not be null");
  }
  impl_->refine_bounds_bridge.criterion = std::move(criterion);
  impl_->refine_bounds_bridge.max_level = max_level;
  impl_->refine_bounds_bridge.domain = impl_->domain;
  impl_->refine_bounds_bridge.vtx_ext =
      ConnectivityVertexExtents(impl_->connectivity);
  impl_->forest->user_pointer = &impl_->refine_bounds_bridge;
  p8est_refine_ext(impl_->forest, 1, max_level, RefineBoundsCallback, nullptr,
                   nullptr);
}

void OctreeForest::balance() {
  // Face balance. (Edge balance -- P8EST_CONNECT_EDGE -- was tried for the ②
  // Stage B edge callback, but it breaks multi-rank AMR face-pair symmetry:
  // the face-connected ghost then enumerates cross-rank faces asymmetrically
  // across ranks -> MPI deadlock in GhostSchedule::exchange. The edge callback
  // does NOT need forest-level edge balance -- FacePairList builds a separate
  // P8EST_CONNECT_EDGE ghost for the edge pass, so the callback fires there.)
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
  const LocalQuadrantRef ref = ResolveLocalQuadrant(impl_->forest, id);
  p8est_tree_t* tree =
      p8est_tree_array_index(impl_->forest->trees, ref.treeid);
  p8est_quadrant_t* q =
      p8est_quadrant_array_index(&tree->quadrants, ref.quad_idx);
  const p4est_qcoord_t len = P8EST_QUADRANT_LEN(q->level);
  double v_lo[3];
  double v_hi[3];
  p8est_qcoord_to_vertex(impl_->connectivity, ref.treeid, q->x, q->y, q->z,
                         v_lo);
  p8est_qcoord_to_vertex(impl_->connectivity, ref.treeid,
                         static_cast<p4est_qcoord_t>(q->x + len),
                         static_cast<p4est_qcoord_t>(q->y + len),
                         static_cast<p4est_qcoord_t>(q->z + len), v_hi);
  const BoundingBox vtx_ext = ConnectivityVertexExtents(impl_->connectivity);
  return ScaleVertexBox(impl_->domain, vtx_ext, v_lo, v_hi);
}

int OctreeForest::quadrant_level(OctantId id) const {
  if (id < 0 || id >= local_num_octants()) {
    throw std::out_of_range("OctantId out of range");
  }
  const LocalQuadrantRef ref = ResolveLocalQuadrant(impl_->forest, id);
  p8est_tree_t* tree =
      p8est_tree_array_index(impl_->forest->trees, ref.treeid);
  p8est_quadrant_t* q =
      p8est_quadrant_array_index(&tree->quadrants, ref.quad_idx);
  return static_cast<int>(q->level);
}

p8est_t* MeshForestAccess::Forest(const OctreeForest& forest) {
  return forest.impl_->forest;
}

p8est_ghost_t* MeshForestAccess::Ghost(const OctreeForest& forest) {
  return forest.impl_->ghost;
}

BoundingBox MeshForestAccess::QuadrantBounds(const OctreeForest& forest,
                                            p4est_topidx_t treeid,
                                            const p8est_quadrant_t* q) {
  const p4est_qcoord_t len = P8EST_QUADRANT_LEN(q->level);
  double v_lo[3];
  double v_hi[3];
  p8est_qcoord_to_vertex(forest.impl_->connectivity, treeid, q->x, q->y, q->z,
                         v_lo);
  p8est_qcoord_to_vertex(forest.impl_->connectivity, treeid,
                         static_cast<p4est_qcoord_t>(q->x + len),
                         static_cast<p4est_qcoord_t>(q->y + len),
                         static_cast<p4est_qcoord_t>(q->z + len), v_hi);
  const BoundingBox vtx_ext = ConnectivityVertexExtents(forest.impl_->connectivity);
  return ScaleVertexBox(forest.impl_->domain, vtx_ext, v_lo, v_hi);
}

}  // namespace octlb
