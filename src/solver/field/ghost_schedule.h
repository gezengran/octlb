#ifndef OCTLB_SRC_SOLVER_FIELD_GHOST_SCHEDULE_H_
#define OCTLB_SRC_SOLVER_FIELD_GHOST_SCHEDULE_H_

#include <cstddef>
#include <cstring>
#include <vector>

#include <mpi.h>

#include "src/mesh/topology/face_pair_list.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/field/face_packable.h"

namespace octlb {

template <FacePackable T>
class GhostSchedule {
 public:
  using Value = typename T::face_value_t;

  GhostSchedule(MPI_Comm comm, const FacePairList& faces,
                BlockCollection<T>& blocks, int nx, int ny, int nz)
      : comm_(comm),
        blocks_(blocks),
        nx_(nx),
        ny_(ny),
        nz_(nz),
        max_elems_(MaxFaceBufferCount(
            nx, ny, nz,
            [](int a, int b, int c, FaceDir d) {
              return T::face_buffer_count(a, b, c, d);
            })) {
    int my_rank = 0;
    MPI_Comm_rank(comm_, &my_rank);

    for (const SameLevelFace& face : faces.same_level_faces()) {
      Entry e{};
      e.local_id = face.local_id;
      e.dir = face.dir;
      e.remote_id = face.remote_id;
      e.remote_rank = face.remote_rank;
      // Symmetric geometric tag (same on both ranks of a face pair).
      e.comm_tag = face.comm_tag;
      e.elem_count = T::face_buffer_count(nx_, ny_, nz_, face.dir);
      e.is_local = (face.remote_rank == my_rank);
      entries_.push_back(e);
    }

    send_bufs_.resize(entries_.size() * static_cast<std::size_t>(max_elems_));
    recv_bufs_.resize(entries_.size() * static_cast<std::size_t>(max_elems_));

    // ② edge-ghost: build edge adjacency by composing same-rank same-level
    // face neighbours (Stage A: same-rank only; cross-rank edges need the
    // p4est corner callback, deferred to Stage B). An edge is two orthogonal
    // faces (d1, d2); the diagonal neighbour is face_nbr[face_nbr[A][d1]][d2].
    const label num_blocks = blocks_.size();
    std::vector<OctantId> face_nbr(
        static_cast<std::size_t>(num_blocks) * 6, OctantId{-1});
    for (const SameLevelFace& f : faces.same_level_faces()) {
      if (f.remote_rank == my_rank) {
        face_nbr[static_cast<std::size_t>(f.local_id) * 6 +
                 static_cast<int>(f.dir)] = f.remote_id;
      }
    }
    for (OctantId a = 0; a < num_blocks; ++a) {
      for (int d1 = 0; d1 < 6; ++d1) {
        for (int d2 = 0; d2 < 6; ++d2) {
          if (d1 / 2 == d2 / 2) continue;  // same axis -> not an edge
          const OctantId b =
              face_nbr[static_cast<std::size_t>(a) * 6 + d1];
          if (b < 0 || b >= num_blocks) continue;
          const OctantId dd =
              face_nbr[static_cast<std::size_t>(b) * 6 + d2];
          if (dd < 0 || dd >= num_blocks) continue;
          EdgeEntry ee{};
          ee.local_id = a;
          ee.d1 = static_cast<FaceDir>(d1);
          ee.d2 = static_cast<FaceDir>(d2);
          ee.remote_id = dd;
          ee.remote_rank = my_rank;
          ee.comm_tag = 0;  // unused for same-rank
          ee.elem_count = T::edge_buffer_count(nx_, ny_, nz_, ee.d1, ee.d2);
          ee.is_local = true;
          edge_entries_.push_back(ee);
        }
      }
    }
    max_edge_elems_ = MaxEdgeBufferCount(
        nx_, ny_, nz_, [](int a2, int b2, int c2, FaceDir ed1, FaceDir ed2) {
          return T::edge_buffer_count(a2, b2, c2, ed1, ed2);
        });
    edge_send_bufs_.resize(edge_entries_.size() *
                           static_cast<std::size_t>(max_edge_elems_));
    edge_recv_bufs_.resize(edge_entries_.size() *
                           static_cast<std::size_t>(max_edge_elems_));
  }

  void exchange() {
    std::vector<MPI_Request> requests;
    requests.reserve(entries_.size() * 2);

    for (std::size_t i = 0; i < entries_.size(); ++i) {
      const Entry& e = entries_[i];
      Value* send_ptr = send_buf(i);
      blocks_[e.local_id].pack_face(e.dir, send_ptr, e.elem_count);

      if (e.is_local) {
        blocks_[e.remote_id].unpack_face(OppositeFace(e.dir), send_ptr,
                                         e.elem_count);
      } else {
        MPI_Request req{};
        MPI_Irecv(recv_ptr(i), static_cast<int>(e.elem_count * sizeof(Value)),
                  MPI_BYTE, e.remote_rank, e.comm_tag, comm_, &req);
        requests.push_back(req);
      }
    }

    for (std::size_t i = 0; i < entries_.size(); ++i) {
      const Entry& e = entries_[i];
      if (e.is_local) {
        continue;
      }
      MPI_Request req{};
      MPI_Isend(send_buf(i), static_cast<int>(e.elem_count * sizeof(Value)),
                MPI_BYTE, e.remote_rank, e.comm_tag, comm_, &req);
      requests.push_back(req);
    }

    if (!requests.empty()) {
      MPI_Waitall(static_cast<int>(requests.size()), requests.data(),
                  MPI_STATUSES_IGNORE);
    }

    for (std::size_t i = 0; i < entries_.size(); ++i) {
      const Entry& e = entries_[i];
      if (e.is_local) {
        continue;
      }
      blocks_[e.local_id].unpack_face(e.dir, recv_ptr(i), e.elem_count);
    }

    // ② edge-ghost exchange. Stage A: same-rank local copy (pack local interior
    // edge -> unpack into the diagonal neighbour's opposite edge ghost). The
    // reverse direction is a separate edge entry, so both ghosts get filled.
    for (std::size_t i = 0; i < edge_entries_.size(); ++i) {
      const EdgeEntry& e = edge_entries_[i];
      Value* send = edge_send_buf(i);
      blocks_[e.local_id].pack_edge(e.d1, e.d2, send, e.elem_count);
      if (e.is_local) {
        blocks_[e.remote_id].unpack_edge(OppositeFace(e.d1), OppositeFace(e.d2),
                                         send, e.elem_count);
      }
      // Stage B: cross-rank edge MPI (Isend/Irecv via EdgePairList tags).
    }
  }

 private:
  struct Entry {
    OctantId local_id{};
    FaceDir dir{};
    OctantId remote_id{};
    int remote_rank{};
    int comm_tag{};
    int elem_count{};
    bool is_local{};
  };

  // ② edge-ghost: an edge is the intersection of two orthogonal faces (d1, d2);
  // the diagonal neighbour across the edge is reached by composing two same
  // face neighbours (Stage A: same-rank only). pack reads the local interior
  // edge line, unpack writes the diagonal neighbour's opposite edge ghost.
  struct EdgeEntry {
    OctantId local_id{};
    FaceDir d1{};
    FaceDir d2{};
    OctantId remote_id{};
    int remote_rank{};
    int comm_tag{};
    int elem_count{};
    bool is_local{};
  };

  Value* send_buf(std::size_t i) {
    return reinterpret_cast<Value*>(&send_bufs_[i * max_elems_]);
  }

  Value* recv_ptr(std::size_t i) {
    return reinterpret_cast<Value*>(&recv_bufs_[i * max_elems_]);
  }

  Value* edge_send_buf(std::size_t i) {
    return reinterpret_cast<Value*>(&edge_send_bufs_[i * max_edge_elems_]);
  }
  Value* edge_recv_buf(std::size_t i) {
    return reinterpret_cast<Value*>(&edge_recv_bufs_[i * max_edge_elems_]);
  }

  MPI_Comm comm_;
  BlockCollection<T>& blocks_;
  int nx_, ny_, nz_;
  int max_elems_;
  std::vector<Entry> entries_;
  std::vector<Value> send_bufs_;
  std::vector<Value> recv_bufs_;
  int max_edge_elems_ = 0;
  std::vector<EdgeEntry> edge_entries_;
  std::vector<Value> edge_send_bufs_;
  std::vector<Value> edge_recv_bufs_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_FIELD_GHOST_SCHEDULE_H_
