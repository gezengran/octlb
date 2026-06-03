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

  Value* send_buf(std::size_t i) {
    return reinterpret_cast<Value*>(&send_bufs_[i * max_elems_]);
  }

  Value* recv_ptr(std::size_t i) {
    return reinterpret_cast<Value*>(&recv_bufs_[i * max_elems_]);
  }

  MPI_Comm comm_;
  BlockCollection<T>& blocks_;
  int nx_, ny_, nz_;
  int max_elems_;
  std::vector<Entry> entries_;
  std::vector<Value> send_bufs_;
  std::vector<Value> recv_bufs_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_FIELD_GHOST_SCHEDULE_H_
