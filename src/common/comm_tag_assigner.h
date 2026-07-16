#ifndef OCTLB_SRC_COMMON_COMM_TAG_ASSIGNER_H_
#define OCTLB_SRC_COMMON_COMM_TAG_ASSIGNER_H_

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace octlb {

// CommTagAssigner assigns MPI communication tags to cross-rank face pairs
// deterministically per rank-pair: for each peer rank, the face keys shared
// with that peer are sorted and each face receives tag = (index + 1) in that
// sorted order.
//
// Both endpoints of a cross-rank face observe the same set of faces for a given
// peer (a face between ranks A and B is local on one side and a ghost on the
// other), so both compute the same sorted order and hence the same tag: tags
// are symmetric across ranks. Tags are unique within a peer; MPI matches
// messages on (src, dst, tag), so per-peer uniqueness is sufficient and tags
// may be reused across different rank-pairs without ambiguity.
//
// Pure / dependency-free: MPI_TAG_UB is taken as a constructor int; this class
// does not call MPI or p4est. The caller queries MPI_TAG_UB and passes it in.
class CommTagAssigner {
 public:
  explicit CommTagAssigner(int tag_ub) : tag_ub_(tag_ub > 0 ? tag_ub : 1) {}

  // Register a cross-rank face (identified by a symmetric canonical face key)
  // shared with peer_rank. Must be called before assign().
  void add_face(uint64_t face_key, int peer_rank) {
    pending_by_peer_[peer_rank].push_back(face_key);
  }

  // Compute per-peer sorted tags for every added face. After this returns,
  // tag_for() returns the assigned tag for each added face key.
  void assign() {
    for (auto& [peer, keys] : pending_by_peer_) {
      std::sort(keys.begin(), keys.end());
      keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
      for (std::size_t i = 0; i < keys.size(); ++i) {
        int tag = static_cast<int>(i) + 1;
        if (tag > tag_ub_) {
          tag = tag_ub_;  // Defensive: faces-per-peer exceeds MPI_TAG_UB.
        }
        tag_by_face_key_[keys[i]] = tag;
      }
    }
  }

  // Returns the assigned tag for face_key, or -1 if unknown / pre-assign().
  int tag_for(uint64_t face_key) const {
    const auto it = tag_by_face_key_.find(face_key);
    return it == tag_by_face_key_.end() ? -1 : it->second;
  }

  int tag_ub() const { return tag_ub_; }

 private:
  int tag_ub_;
  std::unordered_map<int, std::vector<uint64_t>> pending_by_peer_;
  std::unordered_map<uint64_t, int> tag_by_face_key_;
};

}  // namespace octlb

#endif  // OCTLB_SRC_COMMON_COMM_TAG_ASSIGNER_H_