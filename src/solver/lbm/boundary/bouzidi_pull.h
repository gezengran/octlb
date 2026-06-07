#ifndef OCTLB_SRC_SOLVER_LBM_BOUNDARY_BOUZIDI_PULL_H_
#define OCTLB_SRC_SOLVER_LBM_BOUNDARY_BOUZIDI_PULL_H_

namespace octlb {
namespace boundary {

/// Bouzidi interpolated bounce-back (post-collision pull from solid direction).
template <typename T>
inline T BouzidiPostCollisionPull(T f_bb_opp_same_cell, T f_q_interior_post,
                                  T f_same_dir_post_same_cell, double q_frac) {
  if (q_frac >= 1.0 - 1.0e-14) {
    return f_bb_opp_same_cell;
  }
  if (q_frac <= 1.0e-14) {
    return f_bb_opp_same_cell;
  }
  if (q_frac >= 0.5 - 1.0e-14) {
    const double inv2q = 1.0 / (2.0 * q_frac);
    return static_cast<T>(inv2q * static_cast<double>(f_bb_opp_same_cell) +
                          (1.0 - inv2q) *
                              static_cast<double>(f_q_interior_post));
  }
  return static_cast<T>(
      2.0 * q_frac * static_cast<double>(f_bb_opp_same_cell) +
      (1.0 - 2.0 * q_frac) * static_cast<double>(f_same_dir_post_same_cell));
}

}  // namespace boundary
}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_BOUNDARY_BOUZIDI_PULL_H_
