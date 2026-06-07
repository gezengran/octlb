#ifndef OCTLB_SRC_SOLVER_LBM_LATTICE_MATERIAL_INIT_H_
#define OCTLB_SRC_SOLVER_LBM_LATTICE_MATERIAL_INIT_H_

#include "src/common/types.h"
#include "src/mesh/geometry/material_field.h"
#include "src/solver/field/block_collection.h"
#include "src/solver/lbm/block_lattice.h"

namespace octlb {

using MaterialInitLattice = BlockLattice<double, olb::descriptors::D3Q19<>>;

void initialize_from_material(OctantId id, MaterialInitLattice& lattice,
                              const MaterialField& material, double rho0,
                              const double* u0, double omega);

void initialize_from_material(BlockCollection<MaterialInitLattice>& blocks,
                              const MaterialField& material, double rho0,
                              const double* u0, double omega);

}  // namespace octlb

#endif  // OCTLB_SRC_SOLVER_LBM_LATTICE_MATERIAL_INIT_H_
