#include "src/solver/lbm/lattice_material_init.h"

#include "core/platform/platform.h"
#include "descriptor/descriptor.h"
#include "dynamics/lbm.h"

#include "src/solver/lbm/bc_installer.h"

namespace octlb {

void initialize_from_material(OctantId id, MaterialInitLattice& lattice,
                              const MaterialField& material, double rho0,
                              const double* u0, double /*omega*/) {
  bc::StampFromMaterial(id, lattice, material);
  for (int k = 0; k < material.nz(); ++k) {
    for (int j = 0; j < material.ny(); ++j) {
      for (int i = 0; i < material.nx(); ++i) {
        const MaterialKind mk = material.at(id, i, j, k);
        if (mk == MaterialKind::kFluid) {
          auto cell = lattice.get(i, j, k);
          double uSqr = 0.0;
          for (int d = 0; d < 3; ++d) {
            uSqr += u0[d] * u0[d];
          }
          for (int iPop = 0; iPop < olb::descriptors::D3Q19<>::q; ++iPop) {
            cell[iPop] = olb::equilibrium<olb::descriptors::D3Q19<>>::secondOrder(
                iPop, rho0, u0, uSqr);
          }
        }
      }
    }
  }
}

void initialize_from_material(BlockCollection<MaterialInitLattice>& blocks,
                              const MaterialField& material, double rho0,
                              const double* u0, double omega) {
  for (label oid = 0; oid < material.num_octants(); ++oid) {
    initialize_from_material(static_cast<OctantId>(oid),
                             blocks[static_cast<OctantId>(oid)], material, rho0,
                             u0, omega);
  }
}

}  // namespace octlb
