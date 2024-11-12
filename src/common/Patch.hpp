/* _____________________________________________________________________ */
//! \file Patch.hpp

//! \brief Header for the patch class definition

/*

  SCHEMATICS: patch decomposition
  ____________________________________________________

         patch primal origin (ix_origin_m) = i_patch_topology_m * nx_cells_per_patch
         |
         |    nx_cells_per_patch : number of cells
         |    nx_p_m : number of points in the primal grid
         |< -------------- >
         v
         (--|--|--|--|--|--)                                       <- patch primal 0
                           (--|--|--|--|--|--)                     <- patch primal 1
                                             (--|--|--|--|--|--)
                           |                 |
    P    [--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--]   <- primal current global grid
                           |                 |
    D   [--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--]    <- dual current global grid

        (--|--|--|--|--|--|--)                                        <- patch dual 0
        ^                 (--|--|--|--|--|--|--)                      <- patch dual 1
        |                                   (--|--|--|--|--|--|--)
        |
        |
        patch dual origin = i_patch_topology_m * nx_cells_per_patch

      Legend :

      |--|  -> inner cell
      [--| ... |--]  -> global grid
      (--| ... |--)  -> patch grid

      * Particles are contained in the primal boundaries in each direction (not the dual one even if
  the direction is dual)

*/

/* _____________________________________________________________________ */

#pragma once

#include <iostream>

#include "ElectroMagn.hpp"
#include "Field.hpp"
#include "Params.hpp"
#include "Particles.hpp"
#include "Tools.hpp"

//! Represents a subspace in the domain
class Patch {
public:
  Patch() {}
  ~Patch() {}

  //! Number of species in the simulation at init
  int n_species_m;
  //! Number of patch in each direction, global
  int nx_patchs_m, ny_patchs_m, nz_patchs_m;
  //! Number of cells in a patch in each direction, local primal grid
  int nx_cells_m, ny_cells_m, nz_cells_m;

  //! Number of nodes in a patch in each direction, local primal grid
  int nx_p_m, ny_p_m, nz_p_m;
  //! Number of nodes in a patch in each direction, local dual grid
  int nx_d_m, ny_d_m, nz_d_m;

  //! Primal position index of the patch origin in the global space
  int ix_origin_m, iy_origin_m, iz_origin_m;

  //! Index in each direction in the patch topology, local
  int i_patch_topology_m, j_patch_topology_m, k_patch_topology_m;
  //! Index global in the patch topology, computed
  int idx_patch_topology_m;
  //! Flag to determine if the patch is on the border
  bool on_border_m;
  //! Store as a code the type of boundary condition to use for particles
  //! 0 - free
  //! 1 - periodic
  //! 2 - reflective
  int type_border_m;

  //! Boundaries box of the patch, computed
  double inf_m[3];
  double sup_m[3];

  //! Current particles in the patch
  std::vector<Particles<mini_float>> particles_m;
  //! Local current field, one by species
  std::vector<Field<mini_float>> vec_Jx_m, vec_Jy_m, vec_Jz_m;

  //! Flags to determine if the species is projected or not
  std::vector<bool> projected_;

  //! Local current field
  Field<mini_float> Jx_m, Jy_m, Jz_m;

  //! Buffers to exchange particles between patchs
  std::vector<std::vector<Particles<mini_float>>> particles_to_move_m;
  //! Tag particles which leave the patch
  //   std::vector<std::vector<int>> masks_m;

  // _____________________________________________________________
  //
  //! \brief Initialize and allocate properties for the patch
  // _____________________________________________________________
  void allocate(Params &param, Backend &backend, const int i, const int j, const int k);

  // _____________________________________________________________
  //
  //! \brief Initialize and allocate properties for the patch
  // _____________________________________________________________
  void initialize_particles(Params &param);

  // _____________________________________________________________
  //
  //! \brief Give the index patch topology from 3D indexes,
  //!       -1 if out of domain
  // _____________________________________________________________
  INLINE int get_idx_patch(int i, int j, int k) {

    int ixp = i;
    int iyp = j;
    int izp = k;

    // Periodic management of the topology
    if (ixp < 0) {
      ixp = nx_patchs_m - 1;
    } else if (ixp >= nx_patchs_m) {
      ixp = 0;
    }

    if (iyp < 0) {
      iyp = ny_patchs_m - 1;
    } else if (iyp >= ny_patchs_m) {
      iyp = 0;
    }

    if (izp < 0) {
      izp = nz_patchs_m - 1;
    } else if (izp >= nz_patchs_m) {
      izp = 0;
    }

    // Reflective management of the topology
    // if(i<0 || i>=nx_patchs_m)
    //   return -1;
    // if(j<0 || j>=ny_patchs_m)
    //   return -1;
    // if(k<0 || k>=nz_patchs_m)
    //   return -1;

    return ixp * nz_patchs_m * ny_patchs_m + iyp * nz_patchs_m + izp;
  }

  // __________________________________________________
  //
  //! \brief performs some checks for debugging
  // __________________________________________________
  void check() {

    // Identify particles to move
    for (int is = 0; is < n_species_m; is++) {

      // Number of particles for this species is
      unsigned int n_particles = particles_m[is].size();

      for (unsigned int ip = 0; ip < n_particles; ip++) {

        if ((particles_m[is].x_h(ip) < inf_m[0]) || (particles_m[is].x_h(ip) >= sup_m[0]) or
            (particles_m[is].y_h(ip) < inf_m[1]) or (particles_m[is].y_h(ip) >= sup_m[1]) or
            (particles_m[is].z_h(ip) < inf_m[2]) or (particles_m[is].z_h(ip) >= sup_m[2])) {
          std::cerr << "Problem:" << std::endl;
          std::cerr << "Patch: " << i_patch_topology_m << " " << j_patch_topology_m << " "
                    << k_patch_topology_m << " xmin: " << inf_m[0] << " xmax: " << sup_m[0]
                    << std::endl;
          std::cerr << " Species: " << is << " Particles " << ip << "/" << n_particles
                    << " - x: " << particles_m[is].x_h(ip) << " y: " << particles_m[is].y_h(ip)
                    << " z: " << particles_m[is].z_h(ip) << std::endl;
          std::raise(SIGABRT);
        }
      }
    }
  }

  // __________________________________________________
  //
  //! Return the total number of particles in the patch
  // __________________________________________________
  unsigned int get_total_number_of_particles() {
    unsigned int sum = 0;
    for (int is = 0; is < n_species_m; is++) {
      sum += particles_m[is].size();
    }
    return sum;
  }
};