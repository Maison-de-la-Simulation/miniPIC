
/* _____________________________________________________________________ */
//! \file Patch.cpp

//! \brief Methods to build and manage the patch decomposition

/* _____________________________________________________________________ */

#include "Patch.hpp"

/*

  SCHEMATICS: current local and global grids explained
  ____________________________________________________

      patch primal origin for current = i_patch_topology_m * nx_cells_per_patch - 1 = ix_origin_m -
  1
      |
      |    nx_cells_per_patch : number of cells without ghost
      |    nx_p_m : number of points without ghost
      |  < --------------- >
      v
      (~~|--|--|--|--|--|--|~~)                                       <- patch primal 0
                        (~~|--|--|--|--|--|--|~~)                     <- patch primal 1
                                          (~~|--|--|--|--|--|--|~~)
                           |
      [~~|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|~~]   <- primal current global grid
                           |
        [--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--]    <- dual current global grid

        (--|--|--|--|--|--|--)                                        <- patch dual 0
        ^                 (--|--|--|--|--|--|--)                      <- patch dual 1
        |                                   (--|--|--|--|--|--|--)
        |
        |
        patch dual origin = i_patch_topology_m * nx_cells_per_patch

      Legend :

      |~~|  -> ghost cell
      |--|  -> inner cell
      [--| ... |--]  -> global grid
      (--| ... |--)  -> patch grid

      *  +2 ghostcells for primal grids, therefore the primal grid becomes larger than the dual one
      *  Ghostcell are neede because of the time shift of -1/2dt in the projection, as a
  consequence, some particles can be out of the grid
      *  Also need ghost cell in the patch local grids
      *  Each primal patch share a frontier of points
      *  The dual patch does not have ghost cells
      *  Dual patches share a cell at the frontier

*/

// ________________________________________________________________________
//
//! \brief Alloc memory for all members and init all at null
// ________________________________________________________________________
void Patch::allocate(Params &param, Backend &backend, const int i, const int j, const int k) {
  // Store input
  n_species_m = param.get_species_number();

  nx_patchs_m = param.nx_patch;
  ny_patchs_m = param.ny_patch;
  nz_patchs_m = param.nz_patch;

  nx_cells_m = param.nx_cells_by_patch;
  ny_cells_m = param.ny_cells_by_patch;
  nz_cells_m = param.nz_cells_by_patch;

  nx_p_m = param.nx_p_by_patch;
  ny_p_m = param.ny_p_by_patch;
  nz_p_m = param.nz_p_by_patch;

  nx_d_m = param.nx_d_by_patch;
  ny_d_m = param.ny_d_by_patch;
  nz_d_m = param.nz_d_by_patch;

  // Patch topology
  i_patch_topology_m   = i;
  j_patch_topology_m   = j;
  k_patch_topology_m   = k;
  idx_patch_topology_m = i * param.nz_patch * param.ny_patch + j * param.nz_patch + k;

  // Primal index of the patch origin in the global grid
  ix_origin_m = i_patch_topology_m * param.nx_cells_by_patch;
  iy_origin_m = j_patch_topology_m * param.ny_cells_by_patch;
  iz_origin_m = k_patch_topology_m * param.nz_cells_by_patch;

  // Boundary conditions
  on_border_m = false;
  if ((i_patch_topology_m == 0) || (i_patch_topology_m >= param.nx_patch - 1) ||
      (j_patch_topology_m == 0) || (j_patch_topology_m >= param.ny_patch - 1) ||
      (k_patch_topology_m == 0) || (k_patch_topology_m >= param.nz_patch - 1)) {
    on_border_m = true;
  }

  // Compute boundaries box
  inf_m[0] = param.inf_x + (param.dx * nx_cells_m) * i;
  inf_m[1] = param.inf_y + (param.dy * ny_cells_m) * j;
  inf_m[2] = param.inf_z + (param.dz * nz_cells_m) * k;

  sup_m[0] = inf_m[0] + (param.dx * nx_cells_m);
  sup_m[1] = inf_m[1] + (param.dy * ny_cells_m);
  sup_m[2] = inf_m[2] + (param.dz * nz_cells_m);

  // Alloc vector for each species
  if (n_species_m > 0) {
    particles_m.resize(n_species_m);
    particles_to_move_m.resize(n_species_m);

    // masks_m.resize(n_species_m);
    vec_Jx_m.resize(n_species_m);
    vec_Jy_m.resize(n_species_m);
    vec_Jz_m.resize(n_species_m);

    // Alloc projected flags
    projected_.resize(n_species_m);
  }

  for (int is = 0; is < n_species_m; is++) {
    int n_particles = param.n_particles_by_species_in_patch[is] + param.particles_to_add_.size();

    // Alloc memory to store current particles in patch, and init species attributes
    particles_m[is].allocate(param.charge_[is],
                             param.mass_[is],
                             param.temp_[is],
                             n_particles,
                             param.inv_cell_volume,
                             backend);

    // Alloc a tag for each particle, to tag particles which leave the patch and init it at false
    // masks_m[is].resize(n_particles, 0);

    // Alloc empty buffers to exchange particles between patch, 26 neighbors in 3D
    particles_to_move_m[is].resize(26);

    for (auto ibuffer = 0; ibuffer < 26; ibuffer++) {
      particles_to_move_m[is][ibuffer].allocate(param.charge_[is],
                                                param.mass_[is],
                                                param.temp_[is],
                                                0,
                                                param.inv_cell_volume,
                                                backend);

      particles_to_move_m[is][ibuffer].with_electromagnetic_fields_ = false;
    }

    // Alloc local fields, one for each species
    // Must have 2 ghost cells in the primal direction for the projection
    vec_Jx_m[is].allocate(nx_d_m, ny_p_m + 2, nz_p_m + 2, backend, 0.0, 1, 0, 0, "Jx");
    vec_Jy_m[is].allocate(nx_p_m + 2, ny_d_m, nz_p_m + 2, backend, 0.0, 0, 1, 0, "Jy");
    vec_Jz_m[is].allocate(nx_p_m + 2, ny_p_m + 2, nz_d_m, backend, 0.0, 0, 0, 1, "Jz");
  }

  // Alloc local fields
  // Jx_m.allocate(nx_d_m, ny_p_m, nz_p_m);
  // Jy_m.allocate(nx_p_m, ny_d_m, nz_p_m);
  // Jz_m.allocate(nx_p_m, ny_p_m, nz_d_m);
}

// ___________________________________________________________
//
//! \brief Init all particles in the patch with
//! Maxwell-Juttner distribution
//! \param[in] params constant global parameters
// ___________________________________________________________
void Patch::initialize_particles(Params &param) {

  param.n_particles = 0;

  const double cell_volume = param.cell_volume;

  const int n_species = n_species_m;

  if (n_species == 0) {
    return;
  }

  const int total_cells = nx_cells_m * ny_cells_m * nz_cells_m;

  // buffer to store the number of particles per cells per species
  // Needed for proper init with duplication
  std::vector<int> particles_per_cell_counter(n_species * total_cells);

  // Loop over all species
  for (int is = 0; is < n_species; ++is) {

    int n_particles    = particles_m[is].size();
    double temperature = param.temp_[is];
    const double mass  = param.mass_[is];

    // Compute weight
    const int particle_per_cell = param.ppc_[is];
    const double weight_coef    = cell_volume / particle_per_cell;

    // global particle counter
    unsigned int total_particles_counter = 0;

    // compute the species index for position init
    int species_index_for_pos_init = 0;
    while (param.species_names_[species_index_for_pos_init] !=
             param.position_initialization_method_[is] &&
           (species_index_for_pos_init < is)) {
      ++species_index_for_pos_init;
    }

    // Coefficients for drift velocity
    const double vx = -param.drift_velocity_[is][0];
    const double vy = -param.drift_velocity_[is][1];
    const double vz = -param.drift_velocity_[is][2];

    const double v_drift = vx * vx + vy * vy + vz * vz;

    const double gamma_drift = 1.0 / sqrt(1.0 - v_drift);
    const double gm1         = gamma_drift - 1.0;

    // compute the different component of the Matrix block
    // of the Lorentz transformation (drift velocity correction)
    const double Lxx = 1.0 + gm1 * vx * vx / v_drift;
    const double Lyy = 1.0 + gm1 * vy * vy / v_drift;
    const double Lzz = 1.0 + gm1 * vz * vz / v_drift;
    const double Lxy = gm1 * vx * vy / v_drift;
    const double Lxz = gm1 * vx * vz / v_drift;
    const double Lyz = gm1 * vy * vz / v_drift;

    // If param.position_initialization_method_[is] is random_per_cell, we init particles randomly
    // with a new seed for each cell.
    if (param.position_initialization_level_[is] == "patch") {

      // patch index
      int ipatch = i_patch_topology_m * param.ny_patch * param.nz_patch +
                   j_patch_topology_m * param.nz_patch + k_patch_topology_m;

      // Random generator at the patch level
      Random random(param.seed + is + ipatch * n_species);

      // Random Position init
      if (param.position_initialization_method_[is] == "random") {

        // total number of cells

        unsigned int total_cells = nx_cells_m * ny_cells_m * nz_cells_m;

        // total number of particles

        unsigned int total_particles = total_cells * particle_per_cell;

        // Loop on all particles
        for (unsigned int ip = 0; ip < total_particles; ++ip) {

          // Random position
          double x = random.draw(inf_m[0], sup_m[0]);
          double y = random.draw(inf_m[1], sup_m[1]);
          double z = random.draw(inf_m[2], sup_m[2]);

          // Get the density
          const double w = param.density_profiles_[is](x / param.Lx, y / param.Ly, z / param.Lz);

          // only initialize particle if the weight is positive
          if (w > 1e-10) {

            particles_m[is].x_.h(ip) = x;
            particles_m[is].y_.h(ip) = y;
            particles_m[is].z_.h(ip) = z;

            particles_m[is].weight_.h(ip) = w * weight_coef;

            // increment the number of particles in this cell
            ++total_particles_counter;
          }
        }

        // Init at particle positions of species species_index_for_pos_init
      } else {

        unsigned int total_particles = particles_m[species_index_for_pos_init].size();

        for (unsigned int ip = 0; ip < total_particles; ++ip) {

          // Position
          particles_m[is].x_.h(ip) = particles_m[species_index_for_pos_init].x_.h(ip);
          particles_m[is].y_.h(ip) = particles_m[species_index_for_pos_init].y_.h(ip);
          particles_m[is].z_.h(ip) = particles_m[species_index_for_pos_init].z_.h(ip);

          // Get the density
          const double w =
            param.density_profiles_[is](static_cast<double>(particles_m[is].x_.h(ip)) / param.Lx,
                                        static_cast<double>(particles_m[is].y_.h(ip)) / param.Ly,
                                        static_cast<double>(particles_m[is].z_.h(ip)) / param.Lz);

          // weight
          particles_m[is].weight_.h(ip) = w * weight_coef;

          // increment the number of particles in this cell
          ++total_particles_counter;
        }

      } // end if param.position_initialization_method_

      // Random Momentum init
      for (auto ip = 0; ip < total_particles_counter; ++ip) {

        const double energy = Maxwell_Juttner_distribution(temperature / mass, random);

        // Sample angles randomly
        double phi   = std::acos(-random.draw(-1, 1));
        double theta = random.draw(0, 2. * M_PI);
        double psm   = std::sqrt(std::pow(1.0 + energy, 2) - 1.0);

        // Calculate the momentum
        double mx    = psm * cos(theta) * sin(phi);
        double my    = psm * sin(theta) * sin(phi);
        double mz    = psm * cos(phi);
        double gamma = std::sqrt(1.0 + mx * mx + my * my + mz * mz);

        // particles_m[is].mx_.h(ip) = mx;
        // particles_m[is].my_.h(ip) = my;
        // particles_m[is].mz_.h(ip) = mz;

        // Add the drift velocity using the Zenitani correction
        // See Zenitani et al. 2015

        if (v_drift > 0) {

          // Compute the gamma factor using momentum
          double inverse_gamma = 1. / gamma;

          const double check_velocity = (vx * mx + vy * my + vz * mz) * inverse_gamma;

          const double volume_acc = random.draw(0, 1);

          if (check_velocity > volume_acc) {

            const double Phi   = std::atan2(sqrt(vx * vx + vy * vy), vz);
            const double Theta = std::atan2(vy, vx);

            double vpx = mx * inverse_gamma;
            double vpy = my * inverse_gamma;
            double vpz = mz * inverse_gamma;

            const double vfl =
              vpx * cos(Theta) * sin(Phi) + vpy * sin(Theta) * sin(Phi) + vpz * cos(Phi);

            const double vflx = vfl * cos(Theta) * sin(Phi);
            const double vfly = vfl * sin(Theta) * sin(Phi);
            const double vflz = vfl * cos(Phi);

            vpx -= 2. * vflx;
            vpy -= 2. * vfly;
            vpz -= 2. * vflz;

            inverse_gamma = sqrt(1.0 - vpx * vpx - vpy * vpy - vpz * vpz);
            gamma         = 1. / inverse_gamma;

            mx = vpx * gamma;
            my = vpy * gamma;
            mz = vpz * gamma;

          } // here ends the corrections by Zenitani

          particles_m[is].mx_.h(ip) = -gamma * gamma_drift * vx + Lxx * mx + Lxy * my + Lxz * mz;
          particles_m[is].my_.h(ip) = -gamma * gamma_drift * vy + Lxy * my + Lyy * my + Lyz * mz;
          particles_m[is].mz_.h(ip) = -gamma * gamma_drift * vz + Lxz * mz + Lyz * my + Lzz * mz;
        } else {
          particles_m[is].mx_.h(ip) = mx;
          particles_m[is].my_.h(ip) = my;
          particles_m[is].mz_.h(ip) = mz;
        }

        // if (ip >= particles_m[is].size()) {
        //   std::cerr << ip << " " << particles_m[is].size() << std::endl;
        //   std::raise(SIGABRT);
        // }

      } // end for particles

    } else if (param.position_initialization_level_[is] == "cell") {

      // Loop over all cells

      for (int i = 0; i < nx_cells_m; i++) {
        for (int j = 0; j < ny_cells_m; j++) {
          for (int k = 0; k < nz_cells_m; k++) {

            const int i_global = i_patch_topology_m * (param.nx_cells_by_patch) + i;
            const int j_global = j_patch_topology_m * (param.ny_cells_by_patch) + j;
            const int k_global = k_patch_topology_m * (param.nz_cells_by_patch) + k;

            // Local cell index
            const int local_cell_index = i * (ny_cells_m * nz_cells_m) + j * nz_cells_m + k;

            // Global 1d cell index
            const int global_cell_index =
              i_global * (param.ny_cells * param.nz_cells) + j_global * param.nz_cells + k_global;

            // // cell min values
            // const double cell_x_min = i_global * param.dx;
            // const double cell_y_min = j_global * param.dy;
            // const double cell_z_min = k_global * param.dz;

            // // cell max values
            // const double cell_x_max = cell_x_min + param.dx;
            // const double cell_y_max = cell_y_min + param.dy;
            // const double cell_z_max = cell_z_min + param.dz;

            Random random(param.seed + global_cell_index + is);

            // counter to compute the number of particles in the current cell
            particles_per_cell_counter[is * total_cells + local_cell_index] = 0;

            // Random Position init
            if (param.position_initialization_method_[is] == "random") {

              for (auto p = 0; p < particle_per_cell; ++p) {

                const auto ip = total_particles_counter +
                                particles_per_cell_counter[is * total_cells + local_cell_index];

                const double x = (random.draw(0, 1) + i_global) * param.dx;
                const double y = (random.draw(0, 1) + j_global) * param.dy;
                const double z = (random.draw(0, 1) + k_global) * param.dz;

                // Get the density
                const double w =
                  param.density_profiles_[is](x / param.Lx, y / param.Ly, z / param.Lz);

                // only initialize particle if the weight is positive
                if (w > 1e-10) {

                  particles_m[is].x_.h(ip) = x;
                  particles_m[is].y_.h(ip) = y;
                  particles_m[is].z_.h(ip) = z;

                  particles_m[is].weight_.h(ip) = w * weight_coef;

                  // increment the number of particles in this cell
                  ++particles_per_cell_counter[is * total_cells + local_cell_index];
                }
              } // end for particles

              // Position initialization at same positions of species species_index_for_pos_init
            } else {

              for (auto p = 0;
                   p < particles_per_cell_counter[species_index_for_pos_init * total_cells +
                                                  local_cell_index];
                   ++p) {

                const auto ip = total_particles_counter + p;

                // Position
                particles_m[is].x_.h(ip) = particles_m[species_index_for_pos_init].x_.h(ip);
                particles_m[is].y_.h(ip) = particles_m[species_index_for_pos_init].y_.h(ip);
                particles_m[is].z_.h(ip) = particles_m[species_index_for_pos_init].z_.h(ip);

                // Get the density
                const double w = param.density_profiles_[is](
                  static_cast<double>(particles_m[is].x_.h(ip)) / param.Lx,
                  static_cast<double>(particles_m[is].y_.h(ip)) / param.Ly,
                  static_cast<double>(particles_m[is].z_.h(ip)) / param.Lz);

                // weight
                particles_m[is].weight_.h(ip) = w * weight_coef;

                // increment the number of particles in this cell
                ++particles_per_cell_counter[is * total_cells + local_cell_index];
              }
            }

            // Random Momentum init

            for (auto p = 0; p < particles_per_cell_counter[is * total_cells + local_cell_index];
                 ++p) {

              const auto ip = total_particles_counter + p;

              const double energy = Maxwell_Juttner_distribution(temperature / mass, random);

              // Sample angles randomly
              double phi   = std::acos(-random.draw(-1, 1));
              double theta = random.draw(0, 2. * M_PI);
              double psm   = std::sqrt(std::pow(1.0 + energy, 2) - 1.0);

              // Calculate the momentum
              double mx    = psm * cos(theta) * sin(phi);
              double my    = psm * sin(theta) * sin(phi);
              double mz    = psm * cos(phi);
              double gamma = std::sqrt(1.0 + mx * mx + my * my + mz * mz);

              // particles_m[is].mx_.h(ip) = mx;
              // particles_m[is].my_.h(ip) = my;
              // particles_m[is].mz_.h(ip) = mz;

              // Add the drift velocity using the Zenitani correction
              // See Zenitani et al. 2015

              if (v_drift > 0) {

                // Compute the gamma factor using momentum
                double inverse_gamma = 1. / gamma;

                const double check_velocity = (vx * mx + vy * my + vz * mz) * inverse_gamma;

                const double volume_acc = random.draw(0, 1);

                if (check_velocity > volume_acc) {

                  const double Phi   = std::atan2(sqrt(vx * vx + vy * vy), vz);
                  const double Theta = std::atan2(vy, vx);

                  double vpx = mx * inverse_gamma;
                  double vpy = my * inverse_gamma;
                  double vpz = mz * inverse_gamma;

                  const double vfl =
                    vpx * cos(Theta) * sin(Phi) + vpy * sin(Theta) * sin(Phi) + vpz * cos(Phi);

                  const double vflx = vfl * cos(Theta) * sin(Phi);
                  const double vfly = vfl * sin(Theta) * sin(Phi);
                  const double vflz = vfl * cos(Phi);

                  vpx -= 2. * vflx;
                  vpy -= 2. * vfly;
                  vpz -= 2. * vflz;

                  inverse_gamma = sqrt(1.0 - vpx * vpx - vpy * vpy - vpz * vpz);
                  gamma         = 1. / inverse_gamma;

                  mx = vpx * gamma;
                  my = vpy * gamma;
                  mz = vpz * gamma;

                } // here ends the corrections by Zenitani

                particles_m[is].mx_.h(ip) =
                  -gamma * gamma_drift * vx + Lxx * mx + Lxy * my + Lxz * mz;
                particles_m[is].my_.h(ip) =
                  -gamma * gamma_drift * vy + Lxy * my + Lyy * my + Lyz * mz;
                particles_m[is].mz_.h(ip) =
                  -gamma * gamma_drift * vz + Lxz * mz + Lyz * my + Lzz * mz;
              } else {
                particles_m[is].mx_.h(ip) = mx;
                particles_m[is].my_.h(ip) = my;
                particles_m[is].mz_.h(ip) = mz;
              }

              // if (ip >= particles_m[is].size()) {
              //   std::cerr << ip << " " << particles_m[is].size() << std::endl;
              //   std::raise(SIGABRT);
              // }

            } // end for particles

            // Add the new particles of this cell to the counter
            total_particles_counter +=
              particles_per_cell_counter[is * total_cells + local_cell_index];

            // } // if zone
          } // end for z cell
        }
      }

    } // pos init level

    // Add single particles
    for (int ip = 0; ip < param.particles_to_add_.size(); ++ip) {
      if (param.particles_to_add_[ip].is_ == is) {

        const double w = param.particles_to_add_[ip].w_;

        const double x = param.particles_to_add_[ip].x_;
        const double y = param.particles_to_add_[ip].y_;
        const double z = param.particles_to_add_[ip].z_;

        const double mx = param.particles_to_add_[ip].mx_;
        const double my = param.particles_to_add_[ip].my_;
        const double mz = param.particles_to_add_[ip].mz_;

        if ((x >= inf_m[0] && x < sup_m[0]) && (y >= inf_m[1] && y < sup_m[1]) &&
            (z >= inf_m[2] && z < sup_m[2])) {
          particles_m[is].set(total_particles_counter, w, x, y, z, mx, my, mz);
          total_particles_counter += 1;
        }
      }
    }

    // We resize the particles according to the real initialized number
    particles_m[is].resize(total_particles_counter, minipic::host);

    // particles_m[is].print();
    // particles_m[is].check_sum();

    // Copy data initialized on host to device (if exist)
    particles_m[is].sync(minipic::host, minipic::device);

    param.n_particles += total_particles_counter;

  } // end for species
}
