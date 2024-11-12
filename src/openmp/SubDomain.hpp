/* _____________________________________________________________________ */
//! \file SubDomain.hpp

//! \brief Management of the full domain

/* _____________________________________________________________________ */

#pragma once

#include <vector>

#include "Backend.hpp"
#include "Diagnotics.hpp"
#include "Operators.hpp"
#include "Params.hpp"
#include "Patch.hpp"
#include "Profiler.hpp"
#include "Timers.hpp"

//! \brief Wrapper class to clean main
class SubDomain {
public:
  //! List of Patch in this subdomain, linearized by Z, Z fast
  std::vector<Patch> patches_;

  // Init global fields
  ElectroMagn em_;

  // ______________________________________________________
  //
  //! \brief Alloc memory to store all the patch
  //! \param[in] params global parameters
  // ______________________________________________________
  void allocate(Params &params, Backend &backend) {

    std::cout << params.seperator(50) << std::endl << std::endl;
    std::cout << " Initialization" << std::endl << std::endl;

    // Allocate global fields
    em_.allocate(params, backend);

    const double memory_consumption =
      (em_.Ex_m.size() + em_.Ey_m.size() + em_.Ez_m.size() + em_.Bx_m.size() + em_.By_m.size() +
       em_.Bz_m.size() +
       (em_.Jx_m.size() + em_.Jy_m.size() + em_.Jz_m.size()) * (params.species_names_.size() + 1)) *
      8. / (1024. * 1024);

    std::cout << " Field grids: " << memory_consumption << " Mb" << std::endl << std::endl;

    // Allocate and initialize particles for each patch on host
    patches_.resize(params.N_patches);
    for (int i = 0; i < params.nx_patch; i++) {
      for (int j = 0; j < params.ny_patch; j++) {
        for (int k = 0; k < params.nz_patch; k++) {
          int idx = i * params.nz_patch * params.ny_patch + j * params.nz_patch + k;
          // Memory allocate for all particles and local fields
          patches_[idx].allocate(params, backend, i, j, k);

          // Particles position and momentum initialization
          patches_[idx].initialize_particles(params);
        }
      }
    }

    // Momentum correction (to respect the leap frog scheme)
    if (params.momentum_correction) {

      std::cout << " > Apply momentum correction "
                << "\n"
                << std::endl;

      for (int ip = 0; ip < patches_.size(); ip++) {
        operators::interpolate(em_, patches_[ip]);
        operators::push_momentum(patches_[ip], -0.5 * params.dt);
      }
    }

    // For each species, print :
    // - total number of particles
    for (auto is = 0; is < params.species_names_.size(); ++is) {
      unsigned int total_number_of_particles = 0;
      mini_float total_particle_energy       = 0;
      for (auto idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {
        total_number_of_particles += patches_[idx_patch].particles_m[is].size();
        total_particle_energy +=
          patches_[idx_patch].particles_m[is].get_kinetic_energy(minipic::host);
      }
      std::cout << " Species " << params.species_names_[is] << std::endl;

      const double memory_consumption = total_number_of_particles * 14. * 8. / (1024. * 1024);

      std::cout << " - Initialized particles: " << total_number_of_particles << std::endl;
      std::cout << " - Total kinetic energy: " << total_particle_energy << std::endl;
      std::cout << " - Memory footprint: " << memory_consumption << " Mb" << std::endl;
    }
    
    // Checksum for field

    auto sum_Ex_on_host = em_.Ex_m.sum(1, minipic::host);
    auto sum_Ey_on_host = em_.Ey_m.sum(1, minipic::host);
    auto sum_Ez_on_host = em_.Ez_m.sum(1, minipic::host);

    auto sum_Bx_on_host = em_.Bx_m.sum(1, minipic::host);
    auto sum_By_on_host = em_.By_m.sum(1, minipic::host);
    auto sum_Bz_on_host = em_.Bz_m.sum(1, minipic::host);

    auto sum_Ex_on_device = em_.Ex_m.sum(1, minipic::device);
    auto sum_Ey_on_device = em_.Ey_m.sum(1, minipic::device);
    auto sum_Ez_on_device = em_.Ez_m.sum(1, minipic::device);

    auto sum_Bx_on_device = em_.Bx_m.sum(1, minipic::device);
    auto sum_By_on_device = em_.By_m.sum(1, minipic::device);
    auto sum_Bz_on_device = em_.Bz_m.sum(1, minipic::device);

    static const int p = 3;

    std::cout << std::endl;
    std::cout << " -------------------------------- |" << std::endl;
    std::cout << " Check sum for fields             |" << std::endl;
    std::cout << " -------------------------------- |" << std::endl;
    std::cout << " Field  | Host       | Device     |" << std::endl;
    std::cout << " -------------------------------- |" << std::endl;
    std::cout << " Ex     | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Ex_on_host << " | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Ex_on_device << " | " << std::endl;
    std::cout << " Ey     | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Ey_on_host << " | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Ey_on_device << " | " << std::endl;
    std::cout << " Ez     | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Ez_on_host << " | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Ez_on_device << " | " << std::endl;
    std::cout << " Bx     | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Bx_on_host << " | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Bx_on_device << " | " << std::endl;
    std::cout << " By     | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_By_on_host << " | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_By_on_device << " | " << std::endl;
    std::cout << " Bz     | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Bz_on_host << " | " << std::setw(10) << std::scientific << std::setprecision(p)
              << sum_Bz_on_device << " | " << std::endl;

    // Checksum for particles

    double sum_device[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double sum_host[13]   = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    static const std::string vector_name[13] =
      {"weight", "x", "y", "z", "mx", "my", "mz", "Ex", "Ey", "Ez", "Bx", "By", "Bz"};
    for (int ip = 0; ip < patches_.size(); ip++) {
      for (auto is = 0; is < params.species_names_.size(); ++is) {

        sum_host[0] += patches_[ip].particles_m[is].weight_.sum(1, minipic::host);
        sum_host[1] += patches_[ip].particles_m[is].x_.sum(1, minipic::host);
        sum_host[2] += patches_[ip].particles_m[is].y_.sum(1, minipic::host);
        sum_host[3] += patches_[ip].particles_m[is].z_.sum(1, minipic::host);
        sum_host[4] += patches_[ip].particles_m[is].mx_.sum(1, minipic::host);
        sum_host[5] += patches_[ip].particles_m[is].my_.sum(1, minipic::host);
        sum_host[6] += patches_[ip].particles_m[is].mz_.sum(1, minipic::host);
        sum_host[7] += patches_[ip].particles_m[is].Ex_.sum(1, minipic::host);
        sum_host[8] += patches_[ip].particles_m[is].Ey_.sum(1, minipic::host);
        sum_host[9] += patches_[ip].particles_m[is].Ez_.sum(1, minipic::host);
        sum_host[10] += patches_[ip].particles_m[is].Bx_.sum(1, minipic::host);
        sum_host[11] += patches_[ip].particles_m[is].By_.sum(1, minipic::host);
        sum_host[12] += patches_[ip].particles_m[is].Bz_.sum(1, minipic::host);

        sum_device[0] += patches_[ip].particles_m[is].weight_.sum(1, minipic::device);
        sum_device[1] += patches_[ip].particles_m[is].x_.sum(1, minipic::device);
        sum_device[2] += patches_[ip].particles_m[is].y_.sum(1, minipic::device);
        sum_device[3] += patches_[ip].particles_m[is].z_.sum(1, minipic::device);
        sum_device[4] += patches_[ip].particles_m[is].mx_.sum(1, minipic::device);
        sum_device[5] += patches_[ip].particles_m[is].my_.sum(1, minipic::device);
        sum_device[6] += patches_[ip].particles_m[is].mz_.sum(1, minipic::device);
        sum_device[7] += patches_[ip].particles_m[is].Ex_.sum(1, minipic::device);
        sum_device[8] += patches_[ip].particles_m[is].Ey_.sum(1, minipic::device);
        sum_device[9] += patches_[ip].particles_m[is].Ez_.sum(1, minipic::device);
        sum_device[10] += patches_[ip].particles_m[is].Bx_.sum(1, minipic::device);
        sum_device[11] += patches_[ip].particles_m[is].By_.sum(1, minipic::device);
        sum_device[12] += patches_[ip].particles_m[is].Bz_.sum(1, minipic::device);
      }
    }

    std::cout << std::endl;
    std::cout << " -------------------------------- |" << std::endl;
    std::cout << " Check sum for particles          |" << std::endl;
    std::cout << " -------------------------------- |" << std::endl;
    std::cout << " vector | Host       | Device     |" << std::endl;
    std::cout << " -------------------------------- |" << std::endl;

    for (int i = 0; i < 13; i++) {
      std::cout << " " << std::setw(6) << vector_name[i] << " | " << std::setw(10)
                << std::scientific << std::setprecision(p) << sum_host[i] << " | " << std::setw(10)
                << std::scientific << std::setprecision(p) << sum_device[i] << " | " << std::endl;
    }

  }


  // ______________________________________________________________________________
  //
  //! \brief Perform a single PIC iteration
  //! \param[in] Params&  global parameters
  //! \param[in] Timers&  timers
  //! \param[in] Profiler& profiler for thread profiling
  //! \param[in] int it iteration number
  // ______________________________________________________________________________
  void iterate(Params &params, Timers &timers, Profiler &profiler, Backend &backend, int it) {

    if (params.current_projection || params.n_particles > 0) {

      DEBUG("start reset current");

#ifdef __MINIPIC_OMP__
#pragma omp single // single nowait
#endif
      {
        timers.start(timers.reset_current);
        profiler.start(RESET);

        em_.reset_currents(minipic::host);

        timers.stop(timers.reset_current);
        profiler.stop();
      }
      DEBUG("stop reset current");
    }

#ifdef __MINIPIC_OMP__
#pragma omp for schedule(runtime)
#endif

    for (int idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {

      // Interpolate from global field to particles in patch
      timers.start(timers.interpolate, idx_patch);
      profiler.start(EVOLVE_PATCH);

      DEBUG("start interpolate");
      operators::interpolate(em_, patches_[idx_patch]);

      DEBUG("stop interpolate");
      timers.stop(timers.interpolate, idx_patch);
      //profiler.stop();

      // Push all particles in patch
      timers.start(timers.push, idx_patch);
      //profiler.start(PUSH);

      operators::push(patches_[idx_patch], params.dt);

      timers.stop(timers.push, idx_patch);
      //profiler.stop();

      // Do boundary conditions on global domain
      timers.start(timers.pushBC, idx_patch);
      //profiler.start(PUSHBC);

      operators::pushBC(params, patches_[idx_patch]);

      timers.stop(timers.pushBC, idx_patch);
      //profiler.stop();

      // Projection in local field
      if (params.current_projection) {

        timers.start(timers.projection, idx_patch);
        //profiler.start(PROJECT);

        DEBUG("start project");
        // Project in buffers local to the patches
        operators::project(params, patches_[idx_patch]);

        DEBUG("stop project");

        // #endif

        timers.stop(timers.projection, idx_patch);
        //profiler.stop();
      }

      // __________________________________________________________________
      // Identify and copy in buffers particles which leave the patch

      timers.start(timers.id_parts_to_move, idx_patch);
      //profiler.start(EXCHANGE);

      DEBUG("Patch " << idx_patch << ": start identify particles to move");

      if (patches_.size() > 1) {
        operators::identify_particles_to_move(params, patches_[idx_patch], backend);
      }

      DEBUG("Patch " << idx_patch << ": stop identify particles to move");

      timers.stop(timers.id_parts_to_move, idx_patch);
      profiler.stop();
    } // end for patches

    // __________________________________________________________________
    // Exchange particles between patches

    if (patches_.size() > 1) {

#ifdef __MINIPIC_OMP__
#pragma omp for schedule(runtime)
#endif

      for (int idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {
        profiler.start(EXCHANGE);
        timers.start(timers.exchange, idx_patch);

        DEBUG("Patch " << idx_patch << ": exchange");

        operators::exchange_particles(params, patches_, idx_patch);

        DEBUG("Patch " << idx_patch << ": end exchange");

        timers.stop(timers.exchange, idx_patch);
        profiler.stop();
      }
    }

    // __________________________________________________________________
    // Sum all species contribution in the local and global current grids

    if (params.current_projection || params.n_particles > 0) {

#ifdef __MINIPIC_OMP__
#pragma omp for schedule(runtime)
#endif
      for (int idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {

        timers.start(timers.current_local_reduc, idx_patch);
        profiler.start(CURRENT_GLOBAL_BORDERS);

        // Projection directly in the global grid
        // subdomain.patches_[idx_patch].project(param, em_);

        // Projection in local field
        // patches_[idx_patch].project(params);

        // Sum all species contribution in the local fields
        operators::reduc_current(patches_[idx_patch]);

        timers.stop(timers.current_local_reduc, idx_patch);
        profiler.stop();
      }

#ifdef __MINIPIC_OMP__
#pragma omp single // single nowait
#endif
      {

        timers.start(timers.current_global_reduc);

        profiler.start(CURRENT_GLOBAL_BORDERS);
        for (int idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {

          DEBUG("start local 2 global")
          // Copy all local fields in the global fields
          operators::local2global(em_, patches_[idx_patch]);
          DEBUG("end local 2 global")
        }
        profiler.stop();

        timers.stop(timers.current_global_reduc);

        timers.start(timers.currentBC);
        profiler.start(CURRENTBC);

        // Perform the boundary conditions for current
        DEBUG("start current BC")
        operators::currentBC(params, em_);
        DEBUG("end current BC")

        profiler.stop();
        timers.stop(timers.currentBC);
      }

    } // end if current projection

    // __________________________________________________________________
    // Imbalance operator in particles

    if (!params.imbalance_function_.empty()) {

#ifdef __MINIPIC_OMP__
#pragma omp for schedule(runtime)
#endif
      for (int idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {

        timers.start(timers.imbalance, idx_patch);
        operators::imbalance_operator(params,
                                      patches_[idx_patch],
                                      it,
                                      params.imbalance_function_[0]);
        timers.stop(timers.imbalance, idx_patch);

      } // end iterate
    }   // end imbalance

    // __________________________________________________________________
    // Maxwell solver

    if (params.maxwell_solver) {

#ifdef __MINIPIC_OMP__
#pragma omp single // single nowait
#endif
      { timers.start(timers.maxwell_solver); }

#ifdef __MINIPIC_OMP__
#pragma omp single // single nowait
#endif
      {

        // profiler.start(MAXWELL);
        // Generate a laser field with an antenna
        for (auto iantenna = 0; iantenna < params.antenna_profiles_.size(); iantenna++) {
          operators::antenna(params,
                             em_,
                             params.antenna_profiles_[iantenna],
                             params.antenna_positions_[iantenna],
                             it * params.dt);
        }
        // profiler.stop();
      }

      // Solve the Maxwell equation
      operators::solve_maxwell(params, em_, profiler);

#ifdef __MINIPIC_OMP__
#pragma omp single // single nowait
#endif
      { timers.stop(timers.maxwell_solver); }

      // __________________________________________________________________
      // Maxwell Boundary conditions

#ifdef __MINIPIC_OMP__
#pragma omp single // single nowait
#endif
      {
        profiler.start(MAXWELL);
        timers.start(timers.maxwellBC);

        DEBUG("start solve BC")

        // Boundary conditions on EM fields
        operators::solveBC(params, em_);

        DEBUG("end solve BC")
        timers.stop(timers.maxwellBC);
        profiler.stop();
      }

    } // end test params.maxwell_solver
  }   // end iterate

  // ________________________________________________________________
  //! \brief Perform all diagnostics
  //! \param[in] Params&  global parameters
  //! \param[in] Timers&  timers
  //! \param[in] Profiler& profiler for detailed time measurement
  //! \param[in] int it iteration number
  // ________________________________________________________________
  void diagnostics(Params &params, Timers &timers, Profiler &profiler, Backend &backend, int it) {

    if (params.no_diagnostics_at_init and it == 0) {
      return;
    }

#ifdef __MINIPIC_OMP__
#pragma omp single
#endif
    { timers.start(timers.diags_binning); }

    // Particle binning
#ifdef __MINIPIC_OMP__
#pragma omp for schedule(runtime) // collapse(2)
#endif
    for (auto particle_binning : params.particle_binning_properties_) {

      // for each species index of this diagnostic
      for (auto is : particle_binning.species_indexes_) {

        if (!(it % particle_binning.period_)) {

          profiler.start(DIAGS);

          // Call the particle binning function using the properties in particle_binning
          Diags::particle_binning(particle_binning.name_,
                                  params,
                                  patches_,
                                  particle_binning.projected_parameter_,
                                  particle_binning.axis_,
                                  particle_binning.n_cells_,
                                  particle_binning.min_,
                                  particle_binning.max_,
                                  is,
                                  it,
                                  particle_binning.format_,
                                  false);

          profiler.stop();
        } // end if test it % period
      }
    } // end loop on particle_binning_properties_

#ifdef __MINIPIC_OMP__
#pragma omp single // single nowait
#endif
    timers.stop(timers.diags_binning);

    // Particle Clouds
    if ((params.particle_cloud_period < params.n_it) &&
        (!(it % params.particle_cloud_period) or (it == 0))) {

#ifdef __MINIPIC_OMP__
#pragma omp single // single nowait
#endif
      { timers.start(timers.diags_cloud); }

#ifdef __MINIPIC_OMP__
#pragma omp for schedule(runtime)
#endif
      for (auto is = 0; is < params.get_species_number(); ++is) {

        profiler.start(DIAGS);

        Diags::particle_cloud("cloud", params, patches_, is, it, params.particle_cloud_format);

        profiler.stop();
      }

#ifdef __MINIPIC_OMP__
#pragma omp single // single nowait
#endif
      timers.stop(timers.diags_cloud);
    }

    // Field diagnostics
    if (!(it % params.field_diagnostics_period)) {
#ifdef __MINIPIC_OMP__
#pragma omp single
#endif
      {

        timers.start(timers.diags_field);

        profiler.start(DIAGS);

        Diags::fields(params, em_, it, params.field_diagnostics_format);

        profiler.stop();

        timers.stop(timers.diags_field);
      }
    }

    // Scalars diagnostics
    if (!(it % params.scalar_diagnostics_period)) {
#ifdef __MINIPIC_OMP__
#pragma omp single
#endif
      {
        timers.start(timers.diags_scalar);
        for (auto is = 0; is < params.get_species_number(); ++is) {
          profiler.start(DIAGS);

          Diags::scalars(params, patches_, is, it);

          profiler.stop();
        }
        timers.stop(timers.diags_scalar);
      }
    }

    if (!(it % params.scalar_diagnostics_period)) {
#ifdef __MINIPIC_OMP__
#pragma omp single
#endif
      {

        timers.start(timers.diags_scalar);

        profiler.start(DIAGS);

        Diags::scalars(params, em_, it);

        profiler.stop();

        timers.stop(timers.diags_scalar);
      }
    }

  } // end diagnostics

  // __________________________________________________________________
  //
  //! \brief get the total number of particles
  // __________________________________________________________________
  unsigned int get_total_number_of_particles() {
    unsigned int total_number_of_particles = 0;
    for (auto idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {
      total_number_of_particles += patches_[idx_patch].get_total_number_of_particles();
    }
    return total_number_of_particles;
  }

}; // end class
