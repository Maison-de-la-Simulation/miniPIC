/* _____________________________________________________________________ */
//! \file SubDomain.hpp

//! \brief Management of the full domain

/* _____________________________________________________________________ */

#pragma once

#include <atomic>
#include <deque>
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

    int *evolve_particles_flags = backend.evolve_particles_flags;
    int *reset_current_flags = backend.reset_current_flags;
    int *maxwell_solver_flags = backend.maxwell_solver_flags;
    int *reduction_internal_flags = backend.reduction_internal_flags;
    int *reduction_external_flags = backend.reduction_external_flags;               
    int *evolve_patch_flags = backend.evolve_patch_flags;

    std::deque<std::atomic<int>> &task_exchange_count_ = backend.task_exchange_count_;

#pragma omp taskgroup
    {

#pragma omp task untied default(none) firstprivate(params) shared(em_, timers, profiler, minipic::host) \
  depend(out : reset_current_flags)
      {

        if (params.current_projection || params.n_particles > 0) {

          DEBUG("start reset current");

          timers.start(timers.reset_current);
          profiler.start(RESET);

          em_.reset_currents(minipic::host);

          timers.stop(timers.reset_current);
          profiler.stop();
          
          DEBUG("stop reset current");

        } // end if
      
      } 

      //  ______________________________________________________
      //  Create the bins and execution functions by bin

      for (int idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {

        ///////////////////////////
        int n_species = patches_[idx_patch].particles_m.size();

        for (int is = 0; is < n_species; is++) {

          int n_particles = patches_[idx_patch].particles_m[is].size();

          if (n_particles > 0) {

            int bin_size   = params.bin_size;
            int bin_number = 1 + ((n_particles - 1) / bin_size);

            for (unsigned int i_bin = 0; i_bin < bin_number; i_bin++) {

              int rest = n_particles - i_bin * bin_size;
              int size = std::min(bin_size, rest);

              int init = i_bin * bin_size;        // start of bin
              int end  = i_bin * bin_size + size; // end of bin

#pragma omp task untied default(none) firstprivate(idx_patch, is, init, end, rest) \
  shared(timers, profiler, patches_, em_, params)                                  \
  depend(out : evolve_particles_flags[idx_patch]) \
  priority(bin_number)
              {
                profiler.start(EVOLVE_BIN);

                if (rest > 0) {

                  ////

                  // functions using bin
                  timers.start(timers.interpolate, idx_patch);

                  operators::interpolate_bin(em_,
                                             patches_[idx_patch].particles_m[is],
                                             is,
                                             init,
                                             end);

                  timers.stop(timers.interpolate, idx_patch);
                  profiler.stop();

                  timers.start(timers.push, idx_patch);

                  operators::push_bin(params.dt,
                                      patches_[idx_patch].particles_m[is],
                                      is,
                                      init,
                                      end);

                  timers.stop(timers.push, idx_patch);

                  timers.start(timers.pushBC, idx_patch);

                  operators::pushBC_bin(params,
                                        patches_[idx_patch],
                                        patches_[idx_patch].particles_m[is],
                                        patches_[idx_patch].on_border_m,
                                        is,
                                        init,
                                        end);

                  timers.stop(timers.pushBC, idx_patch);
                  
                } // end rest particles

                profiler.stop();
              } // end task
            }   // end bin loop
          }     // end if n particles
        }       // end loop species
      }         // end of the loop

      //  ______________________________________________________
      // Continue witn function can only by solved by patch

      for (int idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {

#pragma omp task untied default(none) firstprivate(idx_patch,it) shared(patches_,                \
                                                                       timers,                   \
                                                                       profiler,                 \
                                                                       params,                   \
                                                                       task_exchange_count_,     \
                                                                       reset_current_flags,      \
                                                                       evolve_particles_flags,   \
                                                                       reduction_internal_flags, \
                                                                       reduction_external_flags, \
                                                                       backend                   \
                                                                    ) \
  depend(in : evolve_particles_flags[idx_patch]) depend(in : reset_current_flags)
        {

          profiler.start(EVOLVE_PATCH);

          ////
          
          if (params.current_projection) {

            // __________________________________________________________________
            // Projection in local field

            timers.start(timers.projection, idx_patch);
            DEBUG("start project");
  
            // Project in buffers local to the patches
            operators::project(params, patches_[idx_patch]);

            DEBUG("stop project");
            timers.stop(timers.projection, idx_patch);          
          }

          // __________________________________________________________________
          // Identify and copy in buffers particles which leave the patch

          timers.start(timers.id_parts_to_move, idx_patch);
          DEBUG("Patch " << idx_patch << ": start identify particles to move");

          if (patches_.size() > 1) {
            operators::identify_particles_to_move(params, patches_[idx_patch], backend);
          }

          DEBUG("Patch " << idx_patch << ": stop identify particles to move");
          timers.stop(timers.id_parts_to_move, idx_patch);

          ////
          
          profiler.stop();

          if (params.current_projection || params.n_particles > 0) {

            // __________________________________________________________________
            // Projection in the local grid

            timers.start(timers.current_local_reduc, idx_patch);

            // Sum all species contribution in the local fields
            operators::reduc_current(patches_[idx_patch]);

            timers.stop(timers.current_local_reduc, idx_patch);

            #pragma omp task untied default(none) firstprivate(idx_patch) \
            shared(timers, patches_, profiler, em_, params) \
            //depend(out : reduction_internal_flags[idx_patch])
            {
              // __________________________________________________________________
              // Projection to global grid internal
              profiler.start(CURRENT_GLOBAL_INTERNAL);
              timers.start(timers.current_global_reduc, idx_patch);
              // Copy all local fields in the global fields for inner cell
              operators::local2global_internal(em_, patches_[idx_patch]);
              profiler.stop();
              timers.stop(timers.current_global_reduc, idx_patch);
            } // end task

            #pragma omp task untied default(none) firstprivate(idx_patch) \
            shared(timers, patches_, profiler, em_, params) \
            //depend(out : reduction_external_flags[idx_patch])
            {
              // __________________________________________________________________
              // Projection to global grid internal
              timers.start(timers.current_global_reduc, idx_patch);
              profiler.start(CURRENT_GLOBAL_BORDERS);
              // Copy all local fields in the global fields
              operators::local2global_borders(em_, patches_[idx_patch]);
              profiler.stop();
              timers.stop(timers.current_global_reduc, idx_patch);
            } // end task
          } // end if

          // __________________________________________________________________
          // Exchange particles between patches

          if (patches_.size() > 1) {

            // Task scheduler
            int i_patch = patches_[idx_patch].i_patch_topology_m;
            int j_patch = patches_[idx_patch].j_patch_topology_m;
            int k_patch = patches_[idx_patch].k_patch_topology_m;

            // Collect new particles from neighbours
            for (int i = -1; i < 2; i++) {
              for (int j = -1; j < 2; j++) {
                for (int k = -1; k < 2; k++) {

                  // id of the Neighbor in vec_patch
                  int idx_neighbor =
                    patches_[idx_patch].get_idx_patch(i_patch + i, j_patch + j, k_patch + k);

                  // Task count
                  if (--task_exchange_count_[idx_neighbor] == 0) {

#pragma omp task untied default(none) firstprivate(idx_neighbor, it) \
  shared(patches_, em_, timers, profiler, params)
                    {

                      timers.start(timers.exchange, idx_neighbor);
                      profiler.start(EXCHANGE)
                      DEBUG("Patch " << idx_neighbor << ": exchange");

                      operators::exchange_particles(params, patches_, idx_neighbor);
                      DEBUG("Patch " << idx_neighbor << ": exchange");
                      profiler.stop();                      
                      timers.stop(timers.exchange, idx_neighbor);

                    } // end exchange task

                    task_exchange_count_[idx_neighbor] = 27; // Reset count

                  } // end if conditional patch ready
                }   // end k cycle
              }     // end j cycle
            }       // end i cycle
          }         // end if task patch > 1
        }           // end task
      }             // end for loop
    }               // end of taskgroup

    // __________________________________________________________________
    // Imbalance operator

    if (!params.imbalance_function_.empty()) {

#pragma omp taskgroup
      {
        for (int idx_patch = 0; idx_patch < patches_.size(); idx_patch++) {

          ///////////////////////////
          int n_species = patches_[idx_patch].particles_m.size();

          for (int is = 0; is < n_species; is++) {

            int n_particles = patches_[idx_patch].particles_m[is].size();

            if (n_particles > 0) {

              int bin_size   = params.bin_size;
              int bin_number = 1 + ((n_particles - 1) / bin_size);

              for (unsigned int i_bin = 0; i_bin < bin_number; i_bin++) {

                int rest = n_particles - i_bin * bin_size;
                int size = std::min(bin_size, rest);

                int init = i_bin * bin_size;        // start of bin
                int end  = i_bin * bin_size + size; // end of bin

#pragma omp task untied default(none) firstprivate(idx_patch, is, init, end, rest) \
  shared(timers, profiler, patches_, em_, it, params) if (rest > 0)
                {
                  if (rest > 0) {
                    timers.start(timers.imbalance, idx_patch);
                    operators::imbalance_operator(params,
                                                  patches_[idx_patch].particles_m[is],
                                                  init,
                                                  end,
                                                  it,
                                                  params.imbalance_function_[0]);
                    timers.stop(timers.imbalance, idx_patch);
                  } // end rest particles

                } // end task
              }   // end bin loop
            }     // end if n particles
          }       // end loop species
        }         // end loop patches
      }           // end taskgroup
    }             // end imbalance if

    if (params.maxwell_solver) {

#pragma omp task untied default(none) firstprivate(it) shared(timers, profiler, em_, params) \
  depend(out : maxwell_solver_flags)
      {

        if (params.current_projection || params.n_particles > 0) {

          timers.start(timers.currentBC);
          profiler.start(CURRENTBC);

          // Perform the boundary conditions for current
          DEBUG("start current BC")
          operators::currentBC(params, em_);
          DEBUG("end current BC")

          profiler.stop();
          timers.stop(timers.currentBC);

        } // end if

        timers.start(timers.maxwell_solver);

        // Generate a laser field with an antenna
        for (auto iantenna = 0; iantenna < params.antenna_profiles_.size(); iantenna++) {
          operators::antenna(params,
                             em_,
                             params.antenna_profiles_[iantenna],
                             params.antenna_positions_[iantenna],
                             it * params.dt);
        }

        // Solve the Maxwell equation
        operators::solve_maxwell(params, em_, profiler);

        timers.stop(timers.maxwell_solver);

        // __________________________________________________________________
        // Maxwell Boundary conditions

        profiler.start(MAXWELL);
        timers.start(timers.maxwellBC);

        DEBUG("start solve BC")

        // Boundary conditions on EM fields
        operators::solveBC(params, em_);

        DEBUG("end solve BC")
        timers.stop(timers.maxwellBC);
        profiler.stop();

      } // end task
    }   // end if maxwell
  }     // end iterate

  // ________________________________________________________________
  //! \brief Perform all diagnostics
  //! \param[in] Params&  global parameters
  //! \param[in] Timers&  timers
  //! \param[in] Profiler& profiler for detailed time measurement
  //! \param[in] int it iteration number
  // ________________________________________________________________
  void diagnostics(Params &params, Timers &timers, Profiler &profiler, Backend &backend, int it) {

    int *maxwell_solver_flags = backend.maxwell_solver_flags;

    int ns = params.get_species_number();

#pragma omp taskgroup
    {

      // Particle binning
      for (auto particle_binning : params.particle_binning_properties_) {

        // for each species index of this diagnostic
        for (auto is : particle_binning.species_indexes_) {

#pragma omp task untied default(none) firstprivate(particle_binning, is, it, ns) \
  shared(timers, profiler, params)
          {

            if (!(it % particle_binning.period_)) {

              profiler.start(DIAGS);
              timers.start(timers.diags_binning);

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
              timers.stop(timers.diags_binning);
              profiler.stop();
            } // end if test it % period
          }   // end task

        } // end inner loop
      }   // end loop on particle_binning_properties_

      // Particle Clouds
      if ((params.particle_cloud_period < params.n_it) &&
          (!(it % params.particle_cloud_period) or (it == 0))) {

        timers.start(timers.diags_cloud);

        for (auto is = 0; is < params.get_species_number(); ++is) {

#pragma omp task untied default(none) firstprivate(is, it) shared(timers, profiler, params)
          {
            profiler.start(DIAGS);
            Diags::particle_cloud("cloud", params, patches_, is, it, params.particle_cloud_format);
            profiler.stop();
          } // End task
        }   // end for

        timers.stop(timers.diags_cloud);
      }

      // Scalars diagnostics
      if (!(it % params.scalar_diagnostics_period)) {

        timers.start(timers.diags_scalar);

        for (auto is = 0; is < params.get_species_number(); ++is) {
#pragma omp task untied default(none) firstprivate(it, is) shared(params, timers, profiler)
          {
            profiler.start(DIAGS);
            Diags::scalars(params, patches_, is, it);
            profiler.stop();
          } // End task
        }

        timers.stop(timers.diags_scalar);
      }

    } // End taskgroup

#if (ns < 1)
#pragma omp taskwait
#endif

#pragma omp taskgroup
    {
      // Scalar diagnostics
      if (!(it % params.scalar_diagnostics_period)) {

#pragma omp task untied default(none) firstprivate(it) shared(params, timers, profiler) \
  depend(in : maxwell_solver_flags)
        {
          timers.start(timers.diags_scalar);
          profiler.start(DIAGS);

          Diags::scalars(params, em_, it);

          profiler.stop();
          timers.stop(timers.diags_scalar);
        } // End task
      }   // If end

      // Field diagnostics
      if (!(it % params.field_diagnostics_period)) {

#pragma omp task untied default(none) firstprivate(it) shared(params, timers, profiler) \
  depend(in : maxwell_solver_flags)
        {
          timers.start(timers.diags_field);
          profiler.start(DIAGS);

          Diags::fields(params, em_, it, params.field_diagnostics_format);

          profiler.stop();
          timers.stop(timers.diags_field);
        } // end task
      }   // If end
    }     // End taskgroup

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
