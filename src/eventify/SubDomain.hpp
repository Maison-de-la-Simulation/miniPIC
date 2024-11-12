/* _____________________________________________________________________ */
//! \file SubDomain.hpp

//! \brief Management of the full domain

/* _____________________________________________________________________ */

#pragma once

#include <atomic>
#include <cassert>
#include <deque>
#include <latch>
#include <ranges>
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

  // ________________________________________________________________
  //! \brief Perform all diagnostics
  //! \param[in] Params&  global parameters
  //! \param[in] Timers&  timers
  //! \param[in] Profiler& profiler for detailed time measurement
  //! \param[in] int it iteration number
  // ________________________________________________________________
  void diagnostics(Params &params, Timers &timers, Profiler &profiler, Backend &backend, int it) {

    // Particle binning
    for (auto particle_binning : params.particle_binning_properties_) {
      // for each species index of this diagnostic
      for (auto is : particle_binning.species_indexes_) {
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
      }
    } // end loop on particle_binning_properties_

    // Particle Clouds
    if ((params.particle_cloud_period < params.n_it) &&
        (!(it % params.particle_cloud_period) or (it == 0))) {

      timers.start(timers.diags_cloud);
      for (auto is = 0; is < params.get_species_number(); ++is) {
        profiler.start(DIAGS);

        Diags::particle_cloud("cloud", params, patches_, is, it, params.particle_cloud_format);

        profiler.stop();
      }
      timers.stop(timers.diags_cloud);
    }

    // Field diagnostics
    if (!(it % params.field_diagnostics_period)) {
      timers.start(timers.diags_field);
      profiler.start(DIAGS);

      Diags::fields(params, em_, it, params.field_diagnostics_format);

      profiler.stop();
      timers.stop(timers.diags_field);
    }

    // Scalars diagnostics
    if (!(it % params.scalar_diagnostics_period)) {
      timers.start(timers.diags_scalar);
      for (auto is = 0; is < params.get_species_number(); ++is) {
        profiler.start(DIAGS);

        Diags::scalars(params, patches_, is, it);

        profiler.stop();
      }
      timers.stop(timers.diags_scalar);
    }

    if (!(it % params.scalar_diagnostics_period)) {
      timers.start(timers.diags_scalar);
      profiler.start(DIAGS);

      Diags::scalars(params, em_, it);

      profiler.stop();
      timers.stop(timers.diags_scalar);
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

namespace mini_pic {
using task_system = eventify::task_system;

// note: The task graph is a combination of
//       a communication channel and a synchronization channel,
//
//       i.e. it handles both the flow of data between tasks
//       and the dependencies between them.
//
//////////////////////////////////////////////////////////////////////////////////////////

struct task_graph {
  std::deque<jsc::event_counter> current_reduction{};
  std::deque<jsc::event_counter> exchange{};
  std::deque<jsc::event_counter> evolve{};
  std::deque<jsc::event_counter> currentBC_antenna_imbalance{};

  jsc::event_counter currentBC_antenna;
  jsc::event_counter imbalance_operator;
  jsc::event_counter particles_print;
  jsc::event_counter field_diags{1};

  jsc::event_counter maxwell_ampere_scheduler{1};
  jsc::event_counter maxwell_solverBC;
  jsc::event_counter maxwell_faraday_scheduler;

  std::latch completion{6};

  explicit task_graph(const std::ptrdiff_t patch_count,
                      const std::ptrdiff_t bin_size,
                      const std::vector<Patch> &patches_,
                      const std::ptrdiff_t n_species,
                      const std::ptrdiff_t n_task_solverBC,
                      const std::ptrdiff_t n_task_solver_faraday)
    : currentBC_antenna{2*patch_count}, imbalance_operator{2*patch_count}, particles_print{patch_count},
      maxwell_solverBC{n_task_solverBC}, maxwell_faraday_scheduler{n_task_solver_faraday} {
    assert(patch_count > 0);

    for (auto _ : std::views::iota(std::ptrdiff_t{0}, patch_count)) {
      current_reduction.emplace_back(2);
    }

    for (auto _ : std::views::iota(std::ptrdiff_t{0}, patch_count)) {
      exchange.emplace_back(27);
    }

    int bin_total = 0; // patch_count;

    for (int i = 0; i < patch_count; i++) {
      int bin_number = 0;

      for (int is = 0; is < n_species; is++) {
        int n_particles = patches_[i].particles_m[is].size();
        bin_number += 1 + ((n_particles - 1) / bin_size);
      }
      evolve.emplace_back(bin_number + 1);

      if (bin_number > 0) {
        bin_total += bin_number; // reducing patch count,
                                 //  do for case no particles
      } else {
        ++bin_total;
      }
    }
    currentBC_antenna_imbalance.emplace_back(bin_total);
  }
};

//////////////////////////////////////////////////////////////////////////////////////////

namespace impl {
auto particle_binning(task_graph &task_graph,
                      Params &params,
                      Timers &timers,
                      Profiler &profiler,
                      std::vector<Patch> &patches_,
                      int &it) -> void {
  operators::particle_binning(params, timers, profiler, patches_, it);
  task_graph.completion.count_down();
}

auto particle_cloud(task_graph &task_graph,
                    Params &params,
                    Timers &timers,
                    Profiler &profiler,
                    std::vector<Patch> &patches_,
                    int &it) -> void {
  operators::particle_cloud(params, timers, profiler, patches_, it);
  task_graph.completion.count_down();
}

auto particle_scalars(task_graph &task_graph,
                      Params &params,
                      Timers &timers,
                      Profiler &profiler,
                      std::vector<Patch> &patches_,
                      int &it) -> void {
  operators::particle_scalars(params, timers, profiler, patches_, it);
  task_graph.completion.count_down();
}

auto diags_scalars(task_graph &task_graph,
                   Params &params,
                   Timers &timers,
                   Profiler &profiler,
                   ElectroMagn &em,
                   int it) -> void {
  operators::diags_scalars(params, timers, profiler, em, it);
  task_graph.completion.count_down();
}

auto diags_fields(task_graph &task_graph,
                  Params &params,
                  Timers &timers,
                  Profiler &profiler,
                  ElectroMagn &em,
                  int it) -> void {
  operators::diags_fields(params, timers, profiler, em, it);
  task_graph.completion.count_down();
}

auto terminal_print(task_graph &task_graph,
                    Params &params,
                    Timers &timers,
                    Profiler &profiler,
                    std::vector<Patch> &patches_,
                    int &it) -> void {
  operators::terminal_print(params, timers, profiler, patches_, it);
  task_graph.completion.count_down();
}

auto maxwell_solverBC(task_system &task_system,
                      task_graph &task_graph,
                      Params &params,
                      Timers &timers,
                      Profiler &profiler,
                      ElectroMagn &em,
                      int &it) -> void {

  timers.stop(timers.maxwell_solver);

  // __________________________________________________________________
  // Maxwell Boundary conditions

  profiler.start(MAXWELL);
  timers.start(timers.maxwellBC);

  DEBUG("start solve BC")

  // Boundary conditions on EM fields
  operators::solveBC(params, em);

  DEBUG("end solve BC")
  timers.stop(timers.maxwellBC);
  profiler.stop();

  // __________________________________________________________________
  // Schedule next section of graph

  if (task_graph.field_diags.count_down()) {
    task_system.submit([&] { impl::diags_scalars(task_graph, params, timers, profiler, em, it); });

    task_system.submit([&] { impl::diags_fields(task_graph, params, timers, profiler, em, it); });
  }
}

auto maxwell_faraday_scheduler(task_system &task_system,
                               task_graph &task_graph,
                               Params &params,
                               Timers &timers,
                               Profiler &profiler,
                               ElectroMagn &em,
                               int &it,
                               const double dt,        
                               const double dt_over_dx,
                               const double dt_over_dy,
                               const double dt_over_dz) -> void {

  // Magnetic field Bx (p,d,d)
  for (unsigned int i = 0; i < em.nx_p_m; i++) {

    task_system.submit([&, i, dt, dt_over_dx, dt_over_dy, dt_over_dz] {
      operators::solve_maxwell_faraday_bx_2d(params, timers, profiler, em, i,
                                        dt, dt_over_dy, dt_over_dz);

      if (task_graph.maxwell_solverBC.count_down()) {
        task_system.submit([&] {
          impl::maxwell_solverBC(task_system, task_graph, params, timers, profiler, em, it);
        });
      
      }
    
    });
  
  }

  // Magnetic field By (d,p,d)
  for (unsigned int i = 1; i < em.nx_d_m - 1; i++) {

    task_system.submit([&, i, dt, dt_over_dx, dt_over_dy, dt_over_dz] {
      operators::solve_maxwell_faraday_by_2d(params, timers, profiler, em, i,
                                      dt, dt_over_dx, dt_over_dz);

      if (task_graph.maxwell_solverBC.count_down()) {
        task_system.submit([&] {
          impl::maxwell_solverBC(task_system, task_graph, params, timers, profiler, em, it);
        });
     
      }
    
    });
  
  }

  // Magnetic field Bz (d,d,p)
  for (unsigned int i = 1; i < em.nx_d_m - 1; i++) {

    task_system.submit([&, i, dt, dt_over_dx, dt_over_dy, dt_over_dz] {
      operators::solve_maxwell_faraday_bz_2d(params, timers, profiler, em, i,
                                      dt, dt_over_dx, dt_over_dy);

      if (task_graph.maxwell_solverBC.count_down()) {
        task_system.submit([&] {
          impl::maxwell_solverBC(task_system, task_graph, params, timers, profiler, em, it);
        });
      
      }
    
    });
  
  }

} // End function

auto maxwell_ampere_scheduler(task_system &task_system,
                              task_graph &task_graph,
                              Params &params,
                              Timers &timers,
                              Profiler &profiler,
                              ElectroMagn &em,
                              int &it,
                              const double dt,        
                              const double dt_over_dx,
                              const double dt_over_dy,
                              const double dt_over_dz) -> void {

  // Electric field Ex (p,d,p)
  for (unsigned int i = 0; i < em.nx_d_m; i++) {

    task_system.submit([&, i, dt, dt_over_dx, dt_over_dy, dt_over_dz] {
      operators::solve_maxwell_ampere_ex_2d(params, timers, profiler, em, i,
                                      dt, dt_over_dy, dt_over_dz);

      if (task_graph.maxwell_faraday_scheduler.count_down()) {
        impl::maxwell_faraday_scheduler(task_system, task_graph, params, timers, profiler, em, it,
                                                            dt, dt_over_dx, dt_over_dy, dt_over_dz);
      }
    
    });
  
  }

  // Electric field Ey (p,d,p)
  for (unsigned int i = 0; i < em.nx_p_m; i++) {

    task_system.submit([&, i, dt, dt_over_dx, dt_over_dy, dt_over_dz] {
      operators::solve_maxwell_ampere_ey_2d(params, timers, profiler, em, i,
                                      dt, dt_over_dx, dt_over_dz);

      if (task_graph.maxwell_faraday_scheduler.count_down()) {
        impl::maxwell_faraday_scheduler(task_system, task_graph, params, timers, profiler, em, it,
                                                            dt, dt_over_dx, dt_over_dy, dt_over_dz);
      }

    });
  
  }

  // Electric field Ez (p,p,d)
  for (unsigned int i = 0; i < em.nx_p_m; i++) {
    
    task_system.submit([&, i, dt, dt_over_dx, dt_over_dy, dt_over_dz] {
      operators::solve_maxwell_ampere_ez_2d(params, timers, profiler, em, i,
                                      dt, dt_over_dx, dt_over_dy);

      if (task_graph.maxwell_faraday_scheduler.count_down()) {
        impl::maxwell_faraday_scheduler(task_system, task_graph, params, timers, profiler, em, it,
                                                            dt, dt_over_dx, dt_over_dy, dt_over_dz);
      }
    
    });
  
  }

} // End function

auto currentBC_antenna(task_system &task_system,
                       task_graph &task_graph,
                       Params &params,
                       Timers &timers,
                       Profiler &profiler,
                       ElectroMagn &em,
                       int &it) -> void {

  // ______________________________________________________
  // Start Maxwell

  if (params.maxwell_solver) {

    if (params.current_projection || params.n_particles > 0) {

      timers.start(timers.currentBC);
      profiler.start(CURRENTBC);

      // Perform the boundary conditions for current
      DEBUG("start current BC")
      operators::currentBC(params, em);
      DEBUG("end current BC")

      profiler.stop();
      timers.stop(timers.currentBC);

    } // end if

    timers.start(timers.maxwell_solver);

    // Generate a laser field with an antenna
    for (auto iantenna = 0; iantenna < params.antenna_profiles_.size(); iantenna++) {
      operators::antenna(params,
                         em,
                         params.antenna_profiles_[iantenna],
                         params.antenna_positions_[iantenna],
                         it * params.dt);
    }

    // __________________________________________________________________
    // Data to share

    const double dt         = params.dt;
    const double dt_over_dx = params.dt * params.inv_dx;
    const double dt_over_dy = params.dt * params.inv_dy;
    const double dt_over_dz = params.dt * params.inv_dz;

    // __________________________________________________________________
    // Schedule next section of graph

    if (task_graph.maxwell_ampere_scheduler.count_down()) {
      task_system.submit([&, dt, dt_over_dx, dt_over_dy, dt_over_dz] {
        impl::maxwell_ampere_scheduler(task_system, task_graph, params, timers, profiler, em, it,
                                        dt, dt_over_dx, dt_over_dy, dt_over_dz);
      });
    }
  } //end if
  else
  {
    task_system.submit([&] { impl::diags_scalars(task_graph, params, timers, profiler, em, it); });
    task_system.submit([&] { impl::diags_fields(task_graph, params, timers, profiler, em, it); });
  }

}

auto imbalance_schedule(task_system &task_system,
                        task_graph &task_graph,
                        Params &params,
                        Timers &timers,
                        Profiler &profiler,
                        std::vector<Patch> &patches_,
                        ElectroMagn &em,
                        const std::ptrdiff_t patch_index,
                        int &it) -> void {
  // timers.start(params, "evolve_bin", patch_index);
  int n_species = patches_[patch_index].particles_m.size();

  for (int is = 0; is < n_species; is++) {
    int n_particles = patches_[patch_index].particles_m[is].size();

    if (n_particles > 0) {

      int bin_size   = params.bin_size;
      int bin_number = 1 + ((n_particles - 1) / bin_size);

      for (unsigned int i_bin = 0; i_bin < bin_number; i_bin++) {

        int rest = n_particles - i_bin * bin_size;
        int size = std::min(bin_size, rest);

        int init = i_bin * bin_size;        // start of bin
        int end  = i_bin * bin_size + size; // end of bin

        task_system.submit([&, patch_index, is, init, end, rest, i_bin] {
          if (rest > 0) {
            timers.start(timers.imbalance, patch_index);
            operators::imbalance_operator(params,
                                          patches_[patch_index].particles_m[is],
                                          init,
                                          end,
                                          it,
                                          params.imbalance_function_[0]);
            timers.stop(timers.imbalance, patch_index);
          }

          if (task_graph.currentBC_antenna_imbalance[0].count_down()) {
            task_system.submit([&] {
              impl::currentBC_antenna(task_system, task_graph, params, timers, profiler, em, it);
            });
          }
        });
      }
    } else {
      if (task_graph.currentBC_antenna_imbalance[0].count_down()) {
        task_system.submit([&] {
          impl::currentBC_antenna(task_system, task_graph, params, timers, profiler, em, it);
        });
      }
    }
  }

  if (n_species < 1) // case no particles
  {
    if (task_graph.currentBC_antenna_imbalance[0].count_down()) {
      task_system.submit([&] {
        impl::currentBC_antenna(task_system, task_graph, params, timers, profiler, em, it);
      });
    }
  }
  // timers.stop(params, "evolve_bin", patch_index);
}

auto projection_borders(task_system &task_system,
                       task_graph &task_graph,
                       Params &params,
                       Timers &timers,
                       Profiler &profiler,
                       std::vector<Patch> &patches_,
                       ElectroMagn &em,
                       std::mutex &mutex,
                       int &it,
                       const std::ptrdiff_t patch_index) -> void {

  // __________________________________________________________________
  // Sum all species contribution in the local and global current grids

  if (params.current_projection || params.n_particles > 0) {

    // __________________________________________________________________
    // Projection to global grid internal
    timers.start(timers.current_global_reduc, patch_index);
    profiler.start(CURRENT_GLOBAL_BORDERS);
    // Copy all local fields in the global fields
    operators::local2global_borders(em, patches_[patch_index]);
    profiler.stop();
    timers.stop(timers.current_global_reduc, patch_index);

  }

  // __________________________________________________________________
  // Schedule next section of graph

  if (task_graph.imbalance_operator.count_down() && !params.imbalance_function_.empty()) {

    const auto patch_count = std::ssize(patches_);

    for (auto patch_index : std::views::iota(std::ptrdiff_t{0}, patch_count)) {
      impl::imbalance_schedule(task_system,
                               task_graph,
                               params,
                               timers,
                               profiler,
                               patches_,
                               em,
                               patch_index,
                               it);
    }
  } else if (task_graph.currentBC_antenna.count_down()) {
    task_system.submit(
      [&] { impl::currentBC_antenna(task_system, task_graph, params, timers, profiler, em, it); });
  }
}


auto projection_internal(task_system &task_system,
                       task_graph &task_graph,
                       Params &params,
                       Timers &timers,
                       Profiler &profiler,
                       std::vector<Patch> &patches_,
                       ElectroMagn &em,
                       std::mutex &mutex,
                       int &it,
                       const std::ptrdiff_t patch_index) -> void {

  // __________________________________________________________________
  // Sum all species contribution in the local and global current grids

  if (params.current_projection || params.n_particles > 0) {

    // __________________________________________________________________
    // Projection to global grid internal
    timers.start(timers.current_global_reduc, patch_index);
    profiler.start(CURRENT_GLOBAL_BORDERS);
    // Copy all local fields in the global fields
    operators::local2global_internal(em, patches_[patch_index]);
    profiler.stop();
    timers.stop(timers.current_global_reduc, patch_index);

  }

  // __________________________________________________________________
  // Schedule next section of graph

  if (task_graph.imbalance_operator.count_down() && !params.imbalance_function_.empty()) {

    const auto patch_count = std::ssize(patches_);

    for (auto patch_index : std::views::iota(std::ptrdiff_t{0}, patch_count)) {
      impl::imbalance_schedule(task_system,
                               task_graph,
                               params,
                               timers,
                               profiler,
                               patches_,
                               em,
                               patch_index,
                               it);
    }
  } else if (task_graph.currentBC_antenna.count_down()) {
    task_system.submit(
      [&] { impl::currentBC_antenna(task_system, task_graph, params, timers, profiler, em, it); });
  }
}


auto exchange(task_system &task_system,
              task_graph &task_graph,
              Params &params,
              Timers &timers,
              Profiler &profiler,
              std::vector<Patch> &patches_,
              ElectroMagn &em,
              int &it,
              std::mutex &mutex,
              const std::ptrdiff_t patch_index) -> void {

  // __________________________________________________________________
  // Exchanges of particles in patch

  timers.start(timers.exchange, patch_index);
  DEBUG("Patch " << patch_index << ": exchange");

  operators::exchange_particles(params, patches_, patch_index);

  DEBUG("Patch " << patch_index << ": exchange");
  timers.stop(timers.exchange, patch_index);

  // __________________________________________________________________
  // Schedule diags output to results

  if (task_graph.particles_print.count_down()) {
    task_system.submit(
      [&] { impl::particle_binning(task_graph, params, timers, profiler, patches_, it); });
    task_system.submit(
      [&] { impl::particle_cloud(task_graph, params, timers, profiler, patches_, it); });
    task_system.submit(
      [&] { impl::particle_scalars(task_graph, params, timers, profiler, patches_, it); });
    task_system.submit(
      [&] { impl::terminal_print(task_graph, params, timers, profiler, patches_, it); });
  }
}

// Function to control the launch of exchanges by neighbors
auto exchange_schedule(task_system &task_system,
                       task_graph &task_graph,
                       std::vector<Patch> &patches_,
                       Params &params,
                       Timers &timers,
                       Profiler &profiler,
                       ElectroMagn &em,
                       int &it,
                       std::mutex &mutex,
                       const std::ptrdiff_t patch_index) -> void {
  int i_patch = patches_[patch_index].i_patch_topology_m;
  int j_patch = patches_[patch_index].j_patch_topology_m;
  int k_patch = patches_[patch_index].k_patch_topology_m;

  // Collect new particles from neighbours
  for (int i = -1; i < 2; i++) {
    for (int j = -1; j < 2; j++) {
      for (int k = -1; k < 2; k++) {

        // id of the Neighbor in vec_patch
        int idx_neighbor =
          patches_[patch_index].get_idx_patch(i_patch + i, j_patch + j, k_patch + k);

        if (task_graph.exchange[idx_neighbor].count_down()) {
          task_system.submit([&, idx_neighbor] {
            impl::exchange(task_system,
                           task_graph,
                           params,
                           timers,
                           profiler,
                           patches_,
                           em,
                           it,
                           mutex,
                           idx_neighbor);
          });
        }
      }
    }
  } // end for each neighbors
}

auto evolve_patch(task_system &task_system,
                  task_graph &task_graph,
                  std::vector<Patch> &patches_,
                  Params &params,
                  Backend &backend,
                  Timers &timers,
                  Profiler &profiler,
                  ElectroMagn &em,
                  int &it,
                  std::mutex &mutex,
                  const std::ptrdiff_t patch_index) -> void {

  profiler.start(EVOLVE_PATCH);

  // __________________________________________________________________
  // Projection in local field

  if (params.current_projection) {
    timers.start(timers.projection, patch_index);
    DEBUG("start project");

    // Project in buffers local to the patches
    operators::project(params, patches_[patch_index]);

    DEBUG("stop project");
    timers.stop(timers.projection, patch_index);
  }

  // __________________________________________________________________
  // Identify and copy in buffers particles which leave the patch

  timers.start(timers.id_parts_to_move, patch_index);
  DEBUG("Patch " << patch_index << ": start identify particles to move");

  if (patches_.size() > 1) {
    operators::identify_particles_to_move(params, patches_[patch_index], backend);
  }

  DEBUG("Patch " << patch_index << ": stop identify particles to move");
  timers.stop(timers.id_parts_to_move, patch_index);

  ////

  profiler.stop();

  if (params.current_projection || params.n_particles > 0) {

    // __________________________________________________________________
    // Projection in the local grid
    timers.start(timers.current_local_reduc, patch_index);
    // Sum all species contribution in the local fields
    operators::reduc_current(patches_[patch_index]);
    timers.stop(timers.current_local_reduc, patch_index);

  }

  task_system.submit([&, patch_index] {
    impl::projection_internal(task_system,
                            task_graph,
                            params,
                            timers,
                            profiler,
                            patches_,
                            em,
                            mutex,
                            it,
                            patch_index);
  });

  task_system.submit([&, patch_index] {
    impl::projection_borders(task_system,
                            task_graph,
                            params,
                            timers,
                            profiler,
                            patches_,
                            em,
                            mutex,
                            it,
                            patch_index);
  });

  // __________________________________________________________________
  // Exchange function to avoid synchronization stages

  impl::exchange_schedule(task_system,
                          task_graph,
                          patches_,
                          params,
                          timers,
                          profiler,
                          em,
                          it,
                          mutex,
                          patch_index);

}

auto evolve_schedule(task_system &task_system,
                     task_graph &task_graph,
                     std::vector<Patch> &patches_,
                     Params &params,
                     Backend &backend,
                     Timers &timers,
                     Profiler &profiler,
                     ElectroMagn &em,
                     std::mutex &mutex,
                     int &it,
                     const std::ptrdiff_t patch_index) -> void {

  int n_species = patches_[patch_index].particles_m.size();

  for (int is = 0; is < n_species; is++) {
    int n_particles = patches_[patch_index].particles_m[is].size();

    if (n_particles > 0) {

      int bin_size   = params.bin_size;
      int bin_number = 1 + ((n_particles - 1) / bin_size);

      for (unsigned int i_bin = 0; i_bin < bin_number; i_bin++) {

        int rest = n_particles - i_bin * bin_size;
        int size = std::min(bin_size, rest);

        int init = i_bin * bin_size;        // start of bin
        int end  = i_bin * bin_size + size; // end of bin

        task_system.submit([&, patch_index, n_species, is, init, end, i_bin, rest] {
          if (rest > 0) {

            profiler.start(EVOLVE_BIN);

            ////

            // functions using bin
            timers.start(timers.interpolate, patch_index);
            operators::interpolate_bin(em, patches_[patch_index].particles_m[is], is, init, end);
            timers.stop(timers.interpolate, patch_index);

            timers.start(timers.push, patch_index);
            operators::push_bin(params.dt, patches_[patch_index].particles_m[is], is, init, end);
            timers.stop(timers.push, patch_index);

            timers.start(timers.pushBC, patch_index);
            operators::pushBC_bin(params,
                                  patches_[patch_index],
                                  patches_[patch_index].particles_m[is],
                                  patches_[patch_index].on_border_m,
                                  is,
                                  init,
                                  end);
            timers.stop(timers.pushBC, patch_index);

            ////

            profiler.stop();

          } // end rest particles

          if (task_graph.evolve[patch_index].count_down()) {
            task_system.submit([&, patch_index] {
              impl::evolve_patch(task_system,
                                 task_graph,
                                 patches_,
                                 params,
                                 backend,
                                 timers,
                                 profiler,
                                 em,
                                 it,
                                 mutex,
                                 patch_index);
            });
          } // end next task
        }); // end task
      }     // end bin loop
    }       // end if end if n particles
    else {
      if (task_graph.evolve[patch_index].count_down()) {
        task_system.submit([&, patch_index] {
          impl::evolve_patch(task_system,
                             task_graph,
                             patches_,
                             params,
                             backend,
                             timers,
                             profiler,
                             em,
                             it,
                             mutex,
                             patch_index);
        });
      } // end next task

    } // end if end if n particles
  }   // end loop species

  if (n_species < 1) {
    task_system.submit([&, patch_index] {
      impl::evolve_patch(task_system,
                         task_graph,
                         patches_,
                         params,
                         backend,
                         timers,
                         profiler,
                         em,
                         it,
                         mutex,
                         patch_index);
    });
  }
}

auto reset_current(task_system &task_system,
                   task_graph &task_graph,
                   std::vector<Patch> &patches_,
                   Params &params,
                   Backend &backend,
                   Timers &timers,
                   Profiler &profiler,
                   ElectroMagn &em,
                   std::mutex &mutex,
                   int &it) -> void {

  operators::reset_current(params, timers, profiler, em);

  const auto patch_count = std::ssize(patches_);

  for (auto patch_index : std::views::iota(std::ptrdiff_t{0}, patch_count)) {
    if (task_graph.evolve[patch_index].count_down()) {
      task_system.submit([&, patch_index] {
        impl::evolve_patch(task_system,
                           task_graph,
                           patches_,
                           params,
                           backend,
                           timers,
                           profiler,
                           em,
                           it,
                           mutex,
                           patch_index);
      });
    }
  }
};
} // namespace impl

//////////////////////////////////////////////////////////////////////////////////////////

auto sync_exec(task_graph &task_graph,
               SubDomain &subdomain,
               Params &params,
               Timers &timers,
               Profiler &profiler,
               Backend &backend,
               int &it) -> void {

  auto &task_system = (*backend.task_system_);
  auto &mutex       = (*backend.mutex);

task_system.submit([&] {
  impl::reset_current(task_system,
                      task_graph,
                      subdomain.patches_,
                      params,
                      backend,
                      timers,
                      profiler,
                      subdomain.em_,
                      mutex,
                      it);
});

  const auto patch_count = std::ssize(subdomain.patches_);

  for (auto patch_index : std::views::iota(std::ptrdiff_t{0}, patch_count)) {
    impl::evolve_schedule(task_system,
                          task_graph,
                          subdomain.patches_,
                          params,
                          backend,
                          timers,
                          profiler,
                          subdomain.em_,
                          mutex,
                          it,
                          patch_index);
  }

  task_graph.completion.wait();
}
} // namespace mini_pic
