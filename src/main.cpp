/* _____________________________________________________________________ */
//! \file Main.cpp

//  _ __ ___ (_)_ __ (_)_ __ (_) ___
// | '_ ` _ \| | '_ \| | '_ \| |/ __|
// | | | | | | | | | | | |_) | | (__
// |_| |_| |_|_|_| |_|_| .__/|_|\\___|
//                     |_|

//! \brief Main file for Minipic
//! NOTE: this program is intended for computer science studies
//! and should be not used for physics simulation
/* _____________________________________________________________________ */

#include <sys/time.h>

#include "Backend.hpp"
#include "Params.hpp"
#include "SubDomain.hpp"

#if defined(__SHAMAN__)
#include "Shaman_analysis.hpp"
#endif

// Load a setup
// #include "study_thermal.hpp"
#include "study_sorted.hpp"
//#include "study_random.hpp"

//! Main function
int main(int argc, char *argv[]) {

  // ___________________________________________
  // Setup input parameters in struct

  // Create the gloal parameters
  Params params;

  // Print the Minipic title
  params.title();

  // default parameters
  setup(params);

  // change from command line arguments
  params.read_from_command_line_arguments(argc, argv);

  // Initialize the main parameters for the simulation
  params.compute();

  // Print a summary of input parameters
  params.info();

  // Create the backend parameters
  Backend backend;

  // Initialize the backend
  backend.init(argc, argv, params);

  // Print the backend information
  backend.info();

  {

    // Timers initialization
    Timers timers(params);
    //timers.start(timers.initialization);

    // ______________________________________________________
    //
    // Initialization
    // ______________________________________________________

    // timers.start(timers. "initialization", 0);

    SubDomain subdomain;

    // Creation of the domain
    subdomain.allocate(params, backend);

    // Initialization of the diagnostics
    Diags::initialize(params);

    // Initialize a profiler
    Profiler profiler(params);

    timers.stop(timers.initialization);
    timers.save_initialization(params);
    timers.start(timers.main_loop);

    // ______________________________________________________
    //
    // Initial diagnostics
    // ______________________________________________________

    timers.start(timers.diags);
    subdomain.diagnostics(params, timers, profiler, backend, 0);
    timers.stop(timers.diags);

    // ______________________________________________________
    //
    // main PIC loop
    // ______________________________________________________

    printf("\n #########    START COMPUTE    ############\n");
    std::cout << " -----------------------------------------------------|" << std::endl;
    std::cout << "           |       |           | Elapsed  | Remaining |" << std::endl;
    std::cout << " Iteration |    %  | Particles | Time (s) | Time (s)  |" << std::endl;
    std::cout << " ----------|-------|-----------|----------|-----------|" << std::endl;

    timers.start(timers.main_loop);

#if defined(__MINIPIC_OMP__) || defined(__MINIPIC_OMP_TASK__)
#pragma omp parallel default(none) firstprivate(params) \
  shared(subdomain, timers, profiler, backend, std::cout)
    {
#endif

      for (int it = 1; it <= params.n_it; it++) {

        // _______________________________________________________
        // Main loop for all programming models except eventify

#if !defined(__MINIPIC_EVENTIFY__)

#if defined(__MINIPIC_OMP_TASK__)
#pragma omp single
        {
#endif

#if defined(__MINIPIC_OMP__)
#pragma omp single
#endif
          { timers.start(timers.pic_iteration); }

          // Single PIC iteration
          subdomain.iterate(params, timers, profiler, backend, it);

#if defined(__MINIPIC_OMP__)
#pragma omp single
#endif
          {
            timers.stop(timers.pic_iteration);
            timers.start(timers.diags);
          }

          // Diagnostics
          subdomain.diagnostics(params, timers, profiler, backend, it);
#if defined(__MINIPIC_OMP__)
#pragma omp single
#endif
          { timers.stop(timers.diags); }

#if defined(__MINIPIC_OMP__)
#pragma omp single // single nowait
#elif defined(__MINIPIC_OMP_TASK__)
#pragma omp task untied default(none) shared(params, subdomain, it, timers, std::cout)
#endif
          {
            if (!(it % params.print_period)) {

              const unsigned int total_number_of_particles =
                subdomain.get_total_number_of_particles();

              double elapsed_time   = timers.get_elapsed_time();
              double remaining_time = elapsed_time / it * (params.n_it - it);

              std::cout << " " << std::setw(9) << it;
              std::cout << " | " << std::fixed << std::setprecision(1) << std::setw(5)
                        << static_cast<float>(it) / static_cast<float>(params.n_it) * 100;
              std::cout << " | " << std::scientific << std::setprecision(2) << std::setw(9)
                        << total_number_of_particles;
              std::cout << " | " << std::scientific << std::setprecision(2) << std::setw(8)
                        << elapsed_time;
              std::cout << " | " << std::scientific << std::setprecision(2) << std::setw(9)
                        << remaining_time;
              std::cout << " | " << std::endl;
            }

            timers.save(params, it);
          }

#if defined(__MINIPIC_OMP_TASK__)
#pragma omp taskwait
        } // end single region
#endif

#elif defined(__MINIPIC_EVENTIFY__)

      // _______________________________________________________
      // Main loop for the eventify programming model
      auto task_graph = mini_pic::task_graph{
        params.N_patches,
        params.bin_size,
        subdomain.patches_,
        static_cast<ptrdiff_t>(params.get_species_number()),
        2 * (subdomain.em_.nx_d_m - 2) + subdomain.em_.nx_p_m, // n_task for maxwell bc
        2 * subdomain.em_.nx_p_m + subdomain.em_.nx_d_m        // n_task for maxwell faraday x
      };

      mini_pic::sync_exec(task_graph, subdomain, params, timers, profiler, backend, it);

#endif

      } // end main loop

#if defined(__MINIPIC_OMP__) || defined(__MINIPIC_OMP_TASK__)
    } // end parallel region
#endif

    DEBUG("End of main loop");

    timers.stop(timers.main_loop);

    std::cout << "#########    END COMPUTE    ############\n\n";
    printf("\n");

    // ____________________________________________________
    // Print timers

    timers.print(params);
    timers.save(params, params.n_it + 1);

    // ____________________________________________________
    // Dump the profiler
    profiler.dump(params);

    // ____________________________________________________
    // Shaman analysis

#if defined(__SHAMAN__)
    shaman_summary(subdomain);
#endif
  }
  // Kokkos::finalize();

  backend.finalize();

  return 0;
}
