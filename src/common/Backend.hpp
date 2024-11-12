/* _____________________________________________________________________ */
//! \file Backend.hpp

//! \brief determine the best backend to use

/* _____________________________________________________________________ */

#ifndef BACKEND_H
#define BACKEND_H

#include "Params.hpp"

// _____________________________________________________________________
//
// Backend class
//
//! \brief manage the backend properties
// _____________________________________________________________________

class Backend {
public:
  // _____________________________________________________________________
  // Public parameters

  int number_of_threads;

#if defined(__MINIPIC_OMP_TASK__)

    int *evolve_particles_flags;
    int *reset_current_flags;
    int *maxwell_solver_flags;
    int *reduction_internal_flags;
    int *reduction_external_flags;
    int *evolve_patch_flags;

    // For exchange counters
    std::deque<std::atomic<int>> task_exchange_count_;

#elif defined(__MINIPIC_EVENTIFY__)

  shared_ptr<eventify::task_system> task_system_;
  shared_ptr<std::mutex> mutex;

#endif

  // _____________________________________________________________________
  // Public methods

  Backend() {}

  ~Backend() {}

  // _____________________________________________________________________
  //
  //! \brief Initialize the backend
  //! \param argc number of arguments
  //! \param argv arguments
  //! \param params global parameters
  // _____________________________________________________________________
  void init(int argc, char *argv[], const Params &params) {

#if defined(__MINIPIC_OMP__) || defined(__MINIPIC_OMP_TASK__) || defined(__MINIPIC_EVENTIFY__)
    number_of_threads = omp_get_max_threads();
#else
    number_of_threads = 1;
#endif

#if defined(__MINIPIC_OMP_TASK__)

    // Flag used for task dependency
    evolve_particles_flags = new int[params.N_patches];
    reduction_internal_flags = new int[params.N_patches];
    reduction_external_flags = new int[params.N_patches];
    evolve_patch_flags = new int[params.N_patches];

    for (unsigned int i_patch = 0; i_patch < params.N_patches; i_patch++) {
      evolve_particles_flags[i_patch] = i_patch;
      task_exchange_count_.emplace_back(27);
    }

#elif defined(__MINIPIC_EVENTIFY__)

    // task_system = new eventify::task_system{static_cast<unsigned int>(number_of_threads)};
    // task_system_ = std::make_shared<eventify::task_system{static_cast<unsigned
    // int>(number_of_threads)}>();
    task_system_ =
      std::make_shared<eventify::task_system>(static_cast<unsigned int>(number_of_threads));
    mutex = std::make_shared<std::mutex>();

#endif
  }

  // _____________________________________________________________________
  //
  //! \brief Finalize the backend
  // _____________________________________________________________________
  void finalize() {
#if defined(__MINIPIC_OMP_TASK__)
    delete[] evolve_particles_flags;
#endif
  }

  // _____________________________________________________________________
  //
  //! \brief Print the backend information
  // _____________________________________________________________________
  void info() {
    std::cout << " > Backend: " << std::endl;

#if defined(__MINIPIC_OMP__)
    std::cout << "   - Selected parallel programming model: OpenMP for" << std::endl;
    std::cout << "   - OMP number of threads: " << number_of_threads << std::endl;
#endif

#if defined(__MINIPIC_OMP_TASK__)
    std::cout << "   - Selected parallel programming model: OpenMP Task" << std::endl;
    std::cout << "   - OMP number of threads: " << number_of_threads << std::endl;
#endif

#if defined(__MINIPIC_EVENTIFY__)
    std::cout << "   - Selected parallel programming model: Eventify" << std::endl;
    std::cout << "   - Eventify number of threads: " << number_of_threads << std::endl;
#endif

  }
};

#endif
