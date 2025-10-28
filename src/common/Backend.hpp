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

#elif defined(__MINIPIC_SYCL__)

  sycl::queue *sycl_queue_;

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
  void init([[maybe_unused]] int argc,
            [[maybe_unused]] char *argv[],
            [[maybe_unused]] const Params &params) {

#if defined(__MINIPIC_OMP__) || defined(__MINIPIC_OMP_TASK__) || defined(__MINIPIC_EVENTIFY__)
    number_of_threads = omp_get_max_threads();
#else
    number_of_threads = 1;
#endif

#if defined(__MINIPIC_OMP_TASK__)

    // Flag used for task dependency
    evolve_particles_flags   = new int[params.N_patches];
    reduction_internal_flags = new int[params.N_patches];
    reduction_external_flags = new int[params.N_patches];
    evolve_patch_flags       = new int[params.N_patches];

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

#elif defined(__MINIPIC_KOKKOS_COMMON__)

    Kokkos::initialize(argc, argv);

#if defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__) || defined(__MINIPIC_KOKKOS_UNIFIED__)
    static_assert(Kokkos::has_shared_space, "code only works on backends with SharedSpace");
#endif

#elif defined(__MINIPIC_SYCL__)

    // if (params.on_gpu_) {
    //   sycl_queue_ = new sycl::queue{sycl::gpu_selector_v};
    // } else {
    //   sycl_queue_ = new sycl::queue{sycl::cpu_selector_v};
    // }

#if defined(__MINIPIC_ON_GPU__)
    sycl_queue_ = new sycl::queue{sycl::gpu_selector_v};
#else
    sycl_queue_ = new sycl::queue{sycl::cpu_selector_v};
#endif

#endif
  }

  // _____________________________________________________________________
  //
  //! \brief Finalize the backend
  // _____________________________________________________________________
  void finalize() {
#if defined(__MINIPIC_OMP_TASK__)
    delete[] evolve_particles_flags;
#elif defined(__MINIPIC_KOKKOS_COMMON__)
    Kokkos::finalize();
#elif defined(__MINIPIC_SYCL__)
    delete sycl_queue_;
#endif
  }

  // _____________________________________________________________________
  //
  //! \brief Print the backend information
  // _____________________________________________________________________
  void info() {
    std::cout << " > Backend: " << std::endl;
#ifdef BACKEND
    std::cout << "   - CMake Name: " << BACKEND << std::endl;
#endif

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

#if defined(__MINIPIC_KOKKOS_COMMON__)
    std::cout << "   - Selected parallel programming model: Kokkos" << std::endl;
    std::cout << "   - Device execution Space: " << typeid(Kokkos::DefaultExecutionSpace).name()
              << std::endl;
    std::cout << "   - Host execution Space: " << typeid(Kokkos::DefaultHostExecutionSpace).name()
              << std::endl;
    // std::cout << "   - Number of threads: " << &Kokkos::num_threads << std::endl;
    Kokkos::print_configuration(std::cout);
#endif

#if defined(__MINIPIC_THRUST__)
    std::cout << "   - Selected parallel programming model: Thrust" << std::endl;
#endif

#if defined(__MINIPIC_THRUST_UNIFIED__)

    std::cout << "   - Selected parallel programming model: Thrust using unified memory"
              << std::endl;

#elif defined(__MINIPIC_SYCL__)

    std::cout << "   - Selected parallel programming model: SYCL" << std::endl;
    std::cout << "   - Device: " << sycl_queue_->get_device().get_info<sycl::info::device::name>()
              << std::endl;
    // Returns the number of parallel compute units available to the device
    std::cout << "   - Max compute units = "
              << sycl_queue_->get_device().get_info<sycl::info::device::max_compute_units>()
              << std::endl;
    // Returns the maximum number of work-items that are permitted in each dimension of a work-group
    // for a kernel running in a three-dimensional index space
    sycl::id<3> wi_sizes =
      sycl_queue_->get_device().get_info<sycl::info::device::max_work_item_sizes<3>>();
    std::cout << "   - Max work item in each dimension of a work group kernel 3D = " << wi_sizes[0]
              << "  " << wi_sizes[1] << "  " << wi_sizes[2] << std::endl;
    // Returns the maximum number of work-items that are permitted in a work-group executing a
    // kernel on a single compute unit
    std::cout << "   - Max work group size = "
              << sycl_queue_->get_device().get_info<sycl::info::device::max_work_group_size>()
              << std::endl;
    // Returns the maximum number of sub-groups in a work-group for any kernel executed on the
    // device
    std::cout << "   - Max number of sub-groups in 1 work group = "
              << sycl_queue_->get_device().get_info<sycl::info::device::max_num_sub_groups>()
              << std::endl;
    // Returns a std::vector of size_t containing the set of sub-group sizes supported by the device
    std::vector<size_t> sg_sizes =
      sycl_queue_->get_device().get_info<sycl::info::device::sub_group_sizes>();
    std::cout << "   - Sub group sizes supported : " << sg_sizes[0] << "  " << sg_sizes[1] << "  "
              << sg_sizes[2] << std::endl;
    std::cout << std::endl;

#elif defined(__MINIPIC_OPENACC__)

    std::cout << "   - Selected parallel programming model: OpenACC" << std::endl;

#endif
  }
};

#endif
