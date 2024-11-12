/* _____________________________________________________________________ */
//! \file Timers.cpp

//! \brief Timer class to measure the time spent in each part of the code

/* _____________________________________________________________________ */

#pragma once

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <vector>

#include "Params.hpp"

namespace level {
// Timers's level
enum { global = 0, thread = 1 };
} // namespace level

// _______________________________________________________________
//
//! Structure representing a code section
// _______________________________________________________________
struct Section {
  int id;
  std::string name;
};

// _______________________________________________________________
//
//! \brief Timer class to measure the time spent in each part of the code
//! \details each timer contains N_patches + 1 counters :
//! - the first one is the global counter
//! - the others are the thread counters (one per patch)
// _______________________________________________________________
template <class T_clock> class Struc_Timers {
public:
  // List of available code sections
  const Section initialization = Section{0, "initialization"};
  const Section main_loop      = Section{1, "main loop"};
  const Section pic_iteration  = Section{2, "pic_iteration"};
  const Section diags          = Section{3, "all diags"};

  const Section interpolate          = Section{4, "interpolate"};
  const Section push                 = Section{5, "push"};
  const Section pushBC               = Section{6, "pushBC"};
  const Section id_parts_to_move     = Section{7, "id_parts_to_move"};
  const Section exchange             = Section{8, "exchange"};
  const Section reset_current        = Section{9, "reset_current"};
  const Section projection           = Section{10, "projection"};
  const Section current_local_reduc  = Section{11, "current_local_reduc"};
  const Section current_global_reduc = Section{12, "current_global_reduc"};
  const Section currentBC            = Section{13, "currentBC"};
  const Section maxwell_solver       = Section{14, "maxwell_solver"};
  const Section maxwellBC            = Section{15, "maxwellBC"};
  const Section diags_sync           = Section{16, "diags_sync"};
  const Section diags_binning        = Section{17, "diags_binning"};
  const Section diags_cloud          = Section{18, "diags_cloud"};
  const Section diags_scalar         = Section{19, "diags_scalar"};
  const Section diags_field          = Section{20, "diags_field"};
  const Section imbalance            = Section{21, "imbalance"};

  // Vector of timers
  std::vector<Section> sections = {initialization,
                                   main_loop,
                                   pic_iteration,
                                   diags,
                                   interpolate,
                                   push,
                                   pushBC,
                                   id_parts_to_move,
                                   exchange,
                                   reset_current,
                                   projection,
                                   current_local_reduc,
                                   current_global_reduc,
                                   currentBC,
                                   maxwell_solver,
                                   maxwellBC,
                                   diags_sync,
                                   diags_binning,
                                   diags_cloud,
                                   diags_scalar,
                                   diags_field,
                                   imbalance};

  // String Buffer to store the timers
  std::stringstream timers_buffer;

  // _______________________________________________________________
  //
  //! Constructor
  // _______________________________________________________________
  Struc_Timers(Params params) {
    // By default there are 2 timers :
    // - initilization
    // - main loop

    N_patches = params.N_patches;

    temporary_times.resize(sections.size() * (1 + N_patches));
    accumulated_times.resize(sections.size() * (1 + N_patches));

    // initialize timers
    for (int i = 0; i < sections.size() * (1 + N_patches); i++) {
      accumulated_times[i] = 0;
      temporary_times[i]   = T_clock::now();
    }

    // Create a new timers.json file
    std::ofstream file;
    file.open("timers.json");

    // Add parameters
    file << "{\n";
    file << "  \"parameters\" : {\n";
    file << "    \"number_of_patches\" : " << N_patches << ",\n";
    file << "    \"number_of_threads\" : " << params.number_of_threads << ",\n";
    file << "    \"iterations\" : " << params.n_it << ",\n";
    file << "    \"save_timers_period\" : " << params.save_timers_period << ",\n";
    file << "    \"save_timers_start\" : " << params.save_timers_start << "\n";
    file << "  },\n";

    // close
    file.close();
  }

  //! Destructor
  ~Struc_Timers(){};

  // _______________________________________________________________
  //
  //! Start a global timer
  // _______________________________________________________________
  void start(Section section) {

    // auto time = std::chrono::high_resolution_clock::now();

    auto index = first_index(section);

    temporary_times[index] = T_clock::now();

    // std::cout << "Start timer " << section.name << " at " <<
    // temporary_times[index].time_since_epoch().count() << std::endl;
  }

  // _______________________________________________________________
  //
  //! Stop a timer
  // _______________________________________________________________
  void stop(Section section) {

    auto index = first_index(section);

    // auto time = std::chrono::high_resolution_clock::now();
    auto time = T_clock::now();

    std::chrono::duration<double> diff = time - temporary_times[index];
    accumulated_times[index] += diff.count();

    // std::cout << "Stop timer " << section.name << " at " << time.time_since_epoch().count()
    //           << " with diff " << diff.count()
    //           << std::endl;
  }

  // _______________________________________________________________
  //
  //! Start a thread timer
  // _______________________________________________________________
  void start(Section section, int i_patch) {

    auto index = first_index(section) + i_patch + 1;

    // auto time = std::chrono::high_resolution_clock::now();
    auto time = T_clock::now();

    temporary_times[index] = T_clock::now();
  }

  // _______________________________________________________________
  //
  //! Stop a timer
  // _______________________________________________________________
  void stop(Section section, int i_patch) {

    auto index = first_index(section) + i_patch + 1;

    // auto time = std::chrono::high_resolution_clock::now();
    auto time = T_clock::now();

    std::chrono::duration<double> diff = time - temporary_times[index];
    accumulated_times[index] += diff.count();
  }

  // _______________________________________________________________
  //
  //! Get the elapsed time since the beginning of the time loop
  // _______________________________________________________________
  double get_elapsed_time() {
    std::chrono::duration<double> diff =
      // std::chrono::high_resolution_clock::now() - temporary_times[1];
      T_clock::now() - temporary_times[1];
    return diff.count();
  }

  // _______________________________________________________________
  //
  //! \brief Print all timers
  //! \param  params : parameters of the simulation
  // _______________________________________________________________
  void print(Params params) {
    double percentage;

    const double initialization_time = accumulated_times[first_index(initialization)];
    const double main_loop_time      = accumulated_times[first_index(main_loop)];
    const double diags_time          = accumulated_times[first_index(diags)];
    const auto pic_iteration_time    = accumulated_times[first_index(pic_iteration)];

    printf(" ---------------------------------------------- |\n");
    printf(" Large timers                                   |\n");
    printf(" ---------------------------------------------- |\n");
    printf("            code part |  time (s)  | percentage |\n");
    printf(" ---------------------|------------|----------- |\n");

    percentage = initialization_time / (initialization_time + main_loop_time) * 100;
    printf("%21s |%11.6lf |%9.2lf %% |\n",
           initialization.name.c_str(),
           initialization_time,
           percentage);

    percentage = main_loop_time / (initialization_time + main_loop_time) * 100;
    printf("%21s |%11.6lf |%9.2lf %% |\n", main_loop.name.c_str(), main_loop_time, percentage);

    percentage = pic_iteration_time / (initialization_time + main_loop_time) * 100;
    const string comp_without_diags_name = "PIC iter (w/o diags)";
    printf("%21s |%11.6lf |%9.2lf %% |\n",
           comp_without_diags_name.c_str(),
           pic_iteration_time,
           percentage);

    percentage = diags_time / (initialization_time + main_loop_time) * 100;
    printf("%21s |%11.6lf |%9.2lf %% |\n", diags.name.c_str(), diags_time, percentage);

    printf(" ---------------------------------------------- |\n");
    printf(" Detailed timers                                |\n");
    printf(" ---------------------------------------------- |\n");
    printf("            code part |  time (s)  | percentage |\n");
    printf(" ---------------------|------------|----------- |\n");

    mini_float coverage = 0;

    for (int itimer = 3; itimer < sections.size(); itimer++) {

      double global_time = accumulated_times[first_index(sections[itimer])];
      double thread_time = 0;

      for (unsigned int i_patch = 0; i_patch < params.N_patches; i_patch++) {
        const int idx = first_index(sections[itimer]) + i_patch + 1;
        thread_time += accumulated_times[idx];
      }

      thread_time /= params.number_of_threads;

      const double total_time = global_time + thread_time;

      coverage += total_time;

      percentage = total_time / (main_loop_time) * 100;
      printf("%21s |%11.6lf |%9.2lf %% |\n", sections[itimer].name.c_str(), total_time, percentage);
    }

    printf(" ------------------------------------ \n");
    printf(" Total coverage : %9.2lf %% \n", coverage / (main_loop_time) * 100);
  }

  // _______________________________________________________________
  //
  //! \brief Save the initialization time in the timers file
  //
  // _______________________________________________________________
  void save_initialization(Params params) {
    const double initialization_time = accumulated_times[first_index(initialization)];
    std::ofstream file;
    file.open("timers.json", std::ios::app);
    file << "  \"initialization\" : " << initialization_time << ",\n";
    file.close();
  }

  // _______________________________________________________________
  //
  //! \brief  Write all timers in a file using json format
  //
  //! \details The file is named "timers.json" and has the following format :
  //! {
  //!   "parameters" : {
  //!       "number_of_patches" : N_patches,
  //!       "number_of_threads" : number_of_threads
  //!   },
  //!   "initilization" : [global],
  //!   "0" : {
  //!       "pic iterations" : [global],
  //!       "diags" : [global],
  //!       "interpolate" : [global, thread1, thread2, ...],
  //!       "push" : [global, thread1, thread2, ...],
  //!     ....
  //!   },
  //!   "10" : {
  //!     "pic iterations" : [global],
  //!     "diags" : [global],
  //!     "interpolate" : [global, thread1, thread2, ...],
  //!     "push" : [global, thread1, thread2, ...],
  //!     ....
  //!   }
  //!   "final" : {
  //!     "pic iterations" : [global],
  //!     "diags" : [global],
  //!     "interpolate" : [global, thread1, thread2, ...],
  //!     "push" : [global, thread1, thread2, ...],
  //!     ....
  //!   }
  //!   "main loop" : [global]
  //! ...
  //! }
  //! Use the scientific format with 6 digits after the comma
  //
  //! \param  params : parameters of the simulation
  //! \param  iteration : current iteration
  // _______________________________________________________________
  void save(Params params, int iteration) {

    if (iteration < params.save_timers_start)
      return;

    // iteration since the start of the saving
    unsigned int timer_iteration = iteration - params.save_timers_start;

    // save the timers every save_timers_period iterations
    // or if the simulation is finished (iteration > n_it)
    if (!(timer_iteration % params.save_timers_period) || (iteration > params.n_it)) {

      const double initialization_time = accumulated_times[first_index(initialization)];
      const double main_loop_time      = accumulated_times[first_index(main_loop)];
      const double diags_time          = accumulated_times[first_index(diags)];
      const auto pic_iteration_time    = accumulated_times[first_index(pic_iteration)];

      std::stringstream local_buffer("");

      local_buffer << std::scientific;
      if (iteration <= params.n_it) {
        local_buffer << "  \"" << iteration << "\" : {\n";
      } else {
        local_buffer << "  \"final\" : {\n";
      }
      local_buffer << "    \"pic iteration\" : " << pic_iteration_time << ",\n";
      local_buffer << "    \"diags\" : " << diags_time << ",\n";

      for (int itimer = 3; itimer < sections.size(); itimer++) {
        local_buffer << "    \"" << sections[itimer].name << "\" : [";
        for (int i = 0; i < N_patches + 1; i++) {
          local_buffer << accumulated_times[first_index(sections[itimer]) + i];
          if (i < N_patches) {
            local_buffer << ", ";
          }
        }
        if (itimer < sections.size() - 1) {
          local_buffer << "],\n";
        } else {
          local_buffer << "]\n";
        }
      }
      local_buffer << "  },\n";

      if (iteration > params.n_it) {
        local_buffer << "  \"main loop\" : " << main_loop_time << "\n";
        local_buffer << "}\n";
      }

      // std::cout << " -> Save timers at iteration " << iteration << std::endl;

      if (params.bufferize_timers_output) {

        timers_buffer << local_buffer.str();

        if (iteration > params.n_it) {
          std::ofstream file;
          file.open("timers.json", std::ios::app);
          file << timers_buffer.str();
          file.close();
        }

      } else {
        std::ofstream file;
        file.open("timers.json", std::ios::app);
        file << local_buffer.str();
        file.close();
      }

    } // end if save_timers_period
  }   // end save

private:
  // Array to store the timers
  std::vector<double> accumulated_times;

  // Array to store temporary values
  std::vector<std::chrono::time_point<T_clock>> temporary_times;

  // Number of patches
  int N_patches;

  //! Get first index from section id
  int first_index(Section section) { return section.id * (N_patches + 1); }

  //! Select a steady clock
  // static constexpr auto clock() {
  //   if constexpr (std::chrono::high_resolution_clock::is_steady) {
  //     return std::chrono::high_resolution_clock();
  //   } else {
  //     return std::chrono::steady_clock();
  //   }
  // }
};

// Timers shortcut
// using Timers = Struc_Timers<std::chrono::high_resolution_clock>;
using Timers = Struc_Timers<std::chrono::steady_clock>;
