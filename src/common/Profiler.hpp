
/* _____________________________________________________________________ */
//! \file Profiler.cpp

//! \brief Profiler class: detailed timing of the code

/* _____________________________________________________________________ */

#pragma once

#include <vector>

#include "Params.hpp"
#include "Tools.hpp"

// Code for each occurrence's type
#define RESET 0
#define EVOLVE_BIN 1
#define EVOLVE_PATCH 2
#define EVOLVE_PATCH 2
#define EXCHANGE 3
#define CURRENT_GLOBAL_INTERNAL 4
#define CURRENT_GLOBAL_BORDERS 5
#define CURRENTBC 6
#define MAXWELL 7
#define DIAGS 8

#define INTERPOLATE 9
#define PUSH 10
#define PUSHBC 11
#define PROJECT 12
#define DIAGNOSTICS 13

//! Structure to store an event with a start, stop and an occurence description
struct Event {
  double start_;
  double stop_;
  int occurrence_code_;
};

//! Profiler class: detailed timing of the code
class Profiler {

private:
  //! thread_id_ is the id of the thread
  int number_of_threads_;

  //! events_ is a vector of Event
  //! Each event is stored in the order of its occurrence
  std::vector<std::vector<Event>> events_;

  //! Reference time
  std::chrono::time_point<std::chrono::high_resolution_clock> reference_time_;

public:
  //! Constructor
  Profiler(Params &params) {
#if defined(__MINIPIC_PROFILER__)
    // Get the number of threads to initialize the events_ vector
    number_of_threads_ = params.number_of_threads;

    events_.resize(number_of_threads_);
    for (auto ithread = 0; ithread < number_of_threads_; ++ithread) {
      events_[ithread].reserve(params.n_it * 4);
    }

    // Get a common reference time
    reference_time_ = std::chrono::high_resolution_clock::now();
#endif
  }

  //! Destructor
  ~Profiler(){};

  //! Start to time the occurrence of code `occurrence_code`
  //! Store the time in seconds since the epoch
  //! \param occurrence_code is the code of the occurrence
  void start(unsigned int occurrence_code) {
#if defined(__MINIPIC_PROFILER__)
    const int thread_id = get_thread_id();
    events_[thread_id].push_back(Event());
    std::chrono::duration<double> diff =
      std::chrono::high_resolution_clock::now() - reference_time_;
    events_[thread_id].back().start_           = diff.count();
    events_[thread_id].back().occurrence_code_ = occurrence_code;
#endif
  }

  //! Stop to the current occurrence
  void stop() {
#if defined(__MINIPIC_PROFILER__)
    const int thread_id = get_thread_id();
    std::chrono::duration<double> diff =
      std::chrono::high_resolution_clock::now() - reference_time_;
    events_[thread_id].back().stop_ = diff.count();
#endif
  }

  //! Output the list of events on disk in a binary format
  //! The file starts by the number of threads
  //! Then each event is stored in the order of its occurrence, starting with the number of events
  void dump(Params &params) {
#if defined(__MINIPIC_PROFILER__)

    std::string file_name(params.output_directory + "/profiler.bin");

    std::ofstream file(file_name, std::ios::binary);
    file.write(reinterpret_cast<char *>(&number_of_threads_), sizeof(int));
    for (auto ithread = 0; ithread < number_of_threads_; ithread++) {
      auto size = events_[ithread].size();
      file.write(reinterpret_cast<char *>(&size), sizeof(int));
      for (auto even : events_[ithread]) {
        file.write(reinterpret_cast<char *>(&even.start_), sizeof(double));
        file.write(reinterpret_cast<char *>(&even.stop_), sizeof(double));
        file.write(reinterpret_cast<char *>(&even.occurrence_code_), sizeof(int));
      }
    }
    file.close();

#endif
  }
};