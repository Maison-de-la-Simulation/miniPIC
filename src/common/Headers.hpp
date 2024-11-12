/* _____________________________________________________________________ */
//! \file Backend.hpp

//! \brief determine the best backend to use

/* _____________________________________________________________________ */

#ifndef HEADERS_H
#define HEADERS_H

// #include "Params.hpp"

// _____________________________________________________________________
//
// Backends
// _____________________________________________________________________

// ____________________________________________________________
// OMP and OMP task

#if defined(__MINIPIC_OMP__) || defined(__MINIPIC_OMP_TASK__)

#include "omp.h"
#include <atomic>
#include <deque>
#include <memory>
#include <vector>

#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE inline __attribute__((always_inline))

#elif defined(__MINIPIC_EVENTIFY__)

#include "omp.h"
#include <atomic>
#include <eventify/task_system.hxx>
#include <jsc/event_counter.hpp>
#include <memory>
#include <vector>

#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE inline __attribute__((always_inline))

#else

#include <memory>
#include <vector>
#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE inline __attribute__((always_inline))

#endif

// _____________________________________________________________________
// Types

#if defined(__SHAMAN__)

#include <shaman.h>

using mini_float = Sdouble;
using namespace Sstd;

#else

// using mini_float = double;
#define mini_float double

using namespace std;

#endif

// _____________________________________________________________________
// Space class

namespace minipic {

class Host {
public:
  static const int value = 1;
};

class Device {
public:
  static const int value = 2;
};

const Host host;
const Device device;

template <typename T> inline void atomicAdd(T *address, T value) {
#if defined(__MINIPIC_OMP__) || defined(__MINIPIC_OMP_TASK__) || defined(__MINIPIC_OMP_TARGET__)
#pragma omp atomic update
  *address += value;
#else
  *address += value;
#endif
}

} // namespace minipic

// onHost  on_host;
// onDevice on_device;

#endif
