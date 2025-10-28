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

// ____________________________________________________________
// Kokkos

#elif defined(__MINIPIC_KOKKOS_COMMON__)

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_ScatterView.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE KOKKOS_INLINE_FUNCTION

// ____________________________________________________________
// Thrust and OMP target

#elif defined(__MINIPIC_THRUST_COMMON__)

#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/host_vector.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>

#if defined(__MINIPIC_THRUST_UNIFIED__)
#include <thrust/universal_vector.h>
#endif

#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE inline __attribute__((always_inline)) __device__

#elif defined(__MINIPIC_HIPTHRUST__)

#include <hip/thrust/device_vector.h>
#include <hip/thrust/host_vector.h>
#include <hip/thrust/transform.h>
#include <hip/thrust/transform_reduce.h>

#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE inline __attribute__((always_inline)) __device__

#elif defined(__MINIPIC_OMP_TARGET__)

#include "omp.h"
#include <memory>
#include <vector>
#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE inline __attribute__((always_inline))

#elif defined(__MINIPIC_OPENACC__)

#include "openacc.h"
#include <memory>
#include <vector>
#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE inline __attribute__((always_inline))

#elif defined(__MINIPIC_SYCL__)

#include <sycl/sycl.hpp>
// #include <sycl/ext/intel/fpga_extensions.hpp>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>

#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE inline __attribute__((always_inline))

#elif defined(__MINIPIC_STDPAR__) || (__MINIPIC_STDPAR_CPU__)
#include <algorithm>
#include <atomic>
#include <cmath>
#include <execution>
#include <memory>
#include <vector>

// #include <cuda/std/atomic> //utiliser la fonction fetch_add dans project de la librairie libcu++
// sur GPU au lieu de addAtomic

#define INLINE inline __attribute__((always_inline))
#define DEVICE_INLINE inline __attribute__((always_inline))

template <typename T>
struct counting_iterator {
private:
  using self = counting_iterator;

public:
  using value_type        = T;
  using difference_type   = typename std::make_signed<T>::type;

  counting_iterator() : value(0) {}
  explicit counting_iterator(value_type v) : value(v) {}

  value_type operator*() const { return value; }
  value_type operator[](difference_type n) const { return value + n; }

  self &operator++() {
    ++value;
    return *this;
  }
  self operator++(int) {
    self result{value};
    ++value;
    return result;
  }
  self &operator--() {
    --value;
    return *this;
  }
  self operator--(int) {
    self result{value};
    --value;
    return result;
  }
  self &operator+=(difference_type n) {
    value += n;
    return *this;
  }
  self &operator-=(difference_type n) {
    value -= n;
    return *this;
  }

  friend self operator+(self const &i, difference_type n) { return self(i.value + n); }
  friend self operator+(difference_type n, self const &i) { return self(i.value + n); }
  friend difference_type operator-(self const &x, self const &y) { return x.value - y.value; }
  friend self operator-(self const &i, difference_type n) { return self(i.value - n); }

  friend bool operator==(self const &x, self const &y) { return x.value == y.value; }
  friend bool operator!=(self const &x, self const &y) { return x.value != y.value; }
  friend bool operator<(self const &x, self const &y) { return x.value < y.value; }
  friend bool operator<=(self const &x, self const &y) { return x.value <= y.value; }
  friend bool operator>(self const &x, self const &y) { return x.value > y.value; }
  friend bool operator>=(self const &x, self const &y) { return x.value >= y.value; }

private:
  value_type value;
};

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

#if defined(__NVIDIA_PROFILER__)
#include <nvToolsExt.h>
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
#elif defined(__MINIPIC_OPENACC__)
#pragma acc atomic update
  *address += value;
#else
  *address += value;
#endif
}

} // namespace minipic

// onHost  on_host;
// onDevice on_device;

#endif
