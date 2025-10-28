
/* _____________________________________________________________________ */
//! \file vector.h

//! \brief Minipic vector class for backend abstraction

/* _____________________________________________________________________ */

// #pragma once
#ifndef VECTOR_H
#define VECTOR_H

#include "Backend.hpp"
#include "Headers.hpp"

// ______________________________________________________________________
//
//! \brief Class Vector for MiniPIC
// ______________________________________________________________________
template <typename T> class Vector {

public:
  // Number of elements
  size_t size_;

  // Main data
#if defined(__MINIPIC_KOKKOS__)
  Kokkos::DualView<T *> data_;
#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
  Kokkos::DualView<T *, Kokkos::SharedSpace> data_;
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
  Kokkos::View<T *, Kokkos::SharedSpace> data_;
#elif defined(__MINIPIC_THRUST__)
  // GPU data vector
  thrust::device_vector<T> device_data_;
  // CPU data vector
  thrust::host_vector<T> host_data_;
  // std::vector<T> host_data_;
#elif defined(__MINIPIC_THRUST_UNIFIED__)
  thrust::universal_vector<T> data_;
#elif defined(__MINIPIC_SYCL__)
  T *host_data_;
  T *device_data_;
  // pointer to the sycl queue
  sycl::queue *sycl_queue_ptr;
#elif defined(__MINIPIC_STDPAR__)
  std::vector<T> data_;
#else
  std::vector<T> data_;
#endif

  // ______________________________________________________________________
  //
  //! \brief constructors
  // ______________________________________________________________________
  Vector() : size_(0) {}
  Vector(size_t size, Backend &backend) { allocate("", backend, size); }

  // ______________________________________________________________________
  //
  //! \brief constructor with allocation
  //! \param[in] size allocation size
  //! \param[in] v default value
  // ______________________________________________________________________
  Vector(size_t size, T v, Backend &backend) {
    allocate("", size, backend);
    fill(v, minipic::host);
    fill(v, minipic::device);
  }

  // ______________________________________________________________________
  //
  //! \brief destructor
  // ______________________________________________________________________
  ~Vector() {
#if defined(__MINIPIC_OPENMP_TARGET__)
    // remove data on device
    T *ptr = get_raw_pointer(minipic::host);
#pragma omp target exit data map(delete : ptr[0 : size_])
#elif defined(__MINIPIC_OPENACC__)
    // remove data on device
    T *ptr = get_raw_pointer(minipic::host);
#pragma acc exit data delete (ptr[0 : size_])
#elif defined(__MINIPIC_SYCL__)
    sycl::free(host_data_, *sycl_queue_ptr);
    sycl::free(device_data_, *sycl_queue_ptr);
#endif
  };

  // ______________________________________________________________________
  //
  //! \brief allocate the data_ object
  // ______________________________________________________________________
  void allocate(std::string name, size_t size, Backend &backend) {
    size_ = size;
#if defined(__MINIPIC_KOKKOS__)
    data_ = Kokkos::DualView<T *>(name, size);
#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
    data_ = Kokkos::DualView<T *, Kokkos::SharedSpace>(name, size);
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    data_ = Kokkos::View<T *, Kokkos::SharedSpace>(name, size);
#elif defined(__MINIPIC_THRUST__)
    device_data_.resize(size);
    host_data_.resize(size);
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    data_.resize(size);
#elif defined(__MINIPIC_OMP_TARGET__)
    data_.resize(size);
    T *ptr = get_raw_pointer(minipic::host);
#pragma omp target enter data map(to : ptr[0 : size])
#elif defined(__MINIPIC_OPENACC__)
    data_.resize(size);
    T *ptr = get_raw_pointer(minipic::host);
#pragma acc enter data create(ptr[0 : size])
#elif defined(__MINIPIC_SYCL__)
    sycl_queue_ptr = backend.sycl_queue_;
    host_data_     = sycl::malloc_host<T>(size, *sycl_queue_ptr);
    device_data_   = sycl::malloc_device<T>(size, *sycl_queue_ptr);
#elif defined(__MINIPIC_STDPAR__)
      data_.resize(size); // since for now there is only one patch, resize is not an issue
#else
    data_.resize(size);
#endif
  }

  // ______________________________________________________________________
  //
  //! \brief [] operator
  //! \return Host data accessor (if not device, point to the host data)
  // ______________________________________________________________________
  INLINE T &operator[](const size_t i) {
#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)
    return data_.h_view(i);
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    return data_(i);
#elif defined(__MINIPIC_THRUST__)
    return host_data_[i];
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    return data_[i];
#elif defined(__MINIPIC_SYCL__)
    return host_data_[i];
#elif defined(__MINIPIC_STDPAR__)
    return data_[i];
#else
    return data_[i];
#endif
  }

  // ______________________________________________________________________
  //
  //! \brief () operator
  //! \return Host data accessor (if not device, point to the host data)
  // ______________________________________________________________________
  INLINE T &operator()(const size_t i) {
#if defined(__MINIPIC_KOKKOS__) || defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
    return data_.h_view(i);
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    return data_(i);
#elif defined(__MINIPIC_THRUST__)
    return host_data_[i];
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    return data_[i];
#elif defined(__MINIPIC_SYCL__)
    return host_data_[i];
#elif defined(__MINIPIC_STDPAR__)
    return data_[i];
#else
    return data_[i];
#endif
  }

  // ______________________________________________________________________
  //
  //! \brief Explicit Device data accessor (if not device, point to the host data)
  //! \return device pointer at index i
  // ______________________________________________________________________
  //   DEVICE_INLINE T &d(const size_t i) {
  // #if defined(__MINIPIC_KOKKOS__)
  //     return data_.d_view(i);
  // #else
  //     return data_[i];
  // #endif
  //   }

  // ______________________________________________________________________
  //
  //! \brief Explicit Host data accessor
  //! \return host pointer at index i
  // ______________________________________________________________________
  INLINE T &h(const size_t i) {
#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)
    return data_.h_view(i);
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    return data_(i);
#elif defined(__MINIPIC_THRUST__)
    return host_data_[i];
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    return data_[i];
#elif defined(__MINIPIC_SYCL__)
    return host_data_[i];
#elif defined(__MINIPIC_STDPAR__)
    return data_[i];
#else
    return data_[i];
#endif
  }

  // ______________________________________________________________________
  //
  //! \brief Get the data pointer
  //! \param[in] space where to keep the data when resizing (must be minipic::host or
  //! minipic::device)
  // ______________________________________________________________________
  template <class T_space> T *get_raw_pointer(const T_space space) {

    // Check that T_Space of Class Host or Device
    static_assert(std::is_same<T_space, minipic::Host>::value ||
                    std::is_same<T_space, minipic::Device>::value,
                  "Must be minipic::host or minipic::device");

    // Host
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)
      return data_.h_view.data();
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
      return data_.data();
#elif defined(__MINIPIC_THRUST__)
      return host_data_.data();
#elif defined(__MINIPIC_THRUST_UNIFIED__)
      return data_.data();
#elif defined(__MINIPIC_SYCL__)
      return host_data_;
#elif defined(__MINIPIC_STDPAR__)
      return data_.data();
#else
      return data_.data();
#endif

      // Device
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)
      return data_.d_view.data();
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
      return data_.data();
#elif defined(__MINIPIC_THRUST__)
      return thrust::raw_pointer_cast(device_data_.data());
#elif defined(__MINIPIC_THRUST_UNIFIED__)
      return thrust::raw_pointer_cast(data_.data());
#elif defined(__MINIPIC_SYCL__)
      return device_data_;
#elif defined(__MINIPIC_STDPAR__)
      return data_.data();
#else
      return data_.data();
#endif

    } else {
      return nullptr;
    }
  }

  // ______________________________________________________________________
  //
  //! \brief return the size
  //! \return size of the vector
  // ______________________________________________________________________
  INLINE T size() { return size_; }

  // ______________________________________________________________________
  //
  //! \brief resize the vector to the new size
  //! \param[in] new_size new vector size
  //! \param[in] space where to keep the data when resizing (must be minipic::host or
  //! minipic::device)
  // ______________________________________________________________________
  template <class T_space> void resize(const size_t new_size, const T_space space) {

    // Check that T_Space of Class Host or Device
    static_assert(std::is_same<T_space, minipic::Host>::value ||
                    std::is_same<T_space, minipic::Device>::value,
                  "Must be minipic::host or minipic::device");

#if defined(__MINIPIC_KOKKOS__) || defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
    // Check that the space is either host or device
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
      data_.modify_host();
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
      data_.modify_device();
    }
    data_.resize(new_size);
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    Kokkos::resize(data_, new_size);
#elif defined(__MINIPIC_THRUST__)
    host_data_.resize(new_size);
    device_data_.resize(new_size);
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    data_.resize(new_size);
#elif defined(__MINIPIC_OMP_TARGET__)
    // bring back data on the host
    sync(minipic::device, minipic::host);
    // Remove data on device
    T *ptr = get_raw_pointer(minipic::host);
#pragma omp target exit data map(delete : ptr[0 : size_])
    data_.resize(new_size);
    ptr = get_raw_pointer(minipic::host);
#pragma omp target enter data map(to : ptr[0 : new_size])
#elif defined(__MINIPIC_OPENACC__)
    // bring back data on the host
    sync(minipic::device, minipic::host);
    // Remove data on device
    T *ptr = get_raw_pointer(minipic::host);
#pragma acc exit data delete (ptr[0 : size_])
    data_.resize(new_size);
    ptr = get_raw_pointer(minipic::host);
#pragma acc enter data create(ptr[0 : new_size])
#elif defined(__MINIPIC_SYCL__)
    // bring back data on the host
    sync(minipic::device, minipic::host);
    // Remove data on device
    sycl::free(device_data_, *sycl_queue_ptr);
    // Allocate new data on device
    device_data_ = sycl::malloc_device<T>(new_size, *sycl_queue_ptr);
    // Allocate new data on host
    sycl::free(host_data_, *sycl_queue_ptr);
    host_data_ = sycl::malloc_host<T>(new_size, *sycl_queue_ptr);
#elif defined(__MINIPIC_STDPAR__)
    data_.resize(new_size);
#else
    data_.resize(new_size);
#endif
    size_ = new_size;
  }

  // ______________________________________________________________________
  //
  //! \brief resize the vector to the new size
  //! \param[in] new_size new vector size
  //! \param[in] value value used to initialize the new elements
  //! \param[in] space where to preserve data (minipic::host or minipic::device)
  //! \tparam T_space class of the space
  // ______________________________________________________________________
  template <class T_space> void resize(const size_t new_size, T value, const T_space space) {

      resize(new_size, space); // calling the previous overload

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)
    for (auto ip = size_; ip < new_size; ++ip) {
      data_.h_view(ip) = value;
    }
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    for (auto ip = size_; ip < new_size; ++ip) {
      data_(ip) = value;
    }
#elif defined(__MINIPIC_THRUST__)
    for (auto ip = size_; ip < new_size; ++ip) {
      host_data_[ip] = value;
    }
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    for (auto ip = size_; ip < new_size; ++ip) {
      data_[ip] = value;
    }
#elif defined(__MINIPIC_SYCL__)
    for (auto ip = size_; ip < new_size; ++ip) {
      host_data_[ip] = value;
    }
#elif defined(__MINIPIC_STDPAR__)
      // initialize to the correct value
    for (auto ip = size_; ip < new_size; ++ip) {
      data_[ip] = value;
    }
#else
    for (auto ip = size_; ip < new_size; ++ip) {
      data_[ip] = value;
    }
#endif
  }

  // ______________________________________________________________________
  //
  //! \brief clear the content, equivalent to size_ = 0
  //! If the raw object has a clear method, we call it
  // ______________________________________________________________________
  void clear() {
    size_ = 0;
#if defined(__MINIPIC_KOKKOS_COMMON__)
    // nothing
#elif defined(__MINIPIC_THRUST__)
    host_data_.clear();
    device_data_.clear();
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    data_.clear();
#elif defined(__MINIPIC_OMP_TARGET__)
    // remove data on device
    T *ptr = get_raw_pointer(minipic::host);
#pragma omp target exit data map(delete : ptr[0 : size_])
    // then on host
    data_.clear();
#elif defined(__MINIPIC_OPENACC__)
    // remove data on device
    T *ptr = get_raw_pointer(minipic::host);
#pragma acc exit data delete (ptr[0 : size_])
    // then on host
    data_.clear();
#elif defined(__MINIPIC_SYCL__)
    // We just put the size to 0 (data still allocated)
#elif defined(__MINIPIC_STDPAR__)
    data_.clear();
#else
    data_.clear();
#endif
  }

  // ______________________________________________________________________
  //
  //! \brief fill the vector with the given value
  // ______________________________________________________________________
  template <class T_space> void fill(const T v, const T_space space) {
    // Check that T_Space of Class Host or Device
    static_assert(std::is_same<T_space, minipic::Host>::value ||
                    std::is_same<T_space, minipic::Device>::value,
                  "Must be minipic::host or minipic::device");

    // Host
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
#if defined(__MINIPIC_KOKKOS__)
      // Fill on host
      for (auto i = 0; i < size_; ++i) {
        data_.h_view(i) = v;
      }

      // quest
      // Boucle non async
      // Kokkos::fence();

#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
      // quest
      //      // Fill only on device (uvm)
      //      Kokkos::parallel_for(
      //        size_,
      //        KOKKOS_CLASS_LAMBDA(const size_t ip) { data_.d_view(ip) = v; });

      //     Kokkos::fence();

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
      Kokkos::Experimental::fill(Kokkos::DefaultHostExecutionSpace(),
                                 Kokkos::Experimental::begin(data_),
                                 Kokkos::Experimental::end(data_),
                                 v);
      Kokkos::fence();

#elif defined(__MINIPIC_THRUST__)
      thrust::fill(host_data_.begin(), host_data_.end(), v);

#elif defined(__MINIPIC_THRUST_UNIFIED__)
      for (size_t i = 0; i < size_; ++i) {
        data_[i] = v;
      }

#elif defined(__MINIPIC_OMP_TARGET__)
      std::fill(data_.begin(), data_.end(), v);

#elif defined(__MINIPIC_OPENACC__)
      std::fill(data_.begin(), data_.end(), v);

#elif defined(__MINIPIC_SYCL__)
      // fill on host
      for (auto i = 0; i < size_; ++i) {
        host_data_[i] = v;
      }

      sycl_queue_ptr->wait();

#elif defined(__MINIPIC_STDPAR__)
      std::fill(std::execution::seq, data_.begin(), data_.end(), v);
#elif defined(__MINIPIC_CUDA__)
      // fill on host
      for (auto i = 0; i < size_; ++i) {
        host_data_[i] = v;
      }

#else
      std::fill(data_.begin(), data_.end(), v);
#endif

    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {

#if defined(__MINIPIC_KOKKOS__)
      // Kokkos::Experimental::fill(Kokkos::DefaultHostExecutionSpace(), data_m.h_view, 0.);

      // Fill on device
      Kokkos::parallel_for(size_, KOKKOS_CLASS_LAMBDA(const size_t ip) { data_.d_view(ip) = v; });

      Kokkos::fence();

#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)

      // Fill only on device (uvm)
      Kokkos::parallel_for(size_, KOKKOS_CLASS_LAMBDA(const size_t ip) { data_.d_view(ip) = v; });

      Kokkos::fence();

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

      Kokkos::Experimental::fill(Kokkos::DefaultExecutionSpace(),
                                 Kokkos::Experimental::begin(data_),
                                 Kokkos::Experimental::end(data_),
                                 v);

      // Kokkos::Experimental::fill(Kokkos::DefaultHostExecutionSpace(),
      // Kokkos::Experimental::begin(data_), Kokkos::Experimental::end(data_), v);

      // // Fill on device
      // Kokkos::parallel_for(
      //   size_,
      //   KOKKOS_CLASS_LAMBDA(const size_t ip) { data_(ip) = v; });

      // // Fill on host
      // for (auto i = 0; i < size_; ++i) {
      //   data_(i) = v;
      // }

      Kokkos::fence();

#elif defined(__MINIPIC_THRUST__)
      thrust::fill(device_data_.begin(), device_data_.end(), v);

#elif defined(__MINIPIC_THRUST_UNIFIED__)
      thrust::fill(data_.begin(), data_.end(), v);

#elif defined(__MINIPIC_OMP_TARGET__)

      // fill on device
      T *ptr = get_raw_pointer(minipic::host);
#pragma omp target teams distribute parallel for
      for (auto ip = 0; ip < size_; ++ip) {
        ptr[ip] = v;
      }

#elif defined(__MINIPIC_OPENACC__)

      // fill on device
      T *ptr = get_raw_pointer(minipic::host);
#pragma acc parallel loop present(ptr[0 : size_])
      for (auto ip = 0; ip < size_; ++ip) {
        ptr[ip] = v;
      }

#elif defined(__MINIPIC_SYCL__)

      T *const device_data = device_data_;

      // fill on device
      sycl_queue_ptr->parallel_for(sycl::range{size_}, [=](sycl::id<1> i) { device_data[i] = v; });
      sycl_queue_ptr->wait();

#elif defined(__MINIPIC_STDPAR__)
      // void fill( ExecutionPolicy&& policy, ForwardIt first, ForwardIt last, const T& value );
      std::fill(std::execution::par_unseq, data_.begin(), data_.end(), v);
      // the CPU copy is not updated directly
#elif defined(__MINIPIC_CUDA__)
      // fill on device
      fillDeviceMemory(device_data_, v, size_);

#else
      std::fill(data_.begin(), data_.end(), v);
#endif

    } else {
      std::cerr << "Vector::sum: Invalid space" << std::endl;
    }
  }

  // _________________________________________________________
  //
  //! \brief sum of the vector
  //! \param[in] power power of the sum
  //! \param[in] space where to perform the reduction (host or device)
  // _________________________________________________________
  template <class T_space> T sum(const int power, T_space space) {
    T sum = 0;

    // ---> Host case
    if constexpr (std::is_same<T_space, minipic::Host>::value) {

#if defined(__MINIPIC_KOKKOS__) || defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)

      typedef Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace> range_policy;

      Kokkos::parallel_reduce(
        "sum",
        range_policy(0, size_),
        KOKKOS_CLASS_LAMBDA(const size_t i, T &lsum) { lsum += Kokkos::pow(data_.h_view(i), power); },
        sum);

      Kokkos::fence();

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

      typedef Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace> range_policy;

      Kokkos::parallel_reduce(
        "sum",
        range_policy(0, size_),
        KOKKOS_CLASS_LAMBDA(const size_t i, T &lsum) { lsum += Kokkos::pow(data_(i), power); },
        sum);

      Kokkos::fence();

#elif defined(__MINIPIC_THRUST__)

      for (size_t i = 0; i < size_; i++) {
        sum += pow(host_data_[i], power);
      }

#elif defined(__MINIPIC_THRUST_UNIFIED__)
      for (size_t i = 0; i < size_; i++) {
        sum += pow(data_[i], power);
      }

#elif defined(__MINIPIC_OMP_TARGET__) || defined(__MINIPIC_OPENACC__)
      for (size_t i = 0; i < size_; i++) {
        sum += pow(data_[i], power);
      }
#elif defined(__MINIPIC_SYCL__)
      for (size_t i = 0; i < size_; i++) {
        sum += sycl::pown(host_data_[i], power);
      }
#elif defined(__MINIPIC_STDPAR__)

      // T transform_reduce( InputIt first, InputIt last, T init, BinaryOp reduce, UnaryOp transform
      // ); Provided that data_ is not of type int. If data_ is int: see the stage_ester directory
      // stage_ester/mini-pic/mini-pic.md

      // sum = std::transform_reduce(data_.begin(),
      //                             data_.end(),
      //                             static_cast<T>(0),
      //                             std::plus<T>(),
      //                             [power](T x) { return std::pow(x, power); });

      for (size_t i = 0; i < size_; i++) {
        sum += pow(data_[i], power);
      }

#else
      for (size_t i = 0; i < size_; i++) {
        sum += pow(data_[i], power);
      }
#endif

      // ---> Device case
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {

#if defined(__MINIPIC_KOKKOS_COMMON__)

#if defined(__MINIPIC_KOKKOS__)
      typename Kokkos::DualView<T *>::t_dev view = data_.d_view;
#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
      typename Kokkos::DualView<T *, Kokkos::SharedSpace>::t_dev view = data_.d_view;
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
      typename Kokkos::View<T *, Kokkos::SharedSpace> view = data_;
#endif

      Kokkos::parallel_reduce(
        "sum",
        size_,
        KOKKOS_CLASS_LAMBDA(const size_t i, T &lsum) { lsum += Kokkos::pow(view(i), power); },
        sum);

      Kokkos::fence();

#elif defined(__MINIPIC_THRUST__)

      sum = thrust::transform_reduce(
        device_data_.begin(),
        device_data_.end(),
        [=] __host__ __device__(T x) { return pow(x, power); },
        static_cast<T>(0.0),
        thrust::plus<T>());

#elif defined(__MINIPIC_THRUST_UNIFIED__)

      sum = thrust::transform_reduce(
        data_.begin(),
        data_.end(),
        [=] __host__ __device__(T x) { return pow(x, power); },
        static_cast<T>(0.0),
        thrust::plus<T>());

#elif defined(__MINIPIC_OMP_TARGET__)
      T *ptr = get_raw_pointer(space);
#pragma acc parallel loop present(ptr[0 : size_]) reduction(+ : sum)
      for (size_t i = 0; i < size_; i++) {
        sum += pow(ptr[i], power);
      }

#elif defined(__MINIPIC_OPENACC__)
      T *ptr = get_raw_pointer(space);
#pragma omp target teams distribute parallel for reduction(+ : sum)
      for (size_t i = 0; i < size_; i++) {
        sum += pow(ptr[i], power);
      }

#elif defined(__MINIPIC_SYCL__)

      // buffer on device for sum
      sycl::buffer<T, 1> d_sum(&sum, 1);

      sycl::buffer<T, 1> data_buf(device_data_, sycl::range<1>(size_));

      sycl_queue_ptr->submit([&](sycl::handler &cgh) {
        sycl::accessor ksum{d_sum, cgh, sycl::write_only};
        sycl::accessor data_acc{data_buf, cgh, sycl::read_only};

        cgh.parallel_for(
          sycl::range<1>{size_},
          // Reduction object, to perform summation - initialises the result to zero
          sycl::reduction(d_sum,
                          cgh,
                          std::plus<T>(),
                          sycl::property::reduction::initialize_to_identity{}),
          [=](sycl::id<1> idx, auto &sum) { sum += sycl::pown(data_acc[idx], power); });
      });

      sycl::host_accessor result{d_sum, sycl::read_only};
      sum = result[0];

#elif defined(__MINIPIC_STDPAR__)
      // T transform_reduce( InputIt first, InputIt last, T init, BinaryOp reduce, UnaryOp transform
      // ); A condition que data_ ne soit pas de type size_t
      sum = std::transform_reduce(std::execution::par_unseq,
                                  data_.begin(),
                                  data_.end(),
                                  static_cast<T>(0),
                                  std::plus<T>(),
                                  [power](T x) { return std::pow(x, power); });

#else
      for (size_t i = 0; i < size_; i++) {
        sum += pow(data_[i], power);
      }
#endif

    } else {
      std::cerr << "Vector::sum: Invalid space" << std::endl;
    }

    return sum;
  }

  // _________________________________________________________
  //
  //! \brief sync host <-> device
  // _________________________________________________________
  template <class from, class to> void sync(const from, const to) {

    // Check the combination of from and to:
    // - from is minipic::Host then to is minipic::Device
    // - from is minipic::Device then to is minipic::Host
    static_assert(
      (std::is_same<from, minipic::Host>::value && std::is_same<to, minipic::Device>::value) ||
        (std::is_same<from, minipic::Device>::value && std::is_same<to, minipic::Host>::value),
      "Vector::sync: Invalid combination of from and to");

    // Host -> Device
    if constexpr (std::is_same<from, minipic::Host>::value &&
                  std::is_same<to, minipic::Device>::value) {

#if defined(__MINIPIC_KOKKOS__)
      data_.modify_host();
      data_.template sync<typename Kokkos::DualView<T *>::execution_space>();
#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
      // nothing (UVM)
#elif defined(__MINIPIC_THRUST__)
      thrust::copy(host_data_.begin(), host_data_.begin() + size_, device_data_.begin());
#elif defined(__MINIPIC_THRUST_UNIFIED__)
      // nothing (UVM)
#elif defined(__MINIPIC_OMP_TARGET__)
      T *ptr = get_raw_pointer(minipic::host);
#pragma omp target update to(ptr[0 : size_])
#elif defined(__MINIPIC_OPENACC__)
      T *ptr = get_raw_pointer(minipic::host);
#pragma acc update device(ptr[0 : size_])
#elif defined(__MINIPIC_SYCL__)
      sycl_queue_ptr->memcpy(device_data_, host_data_, size_ * sizeof(T));
      sycl_queue_ptr->wait();
#elif defined(__MINIPIC_STDPAR__)

#endif

      // Device -> Host
    } else if constexpr (std::is_same<from, minipic::Device>::value &&
                         std::is_same<to, minipic::Host>::value) {

#if defined(__MINIPIC_KOKKOS__)
      data_.modify_device();
      data_.template sync<typename Kokkos::DualView<T *>::host_mirror_space>();
#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
      // nothing (UVM)
#elif defined(__MINIPIC_THRUST__)
      thrust::copy(device_data_.begin(), device_data_.begin() + size_, host_data_.begin());
#elif defined(__MINIPIC_THRUST_UNIFIED__)
      // nothing (UVM)
#elif defined(__MINIPIC_OMP_TARGET__)
      T *ptr = get_raw_pointer(minipic::host);
#pragma omp target update from(ptr[0 : size_])
#elif defined(__MINIPIC_OPENACC__)
      T *ptr = get_raw_pointer(minipic::host);
#pragma acc update host(ptr[0 : size_])
#elif defined(__MINIPIC_SYCL__)
      sycl_queue_ptr->memcpy(host_data_, device_data_, size_ * sizeof(T));
      sycl_queue_ptr->wait();
#elif defined(__MINIPIC_STDPAR__)
#endif
    }
  }
};

// _________________________________________________________________________
// Shortcuts

#if defined(__MINIPIC_OMP__)

using vector_t        = Vector<mini_float>;
using device_vector_t = Vector<mini_float>;

#elif defined(__MINIPIC_KOKKOS__)

using vector_t        = Kokkos::DualView<mini_float *>::t_host;
using device_vector_t = Kokkos::DualView<mini_float *>::t_dev;

#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)

using vector_t        = Kokkos::DualView<mini_float *, Kokkos::SharedSpace>::t_host;
using device_vector_t = Kokkos::DualView<mini_float *, Kokkos::SharedSpace>::t_dev;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

using vector_t        = Kokkos::View<mini_float *, Kokkos::SharedSpace>;
using device_vector_t = Kokkos::View<mini_float *, Kokkos::SharedSpace>;

#elif defined(__MINIPIC_THRUST__)

using vector_t        = thrust::host_vector<mini_float>;
using device_vector_t = thrust::device_vector<mini_float>;
// using device_vector_t = thrust::host_vector<mini_float>;

#elif defined(__MINIPIC_THRUST_UNIFIED__)

using vector_t        = thrust::universal_vector<mini_float>;
using device_vector_t = thrust::universal_vector<mini_float>;

#else

using vector_t        = Vector<mini_float>;
using device_vector_t = Vector<mini_float>;

#endif

#endif // VECTOR_H
