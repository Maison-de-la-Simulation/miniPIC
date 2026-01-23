
/* _____________________________________________________________________ */
//! \file Field.hpp

//! \brief class representing a 3D Field array

/* _____________________________________________________________________ */

// #pragma once
#ifndef FIELD_H
#define FIELD_H

#include "Backend.hpp"
#include "Headers.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

// _________________________________________________________________________________________
//! \brief Data structure, store a 3D field

template <typename T> class Field {
public:
  // _________________________________________________________________________________________
  // Variable members

  //! Name of the field
  std::string name_m;

  //! Sizes in each dimension
  int nx_m, ny_m, nz_m;

  //! Factorized sizes
  int nynz_;

  //! Primal 0 / dual 1
  int dual_x_m, dual_y_m, dual_z_m;

  //! offset to global grid origin
  //! if the direction is dual, the offset is relative to the dual grid
  int ix_offset_;
  int iy_offset_;
  int iz_offset_;

  //! Data linearized, 3rd dimension faster

#if defined(__MINIPIC_KOKKOS__)

  Kokkos::DualView<T ***> data_m;

#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)

  Kokkos::DualView<T ***, Kokkos::SharedSpace> data_m;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

  Kokkos::View<T ***, Kokkos::SharedSpace> data_m;

#elif defined(__MINIPIC_KOKKOS_VIEWS__)

  Kokkos::View<T ***> data_m;
  typename decltype(data_m)::host_mirror_type data_m_h;

#elif defined(__MINIPIC_THRUST__)

  thrust::device_vector<T> device_data_;
  thrust::host_vector<T> host_data_;
  // std::vector<mini_float> host_data_;

#elif defined(__MINIPIC_THRUST_UNIFIED__)

  thrust::universal_vector<T> data_;

#elif defined(__MINIPIC_OPENACC__) || defined(__MINIPIC_OMP_TARGET__)

  std::shared_ptr<std::vector<T>> data_m;

  // Pointer to the data for device management
  T *raw_data_pointer_;

#elif defined(__MINIPIC_SYCL__)

  T *host_data_;
  T *device_data_;

  // pointer to the sycl queue
  sycl::queue *sycl_queue_ptr;

#elif defined(__MINIPIC_STDPAR__)
  std::shared_ptr<std::vector<T>> data_m;

#else

  std::shared_ptr<std::vector<T>> data_m;

#endif

  // _________________________________________________________________________________________
  // Methods

  // _________________________________________________________________________________________
  //! \brief Default constructor - create an empty field
  // _________________________________________________________________________________________
  Field() : name_m("empty"), nx_m(0), ny_m(0), nz_m(0), 
  dual_x_m(0), dual_y_m(0), dual_z_m(0),
  ix_offset_(0), iy_offset_(0), iz_offset_(0) {

    nynz_ = 0;

#if defined(__MINIPIC_KOKKOS_COMMON__)
    // nothing
#elif defined(__MINIPIC_OPENACC__) || defined(__MINIPIC_OMP_TARGET__)
    raw_data_pointer_ = nullptr;
#elif defined(__MINIPIC_SYCL__)
    host_data_     = nullptr;
    device_data_   = nullptr;
    sycl_queue_ptr = nullptr;
#elif defined(__MINIPIC_THRUST__)
    // nothing
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    // nothing
#elif defined(__MINIPIC_STDPAR__)
    data_m = nullptr;
#else
    data_m = nullptr;
#endif
  }

  // _________________________________________________________________________________________
  //! \brief Constructor - allocate memory for the 3D field with 0 values
  //! \param nx number of grid points in the x direction
  //! \param ny number of grid points in the y direction
  //! \param nz number of grid points in the z direction
  //! \param backend backend to use for initialization and memory allocation
  //! \param v default value to fill the field
  //! \param dual_x primal or dual in the x direction
  //! \param dual_y primal or dual in the y direction
  //! \param dual_z primal or dual in the z direction
  //! \param name name of the field
  // _________________________________________________________________________________________
  Field(const int nx,
        const int ny,
        const int nz,
        Backend &backend,
        const T v,
        const int dual_x,
        const int dual_y,
        const int dual_z,
        const std::string name) {
    allocate(nx, ny, nz, backend, v, dual_x, dual_y, dual_z, name);
  }

  // _________________________________________________________________________________________
  //! \brief destructor
  // _________________________________________________________________________________________
  ~Field() {
#if defined(__MINIPIC_OPENACC__)
    if (raw_data_pointer_ != nullptr) {
#pragma acc exit data delete (raw_data_pointer_)
    }
#pragma acc exit data delete (this)
#elif defined(__MINIPIC_OMP_TARGET__)

    if (raw_data_pointer_ != nullptr) {
#pragma omp target exit data map(delete : raw_data_pointer_[0 : nx_m * ny_m * nz_m])
    }
#pragma omp target exit data map(delete : this[0 : 1])
#elif defined(__MINIPIC_SYCL__)
    sycl::free(host_data_, (*sycl_queue_ptr));
    sycl::free(device_data_, (*sycl_queue_ptr));
    sycl_queue_ptr = nullptr;
#elif defined(__MINIPIC_STDPAR__)
// one advantage of shared_ptr is that it can automatically manage the lifetime of objects
#endif
  }

  // _________________________________________________________________________________________
  //
  //! \brief deep copy constructor
  // _________________________________________________________________________________________
  Field(const Field &f) {
    nx_m     = f.nx_m;
    ny_m     = f.ny_m;
    nz_m     = f.nz_m;
    nynz_    = f.nynz_;
    name_m   = f.name_m;
    dual_x_m = f.dual_x_m;
    dual_y_m = f.dual_y_m;
    dual_z_m = f.dual_z_m;

#if defined(__MINIPIC_THRUST__)
    host_data_   = f.host_data_;
    device_data_ = f.device_data_;
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    data_ = f.data_;
#elif defined(__MINIPIC_OMP_TARGET__)
    data_m            = f.data_m;
    raw_data_pointer_ = data_m->data();
#pragma omp target enter data map(alloc : this[0 : 1])
#elif defined(__MINIPIC_OPENACC__)
    data_m            = f.data_m;
    raw_data_pointer_ = data_m->data();
#pragma acc enter data copyin(this)
    acc_attach((void **)&raw_data_pointer_);
#elif defined(__MINIPIC_SYCL__)
    host_data_     = f.host_data_;
    device_data_   = f.device_data_;
    sycl_queue_ptr = f.sycl_queue_ptr;
#elif defined(__MINIPIC_STDPAR__)
    data_m = f.data_m;
#else
    data_m = f.data_m;
#endif
  }

  // _________________________________________________________________________________________
  //
  //! \brief Get 1d index from 3d indexes
  //! \param i index in the x direction
  //! \param j index in the y direction
  //! \param k index in the z direction
  //! \return the 1d index
  // _________________________________________________________________________________________
  inline __attribute__((always_inline)) int index(const int i, const int j, const int k) const {
    return i * (nz_m * ny_m) + j * (nz_m) + k;
  }

  // _________________________________________________________________________________________
  //
  //! \brief Give the total number of points in the grid
  //! \return the total number of points in the grid
  // _________________________________________________________________________________________
  int size() const { return nx_m * ny_m * nz_m; }

  // _________________________________________________________________________________________
  //
  //! \brief Easiest data accessors using 3D indexes
  //! \param i index in the x direction
  //! \param j index in the y direction
  //! \param k index in the z direction
  //! \return the value of the field at the given indexes
  // _________________________________________________________________________________________
  inline __attribute__((always_inline)) T &
  operator()(const int i, const int j, const int k) noexcept {
#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)
    return data_m.h_view(i, j, k);
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    return data_m(i, j, k);
#elif defined(__MINIPIC_KOKKOS_VIEWS__)
    return data_m_h(i, j, k);
#elif defined(__MINIPIC_THRUST__)
    return host_data_[i * (nz_m * ny_m) + j * (nz_m) + k];
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    return data_[i * (nz_m * ny_m) + j * (nz_m) + k];
#elif defined(__MINIPIC_OMP_TARGET__)
    return raw_data_pointer_[i * (nz_m * ny_m) + j * (nz_m) + k];
#elif defined(__MINIPIC_OPENACC__)
    return raw_data_pointer_[i * (nz_m * ny_m) + j * (nz_m) + k];
#elif defined(__MINIPIC_SYCL__)
    return host_data_[i * (nz_m * ny_m) + j * (nz_m) + k];
#elif defined(__MINIPIC_STDPAR__)
  return (*data_m)[i * (nz_m * ny_m) + j * (nz_m) + k];
  // Dereference the shared pointer then access vector elements to retrieve the value
#else
    // return data_m->operator[](i * (nz_m * ny_m) + j * (nz_m) + k);
    return (*data_m)[i * (nz_m * ny_m) + j * (nz_m) + k];
#endif
  }

  // _________________________________________________________________________________________
  //
  //! \brief 1d data accessors
  //! \param idx index in the 1d array
  //! \return the value of the field at the given index
  // _________________________________________________________________________________________
#if defined(__MINIPIC_KOKKOS_COMMON__)
  // no relevant since we do not assume a specific layout
  // Should be computed using the layout properties : stride, etc...
#elif defined(__MINIPIC_THRUST__)
  inline __attribute__((always_inline)) T &operator[](const int idx) { return host_data_[idx]; }
#elif defined(__MINIPIC_THRUST_UNIFIED__)
  inline __attribute__((always_inline)) T &operator[](const int idx) { return data_[idx]; }
#elif defined(__MINIPIC_KOKKOS_VIEWS__)
  inline __attribute__((always_inline)) T &operator[](const int idx) { return data_m_h.data()[idx]; }
#elif defined(__MINIPIC_OMP_TARGET__)
  inline __attribute__((always_inline)) T &operator[](const int idx) {
    return raw_data_pointer_[idx];
  }
#elif defined(__MINIPIC_OPENACC__)
  inline __attribute__((always_inline)) T &operator[](const int idx) {
    return raw_data_pointer_[idx];
  }
#elif defined(__MINIPIC_SYCL__)
  inline __attribute__((always_inline)) T &operator[](const int idx) { return host_data_[idx]; }
#elif defined(__MINIPIC_STDPAR__)
  inline __attribute__((always_inline)) T &operator[](const int idx) { return data_m[idx]; }
#else
  inline __attribute__((always_inline)) T &operator[](const int idx) { return data_m[idx]; }
#endif

  //! \brief return the number of grid points in the x direction
  //! \return return the number of grid points in the x direction
  INLINE int nx() const { return nx_m; }

  //! \brief return the number of grid points in the y direction
  //! \return return the number of grid points in the y direction
  INLINE int ny() const { return ny_m; }

  //! \brief return the number of grid points in the z direction
  //! \return return the number of grid points in the z direction
  INLINE int nz() const { return nz_m; }

  //! \brief return the number of grid points in the y*z direction
  //! \return return the number of grid points in the y*z direction
  INLINE int nynz() const { return nynz_; }

  // _________________________________________________________________________________________
  //
  //! \brief Alloc memory for the 3D field
  //! \param nx number of grid points in the x direction
  //! \param ny number of grid points in the y direction
  //! \param nz number of grid points in the z direction
  //! \param v default value
  //! \param dual_x dual in the x direction
  //! \param dual_y dual in the y direction
  //! \param dual_z dual in the z direction
  //! \param name name of the field
  // _________________________________________________________________________________________
  void allocate(const int nx,
                const int ny,
                const int nz,
                Backend &backend,
                const T v        = 0,
                const int dual_x = 0,
                const int dual_y = 0,
                const int dual_z = 0,
                const int ix_offset = 0,
                const int iy_offset = 0,
                const int iz_offset = 0,
                std::string name = "") {

    nx_m = nx;
    ny_m = ny;
    nz_m = nz;

    nynz_ = ny * nz;

    dual_x_m = dual_x;
    dual_y_m = dual_y;
    dual_z_m = dual_z;

    ix_offset_ = ix_offset;
    iy_offset_ = iy_offset;
    iz_offset_ = iz_offset;

    name_m   = name;

    if (nx_m * ny_m * nz_m == 0) {
      return;
    }

#if defined(__MINIPIC_KOKKOS__)
    data_m = Kokkos::DualView<T ***>(name, nx, ny, nz);
#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
    data_m = Kokkos::DualView<T ***, Kokkos::SharedSpace>(name, nx, ny, nz);
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    data_m = Kokkos::View<T ***, Kokkos::SharedSpace>(name, nx, ny, nz);
#elif defined(__MINIPIC_KOKKOS_VIEWS__)
    data_m   = Kokkos::View<T ***>(name, nx, ny, nz);
    data_m_h = Kokkos::create_mirror_view(data_m);
#elif defined(__MINIPIC_THRUST__)
    resize(nx, ny, nz, v);
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    resize(nx, ny, nz, v);
#elif defined(__MINIPIC_OMP_TARGET__)
    data_m            = std::make_shared<std::vector<T>>(nx * ny * nz, v);
    raw_data_pointer_ = data_m->data();
#elif defined(__MINIPIC_OPENACC__)
    data_m            = std::make_shared<std::vector<T>>(nx * ny * nz, v);
    raw_data_pointer_ = data_m->data();
#elif defined(__MINIPIC_SYCL__)
    sycl_queue_ptr = backend.sycl_queue_;
    host_data_     = sycl::malloc_host<T>(nx * ny * nz, (*sycl_queue_ptr));
    device_data_   = sycl::malloc_device<T>(nx * ny * nz, (*sycl_queue_ptr));
#elif defined(__MINIPIC_STDPAR__)
    data_m = std::make_shared<std::vector<T>>(nx * ny * nz, v);
#else
    // data_m = new std::vector<T>(nx * ny * nz, v);
    data_m = std::make_shared<std::vector<T>>(nx * ny * nz, v);
    // resize(nx, ny, nz, v);
#endif

    fill(v, minipic::host);

#if defined(__MINIPIC_OPENACC__)
#pragma acc enter data copyin(this)
#pragma acc enter data create(raw_data_pointer_[0 : nx * ny * nz])
#elif defined(__MINIPIC_OMP_TARGET__)
#pragma omp target enter data map(alloc : this[0 : 1])
#pragma omp target enter data map(to : raw_data_pointer_[0 : nx * ny * nz])
#endif
  }

  // _________________________________________________________________________________________
  //
  //! \brief Resize field
  //! \warning This function only preserves the data on the host
  //! \warning Data is not updated on device after resizing
  //! \param nx number of grid points in the x direction
  //! \param ny number of grid points in the y direction
  //! \param nz number of grid points in the z direction
  //! \param v default value
  // _________________________________________________________________________________________
  void resize(const int nx, const int ny, const int nz, const T v = 0) {
    nx_m  = nx;
    ny_m  = ny;
    nz_m  = nz;
    nynz_ = ny * nz;
#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)
    data_m.resize(nx, ny, nz);
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    data_m.resize(nx, ny, nz);
#elif defined(__MINIPIC_KOKKOS_VIEWS__)
    Kokkos::resize(data_m, nx, ny, nz);
    Kokkos::resize(data_m_h, nx, ny, nz);
#elif defined(__MINIPIC_THRUST__)
    host_data_.resize(nx * ny * nz, v);
    device_data_.resize(nx * ny * nz, v);
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    data_.resize(nx * ny * nz, v);
#elif defined(__MINIPIC_OMP_TARGET__)
    // We do not resize on device, first bring back data to host
    sync(minipic::device, minipic::host);
// Destroy the device data
#pragma omp target exit data map(delete : raw_data_pointer_[0 : nx * ny * nz])
    data_m->resize(nx * ny * nz, v);
    // get the new pointer to recreate the device data
    ptr = get_raw_pointer(minipic::host);
#pragma omp target enter data map(to : raw_data_pointer_[0 : nx * ny * nz])
#elif defined(__MINIPIC_OPENACC__)
    // We do not resize on device, first bring back data to host
    sync(minipic::device, minipic::host);
// Destroy the device data
#pragma acc exit data delete (raw_data_pointer_)
    data_m->resize(nx * ny * nz, v);
    raw_data_pointer_ = data_m->data();
#pragma acc enter data create(raw_data_pointer_[0 : nx * ny * nz])
#elif defined(__MINIPIC_SYCL__)
    // sycl::free(host_data_, (*sycl_queue_ptr));
    // sycl::free(device_data_, (*sycl_queue_ptr));
    // host_data_ = sycl::malloc_host<T>(nx * ny * nz, (*sycl_queue_ptr));
    // device_data_ = sycl::malloc_device<T>(nx * ny * nz, (*sycl_queue_ptr));
#elif defined(__MINIPIC_STDPAR__)
    data_m->resize(nx * ny * nz, v);
#else
    data_m->resize(nx * ny * nz, v);
#endif
  }

  // _________________________________________________________________________________________
  //
  //! Set the name of the field
  //! \param name name of the field
  // _________________________________________________________________________________________
  void set_name(std::string name) { name_m = name; }

  // _________________________________________________________________________________________
  //
  //! \brief Set all the field at value v
  //! \param v value to set
  //! \param space space where to set the value
  // _________________________________________________________________________________________
  template <class T_space> void fill(const mini_float v, const T_space space) {

    // ---> Host case
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

      // Kokkos::Experimental::fill(Kokkos::DefaultHostExecutionSpace(), data_m.h_view, 0.);

      // for (auto ix = 0; ix < nx_m; ++ix) {
      //   for (auto iy = 0; iy < ny_m; ++iy) {
      //     for (auto iz = 0; iz < nz_m; ++iz) {
      //       data_m.h_view(ix, iy, iz) = v;
      //     }
      //   }
      // }

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<3>>
        mdrange_policy;

      Kokkos::parallel_for(
        mdrange_policy({0, 0, 0}, {nx_m, ny_m, nz_m}),
        KOKKOS_CLASS_LAMBDA(const int ix, const int iy, const int iz) {
          data_m.h_view(ix, iy, iz) = v;
        });

      Kokkos::fence();

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<3>>
        mdrange_policy;

      Kokkos::parallel_for(
        mdrange_policy({0, 0, 0}, {nx_m, ny_m, nz_m}),
        KOKKOS_CLASS_LAMBDA(const int ix, const int iy, const int iz) { data_m(ix, iy, iz) = v; });

      Kokkos::fence();

#elif defined(__MINIPIC_KOKKOS_VIEWS__)

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<3>>
        mdrange_policy;

      auto &data_ref = data_m_h;

      Kokkos::parallel_for(
        mdrange_policy({0, 0, 0}, {nx_m, ny_m, nz_m}),
        KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
          data_ref(ix, iy, iz) = v;
        });

      Kokkos::fence();

#elif defined(__MINIPIC_THRUST__)
      // thrust::fill(host_data_.begin(), host_data_.end(), v);
      for (auto i = 0; i < size(); ++i) {
        host_data_[i] = v;
      }

#elif defined(__MINIPIC_THRUST_UNIFIED__)
      for (auto i = 0; i < size(); ++i) {
        data_[i] = v;
      }

#elif defined(__MINIPIC_SYCL__)

      for (auto i = 0; i < size(); ++i) {
        host_data_[i] = v;
      }

// #elif defined(__MINIPIC__EVENTIFY__)
//       for (auto i = 0; i < size(); ++i) {
//         data_m[i] = v;
//       }
#elif defined(__MINIPIC_STDPAR__)
      std::fill(std::execution::seq, data_m->begin(), data_m->end(), v);
      // for (auto i = 0; i < size(); ++i) {
      //   data_m[i] = v;
      // }
#else
      std::fill(data_m->begin(), data_m->end(), v);
#endif

      // ---> Device case
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
#if defined(__MINIPIC_KOKKOS_COMMON__)

#if defined(__MINIPIC_KOKKOS__)
      typename Kokkos::DualView<T ***>::t_dev &F = data_m.d_view;
#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
      typename Kokkos::DualView<T ***, Kokkos::SharedSpace>::t_dev &F = data_m.d_view;
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
      typename Kokkos::View<T ***, Kokkos::SharedSpace> &F = data_m;
#elif defined(__MINIPIC_KOKKOS_VIEWS__)
      typename Kokkos::View<T ***> &F = data_m;
#endif

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
      Kokkos::parallel_for(
        mdrange_policy({0, 0, 0}, {nx_m, ny_m, nz_m}),
        KOKKOS_CLASS_LAMBDA(const int ix, const int iy, const int iz) { F(ix, iy, iz) = v; });

      Kokkos::fence();

#elif defined(__MINIPIC_THRUST__)
      thrust::fill(device_data_.begin(), device_data_.end(), v);
      // thrust::fill(host_data_.begin(), host_data_.end(), v);
      // for (auto i = 0; i < size(); ++i) {
      //   host_data_[i] = v;
      // }
#elif defined(__MINIPIC_THRUST_UNIFIED__)
      thrust::fill(data_.begin(), data_.end(), v);
      // for (auto i = 0; i < size(); ++i) {
      //   data_[i] = v;
      // }
#elif defined(__MINIPIC_OMP_TARGET__)

#pragma omp target teams distribute parallel for
      for (auto i = 0; i < size(); ++i) {
        raw_data_pointer_[i] = v;
      }

#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel loop gang worker vector present(raw_data_pointer_[0 : nx_m * ny_m * nz_m])
      for (auto i = 0; i < size(); ++i) {
        raw_data_pointer_[i] = v;
      }

#elif defined(__MINIPIC_SYCL__)

      T *device_data = device_data_;

      sycl_queue_ptr->parallel_for(sycl::range<1>{static_cast<size_t>(nx_m * ny_m * nz_m)},
                                   [=](sycl::item<1> i) { device_data[i] = v; });
      sycl_queue_ptr->wait();

#elif defined(__MINIPIC_STDPAR__)
      std::fill(std::execution::par_unseq, data_m->begin(), data_m->end(), v);

#else
      std::fill(data_m->begin(), data_m->end(), v);
#endif
    }
  }

  // _________________________________________________________________________________________
  //
  //! \brief Set all field values to 0
  // _________________________________________________________________________________________
  template <class T_space> void reset(const T_space space) { fill(0, space); }

  // _________________________________________________________________________________________
  //
  //! \brief return the pointer to the data
  //! \param space space where to get the pointer
  //! \tparam T_space execution space
  //! \return return pointer to the first element of the data
  // _________________________________________________________________________________________
  template <class T_space = minipic::Host> T *get_raw_pointer(const T_space space) {

    static_assert(std::is_same<T_space, minipic::Host>::value ||
                    std::is_same<T_space, minipic::Device>::value,
                  "T_space must be either minipic::Host or minipic::Device");

#if defined(__MINIPIC_KOKKOS__) || defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
      return data_m.h_view.data();
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
      return data_m.d_view.data();
    } else {
      return data_m.h_view.data();
    }
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
    return data_m.data();
#elif defined(__MINIPIC_KOKKOS_VIEWS__)
    if constexpr (std::is_same<T_space, minipic::Host>::value)
    {
      return data_m_h.data();
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
      return data_m.data();
    } else {
      return data_m_h.data();
    }
#elif defined(__MINIPIC_THRUST__)
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
      return thrust::raw_pointer_cast(host_data_.data());
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
      return thrust::raw_pointer_cast(device_data_.data());
    } else {
      return thrust::raw_pointer_cast(host_data_.data());
    }
    // required by nvhpc to avaid a warning
    return thrust::raw_pointer_cast(host_data_.data());
#elif defined(__MINIPIC_THRUST_UNIFIED__)
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
      return thrust::raw_pointer_cast(data_.data());
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
      return thrust::raw_pointer_cast(data_.data());
    } else {
      return thrust::raw_pointer_cast(data_.data());
    }
#elif defined(__MINIPIC_OMP_TARGET__)
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
      return data_m->data();
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
      return data_m->data();
    } else {
      return data_m->data();
    }
#elif defined(__MINIPIC_OPENACC__)
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
      return raw_data_pointer_;
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
      return raw_data_pointer_;
    } else {
      return raw_data_pointer_;
    }
#elif defined(__MINIPIC_SYCL__)
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
      return host_data_;
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
      return device_data_;
    } else {
      return host_data_;
    }
#elif defined(__MINIPIC_STDPAR__)
  return data_m->data();
  // return (*data_m).data();
  // returns a pointer to the first element of the vector (data buffer address)
#else
    return data_m->data();
#endif
  }

  // ____________________________________________________________
  //
  //! \brief output the sum of data with power power
  // ____________________________________________________________
  template <class T_space> T sum(const int power, T_space space) const {
    T sum = 0;

    // ---> Host case
    if constexpr (std::is_same<T_space, minipic::Host>::value) {

#if defined(__MINIPIC_KOKKOS__) || defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<3>>
        mdrange_policy;
      Kokkos::parallel_reduce(
        "sum_field_on_host",
        mdrange_policy({0, 0, 0}, {nx_m, ny_m, nz_m}),
        KOKKOS_CLASS_LAMBDA(const int ix, const int iy, const int iz, T &local_sum) {
          local_sum += Kokkos::pow(data_m.h_view(ix, iy, iz), power);
        },
        sum);

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<3>>
        mdrange_policy;

      Kokkos::parallel_reduce(
        "sum_field_on_host",
        mdrange_policy({0, 0, 0}, {nx_m, ny_m, nz_m}),
        KOKKOS_CLASS_LAMBDA(const int ix, const int iy, const int iz, T &local_sum) {
          local_sum += Kokkos::pow(data_m(ix, iy, iz), power);
        },
        sum);

#elif defined(__MINIPIC_KOKKOS_VIEWS__)

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<3>>
        mdrange_policy;

      auto &data_ref = data_m_h;

      Kokkos::parallel_reduce(
        "sum_field_on_host",
        mdrange_policy({0, 0, 0}, {nx_m, ny_m, nz_m}),
        KOKKOS_LAMBDA(const int ix, const int iy, const int iz, T &local_sum) {
          local_sum += Kokkos::pow(data_ref(ix, iy, iz), power);
        },
        sum);

#elif defined(__MINIPIC_THRUST__)
      // sum = thrust::transform_reduce(
      //   host_data_.begin(), host_data_.end(), [] __host__ (T x) { return x * x; }, 0.0,
      //   thrust::plus<T>());
      for (int i = 0; i < size(); i++) {
        sum += pow(host_data_[i], power);
      }
#elif defined(__MINIPIC_THRUST_UNIFIED__)
      for (int i = 0; i < size(); i++) {
        sum += pow(data_[i], power);
      }
#elif defined(__MINIPIC_OMP_TARGET__) || defined(__MINIPIC_OPENACC__)
      for (int i = 0; i < size(); i++) {
        sum += pow(raw_data_pointer_[i], power);
      }
#elif defined(__MINIPIC_SYCL__)
      for (int i = 0; i < size(); i++) {
        sum += sycl::pown(host_data_[i], power);
      }

#elif defined(__MINIPIC_STDPAR__)

      for (int i = 0; i < size(); i++) {
        sum += pow((*data_m)[i], power);
      }

      // sum = std::transform_reduce((*data_m).begin(),
      //                             (*data_m).end(),
      //                             static_cast<T>(0),
      //                             std::plus<T>(),
      //                             [power](T x) { return std::pow(x, power); });

#else
      for (int i = 0; i < size(); i++) {
        sum += pow((*data_m)[i], power);
      }
#endif

      // ---> Device case
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
#if defined(__MINIPIC_KOKKOS_COMMON__)

#if defined(__MINIPIC_KOKKOS__)
      typename Kokkos::DualView<T ***>::t_dev F = data_m.d_view;
#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)
      typename Kokkos::DualView<T ***, Kokkos::SharedSpace>::t_dev F = data_m.d_view;
#elif defined(__MINIPIC_KOKKOS_UNIFIED__)
      typename Kokkos::View<T ***, Kokkos::SharedSpace> F = data_m;
#elif defined(__MINIPIC_KOKKOS_VIEWS__)
      typename Kokkos::View<T ***> F = data_m;
#endif

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
      Kokkos::parallel_reduce(
        "sum_field_on_device",
        mdrange_policy({0, 0, 0}, {nx_m, ny_m, nz_m}),
        KOKKOS_CLASS_LAMBDA(const int ix, const int iy, const int iz, T &local_sum) {
          local_sum += Kokkos::pow(F(ix, iy, iz), power);
        },
        sum);

      Kokkos::fence();

#elif defined(__MINIPIC_THRUST__)
      sum = thrust::transform_reduce(
        device_data_.begin(),
        device_data_.end(),
        [=] __host__ __device__(T x) { return pow(x, power); },
        //[=] __host__ __device__(T x) { return x; },
        0.0,
        thrust::plus<T>());
      // for (int i = 0; i < size(); i++) {
      //   sum += pow(host_data_[i],power);
      // }

#elif defined(__MINIPIC_THRUST_UNIFIED__)

      sum = thrust::transform_reduce(
        data_.begin(),
        data_.end(),
        [=] __host__ __device__(T x) { return pow(x, power); },
        //[=] __host__ __device__(T x) { return x; },
        0.0,
        thrust::plus<T>());
      // for (int i = 0; i < size(); i++) {
      //   sum += pow(host_data_[i],power);
      // }

#elif defined(__MINIPIC_OMP_TARGET__)

#pragma omp target teams distribute parallel for reduction(+ : sum)
      for (int i = 0; i < size(); i++) {
        sum += pow(raw_data_pointer_[i], power);
      }

#elif defined(__MINIPIC_OPENACC__)

#pragma acc parallel loop gang worker vector present(raw_data_pointer_[0 : nx_m * ny_m * nz_m]) \
  reduction(+ : sum)
      for (int i = 0; i < size(); i++) {
        sum += pow(raw_data_pointer_[i], power);
      }
#elif defined(__MINIPIC_SYCL__)

      // buffer on device for sum
      sycl::buffer<T, 1> d_sum(&sum, 1);

      sycl::buffer<T, 1> data_buf(device_data_, nx_m * ny_m * nz_m);

  // Execute a reduction on the device
      // sum = sycl_queue_ptr->submit([&](sycl::handler& cgh) {
      //     auto data_acc = data_buf.get_access(cgh);
      //     return sycl::reduce(cgh, data_acc, 0.0f, sycl::plus<>());
      // }).get();

      sycl_queue_ptr->submit([&](sycl::handler &cgh) {
        sycl::accessor ksum{d_sum, cgh, sycl::write_only};
        sycl::accessor data_acc{data_buf, cgh, sycl::read_only};

        cgh.parallel_for(
          sycl::range<1>{nx_m * ny_m * nz_m},
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
      // ); A condition que data_ ne soit pas de type int !

      sum = std::transform_reduce(std::execution::par_unseq,
                                  data_m->begin(),
                                  data_m->end(),
                                  static_cast<mini_float>(0),
                                  std::plus<T>(),
                                  [power](T x) { return std::pow(x, power); });

      // sum = std::reduce(std::execution::par_unseq, data_m->begin(), data_m->end());

      // sum = std::transform_reduce( std::execution::par_unseq, (*data_m).begin(), (*data_m).end(),
      // static_cast<T>(0), std::plus<T>(), [power](T x) { return std::pow(x, power); } );
#else
      for (int i = 0; i < size(); i++) {
        sum += pow((*data_m)[i], power);
      }
#endif
    }

    return sum;
  }

  // _________________________________________________________________________________________
  //! \brief output the field as a string
  //! \return std::string
  // _________________________________________________________________________________________
  std::string to_string() {
    std::string buffer = "Field " + name_m + "\n";
    buffer += "__________________________________ \n";
    for (auto ix = 0; ix < nx_m; ++ix) {
      buffer += "\n";
      for (auto iy = 0; iy < ny_m; ++iy) {
        buffer += "\n";
        for (auto iz = 0; iz < nz_m; ++iz) {
          // const T field = h(ix, iy, iz);
          // to string with scientific notation
          std::ostringstream out;
          out << std::scientific << this->operator()(ix, iy, iz);
          std::string s = out.str();
          buffer += s + " ";

          // buffer += std::to_string(static_cast<T>(this->operator()(ix, iy, iz))) + " ";
        }
      }
    }
    buffer += "\n __________________________________ \n";
    return buffer;
  }

  // _________________________________________________________________________________________
  //
  //! \brief print all values of the field on host
  // _________________________________________________________________________________________
  void print() {
    std::string buffer = to_string();
    std::cout << buffer << std::endl;
  }

  // _________________________________________________________________________________________
  //
  //! \brief print the sum of the field on host
  // _________________________________________________________________________________________
  void check_sum() {
    T sum = sum();
    std::cout << name_m << " sum: " << sum << std::endl;
  }

  // _________________________________________________________________________________________
  //
  //! \brief Sync Host <-> Device
  // _________________________________________________________________________________________
  template <class T_from, class T_to> void sync(const T_from from, const T_to to) {
    // ---> Host to Device
    if constexpr (std::is_same<T_from, minipic::Host>::value) {
#if defined(__MINIPIC_KOKKOS__)
      data_m.modify_host();
      data_m.template sync<typename Kokkos::DualView<T ***>::execution_space>();
#elif defined(__MINIPIC_KOKKOS_VIEWS__)
      Kokkos::deep_copy(data_m, data_m_h);
#elif defined(__MINIPIC_THRUST__)
      thrust::copy(host_data_.begin(), host_data_.begin() + size(), device_data_.begin());
#elif defined(__MINIPIC_THRUST_UNIFIED__)
      // Nothing to do, data_ is a unified memory pointer
#elif defined(__MINIPIC_OMP_TARGET__)
#pragma omp target update to(raw_data_pointer_[ : size()])
#elif defined(__MINIPIC_OPENACC__)
#pragma acc update device(raw_data_pointer_[ : size()])
#elif defined(__MINIPIC_SYCL__)
      sycl_queue_ptr->memcpy(device_data_, host_data_, nx_m * ny_m * nz_m * sizeof(T));
      sycl_queue_ptr->wait();
#elif defined(__MINIPIC_STDPAR__)
      // nothing (UVM)
#endif
      // ---> Device to Host
    } else if constexpr (std::is_same<T_from, minipic::Device>::value) {
#if defined(__MINIPIC_KOKKOS__)
      data_m.modify_device();
      data_m.template sync<typename Kokkos::DualView<T ***>::host_mirror_space>();
#elif defined(__MINIPIC_KOKKOS_VIEWS__)
      Kokkos::deep_copy(data_m_h, data_m);
#elif defined(__MINIPIC_THRUST__)
      thrust::copy(device_data_.begin(), device_data_.begin() + size(), host_data_.begin());
#elif defined(__MINIPIC_THRUST_UNIFIED__)
      // Nothing to do, data_ is a unified memory pointer
#elif defined(__MINIPIC_OMP_TARGET__)
#pragma omp target update from(raw_data_pointer_[ : size()])
#elif defined(__MINIPIC_OPENACC__)
#pragma acc update host(raw_data_pointer_[ : size()])
#elif defined(__MINIPIC_SYCL__)
      sycl_queue_ptr->memcpy(host_data_, device_data_, nx_m * ny_m * nz_m * sizeof(T));
      sycl_queue_ptr->wait();
#elif defined(__MINIPIC_STDPAR__)
      // nothing (UVM)
#endif
    }
  }
};

// _________________________________________________________________________________________
// Shortucts for the different backends

#if defined(__MINIPIC_KOKKOS__)

using device_field_t = Kokkos::DualView<mini_float ***>::t_dev;
using field_t        = Kokkos::DualView<mini_float ***>::t_host;

#elif defined(__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__)

using device_field_t = Kokkos::DualView<mini_float ***, Kokkos::SharedSpace>::t_dev;
using field_t        = Kokkos::DualView<mini_float ***, Kokkos::SharedSpace>::t_host;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

using device_field_t = Kokkos::View<mini_float ***, Kokkos::SharedSpace>;
using field_t        = Kokkos::View<mini_float ***, Kokkos::SharedSpace>;

#elif defined(__MINIPIC_KOKKOS_VIEWS__)

using device_field_t = Kokkos::View<mini_float ***>;
using field_t        = typename device_field_t::host_mirror_type;

#elif defined(__MINIPIC_THRUST__)

using device_field_t = thrust::device_vector<mini_float>;
// using device_field_t = thrust::host_vector<mini_float>;
using field_t = thrust::host_vector<mini_float>;

#elif defined(__MINIPIC_THRUST_UNIFIED__)
using device_field_t = thrust::universal_vector<mini_float>;
using field_t        = thrust::universal_vector<mini_float>;

// using device_field_t = std::vector<mini_float>;
// using grid_t = std::vector<mini_float>;

#elif defined(__MINIPIC_STDPAR__)
// A voir
using device_field_t = std::vector<mini_float>;
using grid_t         = std::vector<mini_float>;

#else

using device_field_t = std::vector<mini_float>;
using grid_t         = std::vector<mini_float>;

#endif

#endif // end FIELD_H
