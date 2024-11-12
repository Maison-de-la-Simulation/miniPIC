
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
  unsigned int size_;

  // Main data
  std::vector<T> data_;

  // ______________________________________________________________________
  //
  //! \brief constructors
  // ______________________________________________________________________
  Vector() : size_(0) {}
  Vector(unsigned int size, Backend &backend) { allocate("", backend, size); }

  // ______________________________________________________________________
  //
  //! \brief constructor with allocation
  //! \param[in] size allocation size
  //! \param[in] v default value
  // ______________________________________________________________________
  Vector(unsigned int size, T v, Backend &backend) {
    allocate("", size, backend);
    fill(v);
  }

  // ______________________________________________________________________
  //
  //! \brief destructor
  // ______________________________________________________________________
  ~Vector() {
  };

  // ______________________________________________________________________
  //
  //! \brief allocate the data_ object
  // ______________________________________________________________________
  void allocate(std::string name, unsigned int size, Backend &backend) {
    size_ = size;
    data_.resize(size);
  }

  // ______________________________________________________________________
  //
  //! \brief [] operator
  //! \return Host data accessor (if not device, point to the host data)
  // ______________________________________________________________________
  INLINE T &operator[](const int i) {
    return data_[i];
  }

  // ______________________________________________________________________
  //
  //! \brief () operator
  //! \return Host data accessor (if not device, point to the host data)
  // ______________________________________________________________________
  INLINE T &operator()(const int i) {
    return data_[i];
  }

  // ______________________________________________________________________
  //
  //! \brief Explicit Host data accessor
  //! \return host pointer at index i
  // ______________________________________________________________________
  INLINE T &h(const int i) {
    return data_[i];
  }

  // ______________________________________________________________________
  //
  //! \brief Get the data pointer
  //! \param[in] space where to keep the data when resizing (must be minipic::host or
  //! minipic::device)
  // ______________________________________________________________________
  template <class T_space> T * get_raw_pointer(const T_space space) {

    // Check that T_Space of Class Host or Device
    static_assert(std::is_same<T_space, minipic::Host>::value ||
                    std::is_same<T_space, minipic::Device>::value,
                  "Must be minipic::host or minipic::device");

    // Host
    if constexpr (std::is_same<T_space, minipic::Host>::value) {
      return data_.data();

      // Device
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {
      return data_.data();

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
  template <class T_space> void resize(const unsigned int new_size, const T_space space) {

    // Check that T_Space of Class Host or Device
    static_assert(std::is_same<T_space, minipic::Host>::value ||
                    std::is_same<T_space, minipic::Device>::value,
                  "Must be minipic::host or minipic::device");
    data_.resize(new_size);
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
  template <class T_space> void resize(const unsigned int new_size, T value, const T_space space) {

    resize(new_size, space); // j'appele la méthode d'avant
    for (auto ip = size_; ip < new_size; ++ip) {
      data_[ip] = value;
    }
  }

  // ______________________________________________________________________
  //
  //! \brief clear the content, equivalent to size_ = 0
  //! If the raw object has a clear method, we call it
  // ______________________________________________________________________
  void clear() {
    size_ = 0;
    data_.clear();
  }

  // ______________________________________________________________________
  //
  //! \brief fill the vector with the given value
  // ______________________________________________________________________
  void fill(const T v) {

    // void fill( ForwardIt first, ForwardIt last, const T& value );
    std::fill(data_.begin(), data_.end(), v);

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

      for (int i = 0; i < size_; i++) {
        sum += pow(data_[i], power);
      }

      // ---> Device case
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {

      for (int i = 0; i < size_; i++) {
        sum += pow(data_[i], power);
      }

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



      // Device -> Host
    } else if constexpr (std::is_same<from, minipic::Device>::value &&
                         std::is_same<to, minipic::Host>::value) {


    }
  }
};

// _________________________________________________________________________
// Shortcuts

using vector_t        = Vector<mini_float>;
using device_vector_t = Vector<mini_float>;

#endif // VECTOR_H
