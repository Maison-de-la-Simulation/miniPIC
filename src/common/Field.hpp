
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

  //! Primal 0 / dual 1
  int dual_x_m, dual_y_m, dual_z_m;

  //! Data linearized, 3rd dimension faster
  std::shared_ptr<std::vector<T>> data_m;

  // _________________________________________________________________________________________
  // Methods

  // _________________________________________________________________________________________
  //! \brief Default constructor - create an empty field
  // _________________________________________________________________________________________
  Field() : nx_m(0), ny_m(0), nz_m(0), name_m("empty"), dual_x_m(0), dual_y_m(0), dual_z_m(0) {
    data_m = nullptr;
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
  }

  // _________________________________________________________________________________________
  //
  //! \brief shallow copy constructor
  // _________________________________________________________________________________________
  Field(const Field &f) {
    nx_m     = f.nx_m;
    ny_m     = f.ny_m;
    nz_m     = f.nz_m;
    name_m   = f.name_m;
    dual_x_m = f.dual_x_m;
    dual_y_m = f.dual_y_m;
    dual_z_m = f.dual_z_m;
    data_m = f.data_m;
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
    // return data_m->operator[](i * (nz_m * ny_m) + j * (nz_m) + k);
    return (*data_m)[i * (nz_m * ny_m) + j * (nz_m) + k];
  }

  // _________________________________________________________________________________________
  //
  //! \brief 1d data accessors
  //! \param idx index in the 1d array
  //! \return the value of the field at the given index
  // _________________________________________________________________________________________
  inline __attribute__((always_inline)) T &operator[](const int idx) { return data_m[idx]; }

  //! \brief return the number of grid points in the x direction
  //! \return return the number of grid points in the x direction
  INLINE int nx() const { return nx_m; }

  //! \brief return the number of grid points in the y direction
  //! \return return the number of grid points in the y direction
  INLINE int ny() const { return ny_m; }

  //! \brief return the number of grid points in the z direction
  //! \return return the number of grid points in the z direction
  INLINE int nz() const { return nz_m; }

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
                std::string name = "") {

    nx_m = nx;
    ny_m = ny;
    nz_m = nz;

    dual_x_m = dual_x;
    dual_y_m = dual_y;
    dual_z_m = dual_z;
    name_m   = name;

    if (nx_m * ny_m * nz_m == 0) {
      return;
    }

    // data_m = new std::vector<T>(nx * ny * nz, v);
    data_m = std::make_shared<std::vector<T>>(nx * ny * nz, v);
    // resize(nx, ny, nz, v);

    fill(v, minipic::host);

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
    nx_m = nx;
    ny_m = ny;
    nz_m = nz;
    data_m->resize(nx * ny * nz, v);
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

      std::fill(data_m->begin(), data_m->end(), v);

      // ---> Device case
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {

      std::fill(data_m->begin(), data_m->end(), v);

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
  //! \return return pointer to the first element of the data
  // _________________________________________________________________________________________
  template <class T_space = minipic::Host> T *get_raw_pointer(const T_space space) {
    return data_m->data();
  }

  // ____________________________________________________________
  //
  //! \brief output the sum of data with power power
  // ____________________________________________________________
  template <class T_space> T sum(const int power, T_space space) const {
    T sum = 0;

    // ---> Host case
    if constexpr (std::is_same<T_space, minipic::Host>::value) {

      for (int i = 0; i < size(); i++) {
        sum += pow((*data_m)[i], power);
      }

      // ---> Device case
    } else if constexpr (std::is_same<T_space, minipic::Device>::value) {

      for (int i = 0; i < size(); i++) {
        sum += pow((*data_m)[i], power);
      }
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

      // ---> Device to Host
    } else if constexpr (std::is_same<T_from, minipic::Device>::value) {

    }
  }
};

// _________________________________________________________________________________________
// Shortucts for the different backends

using device_field_t = std::vector<mini_float>;
using grid_t         = std::vector<mini_float>;

#endif // FIELD_H
