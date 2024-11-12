
/* _____________________________________________________________________ */
//! \file Particles_SoA.hpp

//! \brief Particles class to store the particles data as contiguous arrays
//!        for each parameter

//! The parameter `n_particles_m` must be used to store the number of particles,
//! the `std::size` method should not be used for this purpose

/* _____________________________________________________________________ */

#pragma once

#include <cstdio>

#include <iomanip>
#include <math.h>
#include <random>

#include "Backend.hpp"
#include "Vector.hpp"

// ________________________________________________________
//
//! \brief Represent an array of particles for 1 species
// ________________________________________________________
template <typename T> class Particles {
public:
  Particles() : n_particles_m(0) {}
  ~Particles() {}

  //! Inverse of the cell volume, usefull to compute the density
  T inv_cell_volume_m;

  //! Number of particles at init
  int n_particles_m;

  //! Species electric charge
  T charge_m;
  //! Species mass
  T mass_m;
  //! Species temperature
  T temperature_m;

  //! Particles positions in 3D
  Vector<T> x_;
  Vector<T> y_;
  Vector<T> z_;
  //! Particles momentums in 3D
  Vector<T> mx_;
  Vector<T> my_;
  Vector<T> mz_;

  //! Weights | Charge density (scalar)
  //! w0,w1,w2,...
  Vector<T> weight_;

  //! Electric field interpolate
  Vector<T> Ex_;
  Vector<T> Ey_;
  Vector<T> Ez_;
  //! Magnetic field interpolate
  Vector<T> Bx_;
  Vector<T> By_;
  Vector<T> Bz_;

  //! Inverse of Lorentz factor: to avoid recompute gamma between PIC steps
  // Vector<T> gamma_inv_;

  //! This flag when false prevents the allocation of E and B fields
  bool with_electromagnetic_fields_ = true;

  // ________________________________________________________
  // Host data accessors

  //! \brief x accessor
  //! \param[in] ip particle index
  INLINE T &x_h(unsigned int ip) { return x_.h(ip); }

  //! \brief y accessor
  //! \param[in] ip particle index
  INLINE T &y_h(unsigned int ip) { return y_.h(ip); }

  //! \brief z accessor
  //! \param[in] ip particle index
  INLINE T &z_h(unsigned int ip) { return z_.h(ip); }

  //! \brief mx accessor
  //! \param[in] ip particle index
  INLINE T &mx_h(unsigned int ip) { return mx_.h(ip); }

  //! \brief my accessor
  //! \param[in] ip particle index
  INLINE T &my_h(unsigned int ip) { return my_.h(ip); }

  //! \brief mz accessor
  //! \param[in] ip particle index
  INLINE T &mz_h(unsigned int ip) { return mz_.h(ip); }

  //! \brief w accessor
  //! \param[in] ip particle index
  INLINE T &w_h(unsigned int ip) { return weight_.h(ip); }

  //! \brief gamma_inv accessor
  //! \param[in] ip particle index
  // INLINE T &gamma_inv_h(unsigned int ip) { return gamma_inv_.h(ip); }

  //! \brief Ex accessor
  //! \param[in] ip particle index
  INLINE T &Ex_h(unsigned int ip) { return Ex_.h(ip); }

  //! \brief Ey accessor
  //! \param[in] ip particle index
  INLINE T &Ey_h(unsigned int ip) { return Ey_.h(ip); }

  //! \brief Ez accessor
  //! \param[in] ip particle index
  INLINE T &Ez_h(unsigned int ip) { return Ez_.h(ip); }

  //! \brief Bx accessor
  //! \param[in] ip particle index
  INLINE T &Bx_h(unsigned int ip) { return Bx_.h(ip); }

  //! \brief By accessor
  //! \param[in] ip particle index
  INLINE T &By_h(unsigned int ip) { return By_.h(ip); }

  //! \brief Bz accessor
  //! \param[in] ip particle index
  INLINE T &Bz_h(unsigned int ip) { return Bz_.h(ip); }

  //! \brief Gamma accessor using the momentum
  //! \param[in] ip particle index
  INLINE T gamma(unsigned int ip) {
    return sqrt(1 + mx_.h(ip) * mx_.h(ip) + my_.h(ip) * my_.h(ip) + mz_.h(ip) * mz_.h(ip));
  }

  // __________________________________________________________________________
  //
  //! \brief Alloc memory for a new species
  // __________________________________________________________________________
  void allocate(T q, T m, T t, int n_particles, T icv, Backend &backend) {
    inv_cell_volume_m = icv;

    n_particles_m = n_particles;

    // Species properties
    charge_m      = q;
    mass_m        = m;
    temperature_m = t;

    x_.allocate("x", n_particles, backend);
    y_.allocate("y", n_particles, backend);
    z_.allocate("z", n_particles, backend);

    mx_.allocate("mx", n_particles, backend);
    my_.allocate("my", n_particles, backend);
    mz_.allocate("mz", n_particles, backend);

    weight_.allocate("w", n_particles, backend);

    // gamma_inv_.allocate("gamma_inv", n_particles);

    Ex_.allocate("Ex", n_particles, backend);
    Ey_.allocate("Ey", n_particles, backend);
    Ez_.allocate("Ez", n_particles, backend);

    Bx_.allocate("Bx", n_particles, backend);
    By_.allocate("By", n_particles, backend);
    Bz_.allocate("Bz", n_particles, backend);
  }

  // __________________________________________________________________________
  //
  //! \brief Give the number of particles, use std::size
  // __________________________________________________________________________
  unsigned int size() const { return n_particles_m; }

  // __________________________________________________________________________
  //
  //! \brief Delete all particles properties, keep species properties
  // __________________________________________________________________________
  void clear() {

    x_.clear();
    y_.clear();
    z_.clear();

    mx_.clear();
    my_.clear();
    mz_.clear();

    weight_.clear();

    if (with_electromagnetic_fields_) {
      Ex_.clear();
      Ey_.clear();
      Ez_.clear();

      Bx_.clear();
      By_.clear();
      Bz_.clear();
    }

    // if (with_gamma_) {
    //   gamma_inv_.clear();
    // }

    n_particles_m = 0;
  }

  // __________________________________________________________________________
  //
  //! \brief Realloc memory to store particles
  // __________________________________________________________________________
  template <class T_space> void resize(int n_particles, const T_space space) {

    // We resize the vectors only if we can gain substantial memory
    // or need more space
    // A particle costs 112 octets

    // This corresponds to a gain of `min_threshold * 112` octets
    const int min_threshold = 500000;

    if (n_particles > n_particles_m || (n_particles_m - n_particles) > min_threshold) {

      x_.resize(n_particles, 0., space);
      y_.resize(n_particles, 0., space);
      z_.resize(n_particles, 0., space);

      mx_.resize(n_particles, 0., space);
      my_.resize(n_particles, 0., space);
      mz_.resize(n_particles, 0., space);

      weight_.resize(n_particles, 0., space);

      if (with_electromagnetic_fields_) {
        Ex_.resize(n_particles, 0., space);
        Ey_.resize(n_particles, 0., space);
        Ez_.resize(n_particles, 0., space);

        Bx_.resize(n_particles, 0., space);
        By_.resize(n_particles, 0., space);
        Bz_.resize(n_particles, 0., space);
      }

      // if (with_gamma_) {
      //   gamma_inv_.resize(n_particles, 0., space);
      // }
    }

    n_particles_m = n_particles;
  }

  // __________________________________________________________________________
  //
  //! \brief Copy particle at index ip in object `particles` at index i of this
  //! \param[in] i index where to put the particles
  // __________________________________________________________________________
  // void set(int i, Particles &particles, int ip) {

  //   x_h(i) = particles.x_h(ip);
  //   y_h(i) = particles.y_h(ip);
  //   z_h(i) = particles.z_h(ip);

  //   mx_h(i) = particles.mx_h(ip);
  //   my_h(i) = particles.my_h(ip);
  //   mz_h(i) = particles.mz_h(ip);
  //   w_h(i)  = particles.w_h(ip);

  //   if (with_electromagnetic_fields_ && particles.with_electromagnetic_fields_) {
  //     Ex_h(i) = particles.Ex_h(ip);
  //     Ey_h(i) = particles.Ey_h(ip);
  //     Ez_h(i) = particles.Ez_h(ip);

  //     Bx_h(i) = particles.Bx_h(ip);
  //     By_h(i) = particles.By_h(ip);
  //     Bz_h(i) = particles.Bz_h(ip);
  //   }
  // }

  // __________________________________________________________________________
  //
  //! \brief Copy particle at index ip in object `particles` at index i of this
  //! \param[in] i index where to put the particles
  //! \param[in] w weight of the particle to add
  //! \param[in] x position of the particle to add
  //! \param[in] y position of the particle to add
  //! \param[in] z position of the particle to add
  //! \param[in] mx momentum of the particle to add
  //! \param[in] my momentum of the particle to add
  //! \param[in] mz momentum of the particle to add
  // __________________________________________________________________________
  void set(int i, T w, T x, T y, T z, T mx, T my, T mz) {

    weight_[i] = w;

    x_[i] = x;
    y_[i] = y;
    z_[i] = z;

    mx_[i] = mx;
    my_[i] = my;
    mz_[i] = mz;

    // gamma_inv_[i] = 1 / sqrt(1 + mx * mx + my * my + mz * mz);

    if (with_electromagnetic_fields_) {
      Ex_[i] = 0;
      Ey_[i] = 0;
      Ez_[i] = 0;

      Bx_[i] = 0;
      By_[i] = 0;
      Bz_[i] = 0;
    }
  }

  // __________________________________________________________________________
  //
  //! \brief Add particles from the buffer in the structure
  // __________________________________________________________________________
  // void add(Particles &buffer) {

  //   if (buffer.size() > 0) {
  //     int last_ip = size();

  //     resize(last_ip + buffer.size(), minipic::host);

  //     for (int ip = 0; ip < buffer.size(); ip++) {

  //       x_h(last_ip) = buffer.x_h(ip);
  //       y_h(last_ip) = buffer.y_h(ip);
  //       z_h(last_ip) = buffer.z_h(ip);

  //       mx_h(last_ip) = buffer.mx_h(ip);
  //       my_h(last_ip) = buffer.my_h(ip);
  //       mz_h(last_ip) = buffer.mz_h(ip);

  //       w_h(last_ip) = buffer.w_h(ip);

  //       if (with_electromagnetic_fields_ && buffer.with_electromagnetic_fields_) {
  //         Ex_h(last_ip) = buffer.Ex_h(ip);
  //         Ey_h(last_ip) = buffer.Ey_h(ip);
  //         Ez_h(last_ip) = buffer.Ez_h(ip);

  //         Bx_h(last_ip) = buffer.Bx_h(ip);
  //         By_h(last_ip) = buffer.By_h(ip);
  //         Bz_h(last_ip) = buffer.Bz_h(ip);
  //       }

  //       // if (with_gamma_ && buffer.with_gamma_) {
  //       //   gamma_inv_h(last_ip) = buffer.gamma_inv_h(ip);
  //       // }

  //       last_ip++;
  //     }
  //   }
  // }

  // __________________________________________________________________________
  //
  //! \brief Add a single particle
  //! \param[in] w weight of the particle to add
  //! \param[in] x position of the particle to add
  //! \param[in] y position of the particle to add
  //! \param[in] z position of the particle to add
  //! \param[in] px momentum of the particle to add
  //! \param[in] py momentum of the particle to add
  //! \param[in] pz momentum of the particle to add
  // __________________________________________________________________________
  // void add(T w, T x, T y, T z, T px, T py, T pz) {

  //   int last_ip = size();

  //   resize(last_ip + 1, minipic::host);

  //   x_h(last_ip) = x;
  //   y_h(last_ip) = y;
  //   z_h(last_ip) = z;

  //   mx_h(last_ip) = px;
  //   my_h(last_ip) = py;
  //   mz_h(last_ip) = pz;

  //   w_h(last_ip) = w;

  //   if (with_electromagnetic_fields_) {
  //     Ex_h(last_ip) = 0;
  //     Ey_h(last_ip) = 0;
  //     Ez_h(last_ip) = 0;

  //     Bx_h(last_ip) = 0;
  //     By_h(last_ip) = 0;
  //     Bz_h(last_ip) = 0;
  //   }
  // }



  // __________________________________________________________________________
  //
  //! \brief Return the total kinetic energy for this particle species
  // __________________________________________________________________________
  template <class T_space> T get_kinetic_energy(T_space space) {

    T kinetic_energy = 0;

    kinetic_energy = get_kinetic_energy_on_host();

    return kinetic_energy * mass_m;
  }

  // __________________________________________________________________________
  //
  //! \brief data transfer host <-> device
  // __________________________________________________________________________
  template <class T_from, class T_to> void sync(const T_from from, const T_to to) {

    weight_.sync(from, to);

    x_.sync(from, to);
    y_.sync(from, to);
    z_.sync(from, to);

    mx_.sync(from, to);
    my_.sync(from, to);
    mz_.sync(from, to);

    // gamma_inv_.sync(from, to);

    Ex_.sync(from, to);
    Ey_.sync(from, to);
    Ez_.sync(from, to);

    Bx_.sync(from, to);
    By_.sync(from, to);
    Bz_.sync(from, to);
  };

  // __________________________________________________________________________
  //
  //! \brief Print all particles properties
  // __________________________________________________________________________
  void print() {
    for (int ip = 0; ip < n_particles_m; ++ip) {
      std::cerr << "" << ip << " - " << x_h(ip) << " " << y_h(ip) << " " << z_h(ip)
                << " mx: " << mx_h(ip) << " my: " << my_h(ip) << " mz: " << mz_h(ip) << std::endl;
    }
  }

  // __________________________________________________________________________
  //
  //! \brief Check all particles properties
  // __________________________________________________________________________
  void check(T xmin, T xmax, T ymin, T ymax, T zmin, T zmax) {

    for (int ip = 0; ip < n_particles_m; ++ip) {

      if ((x_h(ip) <= xmin) || (x_h(ip) >= xmax) || (y_h(ip) <= ymin) || (y_h(ip) >= ymax) ||
          (z_h(ip) <= zmin) || (z_h(ip) >= zmax)) {
        std::cerr << "Particle: " << ip << "/" << n_particles_m << std::endl;
        std::cerr << " x: " << x_h(ip) << " [" << xmin << " " << xmax << "]" << std::endl;
        std::cerr << " y: " << y_h(ip) << " [" << ymin << " " << ymax << "]" << std::endl;
        std::cerr << " z: " << z_h(ip) << " [" << zmin << " " << zmax << "]" << std::endl;
        std::cerr << " mx: " << mx_h(ip) << " my: " << my_h(ip) << " mz: " << mz_h(ip) << std::endl;
      }
    }
  }

  // __________________________________________________________________________
  //
  //! \brief Print all sums
  // __________________________________________________________________________
  void check_sum() {

    T x_sum = 0;
    T y_sum = 0;
    T z_sum = 0;

    T mx_sum = 0;
    T my_sum = 0;
    T mz_sum = 0;

    // T gamma_inv_sum = 0;

    T Ex_sum = 0;
    T Ey_sum = 0;
    T Ez_sum = 0;

    T Bx_sum = 0;
    T By_sum = 0;
    T Bz_sum = 0;

    for (int ip = 0; ip < n_particles_m; ++ip) {

      x_sum += std::abs(x_h(ip));
      y_sum += std::abs(y_h(ip));
      z_sum += std::abs(z_h(ip));

      mx_sum += std::abs(mx_h(ip));
      my_sum += std::abs(my_h(ip));
      mz_sum += std::abs(mz_h(ip));

      // gamma_inv_sum += std::abs(gamma_inv_h(ip));

      Ex_sum += std::abs(Ex_h(ip));
      Ey_sum += std::abs(Ey_h(ip));
      Ez_sum += std::abs(Ez_h(ip));

      Bx_sum += std::abs(Bx_h(ip));
      By_sum += std::abs(By_h(ip));
      Bz_sum += std::abs(Bz_h(ip));
    }

    std::cerr << std::scientific << std::setprecision(15) << "x sum: " << x_sum
              << " - y sum: " << x_sum << " - z sum: " << x_sum << " - mx sum: " << mx_sum
              << " - my sum: " << my_sum << " - mz sum: "
              << mz_sum
              // << " - gamma inv sum: " << gamma_inv_sum
              << " - Ex: " << Ex_sum << " - Ey: " << Ey_sum << " - Ez: " << Ez_sum
              << " - Bx: " << Bx_sum << " - By: " << By_sum << " - Bz: " << Bz_sum << std::endl;
  }

private:
  // __________________________________________________________________________
  //
  //! \brief Return the total kinetic energy for this particle species
  // __________________________________________________________________________
  T get_kinetic_energy_on_host() {
    T kinetic_energy = 0;

    for (auto ip = 0; ip < size(); ++ip) {
      const T gamma = sqrt(1. + mx_h(ip) * mx_h(ip) + my_h(ip) * my_h(ip) + mz_h(ip) * mz_h(ip));
      kinetic_energy += w_h(ip) * (gamma - 1.);
    }

    return kinetic_energy;
  }
};
