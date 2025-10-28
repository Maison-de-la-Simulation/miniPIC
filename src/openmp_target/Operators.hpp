/* _____________________________________________________________________ */
//! \file Operators.hpp

//! \brief contains generic kernels for the particle pusher

/* _____________________________________________________________________ */

#include "ElectroMagn.hpp"
#include "Patch.hpp"
#include "Profiler.hpp"

#ifndef OPERATORS_H
#define OPERATORS_H

namespace operators {

// ______________________________________________________________________________
//
//! \brief Interpolation operator at the patch level :
//! interpolate EM fields from global grid for each particle
//! \param[in] em  global electromagnetic fields
//! \param[in] patch  patch data structure
// ______________________________________________________________________________
auto interpolate(ElectroMagn &em, Patch &patch) -> void {

  const auto inv_dx_m = em.inv_dx_m;
  const auto inv_dy_m = em.inv_dy_m;
  const auto inv_dz_m = em.inv_dz_m;

  Field<mini_float> &Ex = em.Ex_m;
  Field<mini_float> &Ey = em.Ey_m;
  Field<mini_float> &Ez = em.Ez_m;

  Field<mini_float> &Bx = em.Bx_m;
  Field<mini_float> &By = em.By_m;
  Field<mini_float> &Bz = em.Bz_m;

  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

    const mini_float *const __restrict__ x =
      patch.particles_m[is].x_.get_raw_pointer(minipic::host);
    const mini_float *const __restrict__ y =
      patch.particles_m[is].y_.get_raw_pointer(minipic::host);
    const mini_float *const __restrict__ z =
      patch.particles_m[is].z_.get_raw_pointer(minipic::host);

    mini_float *const __restrict__ Exp = patch.particles_m[is].Ex_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Eyp = patch.particles_m[is].Ey_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Ezp = patch.particles_m[is].Ez_.get_raw_pointer(minipic::host);

    mini_float *const __restrict__ Bxp = patch.particles_m[is].Bx_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Byp = patch.particles_m[is].By_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Bzp = patch.particles_m[is].Bz_.get_raw_pointer(minipic::host);

    // For each particle

#if defined(__MINIPIC_OMP_TARGET__)

#pragma omp target
#pragma omp teams distribute parallel for

#elif defined(__MINIPIC_OPENACC__)

#pragma acc parallel present(x, y, z, Exp, Eyp, Ezp, Bxp, Byp, Bzp, Ex, Ey, Ez, Bx, By, Bz)
#pragma acc loop gang worker vector

#endif

    for (size_t ip = 0; ip < n_particles; ++ip) {

      // printf("Particle %d\n", ip);

      // Calculate normalized positions
      const mini_float ixn = x[ip] * inv_dx_m;
      const mini_float iyn = y[ip] * inv_dy_m;
      const mini_float izn = z[ip] * inv_dz_m;

      // Compute indexes in global primal grid
      const unsigned int ixp = static_cast<unsigned int>(floor(ixn));
      const unsigned int iyp = static_cast<unsigned int>(floor(iyn));
      const unsigned int izp = static_cast<unsigned int>(floor(izn));

      // Compute indexes in global dual grid
      const unsigned int ixd = static_cast<unsigned int>(floor(ixn + 0.5));
      const unsigned int iyd = static_cast<unsigned int>(floor(iyn + 0.5));
      const unsigned int izd = static_cast<unsigned int>(floor(izn + 0.5));

      // Compute interpolation coeff, p = primal, d = dual
      const mini_float coeffs[3] = {ixn + 0.5, iyn, izn};

      // interpolation electric field
      // Ex (d, p , p)
      const auto v00 = Ex(ixd, iyp, izp) * (1 - coeffs[0]) + Ex(ixd + 1, iyp, izp) * coeffs[0];
      const auto v01 =
        Ex(ixd, iyp, izp + 1) * (1 - coeffs[0]) + Ex(ixd + 1, iyp, izp + 1) * coeffs[0];
      const auto v10 =
        Ex(ixd, iyp + 1, izp) * (1 - coeffs[0]) + Ex(ixd + 1, iyp + 1, izp) * coeffs[0];
      const auto v11 =
        Ex(ixd, iyp + 1, izp + 1) * (1 - coeffs[0]) + Ex(ixd + 1, iyp + 1, izp + 1) * coeffs[0];

      Exp[ip] = (v00 * (1 - coeffs[1]) + v10 * coeffs[1]) * (1 - coeffs[2]) +
                (v01 * (1 - coeffs[1]) + v11 * coeffs[1]) * coeffs[2];

      // Ey (p, d, p)
      {
        const mini_float coeffs[3] = {ixn, iyn + 0.5, izn};

        const auto v00 = Ey(ixp, iyd, izp) * (1 - coeffs[0]) + Ey(ixp + 1, iyd, izp) * coeffs[0];
        const auto v01 =
          Ey(ixp, iyd, izp + 1) * (1 - coeffs[0]) + Ey(ixp + 1, iyd, izp + 1) * coeffs[0];
        const auto v10 =
          Ey(ixp, iyd + 1, izp) * (1 - coeffs[0]) + Ey(ixp + 1, iyd + 1, izp) * coeffs[0];
        const auto v11 =
          Ey(ixp, iyd + 1, izp + 1) * (1 - coeffs[0]) + Ey(ixp + 1, iyd + 1, izp + 1) * coeffs[0];
        const auto v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const auto v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Eyp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // Ez (p, p, d)
      {
        const mini_float coeffs[3] = {ixn, iyn, izn + 0.5};

        const auto v00 = Ez(ixp, iyp, izd) * (1 - coeffs[0]) + Ez(ixp + 1, iyp, izd) * coeffs[0];
        const auto v01 =
          Ez(ixp, iyp, izd + 1) * (1 - coeffs[0]) + Ez(ixp + 1, iyp, izd + 1) * coeffs[0];
        const auto v10 =
          Ez(ixp, iyp + 1, izd) * (1 - coeffs[0]) + Ez(ixp + 1, iyp + 1, izd) * coeffs[0];
        const auto v11 =
          Ez(ixp, iyp + 1, izd + 1) * (1 - coeffs[0]) + Ez(ixp + 1, iyp + 1, izd + 1) * coeffs[0];
        const auto v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const auto v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Ezp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // interpolation magnetic field
      // Bx (p, d, d)
      {
        const mini_float coeffs[3] = {ixn, iyn + 0.5, izn + 0.5};

        const auto v00 = Bx(ixp, iyd, izd) * (1 - coeffs[0]) + Bx(ixp + 1, iyd, izd) * coeffs[0];
        const auto v01 =
          Bx(ixp, iyd, izd + 1) * (1 - coeffs[0]) + Bx(ixp + 1, iyd, izd + 1) * coeffs[0];
        const auto v10 =
          Bx(ixp, iyd + 1, izd) * (1 - coeffs[0]) + Bx(ixp + 1, iyd + 1, izd) * coeffs[0];
        const auto v11 =
          Bx(ixp, iyd + 1, izd + 1) * (1 - coeffs[0]) + Bx(ixp + 1, iyd + 1, izd + 1) * coeffs[0];
        const auto v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const auto v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Bxp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // By (d, p, d)
      {
        const mini_float coeffs[3] = {ixn + 0.5, iyn, izn + 0.5};

        const auto v00 = By(ixd, iyp, izd) * (1 - coeffs[0]) + By(ixd + 1, iyp, izd) * coeffs[0];
        const auto v01 =
          By(ixd, iyp, izd + 1) * (1 - coeffs[0]) + By(ixd + 1, iyp, izd + 1) * coeffs[0];
        const auto v10 =
          By(ixd, iyp + 1, izd) * (1 - coeffs[0]) + By(ixd + 1, iyp + 1, izd) * coeffs[0];
        const auto v11 =
          By(ixd, iyp + 1, izd + 1) * (1 - coeffs[0]) + By(ixd + 1, iyp + 1, izd + 1) * coeffs[0];
        const auto v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const auto v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Byp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // Bz (d, d, p)
      {
        const mini_float coeffs[3] = {ixn + 0.5, iyn + 0.5, izn};

        const auto v00 = Bz(ixd, iyd, izp) * (1 - coeffs[0]) + Bz(ixd + 1, iyd, izp) * coeffs[0];
        const auto v01 =
          Bz(ixd, iyd, izp + 1) * (1 - coeffs[0]) + Bz(ixd + 1, iyd, izp + 1) * coeffs[0];
        const auto v10 =
          Bz(ixd, iyd + 1, izp) * (1 - coeffs[0]) + Bz(ixd + 1, iyd + 1, izp) * coeffs[0];
        const auto v11 =
          Bz(ixd, iyd + 1, izp + 1) * (1 - coeffs[0]) + Bz(ixd + 1, iyd + 1, izp + 1) * coeffs[0];
        const auto v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const auto v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Bzp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

    } // end for each particles

  } // Species loop
}

// ______________________________________________________________________________
//
//! \brief Move the particle in the space, compute with EM fields interpolate
//! \param[in] patch  patch data structure
//! \param[in] dt time step to use for the pusher
// ______________________________________________________________________________
auto push(Patch &patch, double dt) -> void {

  // For each species
  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

    // q' = dt * (q/2m)
    const mini_float qp = patch.particles_m[is].charge_m * dt * 0.5 / patch.particles_m[is].mass_m;

    mini_float *const __restrict__ x = patch.particles_m[is].x_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ y = patch.particles_m[is].y_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ z = patch.particles_m[is].z_.get_raw_pointer(minipic::host);

    mini_float *const __restrict__ mx = patch.particles_m[is].mx_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ my = patch.particles_m[is].my_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ mz = patch.particles_m[is].mz_.get_raw_pointer(minipic::host);

    mini_float *const __restrict__ Exp = patch.particles_m[is].Ex_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Eyp = patch.particles_m[is].Ey_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Ezp = patch.particles_m[is].Ez_.get_raw_pointer(minipic::host);

    mini_float *const __restrict__ Bxp = patch.particles_m[is].Bx_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Byp = patch.particles_m[is].By_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Bzp = patch.particles_m[is].Bz_.get_raw_pointer(minipic::host);

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(x, y, z, mx, my, mz, Exp, Eyp, Ezp, Bxp, Byp, Bzp)
#pragma acc loop gang worker vector
#endif
    for (size_t ip = 0; ip < n_particles; ++ip) {

      // 1/2 E
      mini_float px = qp * Exp[ip];
      mini_float py = qp * Eyp[ip];
      mini_float pz = qp * Ezp[ip];

      const mini_float ux = mx[ip] + px;
      const mini_float uy = my[ip] + py;
      const mini_float uz = mz[ip] + pz;

      // gamma-factor
      mini_float gamma_inv = qp / sqrt(1 + (ux * ux + uy * uy + uz * uz));

      // B, T = Transform to rotate the particle
      const mini_float tx  = gamma_inv * Bxp[ip];
      const mini_float ty  = gamma_inv * Byp[ip];
      const mini_float tz  = gamma_inv * Bzp[ip];
      const mini_float tsq = 1. + (tx * tx + ty * ty + tz * tz);
      mini_float tsq_inv   = 1. / tsq;

      px += ((1.0 + tx * tx - ty * ty - tz * tz) * ux + 2.0 * (tx * ty + tz) * uy +
             2.0 * (tz * tx - ty) * uz) *
            tsq_inv;

      py += (2.0 * (tx * ty - tz) * ux + (1.0 - tx * tx + ty * ty - tz * tz) * uy +
             2.0 * (ty * tz + tx) * uz) *
            tsq_inv;

      pz += (2.0 * (tz * tx + ty) * ux + 2.0 * (ty * tz - tx) * uy +
             (1.0 - tx * tx - ty * ty + tz * tz) * uz) *
            tsq_inv;

      // gamma-factor
      gamma_inv = 1 / sqrt(1 + (px * px + py * py + pz * pz));

      // Update momentum
      mx[ip] = px;
      my[ip] = py;
      mz[ip] = pz;

      // Update positions
      x[ip] += px * dt * gamma_inv;
      y[ip] += py * dt * gamma_inv;
      z[ip] += pz * dt * gamma_inv;
    }

  } // Loop on species
}

// ______________________________________________________________________________
//
//! \brief Push only the momentum
//! \param[in] patch  patch data structure
//! \param[in] dt time step to use for the pusher
// ______________________________________________________________________________
auto push_momentum(Patch &patch, double dt) -> void {

  // for each species
  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

    // q' = dt * (q/2m)
    const mini_float qp = patch.particles_m[is].charge_m * dt * 0.5 / patch.particles_m[is].mass_m;

    mini_float *const __restrict__ mx = patch.particles_m[is].mx_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ my = patch.particles_m[is].my_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ mz = patch.particles_m[is].mz_.get_raw_pointer(minipic::host);

    mini_float *const __restrict__ Exp = patch.particles_m[is].Ex_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Eyp = patch.particles_m[is].Ey_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Ezp = patch.particles_m[is].Ez_.get_raw_pointer(minipic::host);

    mini_float *const __restrict__ Bxp = patch.particles_m[is].Bx_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Byp = patch.particles_m[is].By_.get_raw_pointer(minipic::host);
    mini_float *const __restrict__ Bzp = patch.particles_m[is].Bz_.get_raw_pointer(minipic::host);

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(mx, my, mz, Exp, Eyp, Ezp, Bxp, Byp, Bzp)
#pragma acc loop gang worker vector
#endif
    for (auto ip = 0; ip < n_particles; ++ip) {

      // 1/2 E
      mini_float px = qp * Exp[ip];
      mini_float py = qp * Eyp[ip];
      mini_float pz = qp * Ezp[ip];

      const mini_float ux = mx[ip] + px;
      const mini_float uy = my[ip] + py;
      const mini_float uz = mz[ip] + pz;

      // gamma-factor
      mini_float gamma_inv = qp / sqrt(1 + (ux * ux + uy * uy + uz * uz));

      // B, T = Transform to rotate the particle
      const mini_float tx  = gamma_inv * Bxp[ip];
      const mini_float ty  = gamma_inv * Byp[ip];
      const mini_float tz  = gamma_inv * Bzp[ip];
      const mini_float tsq = 1. + (tx * tx + ty * ty + tz * tz);
      mini_float tsq_inv   = 1. / tsq;

      px += ((1.0 + tx * tx - ty * ty - tz * tz) * ux + 2.0 * (tx * ty + tz) * uy +
             2.0 * (tz * tx - ty) * uz) *
            tsq_inv;

      py += (2.0 * (tx * ty - tz) * ux + (1.0 - tx * tx + ty * ty - tz * tz) * uy +
             2.0 * (ty * tz + tx) * uz) *
            tsq_inv;

      pz += (2.0 * (tz * tx + ty) * ux + 2.0 * (ty * tz - tx) * uy +
             (1.0 - tx * tx - ty * ty + tz * tz) * uz) *
            tsq_inv;

      // Update momentum
      mx[ip] = px;
      my[ip] = py;
      mz[ip] = pz;

    } // End for each particles
  } // end for species
}

// _____________________________________________________________________
//
//! \brief Boundaries condition on the particles, periodic
//! or reflect the particles which leave the domain
//
//! \param[in] Params & params - constant global simulation parameters
//! \param[in] Patch & patch - current patch
// _____________________________________________________________________
auto pushBC(Params &params, Patch &patch) -> void {

  if (patch.on_border_m) {

    const auto domain_x_min = params.inf_x, domain_y_min = params.inf_y,
               domain_z_min = params.inf_z;
    const auto domain_x_max = params.sup_x, domain_y_max = params.sup_y,
               domain_z_max = params.sup_z;

    // Periodic conditions
    if (params.boundary_condition_code == 1) {

      const int nx_patch = patch.nx_patchs_m, ny_patch = patch.ny_patchs_m,
                nz_patch = patch.nz_patchs_m;
      const auto Lx = params.Lx, Ly = params.Ly, Lz = params.Lz;

      for (int is = 0; is < patch.n_species_m; is++) {

        size_t n_particles = patch.particles_m[is].size();

        mini_float *const __restrict__ x = patch.particles_m[is].x_.get_raw_pointer(minipic::host);
        mini_float *const __restrict__ y = patch.particles_m[is].y_.get_raw_pointer(minipic::host);
        mini_float *const __restrict__ z = patch.particles_m[is].z_.get_raw_pointer(minipic::host);

#ifdef __MINIPIC_OPENMP_TARGET__
#pragma omp target teams distribute parallel for
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(x, y, z)
#pragma acc loop gang worker vector
#endif

        for (auto ip = 0; ip < n_particles; ++ip) {

          // Only relevant if there is just 1 patch in this direction
          // Else the patch exchange with periodicity is managed in the dedicated function
          if (nx_patch == 1) {
            if (x[ip] >= domain_x_max) {
              x[ip] -= Lx;
            } else if (x[ip] < domain_x_min) {
              x[ip] += Lx;
            }
          }
          // y direction
          if (ny_patch == 1) {
            if (y[ip] >= domain_y_max) {
              y[ip] -= Ly;
            } else if (y[ip] < domain_y_min) {
              y[ip] += Ly;
            }
          }
          // z direction
          if (nz_patch == 1) {
            if (z[ip] >= domain_z_max) {
              z[ip] -= Lz;
            } else if (z[ip] < domain_z_min) {
              z[ip] += Lz;
            }
          }
        } // End loop on particles

      } // End loop on species

      // Reflective conditions
    } else if (params.boundary_condition_code == 2) {
      for (int is = 0; is < patch.n_species_m; is++) {

        size_t n_particles = patch.particles_m[is].size();

        mini_float *const __restrict__ x = patch.particles_m[is].x_.get_raw_pointer(minipic::host);
        mini_float *const __restrict__ y = patch.particles_m[is].y_.get_raw_pointer(minipic::host);
        mini_float *const __restrict__ z = patch.particles_m[is].z_.get_raw_pointer(minipic::host);

        mini_float *const __restrict__ mx =
          patch.particles_m[is].mx_.get_raw_pointer(minipic::host);
        mini_float *const __restrict__ my =
          patch.particles_m[is].my_.get_raw_pointer(minipic::host);
        mini_float *const __restrict__ mz =
          patch.particles_m[is].mz_.get_raw_pointer(minipic::host);

#ifdef __MINIPIC_OPENMP_TARGET__
#pragma omp target teams distribute parallel for
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(x, y, z, mx, my, mz)
#pragma acc loop gang worker vector
#endif
        for (auto ip = 0; ip < n_particles; ++ip) {

          if (x[ip] >= domain_x_max) {
            x[ip]  = 2 * domain_x_max - x[ip];
            mx[ip] = -mx[ip];
          } else if (x[ip] < domain_x_min) {
            x[ip]  = 2 * domain_x_min - x[ip];
            mx[ip] = -mx[ip];
          }

          if (y[ip] >= domain_y_max) {
            y[ip]  = 2 * domain_y_max - y[ip];
            my[ip] = -my[ip];
          } else if (y[ip] < domain_y_min) {
            y[ip]  = 2 * domain_y_min - y[ip];
            my[ip] = -my[ip];
          }

          if (z[ip] >= domain_z_max) {
            z[ip]  = 2 * domain_z_max - z[ip];
            mz[ip] = -mz[ip];
          } else if (z[ip] < domain_z_min) {
            z[ip]  = 2 * domain_z_min - z[ip];
            mz[ip] = -mz[ip];
          }

        } // End loop on particles

      } // End loop on species
    } // if type of conditions
  } // if on border
}

// _____________________________________________________________________
//
//! \brief Boundaries condition on the particles, periodic
//! or reflect the particles which leave the domain
//
//! \param[in] Params & params - constant global simulation parameters
//! \param[in] Patch & patch - current patch
//! \param[in] dt time step to use for the pusher
// _____________________________________________________________________
auto imbalance_operator(Params &params,
                        Patch &patch,
                        int it,
                        std::function<double(double, double, double, double)> func_weight) -> void

{

  // std::srand((unsigned) time(NULL));

  // For each species
  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

    for (size_t ip = 0; ip < n_particles; ++ip) {

      double x = patch.particles_m[is].x_h(ip);
      double y = patch.particles_m[is].y_h(ip);
      double z = patch.particles_m[is].z_h(ip);
      double t = it * params.dt;

      double weight = func_weight(x, y, z, t);

      int wlimit = round(weight);

      for (int i_weight = 0; i_weight < wlimit; i_weight++) {
        double result = std::pow(t, 3) + std::pow(t, 2) - t;
      }

    } // End for each particles

  } // end species loop
} // end function

// _______________________________________________________________________
//
//! \brief Current projection from global particles position to local grid
//! \param params  global simulation parameters
//! \param patch  current patch
// _______________________________________________________________________
auto project(Params &params, Patch &patch) -> void {

  for (int is = 0; is < patch.n_species_m; is++) {

    patch.vec_Jx_m[is].reset(minipic::device);
    patch.vec_Jy_m[is].reset(minipic::device);
    patch.vec_Jz_m[is].reset(minipic::device);

#if defined(__MINIPIC_DEBUG__)
    patch.particles_m[is].sync(minipic::device, minipic::host);
    patch.particles_m[is].check(params.inf_x - params.dx,
                                params.sup_x + params.dx,
                                params.inf_y - params.dy,
                                params.sup_y + params.dy,
                                params.inf_z - params.dz,
                                params.sup_z + params.dz);
    // particles_m[is].print();
    // particles_m[is].check_sum();
#endif

    const size_t n_particles = patch.particles_m[is].size();
    if (n_particles > 0) {

      const mini_float inv_cell_volume_x_q =
        params.inv_cell_volume * patch.particles_m[is].charge_m;
      const mini_float dt = params.dt;

      const mini_float inv_dx = params.inv_dx;
      const mini_float inv_dy = params.inv_dy;
      const mini_float inv_dz = params.inv_dz;

      const mini_float xmin = patch.inf_m[0];
      const mini_float ymin = patch.inf_m[1];
      const mini_float zmin = patch.inf_m[2];

      const int nx_Jx = patch.vec_Jx_m[is].nx_m, ny_Jx = patch.vec_Jx_m[is].ny_m,
                nz_Jx = patch.vec_Jx_m[is].nz_m;
      const int nx_Jy = patch.vec_Jy_m[is].nx_m, ny_Jy = patch.vec_Jy_m[is].ny_m,
                nz_Jy = patch.vec_Jy_m[is].nz_m;
      const int nx_Jz = patch.vec_Jz_m[is].nx_m, ny_Jz = patch.vec_Jz_m[is].ny_m,
                nz_Jz = patch.vec_Jz_m[is].nz_m;

      auto &Jx = patch.vec_Jx_m[is];
      auto &Jy = patch.vec_Jy_m[is];
      auto &Jz = patch.vec_Jz_m[is];

      // mini_float*const __restrict__ Jx = patch.vec_Jx_m[is].get_raw_pointer(minipic::host);
      // mini_float*const __restrict__ Jy = patch.vec_Jy_m[is].get_raw_pointer(minipic::host);
      // mini_float*const __restrict__ Jz = patch.vec_Jz_m[is].get_raw_pointer(minipic::host);

      const mini_float *const __restrict__ w =
        patch.particles_m[is].weight_.get_raw_pointer(minipic::host);

      const mini_float *const __restrict__ x =
        patch.particles_m[is].x_.get_raw_pointer(minipic::host);
      const mini_float *const __restrict__ y =
        patch.particles_m[is].y_.get_raw_pointer(minipic::host);
      const mini_float *const __restrict__ z =
        patch.particles_m[is].z_.get_raw_pointer(minipic::host);

      const mini_float *const __restrict__ mx =
        patch.particles_m[is].mx_.get_raw_pointer(minipic::host);
      const mini_float *const __restrict__ my =
        patch.particles_m[is].my_.get_raw_pointer(minipic::host);
      const mini_float *const __restrict__ mz =
        patch.particles_m[is].mz_.get_raw_pointer(minipic::host);

#ifdef __MINIPIC_OPENMP_TARGET__
#pragma omp target teams distribute parallel for
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(x, y, z, mx, my, mz, w, Jx, Jy, Jz)
#pragma acc loop gang worker vector
#endif
      for (auto ip = 0; ip < n_particles; ++ip) {

        const mini_float gamma_inv =
          1 / sqrt(1 + mx[ip] * mx[ip] + my[ip] * my[ip] + mz[ip] * mz[ip]);

        const mini_float charge_weight = inv_cell_volume_x_q * w[ip];

        const mini_float vx = mx[ip] * gamma_inv;
        const mini_float vy = my[ip] * gamma_inv;
        const mini_float vz = mz[ip] * gamma_inv;

        // Current from the particle
        const mini_float Jxp = vx * charge_weight;
        const mini_float Jyp = vy * charge_weight;
        const mini_float Jzp = vz * charge_weight;

        // Calculate normalized position relative to the patch
        // ixn = (particles_m[is].x(part) ) * params.inv_dx;
        // iyn = (particles_m[is].y(part) ) * params.inv_dy;
        // izn = (particles_m[is].z(part) ) * params.inv_dz;
        const mini_float posxn = (x[ip] - 0.5 * dt * vx - xmin) * inv_dx + 1;
        const mini_float posyn = (y[ip] - 0.5 * dt * vy - ymin) * inv_dy + 1;
        const mini_float poszn = (z[ip] - 0.5 * dt * vz - zmin) * inv_dz + 1;

        // Compute indexes in primal grid
        const int ixp = static_cast<int>(floor(posxn));
        const int iyp = static_cast<int>(floor(posyn));
        const int izp = static_cast<int>(floor(poszn));

        // Compute indexes in dual grid
        // For the current, the dual grid is 0.5 * dx shorter on each side of the grid (if dual
        // directions only)
        const int ixd = static_cast<int>(floor(posxn - 0.5));
        const int iyd = static_cast<int>(floor(posyn - 0.5));
        const int izd = static_cast<int>(floor(poszn - 0.5));

        // Projection particle on currant field
        // Compute interpolation coeff, p = primal, d = dual
        // For Jx
        {
          const mini_float coeffs[3] = {posxn - 0.5 - ixd, posyn - iyp, poszn - izp};

          // Project on Jx using atomicAdd

          minipic::atomicAdd(&Jx(ixd, iyp, izp),
                             (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp);
          minipic::atomicAdd(&Jx(ixd, iyp, izp + 1),
                             (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp);
          minipic::atomicAdd(&Jx(ixd, iyp + 1, izp),
                             (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp);
          minipic::atomicAdd(&Jx(ixd, iyp + 1, izp + 1),
                             (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp);
          minipic::atomicAdd(&Jx(ixd + 1, iyp, izp),
                             (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp);
          minipic::atomicAdd(&Jx(ixd + 1, iyp, izp + 1),
                             (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp);
          minipic::atomicAdd(&Jx(ixd + 1, iyp + 1, izp),
                             (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp);
          minipic::atomicAdd(&Jx(ixd + 1, iyp + 1, izp + 1),
                             (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp);

          // //#pragma acc atomic
          // Jx[ixd * (ny_Jx * nz_Jx) + iyp * (nz_Jx) + izp] += (1 - coeffs[0]) * (1 - coeffs[1]) *
          // (1 - coeffs[2]) * Jxp;
          // //#pragma acc atomic
          // Jx[ixd * (ny_Jx * nz_Jx) + iyp * (nz_Jx) + (izp + 1)] += (1 - coeffs[0]) * (1 -
          // coeffs[1]) * (coeffs[2]) * Jxp;
          // //#pragma acc atomic
          // Jx[ixd * (ny_Jx * nz_Jx) + (iyp + 1) * (nz_Jx) + izp] += (1 - coeffs[0]) * (coeffs[1])
          // * (1 - coeffs[2]) * Jxp;
          // //#pragma acc atomic
          // Jx[ixd * (ny_Jx * nz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)] += (1 - coeffs[0]) *
          // (coeffs[1]) * (coeffs[2]) * Jxp;
          // //#pragma acc atomic
          // Jx[(ixd + 1) * (ny_Jx * nz_Jx) + iyp * (nz_Jx) + izp] += (coeffs[0]) * (1 - coeffs[1])
          // * (1 - coeffs[2]) * Jxp;
          // //#pragma acc atomic
          // Jx[(ixd + 1) * (ny_Jx * nz_Jx) + iyp * (nz_Jx) + (izp + 1)] += (coeffs[0]) * (1 -
          // coeffs[1]) * (coeffs[2]) * Jxp;
          // //#pragma acc atomic
          // Jx[(ixd + 1) * (ny_Jx * nz_Jx) + (iyp + 1) * (nz_Jx) + izp] += (coeffs[0]) *
          // (coeffs[1]) * (1 - coeffs[2]) * Jxp;
          // //#pragma acc atomic
          // Jx[(ixd + 1) * (ny_Jx * nz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)] += (coeffs[0]) *
          // (coeffs[1]) * (coeffs[2]) * Jxp;
        }

        // For Jy
        {
          const mini_float coeffs[3] = {posxn - ixp, posyn - 0.5f - iyd, poszn - izp};

          // Project on Jy using atomicAdd
          minipic::atomicAdd(&Jy(ixp, iyd, izp),
                             (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp);
          minipic::atomicAdd(&Jy(ixp, iyd, izp + 1),
                             (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp);
          minipic::atomicAdd(&Jy(ixp, iyd + 1, izp),
                             (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp);
          minipic::atomicAdd(&Jy(ixp, iyd + 1, izp + 1),
                             (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp);
          minipic::atomicAdd(&Jy(ixp + 1, iyd, izp),
                             (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp);
          minipic::atomicAdd(&Jy(ixp + 1, iyd, izp + 1),
                             (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp);
          minipic::atomicAdd(&Jy(ixp + 1, iyd + 1, izp),
                             (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp);
          minipic::atomicAdd(&Jy(ixp + 1, iyd + 1, izp + 1),
                             (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp);

          // #pragma acc atomic
          //  Jy[ixp * (ny_Jy * nz_Jy) + iyd * (nz_Jy) + izp] += (1 - coeffs[0]) * (1 - coeffs[1]) *
          //  (1 - coeffs[2]) * Jyp;
          //  //#pragma acc atomic
          //  Jy[ixp * (ny_Jy * nz_Jy) + iyd * (nz_Jy) + (izp + 1)] += (1 - coeffs[0]) * (1 -
          //  coeffs[1]) * (coeffs[2]) * Jyp;
          //  //#pragma acc atomic
          //  Jy[ixp * (ny_Jy * nz_Jy) + (iyd + 1) * (nz_Jy) + izp] += (1 - coeffs[0]) * (coeffs[1])
          //  * (1 - coeffs[2]) * Jyp;
          //  //#pragma acc atomic
          //  Jy[ixp * (ny_Jy * nz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)] += (1 - coeffs[0]) *
          //  (coeffs[1]) * (coeffs[2]) * Jyp;
          //  //#pragma acc atomic
          //  Jy[(ixp + 1) * (ny_Jy * nz_Jy) + iyd * (nz_Jy) + izp] += (coeffs[0]) * (1 - coeffs[1])
          //  * (1 - coeffs[2]) * Jyp;
          //  //#pragma acc atomic
          //  Jy[(ixp + 1) * (ny_Jy * nz_Jy) + iyd * (nz_Jy) + (izp + 1)] += (coeffs[0]) * (1 -
          //  coeffs[1]) * (coeffs[2]) * Jyp;
          //  //#pragma acc atomic
          //  Jy[(ixp + 1) * (ny_Jy * nz_Jy) + (iyd + 1) * (nz_Jy) + izp] += (coeffs[0]) *
          //  (coeffs[1]) * (1 - coeffs[2]) * Jyp;
          //  //#pragma acc atomic
          //  Jy[(ixp + 1) * (ny_Jy * nz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)] += (coeffs[0]) *
          //  (coeffs[1]) * (coeffs[2]) * Jyp;
        }

        // For Jz
        {

          const mini_float coeffs[3] = {posxn - ixp, posyn - iyp, poszn - 0.5f - izd};

          // Project on Jz using atomicAdd
          minipic::atomicAdd(&Jz(ixp, iyp, izd),
                             (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp);
          minipic::atomicAdd(&Jz(ixp, iyp, izd + 1),
                             (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp);
          minipic::atomicAdd(&Jz(ixp, iyp + 1, izd),
                             (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp);
          minipic::atomicAdd(&Jz(ixp, iyp + 1, izd + 1),
                             (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp);
          minipic::atomicAdd(&Jz(ixp + 1, iyp, izd),
                             (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp);
          minipic::atomicAdd(&Jz(ixp + 1, iyp, izd + 1),
                             (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp);
          minipic::atomicAdd(&Jz(ixp + 1, iyp + 1, izd),
                             (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp);
          minipic::atomicAdd(&Jz(ixp + 1, iyp + 1, izd + 1),
                             (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp);

          // #pragma acc atomic
          //  Jz[ixp * (ny_Jz * nz_Jz) + iyp * (nz_Jz) + izd] += (1 - coeffs[0]) * (1 - coeffs[1]) *
          //  (1 - coeffs[2]) * Jzp;
          //  //#pragma acc atomic
          //  Jz[ixp * (ny_Jz * nz_Jz) + iyp * (nz_Jz) + (izd + 1)] += (1 - coeffs[0]) * (1 -
          //  coeffs[1]) * (coeffs[2]) * Jzp;
          //  //#pragma acc atomic
          //  Jz[ixp * (ny_Jz * nz_Jz) + (iyp + 1) * (nz_Jz) + izd] += (1 - coeffs[0]) * (coeffs[1])
          //  * (1 - coeffs[2]) * Jzp;
          //  //#pragma acc atomic
          //  Jz[ixp * (ny_Jz * nz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)] += (1 - coeffs[0]) *
          //  (coeffs[1]) * (coeffs[2]) * Jzp;
          //  //#pragma acc atomic
          //  Jz[(ixp + 1) * (ny_Jz * nz_Jz) + iyp * (nz_Jz) + izd] += (coeffs[0]) * (1 - coeffs[1])
          //  * (1 - coeffs[2]) * Jzp;
          //  //#pragma acc atomic
          //  Jz[(ixp + 1) * (ny_Jz * nz_Jz) + iyp * (nz_Jz) + (izd + 1)] += (coeffs[0]) * (1 -
          //  coeffs[1]) * (coeffs[2]) * Jzp;
          //  //#pragma acc atomic
          //  Jz[(ixp + 1) * (ny_Jz * nz_Jz) + (iyp + 1) * (nz_Jz) + izd] += (coeffs[0]) *
          //  (coeffs[1]) * (1 - coeffs[2]) * Jzp;
          //  //#pragma acc atomic
          //  Jz[(ixp + 1) * (ny_Jz * nz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)] += (coeffs[0]) *
          //  (coeffs[1]) * (coeffs[2]) * Jzp;
        }
      } // end for each particles

      patch.projected_[is] = true;

    } else {

      patch.projected_[is] = false;

    } // end if n_particles > 0

  } // end loop species
}

// _______________________________________________________
//
//! \brief Solve Maxwell equations to compute EM fields
//! \param params global parameters
// _______________________________________________________
auto solve_maxwell(const Params &params, ElectroMagn &em, Profiler &profiler) -> void {

  const auto dt         = params.dt;
  const auto dt_over_dx = params.dt * params.inv_dx;
  const auto dt_over_dy = params.dt * params.inv_dy;
  const auto dt_over_dz = params.dt * params.inv_dz;

  const auto nx_Ex = em.Ex_m.nx(), ny_Ex = em.Ex_m.ny(), nz_Ex = em.Ex_m.nz();
  const auto nx_Ey = em.Ey_m.nx(), ny_Ey = em.Ey_m.ny(), nz_Ey = em.Ey_m.nz();
  const auto nx_Ez = em.Ez_m.nx(), ny_Ez = em.Ez_m.ny(), nz_Ez = em.Ez_m.nz();

  const auto nx_Bx = em.Bx_m.nx(), ny_Bx = em.Bx_m.ny(), nz_Bx = em.Bx_m.nz();
  const auto nx_By = em.By_m.nx(), ny_By = em.By_m.ny(), nz_By = em.By_m.nz();
  const auto nx_Bz = em.Bz_m.nx(), ny_Bz = em.Bz_m.ny(), nz_Bz = em.Bz_m.nz();

  Field<mini_float> &Ex = em.Ex_m;
  Field<mini_float> &Ey = em.Ey_m;
  Field<mini_float> &Ez = em.Ez_m;

  Field<mini_float> &Bx = em.Bx_m;
  Field<mini_float> &By = em.By_m;
  Field<mini_float> &Bz = em.Bz_m;

  Field<mini_float> &Jx = em.Jx_m;
  Field<mini_float> &Jy = em.Jy_m;
  Field<mini_float> &Jz = em.Jz_m;

  /////     Solve Maxwell Ampere (E)
  // Electric field Ex (d,p,p)

  profiler.start(MAXWELL);

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Ex, Jx, Bz, By)
#pragma acc loop gang worker vector collapse(3)
#endif
  for (unsigned int ix = 0; ix < nx_Ex; ix++) {
    for (unsigned int iy = 0; iy < ny_Ex; iy++) {
      for (unsigned int iz = 0; iz < nz_Ex; iz++) {

        Ex(ix, iy, iz) += -dt * Jx(ix, iy + 1, iz + 1) +
                          dt_over_dy * (Bz(ix, iy + 1, iz) - Bz(ix, iy, iz)) -
                          dt_over_dz * (By(ix, iy, iz + 1) - By(ix, iy, iz));
      }
    }
  }

  // Electric field Ey (p,d,p)

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Ey, Jy, Bx, Bz)
#pragma acc loop gang worker vector collapse(3)
#endif
  for (unsigned int ix = 0; ix < nx_Ey; ix++) {
    for (unsigned int iy = 0; iy < ny_Ey; iy++) {
      for (unsigned int iz = 0; iz < nz_Ey; iz++) {

        Ey(ix, iy, iz) += -dt * Jy(ix + 1, iy, iz + 1) -
                          dt_over_dx * (Bz(ix + 1, iy, iz) - Bz(ix, iy, iz)) +
                          dt_over_dz * (Bx(ix, iy, iz + 1) - Bx(ix, iy, iz));
      }
    }
  }

  // Electric field Ez (p,p,d)

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Ez, Jz, By, Bx)
#pragma acc loop gang worker vector collapse(3)
#endif
  for (unsigned int ix = 0; ix < nx_Ez; ix++) {
    for (unsigned int iy = 0; iy < ny_Ez; iy++) {
      for (unsigned int iz = 0; iz < nz_Ez; iz++) {

        Ez(ix, iy, iz) += -dt * Jz(ix + 1, iy + 1, iz) +
                          dt_over_dx * (By(ix + 1, iy, iz) - By(ix, iy, iz)) -
                          dt_over_dy * (Bx(ix, iy + 1, iz) - Bx(ix, iy, iz));
      }
    }
  }

  /////     Solve Maxwell Faraday (B)

  // Magnetic field Bx (p,d,d)

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bx, Ey, Ez)
#pragma acc loop gang worker vector collapse(3)
#endif
  for (unsigned int ix = 0; ix < nx_Bx; ix++) {
    for (unsigned int iy = 1; iy < ny_Bx - 1; iy++) {
      for (unsigned int iz = 1; iz < nz_Bx - 1; iz++) {

        Bx(ix, iy, iz) += -dt_over_dy * (Ez(ix, iy, iz) - Ez(ix, iy - 1, iz)) +
                          dt_over_dz * (Ey(ix, iy, iz) - Ey(ix, iy, iz - 1));
      }
    }
  }

  // Magnetic field By (d,p,d)

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(By, Ex, Ez)
#pragma acc loop gang worker vector collapse(3)
#endif
  for (unsigned int ix = 1; ix < nx_By - 1; ix++) {
    for (unsigned int iy = 0; iy < ny_By; iy++) {
      for (unsigned int iz = 1; iz < nz_By - 1; iz++) {
        By(ix, iy, iz) += -dt_over_dz * (Ex(ix, iy, iz) - Ex(ix, iy, iz - 1)) +
                          dt_over_dx * (Ez(ix, iy, iz) - Ez(ix - 1, iy, iz));
      }
    }
  }

  // Magnetic field Bz (d,d,p)

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bz, Ey, Ex)
#pragma acc loop gang worker vector collapse(3)
#endif
  for (unsigned int ix = 1; ix < nx_Bz - 1; ix++) {
    for (unsigned int iy = 1; iy < ny_Bz - 1; iy++) {
      for (unsigned int iz = 0; iz < nz_Bz; iz++) {
        Bz(ix, iy, iz) += -dt_over_dx * (Ey(ix, iy, iz) - Ey(ix - 1, iy, iz)) +
                          dt_over_dy * (Ex(ix, iy, iz) - Ex(ix, iy - 1, iz));
      }
    }
  }

  profiler.stop();

} // end solve

// _______________________________________________________________
//
//! \brief Boundaries condition on the global grid
//! \param[in] Params & params - global constant parameters
//! \param[in] ElectroMagn & em - global electromagnetic fields
// _______________________________________________________________
auto currentBC(Params &params, ElectroMagn &em) -> void {

  if (params.boundary_condition == "periodic") {

    Field<mini_float> &Jx = em.Jx_m;
    Field<mini_float> &Jy = em.Jy_m;
    Field<mini_float> &Jz = em.Jz_m;

    const auto nx_Jx = em.Jx_m.nx();
    const auto ny_Jx = em.Jx_m.ny();
    const auto nz_Jx = em.Jx_m.nz();

    const auto nx_Jy = em.Jy_m.nx();
    const auto ny_Jy = em.Jy_m.ny();
    const auto nz_Jy = em.Jy_m.nz();

    const auto nx_Jz = em.Jz_m.nx();
    const auto ny_Jz = em.Jz_m.ny();
    const auto nz_Jz = em.Jz_m.nz();

    // X

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jx)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int iy = 0; iy < ny_Jx; ++iy) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Jx; ++iz) {

        Jx(0, iy, iz) += Jx(nx_Jx - 2, iy, iz);
        Jx(nx_Jx - 2, iy, iz) = Jx(0, iy, iz);

        Jx(1, iy, iz) += Jx(nx_Jx - 1, iy, iz);
        Jx(nx_Jx - 1, iy, iz) = Jx(1, iy, iz);
      }
    }

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jy)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int iy = 0; iy < ny_Jy; ++iy) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Jy; ++iz) {

        Jy(0, iy, iz) += Jy(nx_Jy - 2, iy, iz);
        Jy(nx_Jy - 2, iy, iz) = Jy(0, iy, iz);

        Jy(1, iy, iz) += Jy(nx_Jy - 1, iy, iz);
        Jy(nx_Jy - 1, iy, iz) = Jy(1, iy, iz);
      }
    }

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jz)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int iy = 0; iy < ny_Jz; ++iy) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Jz; ++iz) {
        Jz(0, iy, iz) += Jz(nx_Jz - 2, iy, iz);
        Jz(nx_Jz - 2, iy, iz) = Jz(0, iy, iz);

        Jz(1, iy, iz) += Jz(nx_Jz - 1, iy, iz);
        Jz(nx_Jz - 1, iy, iz) = Jz(1, iy, iz);
      }
    }

    // Y
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jx)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Jx; ++ix) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Jx; ++iz) {

        Jx(ix, 0, iz) += Jx(ix, ny_Jx - 2, iz);
        Jx(ix, ny_Jx - 2, iz) = Jx(ix, 0, iz);

        Jx(ix, 1, iz) += Jx(ix, ny_Jx - 1, iz);
        Jx(ix, ny_Jx - 1, iz) = Jx(ix, 1, iz);
      }
    }

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jy)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Jy; ++ix) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Jy; ++iz) {

        Jy(ix, 0, iz) += Jy(ix, ny_Jy - 2, iz);
        Jy(ix, ny_Jy - 2, iz) = Jy(ix, 0, iz);

        Jy(ix, 1, iz) += Jy(ix, ny_Jy - 1, iz);
        Jy(ix, ny_Jy - 1, iz) = Jy(ix, 1, iz);
      }
    }

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jz)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Jz; ++ix) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Jz; ++iz) {

        Jz(ix, 0, iz) += Jz(ix, ny_Jz - 2, iz);
        Jz(ix, ny_Jz - 2, iz) = Jz(ix, 0, iz);

        Jz(ix, 1, iz) += Jz(ix, ny_Jz - 1, iz);
        Jz(ix, ny_Jz - 1, iz) = Jz(ix, 1, iz);
      }
    }

    // Z
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jx)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Jx; ++ix) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iy = 0; iy < ny_Jx; ++iy) {

        Jx(ix, iy, 0) += Jx(ix, iy, nz_Jx - 2);
        Jx(ix, iy, nz_Jx - 2) = Jx(ix, iy, 0);

        Jx(ix, iy, 1) += Jx(ix, iy, nz_Jx - 1);
        Jx(ix, iy, nz_Jx - 1) = Jx(ix, iy, 1);
      }
    }

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jy)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Jy; ++ix) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iy = 0; iy < ny_Jy; ++iy) {

        Jy(ix, iy, 0) += Jy(ix, iy, nz_Jy - 2);
        Jy(ix, iy, nz_Jy - 2) = Jy(ix, iy, 0);

        Jy(ix, iy, 1) += Jy(ix, iy, nz_Jy - 1);
        Jy(ix, iy, nz_Jy - 1) = Jy(ix, iy, 1);
      }
    }

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jz)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Jz; ++ix) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iy = 0; iy < ny_Jz; ++iy) {

        Jz(ix, iy, 0) += Jz(ix, iy, nz_Jz - 2);
        Jz(ix, iy, nz_Jz - 2) = Jz(ix, iy, 0);

        Jz(ix, iy, 1) += Jz(ix, iy, nz_Jz - 1);
        Jz(ix, iy, nz_Jz - 1) = Jz(ix, iy, 1);
      }
    }

  } // end if periodic
} // end currentBC

// _______________________________________________________________
//
//! \brief Boundaries condition on the global grid
//! \param[in] Params & params - global constant parameters
// _______________________________________________________________
auto solveBC(Params &params, ElectroMagn &em) -> void {

  const auto nx_Bx = em.Bx_m.nx();
  const auto ny_Bx = em.Bx_m.ny();
  const auto nz_Bx = em.Bx_m.nz();

  const auto nx_By = em.By_m.nx();
  const auto ny_By = em.By_m.ny();
  const auto nz_By = em.By_m.nz();

  const auto nx_Bz = em.Bz_m.nx();
  const auto ny_Bz = em.Bz_m.ny();
  const auto nz_Bz = em.Bz_m.nz();

  auto &Bx = em.Bx_m;
  auto &By = em.By_m;
  auto &Bz = em.Bz_m;

  if (params.boundary_condition == "periodic") {

    // X dim
    // By (d,p,d)
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(By)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int iy = 0; iy < ny_By; ++iy) {
      for (unsigned int iz = 0; iz < nz_By; ++iz) {
        // -X
        By(0, iy, iz)         = By(nx_By - 2, iy, iz);
        By(nx_By - 1, iy, iz) = By(1, iy, iz);
      }
    }

    // Bz (d,d,p)
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bz)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int iy = 0; iy < ny_Bz; iy++) {
      for (unsigned int iz = 0; iz < nz_Bz; iz++) {
        // -X
        Bz(0, iy, iz)         = Bz(nx_Bz - 2, iy, iz);
        Bz(nx_Bz - 1, iy, iz) = Bz(1, iy, iz);
      }
    }

    // Y dim
    // Bx (p,d,d)
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bx)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Bx; ix++) {
      for (unsigned int iz = 0; iz < nz_Bx; iz++) {
        // -Y
        Bx(ix, 0, iz)         = Bx(ix, ny_Bx - 2, iz);
        Bx(ix, ny_Bx - 1, iz) = Bx(ix, 1, iz);
      }
    }
    // Bz (d,d,p)
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bz)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Bz; ix++) {
      for (unsigned int iz = 0; iz < nz_Bz; iz++) {
        // -Y
        Bz(ix, 0, iz)         = Bz(ix, ny_Bz - 2, iz);
        Bz(ix, ny_Bz - 1, iz) = Bz(ix, 1, iz);
      }
    }

    // Z dim
    // Bx
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bx)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Bx; ix++) {
      for (unsigned int iy = 0; iy < ny_Bx; iy++) {
        // -Z
        Bx(ix, iy, 0) = Bx(ix, iy, nz_Bx - 2);
        // +Z
        Bx(ix, iy, nz_Bx - 1) = Bx(ix, iy, 1);
      }
    }

    // By
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(By)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_By; ++ix) {
      for (unsigned int iy = 0; iy < ny_By; ++iy) {
        // -Z
        By(ix, iy, 0) = By(ix, iy, nz_By - 2);
        // +Z
        By(ix, iy, nz_By - 1) = By(ix, iy, 1);
      }
    }

  } else if (params.boundary_condition == "reflective") {

    // X dim
    // By (d,p,d)
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(By)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int iy = 0; iy < ny_By; ++iy) {
      for (unsigned int iz = 0; iz < nz_By; ++iz) {
        // -X
        By(0, iy, iz) = By(1, iy, iz);
        // +X
        By(nx_By - 1, iy, iz) = By(nx_By - 2, iy, iz);
      }
    }

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bz)
#pragma acc loop gang worker collapse(2)
#endif
    // Bz (d,d,p)
    for (unsigned int iy = 0; iy < ny_Bz; iy++) {
      for (unsigned int iz = 0; iz < nz_Bz; iz++) {
        // -X
        Bz(0, iy, iz) = Bz(1, iy, iz);
        // +X
        Bz(nx_Bz - 1, iy, iz) = Bz(nx_Bz - 2, iy, iz);
      }
    }

    // Y dim
    // Bx (p,d,d)
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bx)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Bx; ix++) {
      for (unsigned int iz = 0; iz < nz_Bx; iz++) {
        // -Y
        Bx(ix, 0, iz) = Bx(ix, 1, iz);
        // +Y
        Bx(ix, ny_Bx - 1, iz) = Bx(ix, ny_Bx - 2, iz);
      }
    }
    // Bz (-1 to avoid corner)
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bz)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Bz; ix++) {
      for (unsigned int iz = 0; iz < nz_Bz; ++iz) {
        // -Y
        Bz(ix, 0, iz) = Bz(ix, 1, iz);
        // +Y
        Bz(ix, ny_Bz - 1, iz) = Bz(ix, ny_Bz - 2, iz);
      }
    }

    // Z dim
    // Bx
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Bx)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_Bx; ix++) {
      for (unsigned int iy = 0; iy < ny_Bx; iy++) {
        // -Z
        Bx(ix, iy, 0) = Bx(ix, iy, 1);
        // +Z
        Bx(ix, iy, nz_Bx - 1) = Bx(ix, iy, nz_Bx - 2);
      }
    }
    // By
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(2)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(By)
#pragma acc loop gang worker collapse(2)
#endif
    for (unsigned int ix = 0; ix < nx_By; ix++) {
      for (unsigned int iy = 0; iy < ny_By; iy++) {
        // -Z
        By(ix, iy, 0) = By(ix, iy, 1);
        // +Z
        By(ix, iy, nz_By - 1) = By(ix, iy, nz_By - 2);
      }
    }

  } // End if
} // End solveBC

// ______________________________________________________
//
//! \brief This function tags particles which leave the patch,
//!        then, put them in  buffers according to the communication direction
//!        and finally delete them from the patch
//! \param[in] Params & params - global constant parameters
//! \param[in] Patch & patch - current patch
// ______________________________________________________
void identify_particles_to_move(Params &params, Patch &patch, Backend &backend) {

  const mini_float patch_inf[3] = {patch.inf_m[0], patch.inf_m[1], patch.inf_m[2]};
  const mini_float patch_sup[3] = {patch.sup_m[0], patch.sup_m[1], patch.sup_m[2]};

  const mini_float inf[3]    = {params.inf_x, params.inf_y, params.inf_z};
  const mini_float sup[3]    = {params.sup_x, params.sup_y, params.sup_z};
  const mini_float length[3] = {params.Lx, params.Ly, params.Lz};

  const int boundary_condition_code = params.boundary_condition_code;

  // Reset buffers
  for (int is = 0; is < patch.n_species_m; is++) {
    for (int ib = 0; ib < 26; ib++) {
      patch.particles_to_move_m[is][ib].clear();
    }
  }

  // Identify particles to move
  for (int is = 0; is < patch.n_species_m; is++) {

    // Number of particles for this species is
    size_t n_particles = patch.particles_m[is].size();

    if (n_particles == 0)
      continue;

    // Count number of particles to move
    Vector<size_t> n_particle_to_move(26, 0, backend);

    // index for particle copy in the buffers
    Vector<size_t> ip_to_move(26, 0, backend);

    // Mask for buffer direction
    // std::Vector<int> masks(n_particles, -1, backend);
    Vector<int> masks(n_particles, -1, backend);

    Vector<mini_float> &w = patch.particles_m[is].weight_;

    Vector<mini_float> &x = patch.particles_m[is].x_;
    Vector<mini_float> &y = patch.particles_m[is].y_;
    Vector<mini_float> &z = patch.particles_m[is].z_;

    Vector<mini_float> &mx = patch.particles_m[is].mx_;
    Vector<mini_float> &my = patch.particles_m[is].my_;
    Vector<mini_float> &mz = patch.particles_m[is].mz_;

    Vector<int> &masks_accessor = masks;

    // 1 - Compute number of particles to move per buffer and tag them
    for (size_t ip = 0; ip < n_particles; ip++) {

      mini_float shift[3];

      // Compute in which direction the particle goes

      if (x(ip) < patch_inf[0]) {
        shift[0] = -1;
      } else if (x(ip) >= patch_sup[0]) {
        shift[0] = 1;
      } else {
        shift[0] = 0;
      }

      if (y(ip) < patch_inf[1]) {
        shift[1] = -1;
      } else if (y(ip) >= patch_sup[1]) {
        shift[1] = 1;
      } else {
        shift[1] = 0;
      }

      if (z(ip) < patch_inf[2]) {
        shift[2] = -1;
      } else if (z(ip) >= patch_sup[2]) {
        shift[2] = 1;
      } else {
        shift[2] = 0;
      }

      // Tag and count the particles to move
      if (!(shift[0] == 0 && shift[1] == 0 && shift[2] == 0)) {
        int ib = static_cast<int>((shift[0] + 1) * 9 + (shift[1] + 1) * 3 + (shift[2] + 1));
        if (ib > 13)
          ib--;

        // If periodic conditions :
        // We need to update the new position of the particle after identification of the buffers
        if (boundary_condition_code == 1) {

          if (x(ip) >= sup[0]) {
            x(ip) -= length[0];
          } else if (x(ip) < inf[0]) {
            x(ip) += length[0];
          }

          if (y(ip) >= sup[1]) {
            y(ip) -= length[1];
          } else if (y(ip) < inf[1]) {
            y(ip) += length[1];
          }

          if (z(ip) >= sup[2]) {
            z(ip) -= length[2];
          } else if (z(ip) < inf[2]) {
            z(ip) += length[2];
          }
        }

        n_particle_to_move(ib) += 1;

        // we store here the buffer id to use it later
        masks(ip) = ib;
      }
    } // end for particles

    // 2 - Realloc buffers memory

    size_t total_particles_to_remove = 0;

    for (int ib = 0; ib < 26; ib++) {
      patch.particles_to_move_m[is][ib].resize(n_particle_to_move.h(ib), minipic::device);
      total_particles_to_remove += n_particle_to_move.h(ib);

      //      if (n_particle_to_move.h(ib) > 0) {
      //         std::cerr << ib << " " << n_particle_to_move.h(ib) << std::endl;
      //      }
    }

    // 3 -  Move tagged particles in the corresponding buffer

    if (total_particles_to_remove > 0) {

      for (size_t ip = 0; ip < n_particles; ++ip) {

        if (masks.h(ip) >= 0) {

          // if (!(shift[0] == 0 && shift[1] == 0 && shift[2] == 0)) {
          const int ib = masks.h(ip);

          // for (int d = 0; d < 3; d++) {
          //   if (pos[d] < inf_m[d]) {
          //     shift[d] = -1;
          //   } else if (pos[d] >= sup_m[d]) {
          //     shift[d] = 1;
          //   } else {
          //     shift[d] = 0;
          //   }
          // }

          const size_t i = ip_to_move.h(ib);

          patch.particles_to_move_m[is][ib].x_h(i) = patch.particles_m[is].x_h(ip);
          patch.particles_to_move_m[is][ib].y_h(i) = patch.particles_m[is].y_h(ip);
          patch.particles_to_move_m[is][ib].z_h(i) = patch.particles_m[is].z_h(ip);

          patch.particles_to_move_m[is][ib].mx_h(i) = patch.particles_m[is].mx_h(ip);
          patch.particles_to_move_m[is][ib].my_h(i) = patch.particles_m[is].my_h(ip);
          patch.particles_to_move_m[is][ib].mz_h(i) = patch.particles_m[is].mz_h(ip);

          patch.particles_to_move_m[is][ib].w_h(i) = patch.particles_m[is].w_h(ip);

          ip_to_move.h(ib)++;
        }
      } // end for particles

      // for (int ib = 0; ib < 26; ib++) {
      //   if (particles_to_move_m[is][ib].size() > 0) {
      //     particles_to_move_m[is][ib].copy_host_to_device();
      //   }
      // }

    } // end if total_particles_to_remove

    // 4 - Move particles to remove at the end of the vector

    if (total_particles_to_remove > 0) {

      // front particle index
      long ip = 0;

      // last particle index
      long last_ip = n_particles - 1;

      while (ip <= last_ip) {

        // back particle left
        if (masks_accessor(last_ip) >= 0) {
          last_ip--;
          continue;
        }

        // Front particle left :
        // if the mask value > 0,
        // then the corresponding index is available for a particle
        if (masks_accessor(ip) >= 0) {
          // Copy particle last_ip at ip index

          x(ip) = x(last_ip);
          y(ip) = y(last_ip);
          z(ip) = z(last_ip);

          mx(ip) = mx(last_ip);
          my(ip) = my(last_ip);
          mz(ip) = mz(last_ip);

          w(ip) = w(last_ip);

          last_ip--;
          ip++;
          // else front particle stay, check next one
        } else {
          ip++;
        }
      }

      // Delete tagged particles by resizing particles_m[is]
      patch.particles_m[is].resize(n_particles - total_particles_to_remove, minipic::device);

    } // end if total_particles_to_remove > 0

    // std::cerr << "patch: " << idx_patch_topology_m << " sp: " << is << " - after erase: " <<
    // particles_m[is].get_kinetic_energy() << " size: "  << particles_m[is].size() << std::endl;
    // particles_m[is].print();

  } // end for species

  // erase_particles();
}

// ___________________________________________________________
//
//! \brief Get the particles from neighbors
//! \param[in] params constant global parameters
//! \param[in] vec_patch vector of all patches
// ___________________________________________________________
auto exchange_particles(Params &params, std::vector<Patch> &vec_patch, int id_patch) -> void {

  // Current patch to handle
  Patch &patch = vec_patch[id_patch];

  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t number_of_particles = patch.particles_m[is].size();

    // total number of particles coming from other patches
    size_t coming_number_of_particles = 0;

    // Compute the total number of particles that will come from other patches
    for (int i = -1; i < 2; i++) {
      for (int j = -1; j < 2; j++) {
        for (int k = -1; k < 2; k++) {

          // Buffer if where to get the coming particles in my neighbor
          int idx_buffer = (i * -1 + 1) * 9 + (j * -1 + 1) * 3 + (k * -1 + 1);

          // 13 eq. i=j=k=0
          if (idx_buffer == 13) {
            continue;
          }

          // id of the Neighbor in vec_patch
          const int idx_neighbor = params.get_patch_index(patch.i_patch_topology_m + i,
                                                          patch.j_patch_topology_m + j,
                                                          patch.k_patch_topology_m + k);

          // 13 eq. i=j=k=0
          if (idx_buffer > 13) {
            idx_buffer--;
          }

          coming_number_of_particles +=
            vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].size();

        } // end for each neighbors
      } // end for each neighbors
    } // end for each neighbors

    if (coming_number_of_particles > 0) {
      patch.particles_m[is].resize(number_of_particles + coming_number_of_particles,
                                   minipic::device);
    }

    // Index where to start to copy the incoming particles in Particles
    size_t ip_buffer_start = number_of_particles;

    // Collect new particles from neighbours
    for (int i = -1; i < 2; i++) {
      for (int j = -1; j < 2; j++) {
        for (int k = -1; k < 2; k++) {

          // Buffer if where to get the coming particles in my neighbor
          int idx_buffer = (i * -1 + 1) * 9 + (j * -1 + 1) * 3 + (k * -1 + 1);
          // 13 eq. i=j=k=0
          if (idx_buffer == 13) {
            continue;
          }

          // id of the Neighbor in vec_patch
          int idx_neighbor = params.get_patch_index(patch.i_patch_topology_m + i,
                                                    patch.j_patch_topology_m + j,
                                                    patch.k_patch_topology_m + k);

          // 13 eq. i=j=k=0
          if (idx_buffer > 13) {
            idx_buffer--;
          }

          const int buffer_size =
            vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].size();

          if (buffer_size > 0) {

#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
            for (size_t ip = 0; ip < buffer_size; ++ip) {

              patch.particles_m[is].x_.h(ip_buffer_start + ip) =
                vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].x_.h(ip);
              patch.particles_m[is].y_.h(ip_buffer_start + ip) =
                vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].y_.h(ip);
              patch.particles_m[is].z_.h(ip_buffer_start + ip) =
                vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].z_.h(ip);

              patch.particles_m[is].mx_.h(ip_buffer_start + ip) =
                vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].mx_.h(ip);
              patch.particles_m[is].my_.h(ip_buffer_start + ip) =
                vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].my_.h(ip);
              patch.particles_m[is].mz_.h(ip_buffer_start + ip) =
                vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].mz_.h(ip);

              patch.particles_m[is].weight_.h(ip_buffer_start + ip) =
                vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].weight_.h(ip);

              // Ex.h(last_ip) = buffer.Ex_h(ip);
              // Ey.h(last_ip) = buffer.Ey_h(ip);
              // Ez.h(last_ip) = buffer.Ez_h(ip);

              // Bx.h(last_ip) = buffer.Bx_h(ip);
              // By.h(last_ip) = buffer.By_h(ip);
              // Bz/h(last_ip) = buffer.Bz_h(ip);

            } // end ip loop

            // We check that there are particles to move inside the function `add`
            // particles_m[is].add(vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer]);

            ip_buffer_start += buffer_size;

          } // if buffer size > 0

        } // end for each neighbors
      } // end for each neighbors
    } // end for each neighbors

    // std::cerr << "end exchange" << std::endl;

    // std::cerr << "patch: " << idx_patch_topology_m << " sp: " << is << " - after exchange: " <<
    // particles_m[is].get_kinetic_energy() << " size: "  << particles_m[is].size() << std::endl;

#if defined(__MINIPIC_DEBUG__)
    patch.particles_m[is].sync(minipic::device, minipic::host);
    // particles_m[is].check(params.inf_x, params.sup_x, params.inf_y, params.sup_y, params.inf_z,
    // params.sup_z);
    patch.particles_m[is].print();
    // particles_m[is].check_sum();
#endif

  } // end for species
}

// ______________________________________________________
//
//! \brief Sum all species local current grids in local grid
//! \param[in] patch  current patch to handle
// ______________________________________________________
auto reduc_current(Patch &patch) -> void {

  for (int is = 1; is < patch.n_species_m; is++) {

    // Only if particles projected
    if (patch.projected_[is]) {

      // Jx
      {
        auto nx = patch.vec_Jx_m[0].nx();
        auto ny = patch.vec_Jx_m[0].ny();
        auto nz = patch.vec_Jx_m[0].nz();

        auto &Jx0 = patch.vec_Jx_m[0];
        auto &Jxs = patch.vec_Jx_m[is];
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jx0, Jxs)
#pragma acc loop gang worker collapse(3)
#endif
        for (int ix = 0; ix < nx; ix++) {
          for (int iy = 0; iy < ny; iy++) {
            for (int iz = 0; iz < nz; iz++) {
              Jx0(ix, iy, iz) += Jxs(ix, iy, iz);
            }
          }
        }
      }

      // Jy
      {
        auto nx = patch.vec_Jy_m[0].nx();
        auto ny = patch.vec_Jy_m[0].ny();
        auto nz = patch.vec_Jy_m[0].nz();

        auto &Jy0 = patch.vec_Jy_m[0];
        auto &Jys = patch.vec_Jy_m[is];

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jy0, Jys)
#pragma acc loop gang worker collapse(3)
#endif
        for (int ix = 0; ix < nx; ix++) {
          for (int iy = 0; iy < ny; iy++) {
            for (int iz = 0; iz < nz; iz++) {
              Jy0(ix, iy, iz) += Jys(ix, iy, iz);
            }
          }
        }
      }

      // Jz
      {
        auto nx = patch.vec_Jz_m[0].nx();
        auto ny = patch.vec_Jz_m[0].ny();
        auto nz = patch.vec_Jz_m[0].nz();

        auto &Jz0 = patch.vec_Jz_m[0];
        auto &Jzs = patch.vec_Jz_m[is];
#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jz0, Jzs)
#pragma acc loop gang worker collapse(3)
#endif
        for (int ix = 0; ix < nx; ix++) {
          for (int iy = 0; iy < ny; iy++) {
            for (int iz = 0; iz < nz; iz++) {
              Jz0(ix, iy, iz) += Jzs(ix, iy, iz);
            }
          }
        }
      }

    } // end check if particles
  } // end for species
}

// ____________________________________________________________________________
//! \brief Copy all local current grid in the global grid
//! \param[in] ElectroMagn & em - global electromagnetic fields
//! \param[in] Patch & patch - current patch
// ____________________________________________________________________________
auto local2global(ElectroMagn &em, Patch &patch) -> void {

  bool projected = false;

  for (int is = 0; is < patch.n_species_m; is++) {
    projected = projected || patch.projected_[is];
  }

  // projection only if particles in this patch
  if (projected) {

    // for (int is = 0; is < n_species_m; is++) {
    const int i_global_p = patch.ix_origin_m;
    const int j_global_p = patch.iy_origin_m;
    const int k_global_p = patch.iz_origin_m;
    const int i_global_d = patch.ix_origin_m;
    const int j_global_d = patch.iy_origin_m;
    const int k_global_d = patch.iz_origin_m;

    // Jx
    {
      auto nx = patch.vec_Jx_m[0].nx();
      auto ny = patch.vec_Jx_m[0].ny();
      auto nz = patch.vec_Jx_m[0].nz();

      auto &Jx  = em.Jx_m;
      auto &Jxs = patch.vec_Jx_m[0];

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jx, Jxs)
#pragma acc loop gang worker collapse(3)
#endif
      for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
          for (int iz = 0; iz < nz; iz++) {
            Jx(i_global_d + ix, j_global_p + iy, k_global_p + iz) += Jxs(ix, iy, iz);
          }
        }
      }
    }

    // Jy
    {
      auto nx = patch.vec_Jy_m[0].nx();
      auto ny = patch.vec_Jy_m[0].ny();
      auto nz = patch.vec_Jy_m[0].nz();

      auto &Jy  = em.Jy_m;
      auto &Jys = patch.vec_Jy_m[0];

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jy, Jys)
#pragma acc loop gang worker collapse(3)
#endif
      for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
          for (int iz = 0; iz < nz; iz++) {
            Jy(i_global_p + ix, j_global_d + iy, k_global_p + iz) += Jys(ix, iy, iz);
          }
        }
      }
    }

    // Jz
    {
      auto nx = patch.vec_Jz_m[0].nx();
      auto ny = patch.vec_Jz_m[0].ny();
      auto nz = patch.vec_Jz_m[0].nz();

      auto &Jz  = em.Jz_m;
      auto &Jzs = patch.vec_Jz_m[0];

#if defined(__MINIPIC_OPENMP_TARGET__)
#pragma omp target teams distribute parallel for collapse(3)
#elif defined(__MINIPIC_OPENACC__)
#pragma acc parallel present(Jz, Jzs)
#pragma acc loop gang worker collapse(3)
#endif
      for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
          for (int iz = 0; iz < nz; iz++) {
            Jz(i_global_p + ix, j_global_p + iy, k_global_d + iz) += Jzs(ix, iy, iz);
          }
        }
      }
    }

  } // end if total_particles
}

// ____________________________________________________________________________
//
//! \brief Emit a laser field in the x direction using an antenna
//! \param[in] Params & params - global constant parameters
//! \param[in] profile - (std::function<double(double y, double z, double t)>) profile of the
//! antenna
//! \param[in] x - (double) position of the antenna \param[in] double t - (double) current
//! time
// ____________________________________________________________________________
auto antenna(Params &params,
             ElectroMagn &em,
             std::function<double(double, double, double)> profile,
             double x,
             double t) -> void {

  em.Jz_m.sync(minipic::device, minipic::host);

  Field<mini_float> *J = &em.Jz_m;

  const int ix = floor((x - params.inf_x - J->dual_x_m * 0.5 * params.dx) / params.dx);

  const double yfs = 0.5 * params.Ly + params.inf_y;
  const double zfs = 0.5 * params.Lz + params.inf_z;

  for (unsigned int iy = 0; iy < J->ny_m; ++iy) {
    for (unsigned int iz = 0; iz < J->nz_m; ++iz) {

      const double y = (iy - J->dual_y_m * 0.5) * params.dy + params.inf_y - yfs;
      const double z = (iz - J->dual_z_m * 0.5) * params.dz + params.inf_z - zfs;

      (*J)(ix, iy, iz) = profile(y, z, t);
    }
  }

  em.Jz_m.sync(minipic::host, minipic::device);
}

} // end namespace operators

#endif // OPERATORS_H