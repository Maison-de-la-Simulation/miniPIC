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
auto interpolate(ElectroMagn &em, Patch &patch, Backend &backend) -> void {

  const auto inv_dx_m = em.inv_dx_m;
  const auto inv_dy_m = em.inv_dy_m;
  const auto inv_dz_m = em.inv_dz_m;

  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

    if (n_particles == 0)
      continue;

    sycl::queue *q = backend.sycl_queue_;
    sycl::range<1> n_part{n_particles};

    const int nx_Ex = em.Ex_m.nx_m, ny_Ex = em.Ex_m.ny_m, nz_Ex = em.Ex_m.nz_m;
    const int nx_Ey = em.Ey_m.nx_m, ny_Ey = em.Ey_m.ny_m, nz_Ey = em.Ey_m.nz_m;
    const int nx_Ez = em.Ez_m.nx_m, ny_Ez = em.Ez_m.ny_m, nz_Ez = em.Ez_m.nz_m;

    const int nx_Bx = em.Bx_m.nx_m, ny_Bx = em.Bx_m.ny_m, nz_Bx = em.Bx_m.nz_m;
    const int nx_By = em.By_m.nx_m, ny_By = em.By_m.ny_m, nz_By = em.By_m.nz_m;
    const int nx_Bz = em.Bz_m.nx_m, ny_Bz = em.Bz_m.ny_m, nz_Bz = em.Bz_m.nz_m;

    const mini_float *x = patch.particles_m[is].x_.get_raw_pointer(minipic::device);
    const mini_float *y = patch.particles_m[is].y_.get_raw_pointer(minipic::device);
    const mini_float *z = patch.particles_m[is].z_.get_raw_pointer(minipic::device);

    mini_float *const __restrict__ Exp = patch.particles_m[is].Ex_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ Eyp = patch.particles_m[is].Ey_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ Ezp = patch.particles_m[is].Ez_.get_raw_pointer(minipic::device);

    mini_float *const __restrict__ Bxp = patch.particles_m[is].Bx_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ Byp = patch.particles_m[is].By_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ Bzp = patch.particles_m[is].Bz_.get_raw_pointer(minipic::device);

    const mini_float *const __restrict__ Ex = em.Ex_m.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Ey = em.Ey_m.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Ez = em.Ez_m.get_raw_pointer(minipic::device);

    const mini_float *const __restrict__ Bx = em.Bx_m.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ By = em.By_m.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Bz = em.Bz_m.get_raw_pointer(minipic::device);

    q->submit([&](auto &handler) {
      handler.parallel_for(n_part, [=](sycl::id<1> part) {
        // // Calculate normalized positions
        const mini_float ixn = x[part] * inv_dx_m;
        const mini_float iyn = y[part] * inv_dy_m;
        const mini_float izn = z[part] * inv_dz_m;

        // // Compute indexes in global primal grid
        const unsigned int ixp = std::floor(ixn);
        const unsigned int iyp = std::floor(iyn);
        const unsigned int izp = std::floor(izn);

        // Compute indexes in global dual grid
        const unsigned int ixd = std::floor(ixn + 0.5f);
        const unsigned int iyd = std::floor(iyn + 0.5f);
        const unsigned int izd = std::floor(izn + 0.5f);

        // Compute interpolation coeff, p = primal, d = dual

        // interpolation electric field
        // Ex (d, p , p)
        {
          mini_float coeffs[3] = {ixn + 0.5f, iyn, izn};
          const mini_float v00 = Ex[ixd * (ny_Ex * nz_Ex) + iyp * (nz_Ex) + izp] * (1 - coeffs[0]) +
                                 Ex[(ixd + 1) * (ny_Ex * nz_Ex) + iyp * (nz_Ex) + izp] * coeffs[0];
          const mini_float v01 =
            Ex[ixd * (ny_Ex * nz_Ex) + iyp * (nz_Ex) + (izp + 1)] * (1 - coeffs[0]) +
            Ex[(ixd + 1) * (ny_Ex * nz_Ex) + iyp * (nz_Ex) + (izp + 1)] * coeffs[0];
          const mini_float v10 =
            Ex[ixd * (ny_Ex * nz_Ex) + (iyp + 1) * (nz_Ex) + izp] * (1 - coeffs[0]) +
            Ex[(ixd + 1) * (ny_Ex * nz_Ex) + (iyp + 1) * (nz_Ex) + izp] * coeffs[0];
          const mini_float v11 =
            Ex[ixd * (ny_Ex * nz_Ex) + (iyp + 1) * (nz_Ex) + (izp + 1)] * (1 - coeffs[0]) +
            Ex[(ixd + 1) * (ny_Ex * nz_Ex) + (iyp + 1) * (nz_Ex) + (izp + 1)] * coeffs[0];
          const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
          const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];
          Exp[part]           = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
        }

        // Ey (p, d, p)
        {
          const mini_float coeffs[3] = {ixn, iyn + 0.5f, izn};
          const mini_float v00 = Ey[ixp * (ny_Ey * nz_Ey) + iyd * (nz_Ey) + izp] * (1 - coeffs[0]) +
                                 Ey[(ixp + 1) * (ny_Ey * nz_Ey) + iyd * (nz_Ey) + izp] * coeffs[0];
          const mini_float v01 =
            Ey[ixp * (ny_Ey * nz_Ey) + iyd * (nz_Ey) + (izp + 1)] * (1 - coeffs[0]) +
            Ey[(ixp + 1) * (ny_Ey * nz_Ey) + iyd * (nz_Ey) + (izp + 1)] * coeffs[0];
          const mini_float v10 =
            Ey[ixp * (ny_Ey * nz_Ey) + (iyd + 1) * (nz_Ey) + izp] * (1 - coeffs[0]) +
            Ey[(ixp + 1) * (ny_Ey * nz_Ey) + (iyd + 1) * (nz_Ey) + izp] * coeffs[0];
          const mini_float v11 =
            Ey[ixp * (ny_Ey * nz_Ey) + (iyd + 1) * (nz_Ey) + (izp + 1)] * (1 - coeffs[0]) +
            Ey[(ixp + 1) * (ny_Ey * nz_Ey) + (iyd + 1) * (nz_Ey) + (izp + 1)] * coeffs[0];
          const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
          const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

          Eyp[part] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
        }

        // // //particles_m[is].Ey_.d_view(part) = compute_interpolation(ixp, b, izp, coeffs, Ey);
        // Ez (p, p, d)
        {
          const mini_float coeffs[3] = {ixn, iyn, izn + 0.5f};

          const mini_float v00 = Ez[ixp * (ny_Ez * nz_Ez) + iyp * (nz_Ez) + izd] * (1 - coeffs[0]) +
                                 Ez[(ixp + 1) * (ny_Ez * nz_Ez) + iyp * (nz_Ez) + izd] * coeffs[0];
          const mini_float v01 =
            Ez[ixp * (ny_Ez * nz_Ez) + iyp * (nz_Ez) + (izd + 1)] * (1 - coeffs[0]) +
            Ez[(ixp + 1) * (ny_Ez * nz_Ez) + iyp * (nz_Ez) + (izd + 1)] * coeffs[0];
          const mini_float v10 =
            Ez[ixp * (ny_Ez * nz_Ez) + (iyp + 1) * (nz_Ez) + izd] * (1 - coeffs[0]) +
            Ez[(ixp + 1) * (ny_Ez * nz_Ez) + (iyp + 1) * (nz_Ez) + izd] * coeffs[0];
          const mini_float v11 =
            Ez[ixp * (ny_Ez * nz_Ez) + (iyp + 1) * (nz_Ez) + (izd + 1)] * (1 - coeffs[0]) +
            Ez[(ixp + 1) * (ny_Ez * nz_Ez) + (iyp + 1) * (nz_Ez) + (izd + 1)] * coeffs[0];
          const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
          const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

          Ezp[part] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
        }
        // particles_m[is].Ez_.d_view(part) = compute_interpolation(ixp, iyp, g, coeffs, Ez);

        // interpolation magnetic field
        // Bx (p, d, d)
        {
          const mini_float coeffs[3] = {ixn, iyn + 0.5f, izn + 0.5f};

          const mini_float v00 = Bx[ixp * (ny_Bx * nz_Bx) + iyd * (nz_Bx) + izd] * (1 - coeffs[0]) +
                                 Bx[(ixp + 1) * (ny_Bx * nz_Bx) + iyd * (nz_Bx) + izd] * coeffs[0];
          const mini_float v01 =
            Bx[ixp * (ny_Bx * nz_Bx) + iyd * (nz_Bx) + (izd + 1)] * (1 - coeffs[0]) +
            Bx[(ixp + 1) * (ny_Bx * nz_Bx) + iyd * (nz_Bx) + (izd + 1)] * coeffs[0];
          const mini_float v10 =
            Bx[ixp * (ny_Bx * nz_Bx) + (iyd + 1) * (nz_Bx) + izd] * (1 - coeffs[0]) +
            Bx[(ixp + 1) * (ny_Bx * nz_Bx) + (iyd + 1) * (nz_Bx) + izd] * coeffs[0];
          const mini_float v11 =
            Bx[ixp * (ny_Bx * nz_Bx) + (iyd + 1) * (nz_Bx) + (izd + 1)] * (1 - coeffs[0]) +
            Bx[(ixp + 1) * (ny_Bx * nz_Bx) + (iyd + 1) * (nz_Bx) + (izd + 1)] * coeffs[0];
          const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
          const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

          Bxp[part] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
        }
        // particles_m[is].Bx_.d_view(part) = compute_interpolation(ixp, b, g, coeffs, Bx);

        // By (d, p, d)
        {
          const mini_float coeffs[3] = {ixn + 0.5f, iyn, izn + 0.5f};

          const mini_float v00 = By[ixd * (ny_By * nz_By) + iyp * (nz_By) + izd] * (1 - coeffs[0]) +
                                 By[(ixd + 1) * (ny_By * nz_By) + iyp * (nz_By) + izd] * coeffs[0];
          const mini_float v01 =
            By[ixd * (ny_By * nz_By) + iyp * (nz_By) + (izd + 1)] * (1 - coeffs[0]) +
            By[(ixd + 1) * (ny_By * nz_By) + iyp * (nz_By) + (izd + 1)] * coeffs[0];
          const mini_float v10 =
            By[ixd * (ny_By * nz_By) + (iyp + 1) * (nz_By) + izd] * (1 - coeffs[0]) +
            By[(ixd + 1) * (ny_By * nz_By) + (iyp + 1) * (nz_By) + izd] * coeffs[0];
          const mini_float v11 =
            By[ixd * (ny_By * nz_By) + (iyp + 1) * (nz_By) + (izd + 1)] * (1 - coeffs[0]) +
            By[(ixd + 1) * (ny_By * nz_By) + (iyp + 1) * (nz_By) + (izd + 1)] * coeffs[0];
          const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
          const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

          Byp[part] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
        }
        // particles_m[is].By_.d_view(part) = compute_interpolation(a, iyp, g, coeffs, By);

        // Bz (d, d, p)
        {
          const mini_float coeffs[3] = {ixn + 0.5f, iyn + 0.5f, izn};

          const mini_float v00 = Bz[ixd * (ny_Bz * nz_Bz) + iyd * (nz_Bz) + izp] * (1 - coeffs[0]) +
                                 Bz[(ixd + 1) * (ny_Bz * nz_Bz) + iyd * (nz_Bz) + izp] * coeffs[0];
          const mini_float v01 =
            Bz[ixd * (ny_Bz * nz_Bz) + iyd * (nz_Bz) + (izp + 1)] * (1 - coeffs[0]) +
            Bz[(ixd + 1) * (ny_Bz * nz_Bz) + iyd * (nz_Bz) + (izp + 1)] * coeffs[0];
          const mini_float v10 =
            Bz[ixd * (ny_Bz * nz_Bz) + (iyd + 1) * (nz_Bz) + izp] * (1 - coeffs[0]) +
            Bz[(ixd + 1) * (ny_Bz * nz_Bz) + (iyd + 1) * (nz_Bz) + izp] * coeffs[0];
          const mini_float v11 =
            Bz[ixd * (ny_Bz * nz_Bz) + (iyd + 1) * (nz_Bz) + (izp + 1)] * (1 - coeffs[0]) +
            Bz[(ixd + 1) * (ny_Bz * nz_Bz) + (iyd + 1) * (nz_Bz) + (izp + 1)] * coeffs[0];
          const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
          const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

          Bzp[part] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
        }
      }); // end lambda and sycl parallel_for
    });   // end sycl submit

    q->wait();

  } // Species loop
}

// ______________________________________________________________________________
//
//! \brief Move the particle in the space, compute with EM fields interpolate
//! \param[in] patch  patch data structure
//! \param[in] dt time step to use for the pusher
// ______________________________________________________________________________
auto push(Patch &patch, Backend &backend, double dt) -> void {

  // For each species
  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

    // q' = dt * (q/2m)
    const mini_float qp = patch.particles_m[is].charge_m * dt * 0.5 / patch.particles_m[is].mass_m;

    sycl::queue *q = backend.sycl_queue_;
    sycl::range<1> n_part{n_particles};

    mini_float *const __restrict__ x = patch.particles_m[is].x_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ y = patch.particles_m[is].y_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ z = patch.particles_m[is].z_.get_raw_pointer(minipic::device);

    mini_float *const __restrict__ mx = patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ my = patch.particles_m[is].my_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ mz = patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

    const mini_float *const __restrict__ Exp =
      patch.particles_m[is].Ex_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Eyp =
      patch.particles_m[is].Ey_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Ezp =
      patch.particles_m[is].Ez_.get_raw_pointer(minipic::device);

    const mini_float *const __restrict__ Bxp =
      patch.particles_m[is].Bx_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Byp =
      patch.particles_m[is].By_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Bzp =
      patch.particles_m[is].Bz_.get_raw_pointer(minipic::device);

    q->parallel_for(n_part,
                    [=](sycl::id<1> ip) {
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
                    } // end lambda
    ); // end sycl parallel_for

    q->wait(); // end sycl submit

  } // Loop on species
}

// ______________________________________________________________________________
//
//! \brief Push only the momentum
//! \param[in] patch  patch data structure
//! \param[in] dt time step to use for the pusher
// ______________________________________________________________________________
auto push_momentum(Patch &patch, Backend &backend, double dt) -> void {

  // for each species
  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

    // q' = dt * (q/2m)
    const mini_float qp = patch.particles_m[is].charge_m * dt * 0.5 / patch.particles_m[is].mass_m;

    sycl::queue *q = patch.particles_m[is].x_.sycl_queue_ptr;
    sycl::range<1> n_part{n_particles};

    mini_float *const __restrict__ mx = patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ my = patch.particles_m[is].my_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ mz = patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

    const mini_float *const __restrict__ Exp =
      patch.particles_m[is].Ex_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Eyp =
      patch.particles_m[is].Ey_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Ezp =
      patch.particles_m[is].Ez_.get_raw_pointer(minipic::device);

    const mini_float *const __restrict__ Bxp =
      patch.particles_m[is].Bx_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Byp =
      patch.particles_m[is].By_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ Bzp =
      patch.particles_m[is].Bz_.get_raw_pointer(minipic::device);

    q->parallel_for(n_part,
                    [=](sycl::id<1> ip) {
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
                    } // end lambda
    ); // end sycl parallel_for

    q->wait(); // end sycl submit
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
auto pushBC(Params &params, Patch &patch, Backend &backend) -> void {

  if (patch.on_border_m) {

    const mini_float inf_global[3] = {params.inf_x, params.inf_y, params.inf_z};
    const mini_float sup_global[3] = {params.sup_x, params.sup_y, params.sup_z};

    // Periodic conditions
    if (params.boundary_condition_code == 1) {

      const int N_patches[3]     = {patch.nx_patchs_m, patch.ny_patchs_m, patch.nz_patchs_m};
      const mini_float length[3] = {params.Lx, params.Ly, params.Lz};

      for (int is = 0; is < patch.n_species_m; is++) {

        size_t n_particles = patch.particles_m[is].size();

        sycl::queue *q = backend.sycl_queue_;

        mini_float *const __restrict__ x =
          patch.particles_m[is].x_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ y =
          patch.particles_m[is].y_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ z =
          patch.particles_m[is].z_.get_raw_pointer(minipic::device);

        q->parallel_for(n_particles,
                        [=](sycl::id<1> part) {
                          mini_float *pos[3] = {&x[part], &y[part], &z[part]};

                          for (int d = 0; d < 3; d++) {

                            // Only relevant if there is just 1 patch in this direction
                            // Else the patch exchange with periodicity is managed in the dedicated
                            // function
                            if (N_patches[d] == 1) {
                              if (*pos[d] >= sup_global[d]) {

                                *pos[d] -= length[d];

                              } else if (*pos[d] < inf_global[d]) {

                                *pos[d] += length[d];
                              }
                            }
                          }
                        } // End loop on particles
        ); // end sycl parallel_for

        q->wait();

      } // End loop on species

      // Reflective conditions
    } else if (params.boundary_condition_code == 2) {
      for (int is = 0; is < patch.n_species_m; is++) {

        size_t n_particles = patch.particles_m[is].size();

        sycl::queue *q = patch.particles_m[is].x_.sycl_queue_ptr;

        mini_float *const __restrict__ x =
          patch.particles_m[is].x_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ y =
          patch.particles_m[is].y_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ z =
          patch.particles_m[is].z_.get_raw_pointer(minipic::device);

        mini_float *const __restrict__ mx =
          patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ my =
          patch.particles_m[is].my_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ mz =
          patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

        q->parallel_for(n_particles,
                        [=](sycl::id<1> part) {
                          mini_float *pos[3] = {&x[part], &y[part], &z[part]};

                          mini_float *momentum[3] = {&mx[part], &my[part], &mz[part]};

                          for (int d = 0; d < 3; d++) {

                            if (*pos[d] >= sup_global[d]) {

                              *pos[d]      = 2 * sup_global[d] - *pos[d];
                              *momentum[d] = -*momentum[d];

                            } else if (*pos[d] < inf_global[d]) {

                              *pos[d]      = 2 * inf_global[d] - *pos[d];
                              *momentum[d] = -*momentum[d];
                            }
                          }
                        } // End lambda
        ); // end sycl parallel_for

        q->wait();

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
auto project(Params &params, Patch &patch, Backend &backend) -> void {

  for (int is = 0; is < patch.n_species_m; is++) {

    patch.vec_Jx_m[is].reset(minipic::host);
    patch.vec_Jy_m[is].reset(minipic::host);
    patch.vec_Jz_m[is].reset(minipic::host);

    patch.vec_Jx_m[is].reset(minipic::device);
    patch.vec_Jy_m[is].reset(minipic::device);
    patch.vec_Jz_m[is].reset(minipic::device);

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

      sycl::queue *q = backend.sycl_queue_;

      const int nx_Jx = patch.vec_Jx_m[is].nx_m, ny_Jx = patch.vec_Jx_m[is].ny_m,
                nz_Jx = patch.vec_Jx_m[is].nz_m;
      const int nx_Jy = patch.vec_Jy_m[is].nx_m, ny_Jy = patch.vec_Jy_m[is].ny_m,
                nz_Jy = patch.vec_Jy_m[is].nz_m;
      const int nx_Jz = patch.vec_Jz_m[is].nx_m, ny_Jz = patch.vec_Jz_m[is].ny_m,
                nz_Jz = patch.vec_Jz_m[is].nz_m;

      // Pointer toward Sycl views

      const mini_float *const __restrict__ w =
        patch.particles_m[is].weight_.get_raw_pointer(minipic::device);

      const mini_float *const __restrict__ x =
        patch.particles_m[is].x_.get_raw_pointer(minipic::device);
      const mini_float *const __restrict__ y =
        patch.particles_m[is].y_.get_raw_pointer(minipic::device);
      const mini_float *const __restrict__ z =
        patch.particles_m[is].z_.get_raw_pointer(minipic::device);

      const mini_float *const __restrict__ mx =
        patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
      const mini_float *const __restrict__ my =
        patch.particles_m[is].my_.get_raw_pointer(minipic::device);
      const mini_float *const __restrict__ mz =
        patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

      mini_float *const __restrict__ Jx_loc = patch.vec_Jx_m[is].get_raw_pointer(minipic::device);
      mini_float *const __restrict__ Jy_loc = patch.vec_Jy_m[is].get_raw_pointer(minipic::device);
      mini_float *const __restrict__ Jz_loc = patch.vec_Jz_m[is].get_raw_pointer(minipic::device);

      q->parallel_for(
        n_particles,
        [=](sycl::id<1> part) {
          const mini_float gamma_inv =
            1 / sqrt(1 + mx[part] * mx[part] + my[part] * my[part] + mz[part] * mz[part]);

          const mini_float charge_weight = inv_cell_volume_x_q * w[part];

          const mini_float vx = mx[part] * gamma_inv;
          const mini_float vy = my[part] * gamma_inv;
          const mini_float vz = mz[part] * gamma_inv;

          // Current from the particle
          const mini_float Jxp = vx * charge_weight;
          const mini_float Jyp = vy * charge_weight;
          const mini_float Jzp = vz * charge_weight;

          // Calculate normalized position relative to the patch
          // ixn = (particles_m[is].x(part) ) * params.inv_dx;
          // iyn = (particles_m[is].y(part) ) * params.inv_dy;
          // izn = (particles_m[is].z(part) ) * params.inv_dz;
          const mini_float posxn = (x[part] - 0.5 * dt * vx - xmin) * inv_dx + 1;
          const mini_float posyn = (y[part] - 0.5 * dt * vy - ymin) * inv_dy + 1;
          const mini_float poszn = (z[part] - 0.5 * dt * vz - zmin) * inv_dz + 1;

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

          // Project on Jx
          {
            const mini_float coeffs[3] = {posxn - 0.5f - ixd, posyn - iyp, poszn - izp};

            auto Jx000 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jx_loc[ixd * (ny_Jx * nz_Jx) + iyp * (nz_Jx) + izp]);
            auto Jx001 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jx_loc[ixd * (ny_Jx * nz_Jx) + iyp * (nz_Jx) + (izp + 1)]);
            auto Jx010 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jx_loc[ixd * (ny_Jx * nz_Jx) + (iyp + 1) * (nz_Jx) + izp]);
            auto Jx011 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jx_loc[ixd * (ny_Jx * nz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)]);
            auto Jx100 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jx_loc[(ixd + 1) * (ny_Jx * nz_Jx) + iyp * (nz_Jx) + izp]);
            auto Jx101 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jx_loc[(ixd + 1) * (ny_Jx * nz_Jx) + iyp * (nz_Jx) + (izp + 1)]);
            auto Jx110 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jx_loc[(ixd + 1) * (ny_Jx * nz_Jx) + (iyp + 1) * (nz_Jx) + izp]);
            auto Jx111 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jx_loc[(ixd + 1) * (ny_Jx * nz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)]);

            Jx000 += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
            Jx001 += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
            Jx010 += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
            Jx011 += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;
            Jx100 += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
            Jx101 += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
            Jx110 += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
            Jx111 += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;
          }

          // Project Jy
          {
            const mini_float coeffs[3] = {posxn - ixp, posyn - 0.5f - iyd, poszn - izp};

            auto Jy000 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jy_loc[ixp * (ny_Jy * nz_Jy) + iyd * (nz_Jy) + izp]);
            auto Jy001 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jy_loc[ixp * (ny_Jy * nz_Jy) + iyd * (nz_Jy) + (izp + 1)]);
            auto Jy010 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jy_loc[ixp * (ny_Jy * nz_Jy) + (iyd + 1) * (nz_Jy) + izp]);
            auto Jy011 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jy_loc[ixp * (ny_Jy * nz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)]);
            auto Jy100 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jy_loc[(ixp + 1) * (ny_Jy * nz_Jy) + iyd * (nz_Jy) + izp]);
            auto Jy101 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jy_loc[(ixp + 1) * (ny_Jy * nz_Jy) + iyd * (nz_Jy) + (izp + 1)]);
            auto Jy110 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jy_loc[(ixp + 1) * (ny_Jy * nz_Jy) + (iyd + 1) * (nz_Jy) + izp]);
            auto Jy111 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jy_loc[(ixp + 1) * (ny_Jy * nz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)]);

            Jy000 += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
            Jy001 += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
            Jy010 += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
            Jy011 += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;
            Jy100 += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
            Jy101 += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
            Jy110 += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
            Jy111 += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;
          }

          // Project Jz
          {
            mini_float coeffs[3] = {posxn - ixp, posyn - iyp, poszn - 0.5f - izd};

            auto Jz000 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jz_loc[ixp * (ny_Jz * nz_Jz) + iyp * (nz_Jz) + izd]);
            auto Jz001 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jz_loc[ixp * (ny_Jz * nz_Jz) + iyp * (nz_Jz) + (izd + 1)]);
            auto Jz010 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jz_loc[ixp * (ny_Jz * nz_Jz) + (iyp + 1) * (nz_Jz) + izd]);
            auto Jz011 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jz_loc[ixp * (ny_Jz * nz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)]);
            auto Jz100 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jz_loc[(ixp + 1) * (ny_Jz * nz_Jz) + iyp * (nz_Jz) + izd]);
            auto Jz101 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jz_loc[(ixp + 1) * (ny_Jz * nz_Jz) + iyp * (nz_Jz) + (izd + 1)]);
            auto Jz110 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jz_loc[(ixp + 1) * (ny_Jz * nz_Jz) + (iyp + 1) * (nz_Jz) + izd]);
            auto Jz111 = sycl::atomic_ref<mini_float,
                                          sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>(
              Jz_loc[(ixp + 1) * (ny_Jz * nz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)]);

            Jz000 += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
            Jz001 += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
            Jz010 += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
            Jz011 += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
            Jz100 += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
            Jz101 += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
            Jz110 += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
            Jz111 += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
          }
        } // End lambda
      ); // end sycl parallel_for

      q->wait();

      patch.projected_[is] = true;

    } else {

      patch.projected_[is] = false;

    } // end if n_particles > 0

  } // end loop species
}

// _______________________________________________________________________
//
//! \brief Current projection directly in the global array
//! \param[in] params constant global parameters
//! \param[in] em electromagnetic fields
//! \param[in] patch current patch to handle
// _______________________________________________________________________
auto project(Params &params, ElectroMagn &em, Patch &patch, Backend &backend) -> void {}

// _______________________________________________________
//
//! \brief Solve Maxwell equations to compute EM fields
//! \param params global parameters
// _______________________________________________________
auto solve_maxwell(const Params &params,
                   ElectroMagn &em,
                   Backend &backend,
                   Profiler &profiler) -> void {

  const auto dt         = params.dt;
  const auto dt_over_dx = params.dt * params.inv_dx;
  const auto dt_over_dy = params.dt * params.inv_dy;
  const auto dt_over_dz = params.dt * params.inv_dz;

  sycl::queue *q = backend.sycl_queue_;

  const int nx_p = em.nx_p_m, ny_p = em.ny_p_m, nz_p = em.nz_p_m;
  const int nx_d = em.nx_d_m, ny_d = em.ny_d_m, nz_d = em.nz_d_m;

  const int nx_Jx = em.Jx_m.nx_m, ny_Jx = em.Jx_m.ny_m, nz_Jx = em.Jx_m.nz_m;
  const int nx_Jy = em.Jy_m.nx_m, ny_Jy = em.Jy_m.ny_m, nz_Jy = em.Jy_m.nz_m;
  const int nx_Jz = em.Jz_m.nx_m, ny_Jz = em.Jz_m.ny_m, nz_Jz = em.Jz_m.nz_m;

  mini_float *const __restrict__ Jx = em.Jx_m.get_raw_pointer(minipic::device);
  mini_float *const __restrict__ Jy = em.Jy_m.get_raw_pointer(minipic::device);
  mini_float *const __restrict__ Jz = em.Jz_m.get_raw_pointer(minipic::device);

  mini_float *const __restrict__ Ex = em.Ex_m.get_raw_pointer(minipic::device);
  mini_float *const __restrict__ Ey = em.Ey_m.get_raw_pointer(minipic::device);
  mini_float *const __restrict__ Ez = em.Ez_m.get_raw_pointer(minipic::device);

  mini_float *const __restrict__ Bx = em.Bx_m.get_raw_pointer(minipic::device);
  mini_float *const __restrict__ By = em.By_m.get_raw_pointer(minipic::device);
  mini_float *const __restrict__ Bz = em.Bz_m.get_raw_pointer(minipic::device);

  sycl::range<3> global_dpp{static_cast<size_t>(nx_d),
                            static_cast<size_t>(ny_p),
                            static_cast<size_t>(nz_p)};
  sycl::range<3> global_pdp{static_cast<size_t>(nx_p),
                            static_cast<size_t>(ny_d),
                            static_cast<size_t>(nz_p)};
  sycl::range<3> global_ppd{static_cast<size_t>(nx_p),
                            static_cast<size_t>(ny_p),
                            static_cast<size_t>(nz_d)};
  sycl::range<3> local_dpp{1, 1, static_cast<size_t>(nz_p)};
  sycl::range<3> local_pdp{1, 1, static_cast<size_t>(nz_p)};
  sycl::range<3> local_ppd{1, 1, static_cast<size_t>(nz_d)};

  sycl::range<3> global_pdd{static_cast<size_t>(nx_p),
                            static_cast<size_t>(ny_d - 2),
                            static_cast<size_t>(nz_d - 2)};
  sycl::range<3> global_dpd{static_cast<size_t>(nx_d - 2),
                            static_cast<size_t>(ny_p),
                            static_cast<size_t>(nz_d - 2)};
  sycl::range<3> global_ddp{static_cast<size_t>(nx_d - 2),
                            static_cast<size_t>(ny_d - 2),
                            static_cast<size_t>(nz_p)};
  sycl::range<3> local_pdd{1, 1, static_cast<size_t>(nz_d - 2)};
  sycl::range<3> local_dpd{1, 1, static_cast<size_t>(nz_d - 2)};
  sycl::range<3> local_ddp{1, 1, static_cast<size_t>(nz_p)};

  /////     Solve Maxwell Ampere (E)
  // Electric field Ex (d,p,p)

  q->parallel_for(sycl::nd_range{global_dpp, local_dpp}, [=](sycl::nd_item<3> item) {
    size_t ix = item.get_global_id(0);
    size_t iy = item.get_global_id(1);
    size_t iz = item.get_global_id(2);

    Ex[ix * (ny_p * nz_p) + iy * (nz_p) + iz] +=
      -dt * Jx[ix * (ny_Jx * nz_Jx) + (iy + 1) * (nz_Jx) + iz + 1] +
      dt_over_dy * (Bz[ix * (ny_d * nz_p) + (iy + 1) * (nz_p) + iz] -
                    Bz[ix * (ny_d * nz_p) + iy * (nz_p) + iz]) -
      dt_over_dz * (By[ix * (ny_p * nz_d) + iy * (nz_d) + (iz + 1)] -
                    By[ix * (ny_p * nz_d) + iy * (nz_d) + iz]);
  }); // end sycl parallel_for

  // Electric field Ey (p,d,p)
  q->parallel_for(sycl::nd_range{global_pdp, local_pdp}, [=](sycl::nd_item<3> item) {
    size_t ix = item.get_global_id(0);
    size_t iy = item.get_global_id(1);
    size_t iz = item.get_global_id(2);

    Ey[ix * (ny_d * nz_p) + iy * (nz_p) + iz] +=
      -dt * Jy[(ix + 1) * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz + 1] -
      dt_over_dx * (Bz[(ix + 1) * (ny_d * nz_p) + iy * (nz_p) + iz] -
                    Bz[ix * (ny_d * nz_p) + iy * (nz_p) + iz]) +
      dt_over_dz * (Bx[ix * (ny_d * nz_d) + iy * (nz_d) + (iz + 1)] -
                    Bx[ix * (ny_d * nz_d) + iy * (nz_d) + iz]);
  }); // end sycl parallel_for

  // Electric field Ez (p,p,d)
  q->parallel_for(sycl::nd_range{global_ppd, local_ppd}, [=](sycl::nd_item<3> item) {
    size_t ix = item.get_global_id(0);
    size_t iy = item.get_global_id(1);
    size_t iz = item.get_global_id(2);
    Ez[ix * (ny_p * nz_d) + iy * (nz_d) + iz] +=
      -dt * Jz[(ix + 1) * (ny_Jz * nz_Jz) + (iy + 1) * (nz_Jz) + iz] +
      dt_over_dx * (By[(ix + 1) * (ny_p * nz_d) + iy * (nz_d) + iz] -
                    By[ix * (ny_p * nz_d) + iy * (nz_d) + iz]) -
      dt_over_dy * (Bx[ix * (ny_d * nz_d) + (iy + 1) * (nz_d) + iz] -
                    Bx[ix * (ny_d * nz_d) + iy * (nz_d) + iz]);
  }); // end sycl parallel_for

  q->wait();

  /////     Solve Maxwell Faraday (B)

  // Magnetic field Bx (p,d,d)

  q->parallel_for(sycl::nd_range{global_pdd, local_pdd}, [=](sycl::nd_item<3> item) {
    size_t ix = item.get_global_id(0);
    size_t iy = item.get_global_id(1) + 1;
    size_t iz = item.get_global_id(2) + 1;

    Bx[ix * (ny_d * nz_d) + iy * (nz_d) + iz] +=
      -dt_over_dy * (Ez[ix * (ny_p * nz_d) + iy * (nz_d) + iz] -
                     Ez[ix * (ny_p * nz_d) + (iy - 1) * (nz_d) + iz]) +
      dt_over_dz *
        (Ey[ix * (ny_d * nz_p) + iy * (nz_p) + iz] - Ey[ix * (ny_d * nz_p) + iy * (nz_p) + iz - 1]);
  }); // end sycl parallel_for

  // Magnetic field By (d,p,d)

  q->parallel_for(sycl::nd_range{global_dpd, local_dpd}, [=](sycl::nd_item<3> item) {
    size_t ix = item.get_global_id(0) + 1;
    size_t iy = item.get_global_id(1);
    size_t iz = item.get_global_id(2) + 1;

    By[ix * (ny_p * nz_d) + iy * (nz_d) + iz] +=
      -dt_over_dz * (Ex[ix * (ny_p * nz_p) + iy * (nz_p) + iz] -
                     Ex[ix * (ny_p * nz_p) + iy * (nz_p) + iz - 1]) +
      dt_over_dx * (Ez[ix * (ny_p * nz_d) + iy * (nz_d) + iz] -
                    Ez[(ix - 1) * (ny_p * nz_d) + iy * (nz_d) + iz]);
  }); // end sycl parallel_for

  // Magnetic field Bz (d,d,p)

  q->parallel_for(sycl::nd_range{global_ddp, local_ddp}, [=](sycl::nd_item<3> item) {
    size_t ix = item.get_global_id(0) + 1;
    size_t iy = item.get_global_id(1) + 1;
    size_t iz = item.get_global_id(2);

    Bz[ix * (ny_d * nz_p) + iy * (nz_p) + iz] +=
      -dt_over_dx * (Ey[ix * (ny_d * nz_p) + iy * (nz_p) + iz] -
                     Ey[(ix - 1) * (ny_d * nz_p) + iy * (nz_p) + iz]) +
      dt_over_dy * (Ex[ix * (ny_p * nz_p) + iy * (nz_p) + iz] -
                    Ex[ix * (ny_p * nz_p) + (iy - 1) * (nz_p) + iz]);
  }); // end sycl parallel_for

  q->wait();

} // end solve

// _______________________________________________________________
//
//! \brief Boundaries condition on the global grid
//! \param[in] Params & params - global constant parameters
//! \param[in] ElectroMagn & em - global electromagnetic fields
// _______________________________________________________________
auto currentBC(Params &params, ElectroMagn &em, Backend &backend) -> void {

  if (params.boundary_condition == "periodic") {

    sycl::queue *q = backend.sycl_queue_;

    mini_float *const __restrict__ Jx = em.Jx_m.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ Jy = em.Jy_m.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ Jz = em.Jz_m.get_raw_pointer(minipic::device);

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

    q->parallel_for(sycl::range<2>{static_cast<size_t>(ny_Jx), static_cast<size_t>(nz_Jx)},
                    [=](sycl::id<2> idx) {
                      const auto iy = idx[0];
                      const auto iz = idx[1];

                      Jx[0 * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz] +=
                        Jx[(nx_Jx - 2) * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz];
                      Jx[(nx_Jx - 2) * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz] =
                        Jx[0 * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz];

                      Jx[1 * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz] +=
                        Jx[(nx_Jx - 1) * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz];
                      Jx[(nx_Jx - 1) * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz] =
                        Jx[1 * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz];
                    });

    q->parallel_for(sycl::range<2>{static_cast<size_t>(ny_Jy), static_cast<size_t>(nz_Jy)},
                    [=](sycl::id<2> idx) {
                      const auto iy = idx[0];
                      const auto iz = idx[1];

                      Jy[0 * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz] +=
                        Jy[(nx_Jy - 2) * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz];
                      Jy[(nx_Jy - 2) * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz] =
                        Jy[0 * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz];

                      Jy[1 * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz] +=
                        Jy[(nx_Jy - 1) * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz];
                      Jy[(nx_Jy - 1) * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz] =
                        Jy[1 * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz];
                    });

    q->parallel_for(sycl::range<2>{static_cast<size_t>(ny_Jz), static_cast<size_t>(nz_Jz)},
                    [=](sycl::id<2> idx) {
                      const auto iy = idx[0];
                      const auto iz = idx[1];

                      Jz[0 * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz] +=
                        Jz[(nx_Jz - 2) * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz];
                      Jz[(nx_Jz - 2) * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz] =
                        Jz[0 * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz];

                      Jz[1 * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz] +=
                        Jz[(nx_Jz - 1) * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz];
                      Jz[(nx_Jz - 1) * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz] =
                        Jz[1 * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz];
                    });

    q->wait();

    // Y

    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_Jx), static_cast<size_t>(nz_Jx)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iz = idx[1];

                      Jx[ix * (ny_Jx * nz_Jx) + 0 * (nz_Jx) + iz] +=
                        Jx[ix * (ny_Jx * nz_Jx) + (ny_Jx - 2) * (nz_Jx) + iz];
                      Jx[ix * (ny_Jx * nz_Jx) + (ny_Jx - 2) * (nz_Jx) + iz] =
                        Jx[ix * (ny_Jx * nz_Jx) + 0 * (nz_Jx) + iz];

                      Jx[ix * (ny_Jx * nz_Jx) + 1 * (nz_Jx) + iz] +=
                        Jx[ix * (ny_Jx * nz_Jx) + (ny_Jx - 1) * (nz_Jx) + iz];
                      Jx[ix * (ny_Jx * nz_Jx) + (ny_Jx - 1) * (nz_Jx) + iz] =
                        Jx[ix * (ny_Jx * nz_Jx) + 1 * (nz_Jx) + iz];
                    });

    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_Jy), static_cast<size_t>(nz_Jy)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iz = idx[1];

                      Jy[ix * (ny_Jy * nz_Jy) + 0 * (nz_Jy) + iz] +=
                        Jy[ix * (ny_Jy * nz_Jy) + (ny_Jy - 2) * (nz_Jy) + iz];
                      Jy[ix * (ny_Jy * nz_Jy) + (ny_Jy - 2) * (nz_Jy) + iz] =
                        Jy[ix * (ny_Jy * nz_Jy) + 0 * (nz_Jy) + iz];

                      Jy[ix * (ny_Jy * nz_Jy) + 1 * (nz_Jy) + iz] +=
                        Jy[ix * (ny_Jy * nz_Jy) + (ny_Jy - 1) * (nz_Jy) + iz];
                      Jy[ix * (ny_Jy * nz_Jy) + (ny_Jy - 1) * (nz_Jy) + iz] =
                        Jy[ix * (ny_Jy * nz_Jy) + 1 * (nz_Jy) + iz];
                    });

    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_Jz), static_cast<size_t>(nz_Jz)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iz = idx[1];

                      Jz[ix * (ny_Jz * nz_Jz) + 0 * (nz_Jz) + iz] +=
                        Jz[ix * (ny_Jz * nz_Jz) + (ny_Jz - 2) * (nz_Jz) + iz];
                      Jz[ix * (ny_Jz * nz_Jz) + (ny_Jz - 2) * (nz_Jz) + iz] =
                        Jz[ix * (ny_Jz * nz_Jz) + 0 * (nz_Jz) + iz];

                      Jz[ix * (ny_Jz * nz_Jz) + 1 * (nz_Jz) + iz] +=
                        Jz[ix * (ny_Jz * nz_Jz) + (ny_Jz - 1) * (nz_Jz) + iz];
                      Jz[ix * (ny_Jz * nz_Jz) + (ny_Jz - 1) * (nz_Jz) + iz] =
                        Jz[ix * (ny_Jz * nz_Jz) + 1 * (nz_Jz) + iz];
                    });

    q->wait();

    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_Jx), static_cast<size_t>(ny_Jx)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iy = idx[1];

                      Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + 0] +=
                        Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + (nz_Jx - 2)];
                      Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + (nz_Jx - 2)] =
                        Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + 0];

                      Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + 1] +=
                        Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + (nz_Jx - 1)];
                      Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + (nz_Jx - 1)] =
                        Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + 1];
                    });

    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_Jy), static_cast<size_t>(ny_Jy)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iy = idx[1];

                      Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + 0] +=
                        Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + (nz_Jy - 2)];
                      Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + (nz_Jy - 2)] =
                        Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + 0];

                      Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + 1] +=
                        Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + (nz_Jy - 1)];
                      Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + (nz_Jy - 1)] =
                        Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + 1];
                    });

    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_Jz), static_cast<size_t>(ny_Jz)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iy = idx[1];

                      Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + 0] +=
                        Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + (nz_Jz - 2)];
                      Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + (nz_Jz - 2)] =
                        Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + 0];

                      Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + 1] +=
                        Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + (nz_Jz - 1)];
                      Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + (nz_Jz - 1)] =
                        Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + 1];
                    });

    q->wait();

    // q->single_task([=]() {

    // for (unsigned int iy = 0; iy < ny_Jx; ++iy) {
    //   for (unsigned int iz = 0; iz < nz_Jx; ++iz) {

    //     Jx[0 * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz] +=
    //       Jx[(nx_Jx - 2) * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz];
    //     Jx[(nx_Jx - 2) * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz] =
    //       Jx[0 * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz];

    //     Jx[1 * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz] +=
    //       Jx[(nx_Jx - 1) * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz];
    //     Jx[(nx_Jx - 1) * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz] =
    //       Jx[1 * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz];
    //   }
    // }

    // for (unsigned int iy = 0; iy < ny_Jy; ++iy) {
    //   for (unsigned int iz = 0; iz < nz_Jy; ++iz) {

    //     Jy[0 * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz] +=
    //       Jy[(nx_Jy - 2) * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz];
    //     Jy[(nx_Jy - 2) * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz] =
    //       Jy[0 * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz];

    //     Jy[1 * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz] +=
    //       Jy[(nx_Jy - 1) * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz];
    //     Jy[(nx_Jy - 1) * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz] =
    //       Jy[1 * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz];
    //   }
    // }

    // for (unsigned int iy = 0; iy < ny_Jz; ++iy) {
    //   for (unsigned int iz = 0; iz < nz_Jz; ++iz) {
    //     Jz[0 * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz] +=
    //       Jz[(nx_Jz - 2) * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz];
    //     Jz[(nx_Jz - 2) * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz] =
    //       Jz[0 * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz];

    //     Jz[1 * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz] +=
    //       Jz[(nx_Jz - 1) * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz];
    //     Jz[(nx_Jz - 1) * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz] =
    //       Jz[1 * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz];
    //   }
    // }

    // Y

    // for (unsigned int ix = 0; ix < nx_Jx; ++ix) {
    //   for (unsigned int iz = 0; iz < nz_Jx; ++iz) {
    //     Jx[ix * (ny_Jx * nz_Jx) + 0 * (nz_Jx) + iz] +=
    //       Jx[ix * (ny_Jx * nz_Jx) + (ny_Jx - 2) * (nz_Jx) + iz];
    //     Jx[ix * (ny_Jx * nz_Jx) + (ny_Jx - 2) * (nz_Jx) + iz] =
    //       Jx[ix * (ny_Jx * nz_Jx) + 0 * (nz_Jx) + iz];

    //     Jx[ix * (ny_Jx * nz_Jx) + 1 * (nz_Jx) + iz] +=
    //       Jx[ix * (ny_Jx * nz_Jx) + (ny_Jx - 1) * (nz_Jx) + iz];
    //     Jx[ix * (ny_Jx * nz_Jx) + (ny_Jx - 1) * (nz_Jx) + iz] =
    //       Jx[ix * (ny_Jx * nz_Jx) + 1 * (nz_Jx) + iz];
    //   }
    // }

    // for (unsigned int ix = 0; ix < nx_Jy; ++ix) {
    //   for (unsigned int iz = 0; iz < nz_Jy; ++iz) {

    //     Jy[ix * (ny_Jy * nz_Jy) + 0 * (nz_Jy) + iz] +=
    //       Jy[ix * (ny_Jy * nz_Jy) + (ny_Jy - 2) * (nz_Jy) + iz];
    //     Jy[ix * (ny_Jy * nz_Jy) + (ny_Jy - 2) * (nz_Jy) + iz] =
    //       Jy[ix * (ny_Jy * nz_Jy) + 0 * (nz_Jy) + iz];

    //     Jy[ix * (ny_Jy * nz_Jy) + 1 * (nz_Jy) + iz] +=
    //       Jy[ix * (ny_Jy * nz_Jy) + (ny_Jy - 1) * (nz_Jy) + iz];
    //     Jy[ix * (ny_Jy * nz_Jy) + (ny_Jy - 1) * (nz_Jy) + iz] =
    //       Jy[ix * (ny_Jy * nz_Jy) + 1 * (nz_Jy) + iz];
    //   }
    // }

    // for (unsigned int ix = 0; ix < nx_Jz; ++ix) {
    //   for (unsigned int iz = 0; iz < nz_Jz; ++iz) {

    //     Jz[ix * (ny_Jz * nz_Jz) + 0 * (nz_Jz) + iz] +=
    //       Jz[ix * (ny_Jz * nz_Jz) + (ny_Jz - 2) * (nz_Jz) + iz];
    //     Jz[ix * (ny_Jz * nz_Jz) + (ny_Jz - 2) * (nz_Jz) + iz] =
    //       Jz[ix * (ny_Jz * nz_Jz) + 0 * (nz_Jx) + iz];

    //     Jz[ix * (ny_Jz * nz_Jz) + 1 * (nz_Jz) + iz] +=
    //       Jz[ix * (ny_Jz * nz_Jz) + (ny_Jz - 1) * (nz_Jz) + iz];
    //     Jz[ix * (ny_Jz * nz_Jz) + (ny_Jz - 1) * (nz_Jz) + iz] =
    //       Jz[ix * (ny_Jz * nz_Jz) + 1 * (nz_Jz) + iz];
    //   }
    // }

    // Z

    // for (unsigned int ix = 0; ix < nx_Jx; ++ix) {
    //   for (unsigned int iy = 0; iy < ny_Jx; ++iy) {

    //     Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + 0] +=
    //       Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + (nz_Jx - 2)];
    //     Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + (nz_Jx - 2)] =
    //       Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + 0];

    //     Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + 1] +=
    //       Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + (nz_Jx - 1)];
    //     Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + (nz_Jx - 1)] =
    //       Jx[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + 1];
    //   }
    // }

    // for (unsigned int ix = 0; ix < nx_Jy; ++ix) {
    //   for (unsigned int iy = 0; iy < ny_Jy; ++iy) {

    //     Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + 0] +=
    //       Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + (nz_Jy - 2)];
    //     Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + (nz_Jy - 2)] =
    //       Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + 0];

    //     Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + 1] +=
    //       Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + (nz_Jy - 1)];
    //     Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + (nz_Jy - 1)] =
    //       Jy[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + 1];
    //   }
    // }

    // for (unsigned int ix = 0; ix < nx_Jz; ++ix) {
    //   for (unsigned int iy = 0; iy < ny_Jz; ++iy) {

    //     Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + 0] +=
    //       Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + (nz_Jz - 2)];
    //     Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + (nz_Jz - 2)] =
    //       Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + 0];

    //     Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + 1] +=
    //       Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + (nz_Jz - 1)];
    //     Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + (nz_Jz - 1)] =
    //       Jz[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + 1];
    //   }
    // }
    // }); // end sycl parallel_for
    // q->wait();

  } // end if periodic
} // end currentBC

// _______________________________________________________________
//
//! \brief Boundaries condition on the global grid
//! \param[in] Params & params - global constant parameters
// _______________________________________________________________
auto solveBC(Params &params, ElectroMagn &em, Backend &backend) -> void {

  sycl::queue *q = backend.sycl_queue_;

  mini_float *const __restrict__ Bx = em.Bx_m.get_raw_pointer(minipic::device);
  mini_float *const __restrict__ By = em.By_m.get_raw_pointer(minipic::device);
  mini_float *const __restrict__ Bz = em.Bz_m.get_raw_pointer(minipic::device);

  const auto nx_Bx = em.Bx_m.nx();
  const auto ny_Bx = em.Bx_m.ny();
  const auto nz_Bx = em.Bx_m.nz();

  const auto nx_By = em.By_m.nx();
  const auto ny_By = em.By_m.ny();
  const auto nz_By = em.By_m.nz();

  const auto nx_Bz = em.Bz_m.nx();
  const auto ny_Bz = em.Bz_m.ny();
  const auto nz_Bz = em.Bz_m.nz();

  if (params.boundary_condition == "periodic") {

    // X dim

    // By (d,p,d)
    q->parallel_for(sycl::range<2>{static_cast<size_t>(ny_By), static_cast<size_t>(nz_By)},
                    [=](sycl::id<2> idx) {
                      const auto iy = idx[0];
                      const auto iz = idx[1];

                      By[0 * (ny_By * nz_By) + iy * (nz_By) + iz] =
                        By[(nx_By - 2) * (ny_By * nz_By) + iy * (nz_By) + iz];
                      By[(nx_By - 1) * (ny_By * nz_By) + iy * (nz_By) + iz] =
                        By[1 * (ny_By * nz_By) + iy * (nz_By) + iz];
                    });

    // Bz (d,d,p)
    q->parallel_for(sycl::range<2>{static_cast<size_t>(ny_Bz), static_cast<size_t>(nz_Bz)},
                    [=](sycl::id<2> idx) {
                      const auto iy = idx[0];
                      const auto iz = idx[1];

                      Bz[0 * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz] =
                        Bz[(nx_Bz - 2) * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz];
                      Bz[(nx_Bz - 1) * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz] =
                        Bz[1 * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz];
                    });

    q->wait();

    // Y dim
    // Bx (p,d,d)
    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_Bx), static_cast<size_t>(nz_Bx)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iz = idx[1];

                      Bx[ix * (ny_Bx * nz_Bx) + 0 * (nz_Bx) + iz] =
                        Bx[ix * (ny_Bx * nz_Bx) + (ny_Bx - 2) * (nz_Bx) + iz];
                      Bx[ix * (ny_Bx * nz_Bx) + (ny_Bx - 1) * (nz_Bx) + iz] =
                        Bx[ix * (ny_Bx * nz_Bx) + 1 * (nz_Bx) + iz];
                    });

    // Bz (d,d,p)
    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_Bz), static_cast<size_t>(nz_Bz)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iz = idx[1];

                      Bz[ix * (ny_Bz * nz_Bz) + 0 * (nz_Bz) + iz] =
                        Bz[ix * (ny_Bz * nz_Bz) + (ny_Bz - 2) * (nz_Bz) + iz];
                      Bz[ix * (ny_Bz * nz_Bz) + (ny_Bz - 1) * (nz_Bz) + iz] =
                        Bz[ix * (ny_Bz * nz_Bz) + 1 * (nz_Bz) + iz];
                    });

    q->wait();

    // Z dim
    // Bx
    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_Bx), static_cast<size_t>(ny_Bx)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iy = idx[1];

                      Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + 0] =
                        Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + (nz_Bx - 2)];
                      Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + (nz_Bx - 1)] =
                        Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + 1];
                    });

    // By
    q->parallel_for(sycl::range<2>{static_cast<size_t>(nx_By), static_cast<size_t>(ny_By)},
                    [=](sycl::id<2> idx) {
                      const auto ix = idx[0];
                      const auto iy = idx[1];

                      By[ix * (ny_By * nz_By) + iy * (nz_By) + 0] =
                        By[ix * (ny_By * nz_By) + iy * (nz_By) + (nz_By - 2)];
                      By[ix * (ny_By * nz_By) + iy * (nz_By) + (nz_By - 1)] =
                        By[ix * (ny_By * nz_By) + iy * (nz_By) + 1];
                    });

    q->wait();

    // q->single_task([=]() {

    // X dim
    // By (d,p,d)
    // for (unsigned int iy = 0; iy < ny_By; ++iy) {

    //   for (unsigned int iz = 0; iz < nz_By; ++iz) {
    //     // -X
    //     By[0 * (ny_By * nz_By) + iy * (nz_By) + iz] =
    //       By[(nx_By - 2) * (ny_By * nz_By) + iy * (nz_By) + iz];
    //     By[(nx_By - 1) * (ny_By * nz_By) + iy * (nz_By) + iz] =
    //       By[1 * (ny_By * nz_By) + iy * (nz_By) + iz];
    //   }
    // }

    // Bz (d,d,p)
    // for (unsigned int iy = 0; iy < ny_Bz; iy++) {
    //   for (unsigned int iz = 0; iz < nz_Bz; iz++) {
    //     // -X
    //     Bz[0 * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz] =
    //       Bz[(nx_Bz - 2) * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz];
    //     Bz[(nx_Bz - 1) * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz] =
    //       Bz[1 * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz];
    //   }
    // }

    // Y dim
    // Bx (p,d,d)
    // for (unsigned int ix = 0; ix < nx_Bx; ix++) {
    //   for (unsigned int iz = 0; iz < nz_Bx; iz++) {
    //     // -Y
    //     Bx[ix * (ny_Bx * nz_Bx) + 0 * (nz_Bx) + iz] =
    //       Bx[ix * (ny_Bx * nz_Bx) + (ny_Bx - 2) * (nz_Bx) + iz];
    //     Bx[ix * (ny_Bx * nz_Bx) + (ny_Bx - 1) * (nz_Bx) + iz] =
    //       Bx[ix * (ny_Bx * nz_Bx) + 1 * (nz_Bx) + iz];
    //   }
    // }
    // Bz (d,d,p)
    // for (unsigned int ix = 0; ix < nx_Bz; ix++) {
    //   for (unsigned int iz = 0; iz < nz_Bz; iz++) {
    //     // -Y
    //     Bz[ix * (ny_Bz * nz_Bz) + 0 * (nz_Bz) + iz] =
    //       Bz[ix * (ny_Bz * nz_Bz) + (ny_Bz - 2) * (nz_Bz) + iz];
    //     Bz[ix * (ny_Bz * nz_Bz) + (ny_Bz - 1) * (nz_Bz) + iz] =
    //       Bz[ix * (ny_Bz * nz_Bz) + 1 * (nz_Bz) + iz];
    //   }
    // }

    // Z dim
    // Bx
    //   for (unsigned int ix = 0; ix < nx_Bx; ix++) {
    //     for (unsigned int iy = 0; iy < ny_Bx; iy++) {
    //       // -Z
    //       Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + 0] =
    //         Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + (nz_Bx - 2)];
    //       Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + (nz_Bx - 1)] =
    //         Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + 1];
    //     }
    //   }
    //   // By
    //   for (unsigned int ix = 0; ix < nx_By; ++ix) {
    //     for (unsigned int iy = 0; iy < ny_By; ++iy) {
    //       // -Z
    //       By[ix * (ny_By * nz_By) + iy * (nz_By) + 0] =
    //         By[ix * (ny_By * nz_By) + iy * (nz_By) + (nz_By - 2)];
    //       By[ix * (ny_By * nz_By) + iy * (nz_By) + (nz_By - 1)] =
    //         By[ix * (ny_By * nz_By) + iy * (nz_By) + 1];
    //     }
    //   }
    // }); // end sycl parallel_for
    // q->wait();

  } else if (params.boundary_condition == "reflective") {

    q->single_task([=]() {
      // X dim
      // By (d,p,d)
      for (unsigned int iy = 0; iy < ny_By; ++iy) {
        for (unsigned int iz = 0; iz < nz_By; ++iz) {
          // -X
          By[0 + iy * (nz_By) + iz] = By[1 * (ny_By * nz_By) + iy * (nz_By) + iz];
          By[(nx_By - 1) * (ny_By * nz_By) + iy * (nz_By) + iz] =
            By[(nx_By - 2) * (ny_By * nz_By) + iy * (nz_By) + iz];
        }
      }

      // Bz (d,d,p)
      for (unsigned int iy = 0; iy < ny_Bz; iy++) {
        for (unsigned int iz = 0; iz < nz_Bz; iz++) {
          // -X
          Bz[0 + iy * (nz_Bz) + iz] = Bz[1 * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz];
          Bz[(nx_Bz - 1) * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz] =
            Bz[(nx_Bz - 2) * (ny_Bz * nz_Bz) + iy * (nz_Bz) + iz];
        }
      }

      // Y dim
      // Bx (p,d,d)
      for (unsigned int ix = 0; ix < nx_Bx; ix++) {
        for (unsigned int iz = 0; iz < nz_Bx; iz++) {
          // -Y
          Bx[ix * (ny_Bx * nz_Bx) + 0 + iz] = Bx[ix * (ny_Bx * nz_Bx) + 1 * (nz_Bx) + iz];
          Bx[ix * (ny_Bx * nz_Bx) + (ny_Bx - 1) * (nz_Bx) + iz] =
            Bx[ix * (ny_Bx * nz_Bx) + (ny_Bx - 2) * (nz_Bx) + iz];
        }
      }
      // Bz (-1 to avoid corner)
      for (unsigned int ix = 0; ix < nx_Bz; ix++) {
        for (unsigned int iz = 0; iz < nz_Bz; ++iz) {
          // -Y
          Bz[ix * (ny_Bz * nz_Bz) + 0 + iz] = Bz[ix * (ny_Bz * nz_Bz) + 1 * (nz_Bz) + iz];
          Bz[ix * (ny_Bz * nz_Bz) + (ny_Bz - 1) * (nz_Bz) + iz] =
            Bz[ix * (ny_Bz * nz_Bz) + (ny_Bz - 2) * (nz_Bz) + iz];
        }
      }

      // Z dim
      // Bx
      for (unsigned int ix = 0; ix < nx_Bx; ix++) {
        for (unsigned int iy = 0; iy < ny_Bx; iy++) {
          // -Z
          Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + 0] = Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + 1];
          Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + (nz_Bx - 1)] =
            Bx[ix * (ny_Bx * nz_Bx) + iy * (nz_Bx) + (nz_Bx - 2)];
        }
      }
      // By
      for (unsigned int ix = 0; ix < nx_By; ix++) {
        for (unsigned int iy = 0; iy < ny_By; iy++) {
          // -Z
          By[ix * (ny_By * nz_By) + iy * (nz_By) + 0] = By[ix * (ny_By * nz_By) + iy * (nz_By) + 1];
          By[ix * (ny_By * nz_By) + iy * (nz_By) + (nz_By - 1)] =
            By[ix * (ny_By * nz_By) + iy * (nz_By) + (nz_By - 2)];
        }
      }
    }); // end sycl parallel_for
    q->wait();

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
void identify_particles_to_move(Params &params, Patch &patch, Backend &backend) {}

// ___________________________________________________________
//
//! \brief Get the particles from neighbors
//! \param[in] params constant global parameters
//! \param[in] vec_patch vector of all patches
// ___________________________________________________________
auto exchange_particles(Params &params, std::vector<Patch> &vec_patch, int id_patch) -> void {}

// ______________________________________________________
//
//! \brief Sum all species local current grids in local grid
//! \param[in] patch  current patch to handle
// ______________________________________________________
auto reduc_current(Patch &patch, Backend &backend) -> void {

  for (int is = 1; is < patch.n_species_m; is++) {

    // Only if particles projected
    if (patch.projected_[is]) {

      sycl::queue *q = backend.sycl_queue_;

      // Jx

      const int nx_Jx = patch.vec_Jx_m[is].nx_m, ny_Jx = patch.vec_Jx_m[is].ny_m,
                nz_Jx = patch.vec_Jx_m[is].nz_m;

      sycl::range<3> global_Jx{static_cast<size_t>(nx_Jx),
                               static_cast<size_t>(ny_Jx),
                               static_cast<size_t>(nz_Jx)};
      sycl::range<3> local_Jx{1, 1, static_cast<size_t>(nz_Jx)};

      mini_float *const __restrict__ Jx_0 = patch.vec_Jx_m[0].get_raw_pointer(minipic::device);
      const mini_float *const __restrict__ Jx_is =
        patch.vec_Jx_m[is].get_raw_pointer(minipic::device);

      q->parallel_for(sycl::nd_range{global_Jx, local_Jx}, [=](sycl::nd_item<3> item) {
        size_t ix = item.get_global_id(0);
        size_t iy = item.get_global_id(1);
        size_t iz = item.get_global_id(2);

        Jx_0[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz] +=
          Jx_is[ix * (ny_Jx * nz_Jx) + iy * (nz_Jx) + iz];
      });
      q->wait();

      // Jy

      const int nx_Jy = patch.vec_Jy_m[is].nx_m, ny_Jy = patch.vec_Jy_m[is].ny_m,
                nz_Jy = patch.vec_Jy_m[is].nz_m;

      sycl::range<3> global_Jy{static_cast<size_t>(nx_Jy),
                               static_cast<size_t>(ny_Jy),
                               static_cast<size_t>(nz_Jy)};
      sycl::range<3> local_Jy{1, 1, static_cast<size_t>(nz_Jy)};

      mini_float *const __restrict__ Jy_0 = patch.vec_Jy_m[0].get_raw_pointer(minipic::device);
      const mini_float *const __restrict__ Jy_is =
        patch.vec_Jy_m[is].get_raw_pointer(minipic::device);

      q->parallel_for(sycl::nd_range{global_Jy, local_Jy}, [=](sycl::nd_item<3> item) {
        size_t ix = item.get_global_id(0);
        size_t iy = item.get_global_id(1);
        size_t iz = item.get_global_id(2);

        Jy_0[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz] +=
          Jy_is[ix * (ny_Jy * nz_Jy) + iy * (nz_Jy) + iz];
      });
      q->wait();

      // Jz

      const int nx_Jz = patch.vec_Jz_m[is].nx_m, ny_Jz = patch.vec_Jz_m[is].ny_m,
                nz_Jz = patch.vec_Jz_m[is].nz_m;

      sycl::range<3> global_Jz{static_cast<size_t>(nx_Jz),
                               static_cast<size_t>(ny_Jz),
                               static_cast<size_t>(nz_Jz)};
      sycl::range<3> local_Jz{1, 1, static_cast<size_t>(nz_Jz)};

      mini_float *const __restrict__ Jz_0 = patch.vec_Jz_m[0].get_raw_pointer(minipic::device);
      const mini_float *const __restrict__ Jz_is =
        patch.vec_Jz_m[is].get_raw_pointer(minipic::device);

      q->parallel_for(sycl::nd_range{global_Jz, local_Jz}, [=](sycl::nd_item<3> item) {
        size_t ix = item.get_global_id(0);
        size_t iy = item.get_global_id(1);
        size_t iz = item.get_global_id(2);

        Jz_0[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz] +=
          Jz_is[ix * (ny_Jz * nz_Jz) + iy * (nz_Jz) + iz];
      });
      q->wait();

    } // end if projected
  } // end for species
}

// ____________________________________________________________________________
//! \brief Copy all local current grid in the global grid
//! \param[in] ElectroMagn & em - global electromagnetic fields
//! \param[in] Patch & patch - current patch
// ____________________________________________________________________________
auto local2global(ElectroMagn &em, Patch &patch, Backend &backend) -> void {

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

    sycl::queue *q = backend.sycl_queue_;

    // Jx

    const int nx_Jx = em.Jx_m.nx(), ny_Jx = em.Jx_m.ny(), nz_Jx = em.Jx_m.nz();
    const int nx_Jx_0 = patch.vec_Jx_m[0].nx(), ny_Jx_0 = patch.vec_Jx_m[0].ny(),
              nz_Jx_0 = patch.vec_Jx_m[0].nz();

    const mini_float *const __restrict__ Jx_0 = patch.vec_Jx_m[0].get_raw_pointer(minipic::device);
    mini_float *const __restrict__ Jx         = em.Jx_m.get_raw_pointer(minipic::device);

    sycl::range<3> global_Jx{static_cast<size_t>(nx_Jx_0),
                             static_cast<size_t>(ny_Jx_0),
                             static_cast<size_t>(nz_Jx_0)};
    sycl::range<3> local_Jx{1, 1, static_cast<size_t>(nz_Jx_0)};

    q->parallel_for(sycl::nd_range{global_Jx, local_Jx}, [=](sycl::nd_item<3> item) {
      size_t ix = item.get_global_id(0);
      size_t iy = item.get_global_id(1);
      size_t iz = item.get_global_id(2);

      Jx[(i_global_d + ix) * (ny_Jx * nz_Jx) + (j_global_p + iy) * (nz_Jx) + (k_global_p + iz)] +=
        Jx_0[ix * (ny_Jx_0 * nz_Jx_0) + iy * (nz_Jx_0) + iz];
    });
    q->wait();

    // Jy

    const int nx_Jy = em.Jy_m.nx(), ny_Jy = em.Jy_m.ny(), nz_Jy = em.Jy_m.nz();
    const int nx_Jy_0 = patch.vec_Jy_m[0].nx(), ny_Jy_0 = patch.vec_Jy_m[0].ny(),
              nz_Jy_0 = patch.vec_Jy_m[0].nz();

    const mini_float *const __restrict__ Jy_0 = patch.vec_Jy_m[0].get_raw_pointer(minipic::device);
    mini_float *const __restrict__ Jy         = em.Jy_m.get_raw_pointer(minipic::device);

    sycl::range<3> global_Jy{static_cast<size_t>(nx_Jy_0),
                             static_cast<size_t>(ny_Jy_0),
                             static_cast<size_t>(nz_Jy_0)};
    sycl::range<3> local_Jy{1, 1, static_cast<size_t>(nz_Jy_0)};

    q->parallel_for(
      sycl::nd_range{global_Jy, local_Jy},
      [=](sycl::nd_item<3> item) {
        size_t ix = item.get_global_id(0);
        size_t iy = item.get_global_id(1);
        size_t iz = item.get_global_id(2);

        Jy[(i_global_p + ix) * (ny_Jy * nz_Jy) + (j_global_d + iy) * (nz_Jy) + (k_global_p + iz)] +=
          Jy_0[ix * (ny_Jy_0 * nz_Jy_0) + iy * (nz_Jy_0) + iz];
      }

    );
    q->wait();

    // Jz

    const int nx_Jz = em.Jz_m.nx(), ny_Jz = em.Jz_m.ny(), nz_Jz = em.Jz_m.nz();
    const int nx_Jz_0 = patch.vec_Jz_m[0].nx(), ny_Jz_0 = patch.vec_Jz_m[0].ny(),
              nz_Jz_0 = patch.vec_Jz_m[0].nz();

    const mini_float *const __restrict__ Jz_0 = patch.vec_Jz_m[0].get_raw_pointer(minipic::device);
    mini_float *const __restrict__ Jz         = em.Jz_m.get_raw_pointer(minipic::device);

    sycl::range<3> global_Jz{static_cast<size_t>(nx_Jz_0),
                             static_cast<size_t>(ny_Jz_0),
                             static_cast<size_t>(nz_Jz_0)};
    sycl::range<3> local_Jz{1, 1, static_cast<size_t>(nz_Jz_0)};

    q->parallel_for(sycl::nd_range{global_Jz, local_Jz}, [=](sycl::nd_item<3> item) {
      size_t ix = item.get_global_id(0);
      size_t iy = item.get_global_id(1);
      size_t iz = item.get_global_id(2);

      Jz[(i_global_p + ix) * (ny_Jz * nz_Jz) + (j_global_p + iy) * (nz_Jz) + (k_global_d + iz)] +=
        Jz_0[ix * (ny_Jz_0 * nz_Jz_0) + iy * (nz_Jz_0) + iz];
    });
    q->wait();

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