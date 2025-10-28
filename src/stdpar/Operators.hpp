/* _____________________________________________________________________ */
//! \file Operators.hpp

//! \brief contains generic kernels for the particle pusher

/* _____________________________________________________________________ */

#include "ElectroMagn.hpp"
#include "Patch.hpp"
#include "Profiler.hpp"

#ifndef OPERATORS_H
#define OPERATORS_H

#ifdef __MINIPIC_STDPAR_CPU__
const std::execution::parallel_unsequenced_policy &policy = std::execution::par_unseq;
// const std::execution::parallel_policy &policy = std::execution::par;
#elif __MINIPIC_STDPAR__
const std::execution::parallel_unsequenced_policy &policy = std::execution::par_unseq;
#endif

namespace operators {

// ______________________________________________________________________________
//
//! \brief Interpolation operator at the patch level :
//! interpolate EM fields from global grid for each particle
//! \param[in] em  global electromagnetic fields
//! \param[in] patch  patch data structure
// ______________________________________________________________________________
auto interpolate(ElectroMagn &em, Patch &patch) -> void {
  const mini_float inv_dx_m = em.inv_dx_m;
  const mini_float inv_dy_m = em.inv_dy_m;
  const mini_float inv_dz_m = em.inv_dz_m;

  const auto ny_Ex = em.Ex_m.ny(), nz_Ex = em.Ex_m.nz(), nynz_Ex = nz_Ex * ny_Ex;
  const auto ny_Ey = em.Ey_m.ny(), nz_Ey = em.Ey_m.nz(), nynz_Ey = nz_Ey * ny_Ey;
  const auto ny_Ez = em.Ez_m.ny(), nz_Ez = em.Ez_m.nz(), nynz_Ez = nz_Ez * ny_Ez;

  const auto ny_Bx = em.Bx_m.ny(), nz_Bx = em.Bx_m.nz(), nynz_Bx = nz_Bx * ny_Bx;
  const auto ny_By = em.By_m.ny(), nz_By = em.By_m.nz(), nynz_By = nz_By * ny_By;
  const auto ny_Bz = em.Bz_m.ny(), nz_Bz = em.Bz_m.nz(), nynz_Bz = nz_Bz * ny_Bz;

  const mini_float * const __restrict__ Ex = em.Ex_m.get_raw_pointer(minipic::device);
  const mini_float * const __restrict__ Ey = em.Ey_m.get_raw_pointer(minipic::device);
  const mini_float * const __restrict__ Ez = em.Ez_m.get_raw_pointer(minipic::device);

  const mini_float * const __restrict__ Bx = em.Bx_m.get_raw_pointer(minipic::device);
  const mini_float * const __restrict__ By = em.By_m.get_raw_pointer(minipic::device);
  const mini_float * const __restrict__ Bz = em.Bz_m.get_raw_pointer(minipic::device);

  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

    //std::vector<int> res(n_particles);

    const mini_float * const __restrict__ x = patch.particles_m[is].x_.get_raw_pointer(minipic::device);
    const mini_float * const __restrict__ y = patch.particles_m[is].y_.get_raw_pointer(minipic::device);
    const mini_float * const __restrict__ z = patch.particles_m[is].z_.get_raw_pointer(minipic::device);

    mini_float * const __restrict__ Exp = patch.particles_m[is].Ex_.get_raw_pointer(minipic::device);
    mini_float * const __restrict__ Eyp = patch.particles_m[is].Ey_.get_raw_pointer(minipic::device);
    mini_float * const __restrict__ Ezp = patch.particles_m[is].Ez_.get_raw_pointer(minipic::device);

    mini_float * const __restrict__ Bxp = patch.particles_m[is].Bx_.get_raw_pointer(minipic::device);
    mini_float * const __restrict__ Byp = patch.particles_m[is].By_.get_raw_pointer(minipic::device);
    mini_float * const __restrict__ Bzp = patch.particles_m[is].Bz_.get_raw_pointer(minipic::device);

    std::for_each(policy, counting_iterator<size_t>(0), counting_iterator<size_t>(n_particles), [=](size_t ip) {
      // // Calculate normalized positions
      const mini_float ixn = x[ip] * inv_dx_m;
      const mini_float iyn = y[ip] * inv_dy_m;
      const mini_float izn = z[ip] * inv_dz_m;

      // // Compute indexes in global primal grid
      const unsigned int ixp = static_cast<unsigned int>(floor(ixn));
      const unsigned int iyp = static_cast<unsigned int>(floor(iyn));
      const unsigned int izp = static_cast<unsigned int>(floor(izn));

      // Compute indexes in global dual grid
      const unsigned int ixd = static_cast<unsigned int>(floor(ixn + 0.5));
      const unsigned int iyd = static_cast<unsigned int>(floor(iyn + 0.5));
      const unsigned int izd = static_cast<unsigned int>(floor(izn + 0.5));

      // Compute interpolation coeff, p = primal, d = dual
      // interpolation electric field
      // Ex (d, p , p)
      {
        const mini_float coeffs[3] = {ixn + 0.5, iyn, izn};

        const auto v00 = Ex[ixd * nynz_Ex + iyp * nz_Ex + izp] * (1 - coeffs[0]) +
                         Ex[(ixd + 1) * nynz_Ex + iyp * nz_Ex + izp] * coeffs[0];
        const auto v01 = Ex[ixd * nynz_Ex + iyp * nz_Ex + izp + 1] * (1 - coeffs[0]) +
                         Ex[(ixd + 1) * nynz_Ex + iyp * nz_Ex + izp + 1] * coeffs[0];
        const auto v10 = Ex[ixd * nynz_Ex + (iyp + 1) * nz_Ex + izp] * (1 - coeffs[0]) +
                         Ex[(ixd + 1) * nynz_Ex + (iyp + 1) * nz_Ex + izp] * coeffs[0];
        const auto v11 = Ex[ixd * nynz_Ex + (iyp + 1) * nz_Ex + izp + 1] * (1 - coeffs[0]) +
                         Ex[(ixd + 1) * nynz_Ex + (iyp + 1) * nz_Ex + izp + 1] * coeffs[0];

        const auto v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const auto v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Exp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // Ey (p, d, p)
      {
        const mini_float coeffs[3] = {ixn, iyn + 0.5, izn};

        const mini_float v00 = Ey[ixp * nynz_Ey + iyd * nz_Ey + izp] * (1 - coeffs[0]) +
                               Ey[(ixp + 1) * nynz_Ey + iyd * nz_Ey + izp] * coeffs[0];
        const mini_float v01 = Ey[ixp * nynz_Ey + iyd * nz_Ey + izp + 1] * (1 - coeffs[0]) +
                               Ey[(ixp + 1) * nynz_Ey + iyd * nz_Ey + izp + 1] * coeffs[0];
        const mini_float v10 = Ey[ixp * nynz_Ey + (iyd + 1) * nz_Ey + izp] * (1 - coeffs[0]) +
                               Ey[(ixp + 1) * nynz_Ey + (iyd + 1) * nz_Ey + izp] * coeffs[0];
        const mini_float v11 = Ey[ixp * nynz_Ey + (iyd + 1) * nz_Ey + izp + 1] * (1 - coeffs[0]) +
                               Ey[(ixp + 1) * nynz_Ey + (iyd + 1) * nz_Ey + izp + 1] * coeffs[0];

        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Eyp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // Ez (p, p, d)
      {
        const mini_float coeffs[3] = {ixn, iyn, izn + 0.5};

        const mini_float v00 = Ez[ixp * nynz_Ez + iyp * nz_Ez + izd] * (1 - coeffs[0]) +
                               Ez[(ixp + 1) * nynz_Ez + iyp * nz_Ez + izd] * coeffs[0];
        const mini_float v01 = Ez[ixp * nynz_Ez + iyp * nz_Ez + izd + 1] * (1 - coeffs[0]) +
                               Ez[(ixp + 1) * nynz_Ez + iyp * nz_Ez + izd + 1] * coeffs[0];
        const mini_float v10 = Ez[ixp * nynz_Ez + (iyp + 1) * nz_Ez + izd] * (1 - coeffs[0]) +
                               Ez[(ixp + 1) * nynz_Ez + (iyp + 1) * nz_Ez + izd] * coeffs[0];
        const mini_float v11 = Ez[ixp * nynz_Ez + (iyp + 1) * nz_Ez + izd + 1] * (1 - coeffs[0]) +
                               Ez[(ixp + 1) * nynz_Ez + (iyp + 1) * nz_Ez + izd + 1] * coeffs[0];

        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Ezp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // interpolation magnetic field
      // Bx [p, d, d)
      {
        const mini_float coeffs[3] = {ixn, iyn + 0.5, izn + 0.5};

        const mini_float v00 = Bx[ixp * nynz_Bx + iyd * nz_Bx + izd] * (1 - coeffs[0]) +
                               Bx[(ixp + 1) * nynz_Bx + iyd * nz_Bx + izd] * coeffs[0];
        const mini_float v01 = Bx[ixp * nynz_Bx + iyd * nz_Bx + izd + 1] * (1 - coeffs[0]) +
                               Bx[(ixp + 1) * nynz_Bx + iyd * nz_Bx + izd + 1] * coeffs[0];
        const mini_float v10 = Bx[ixp * nynz_Bx + (iyd + 1) * nz_Bx + izd] * (1 - coeffs[0]) +
                               Bx[(ixp + 1) * nynz_Bx + (iyd + 1) * nz_Bx + izd] * coeffs[0];
        const mini_float v11 = Bx[ixp * nynz_Bx + (iyd + 1) * nz_Bx + izd + 1] * (1 - coeffs[0]) +
                               Bx[(ixp + 1) * nynz_Bx + (iyd + 1) * nz_Bx + izd + 1] * coeffs[0];

        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Bxp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // By [d, p, d)
      {
        const mini_float coeffs[3] = {ixn + 0.5, iyn, izn + 0.5};

        const mini_float v00 = By[ixd * nynz_By + iyp * nz_By + izd] * (1 - coeffs[0]) +
                               By[(ixd + 1) * nynz_By + iyp * nz_By + izd] * coeffs[0];
        const mini_float v01 = By[ixd * nynz_By + iyp * nz_By + izd + 1] * (1 - coeffs[0]) +
                               By[(ixd + 1) * nynz_By + iyp * nz_By + izd + 1] * coeffs[0];
        const mini_float v10 = By[ixd * nynz_By + (iyp + 1) * nz_By + izd] * (1 - coeffs[0]) +
                               By[(ixd + 1) * nynz_By + (iyp + 1) * nz_By + izd] * coeffs[0];
        const mini_float v11 = By[ixd * nynz_By + (iyp + 1) * nz_By + izd + 1] * (1 - coeffs[0]) +
                               By[(ixd + 1) * nynz_By + (iyp + 1) * nz_By + izd + 1] * coeffs[0];

        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Byp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // Bz (d, d, p)
      {
        const mini_float coeffs[3] = {ixn + 0.5, iyn + 0.5, izn};

        const mini_float v00 = Bz[ixd * nynz_Bz + iyd * nz_Bz + izp] * (1 - coeffs[0]) +
                               Bz[(ixd + 1) * nynz_Bz + iyd * nz_Bz + izp] * coeffs[0];
        const mini_float v01 = Bz[ixd * nynz_Bz + iyd * nz_Bz + izp + 1] * (1 - coeffs[0]) +
                               Bz[(ixd + 1) * nynz_Bz + iyd * nz_Bz + izp + 1] * coeffs[0];
        const mini_float v10 = Bz[ixd * nynz_Bz + (iyd + 1) * nz_Bz + izp] * (1 - coeffs[0]) +
                               Bz[(ixd + 1) * nynz_Bz + (iyd + 1) * nz_Bz + izp] * coeffs[0];
        const mini_float v11 = Bz[ixd * nynz_Bz + (iyd + 1) * nz_Bz + izp + 1] * (1 - coeffs[0]) +
                               Bz[(ixd + 1) * nynz_Bz + (iyd + 1) * nz_Bz + izp + 1] * coeffs[0];

        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        Bzp[ip] = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }
    });

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
    const mini_float qp = patch.particles_m[is].charge_m * dt * 0.5 / patch.particles_m[is].mass_m;

    mini_float *Exp = patch.particles_m[is].Ex_.get_raw_pointer(minipic::device);
    mini_float *Eyp = patch.particles_m[is].Ey_.get_raw_pointer(minipic::device);
    mini_float *Ezp = patch.particles_m[is].Ez_.get_raw_pointer(minipic::device);

    mini_float *Bxp = patch.particles_m[is].Bx_.get_raw_pointer(minipic::device);
    mini_float *Byp = patch.particles_m[is].By_.get_raw_pointer(minipic::device);
    mini_float *Bzp = patch.particles_m[is].Bz_.get_raw_pointer(minipic::device);

    mini_float *mxp = patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
    mini_float *myp = patch.particles_m[is].my_.get_raw_pointer(minipic::device);
    mini_float *mzp = patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

    mini_float *xp = patch.particles_m[is].x_.get_raw_pointer(minipic::device);
    mini_float *yp = patch.particles_m[is].y_.get_raw_pointer(minipic::device);
    mini_float *zp = patch.particles_m[is].z_.get_raw_pointer(minipic::device);

    std::for_each(policy, counting_iterator<size_t>(0), counting_iterator<size_t>(n_particles), [=](size_t ip) {
      mini_float px = qp * Exp[ip];
      mini_float py = qp * Eyp[ip];
      mini_float pz = qp * Ezp[ip];

      const mini_float ux = mxp[ip] + px;
      const mini_float uy = myp[ip] + py;
      const mini_float uz = mzp[ip] + pz;

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
      mxp[ip] = px;
      myp[ip] = py;
      mzp[ip] = pz;

      // Update positions
      xp[ip] += px * dt * gamma_inv;
      yp[ip] += py * dt * gamma_inv;
      zp[ip] += pz * dt * gamma_inv;
    });

  } // Loop on species
}

// ______________________________________________________________________________
//
//! \brief Push only the momentum
//! \param[in] patch  patch data structure
//! \param[in] dt time step to use for the pusher
// ______________________________________________________________________________
// NOT TESTED because momentum_correction is false by default. In Default_gpu and antenna setups,
// the momentum_correction is false
auto push_momentum(Patch &patch, double dt) -> void {

  // for each species
  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

    // q' = dt * (q/2m)
    const mini_float qp = patch.particles_m[is].charge_m * dt * 0.5 / patch.particles_m[is].mass_m;

    mini_float *Exp = patch.particles_m[is].Ex_.get_raw_pointer(minipic::device);
    mini_float *Eyp = patch.particles_m[is].Ey_.get_raw_pointer(minipic::device);
    mini_float *Ezp = patch.particles_m[is].Ez_.get_raw_pointer(minipic::device);

    mini_float *Bxp = patch.particles_m[is].Bx_.get_raw_pointer(minipic::device);
    mini_float *Byp = patch.particles_m[is].By_.get_raw_pointer(minipic::device);
    mini_float *Bzp = patch.particles_m[is].Bz_.get_raw_pointer(minipic::device);

    mini_float *mxp = patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
    mini_float *myp = patch.particles_m[is].my_.get_raw_pointer(minipic::device);
    mini_float *mzp = patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

    std::for_each(policy, counting_iterator<size_t>(0), counting_iterator<size_t>(n_particles), [=](size_t ip) {
      // 1/2 E
      mini_float px = qp * Exp[ip];
      mini_float py = qp * Eyp[ip];
      mini_float pz = qp * Ezp[ip];

      const mini_float ux = mxp[ip] + px;
      const mini_float uy = myp[ip] + py;
      const mini_float uz = mzp[ip] + pz;

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
      mxp[ip] = px;
      myp[ip] = py;
      mzp[ip] = pz;
    }); // End for each particles
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

    const mini_float inf_global[3] = {params.inf_x, params.inf_y, params.inf_z};
    const mini_float sup_global[3] = {params.sup_x, params.sup_y, params.sup_z};

    // Periodic conditions
    if (params.boundary_condition_code == 1) {

      const int N_patches[3]     = {patch.nx_patchs_m, patch.ny_patchs_m, patch.nz_patchs_m};
      const mini_float length[3] = {params.Lx, params.Ly, params.Lz};

      for (int is = 0; is < patch.n_species_m; is++) {
        size_t n_particles = patch.particles_m[is].size();

        mini_float *xp = patch.particles_m[is].x_.get_raw_pointer(minipic::device);
        mini_float *yp = patch.particles_m[is].y_.get_raw_pointer(minipic::device);
        mini_float *zp = patch.particles_m[is].z_.get_raw_pointer(minipic::device);

        std::for_each(policy, counting_iterator<size_t>(0), counting_iterator<size_t>(n_particles), [=](size_t ip) {
          mini_float *pos[3] = {&xp[ip], &yp[ip], &zp[ip]};

          for (int d = 0; d < 3; d++) {
            // Only relevant if there is just 1 patch in this direction
            // Else the patch exchange with periodicity is managed in the dedicated function
            if (N_patches[d] == 1) {
              if (*pos[d] >= sup_global[d]) {
                *pos[d] -= length[d];
              } else if (*pos[d] < inf_global[d]) {
                *pos[d] += length[d];
              }
            }
          }
        }); // End loop on particles

      } // End loop on species

      // Reflective conditions
      // NOT TESTED because antenna, default and beam are periodic so boundary_condition_code equals
      // to 1 and not 2
    } else if (params.boundary_condition_code == 2) {
      for (int is = 0; is < patch.n_species_m; is++) {

        size_t n_particles = patch.particles_m[is].size();

        mini_float *xp = patch.particles_m[is].x_.get_raw_pointer(minipic::device);
        mini_float *yp = patch.particles_m[is].y_.get_raw_pointer(minipic::device);
        mini_float *zp = patch.particles_m[is].z_.get_raw_pointer(minipic::device);

        mini_float *mxp = patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
        mini_float *myp = patch.particles_m[is].my_.get_raw_pointer(minipic::device);
        mini_float *mzp = patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

        std::for_each(policy, counting_iterator<size_t>(0), counting_iterator<size_t>(n_particles), [=](size_t ip) {
          mini_float *pos[3]      = {&xp[ip], &yp[ip], &zp[ip]};
          mini_float *momentum[3] = {&mxp[ip], &myp[ip], &mzp[ip]};

          for (int d = 0; d < 3; d++) {
            if (*pos[d] >= sup_global[d]) {

              *pos[d]      = 2 * sup_global[d] - *pos[d];
              *momentum[d] = -*momentum[d];

            } else if (*pos[d] < inf_global[d]) {

              *pos[d]      = 2 * inf_global[d] - *pos[d];
              *momentum[d] = -*momentum[d];
            }
          }
        }); // End loop on particles
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

    const size_t n_particles = patch.particles_m[is].size();
    if (n_particles > 0) {
      const int nx_Jx = patch.vec_Jx_m[is].nx_m, ny_Jx = patch.vec_Jx_m[is].ny_m,
                nz_Jx = patch.vec_Jx_m[is].nz_m, nynz_Jx = ny_Jx * nz_Jx;
      const int nx_Jy = patch.vec_Jy_m[is].nx_m, ny_Jy = patch.vec_Jy_m[is].ny_m,
                nz_Jy = patch.vec_Jy_m[is].nz_m, nynz_Jy = ny_Jy * nz_Jy;
      const int nx_Jz = patch.vec_Jz_m[is].nx_m, ny_Jz = patch.vec_Jz_m[is].ny_m,
                nz_Jz = patch.vec_Jz_m[is].nz_m, nynz_Jz = ny_Jz * nz_Jz;

      const mini_float inv_cell_volume_x_q =
        params.inv_cell_volume * patch.particles_m[is].charge_m;
      const mini_float dt = params.dt;

      const mini_float inv_dx = params.inv_dx;
      const mini_float inv_dy = params.inv_dy;
      const mini_float inv_dz = params.inv_dz;

      const mini_float xmin = patch.inf_m[0];
      const mini_float ymin = patch.inf_m[1];
      const mini_float zmin = patch.inf_m[2];

      mini_float * const __restrict__ Jx = patch.vec_Jx_m[is].get_raw_pointer(minipic::device);
      mini_float * const __restrict__ Jy = patch.vec_Jy_m[is].get_raw_pointer(minipic::device);
      mini_float * const __restrict__ Jz = patch.vec_Jz_m[is].get_raw_pointer(minipic::device);

      const mini_float * const __restrict__ w = patch.particles_m[is].weight_.get_raw_pointer(minipic::device);

      const mini_float * const __restrict__ x = patch.particles_m[is].x_.get_raw_pointer(minipic::device);
      const mini_float * const __restrict__ y = patch.particles_m[is].y_.get_raw_pointer(minipic::device);
      const mini_float * const __restrict__ z = patch.particles_m[is].z_.get_raw_pointer(minipic::device);

      const mini_float * const __restrict__ mx = patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
      const mini_float * const __restrict__ my = patch.particles_m[is].my_.get_raw_pointer(minipic::device);
      const mini_float * const __restrict__ mz = patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

      std::for_each(policy, counting_iterator<size_t>(0), counting_iterator<size_t>(n_particles), [=](size_t ip) {
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

#ifdef __MINIPIC_STDPAR__
        // For Jx
        {
          const mini_float coeffs[3] = {posxn - 0.5 - ixd, posyn - iyp, poszn - izp};

          // Project on Jx using atomicAdd

          atomicAdd(&Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) + izp],
                    (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp);
          atomicAdd(&Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) + (izp + 1)],
                    (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp);
          atomicAdd(&Jx[ixd * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp],
                    (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp);
          atomicAdd(&Jx[ixd * (nynz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)],
                    (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp);
          atomicAdd(&Jx[(ixd + 1) * (nynz_Jx) + iyp * (nz_Jx) + izp],
                    (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp);
          atomicAdd(&Jx[(ixd + 1) * (nynz_Jx) + iyp * (nz_Jx) + (izp + 1)],
                    (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp);
          atomicAdd(&Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp],
                    (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp);
          atomicAdd(&Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)],
                    (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp);
        }
        // For Jy
        {
          const mini_float coeffs[3] = {posxn - ixp, posyn - 0.5f - iyd, poszn - izp};

          // Project on Jy using atomicAdd
          atomicAdd(&Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) + izp],
                    (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp);
          atomicAdd(&Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) + (izp + 1)],
                    (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp);
          atomicAdd(&Jy[ixp * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp],
                    (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp);
          atomicAdd(&Jy[ixp * (nynz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)],
                    (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp);
          atomicAdd(&Jy[(ixp + 1) * (nynz_Jy) + iyd * (nz_Jy) + izp],
                    (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp);
          atomicAdd(&Jy[(ixp + 1) * (nynz_Jy) + iyd * (nz_Jy) + (izp + 1)],
                    (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp);
          atomicAdd(&Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp],
                    (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp);
          atomicAdd(&Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)],
                    (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp);
        }

        // For Jz
        {

          const mini_float coeffs[3] = {posxn - ixp, posyn - iyp, poszn - 0.5f - izd};

          // Project on Jz using atomicAdd
          atomicAdd(&Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) + izd],
                    (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp);
          atomicAdd(&Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) + (izd + 1)],
                    (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp);
          atomicAdd(&Jz[ixp * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd],
                    (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp);
          atomicAdd(&Jz[ixp * (nynz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)],
                    (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp);
          atomicAdd(&Jz[(ixp + 1) * (nynz_Jz) + iyp * (nz_Jz) + izd],
                    (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp);
          atomicAdd(&Jz[(ixp + 1) * (nynz_Jz) + iyp * (nz_Jz) + (izd + 1)],
                    (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp);
          atomicAdd(&Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd],
                    (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp);
          atomicAdd(&Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)],
                    (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp);
        }

#elif __MINIPIC_STDPAR_CPU__
      /*
          // Project on Jx
          {
          const mini_float coeffs[3] = {posxn - 0.5 - ixd, posyn - iyp, poszn - izp};
      
          Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) + izp]+=((1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp ) ;
          Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) + (izp + 1)] +=((1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp) ;
          Jx[ixd * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp] +=((1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp) ;
          Jx[ixd * (nynz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)] +=((1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp) ;
          Jx[(ixd + 1) * (nynz_Jx) + iyp * (nz_Jx) + izp] +=((coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp) ;
          Jx[(ixd + 1) * (nynz_Jx) + iyp * (nz_Jx) + (izp + 1)] +=((coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp) ;
          Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp] +=((coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp) ;
          Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)] +=((coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp) ;
        }

        //Project on Jy
        {
          const mini_float coeffs[3] = {posxn - ixp, posyn - 0.5f - iyd, poszn - izp};

          Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) + izp] +=((1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp) ;
          Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) + (izp + 1)] +=((1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp)  ;
          Jy[ixp * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp] +=((1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp) ;
          Jy[ixp * (nynz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)] +=((1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp) ;
          Jy[(ixp + 1) * (nynz_Jy) + iyd * (nz_Jy) + izp] +=((coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp) ;
          Jy[(ixp + 1) * (nynz_Jy) + iyd * (nz_Jy) + (izp + 1)] +=((coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp) ;
          Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp] +=((coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp) ;
          Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)] +=((coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp) ;

          }

          // Project on Jz
          {
            const mini_float coeffs[3] = {posxn - ixp, posyn - iyp, poszn - 0.5f - izd};

            Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) + izd] +=((1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp) ;
            Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) + (izd + 1)] += ((1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp) ;
            Jz[ixp * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd] +=((1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp) ;
            Jz[ixp * (nynz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)] +=((1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp) ;
            Jz[(ixp + 1) * (nynz_Jz) + iyp * (nz_Jz) + izd] +=((coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp) ;
            Jz[(ixp + 1) * (nynz_Jz) + iyp * (nz_Jz) + (izd + 1)] +=((coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp) ;
            Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd] +=((coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp) ;
            Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)] +=((coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp) ;
          } 
        */

          
            // Project on Jx
          {
          const mini_float coeffs[3] = {posxn - 0.5 - ixd, posyn - iyp, poszn - izp};
      
          reinterpret_cast<std::atomic<mini_float>&>(Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) + izp]).fetch_add(((1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp ), std::memory_order_relaxed);
          reinterpret_cast<std::atomic<mini_float>&>(Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) + (izp + 1)]).fetch_add(((1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp), std::memory_order_relaxed);
          reinterpret_cast<std::atomic<mini_float>&>(Jx[ixd * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp]).fetch_add(((1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp), std::memory_order_relaxed);
          reinterpret_cast<std::atomic<mini_float>&>(Jx[ixd * (nynz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)]).fetch_add(((1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp), std::memory_order_relaxed);
          reinterpret_cast<std::atomic<mini_float>&>(Jx[(ixd + 1) * (nynz_Jx) + iyp * (nz_Jx) + izp]).fetch_add(((coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp), std::memory_order_relaxed);
          reinterpret_cast<std::atomic<mini_float>&>(Jx[(ixd + 1) * (nynz_Jx) + iyp * (nz_Jx) + (izp + 1)]).fetch_add(((coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp), std::memory_order_relaxed);
          reinterpret_cast<std::atomic<mini_float>&>(Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp]).fetch_add(((coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp), std::memory_order_relaxed);
          reinterpret_cast<std::atomic<mini_float>&>(Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)]).fetch_add(((coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp), std::memory_order_relaxed);
          }

          //Project on Jy
          {
            const mini_float coeffs[3] = {posxn - ixp, posyn - 0.5f - iyd, poszn - izp};

            reinterpret_cast<std::atomic<mini_float>&>(Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) + izp]).fetch_add(((1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp), std::memory_order_relaxed);
            reinterpret_cast<std::atomic<mini_float>&>(Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) + (izp + 1)]).fetch_add(((1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp) , std::memory_order_relaxed);
            reinterpret_cast<std::atomic<mini_float>&>(Jy[ixp * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp]).fetch_add(((1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp), std::memory_order_relaxed);
            reinterpret_cast<std::atomic<mini_float>&>(Jy[ixp * (nynz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)]).fetch_add(((1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp), std::memory_order_relaxed);
            reinterpret_cast<std::atomic<mini_float>&>(Jy[(ixp + 1) * (nynz_Jy) + iyd * (nz_Jy) + izp]).fetch_add(((coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp), std::memory_order_relaxed);
            reinterpret_cast<std::atomic<mini_float>&>(Jy[(ixp + 1) * (nynz_Jy) + iyd * (nz_Jy) + (izp + 1)]).fetch_add(((coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp), std::memory_order_relaxed);
            reinterpret_cast<std::atomic<mini_float>&>(Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp]).fetch_add(((coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp), std::memory_order_relaxed);
            reinterpret_cast<std::atomic<mini_float>&>(Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)]).fetch_add(((coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp), std::memory_order_relaxed);
   
            }

            // Project on Jz
            {
              const mini_float coeffs[3] = {posxn - ixp, posyn - iyp, poszn - 0.5f - izd};

              reinterpret_cast<std::atomic<mini_float>&>(Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) + izd]).fetch_add(((1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp), std::memory_order_relaxed);
              reinterpret_cast<std::atomic<mini_float>&>(Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) + (izd + 1)]).fetch_add( ((1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp), std::memory_order_relaxed);
              reinterpret_cast<std::atomic<mini_float>&>(Jz[ixp * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd]).fetch_add(((1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp), std::memory_order_relaxed);
              reinterpret_cast<std::atomic<mini_float>&>(Jz[ixp * (nynz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)]).fetch_add(((1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp), std::memory_order_relaxed);
              reinterpret_cast<std::atomic<mini_float>&>(Jz[(ixp + 1) * (nynz_Jz) + iyp * (nz_Jz) + izd]).fetch_add(((coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp), std::memory_order_relaxed);
              reinterpret_cast<std::atomic<mini_float>&>(Jz[(ixp + 1) * (nynz_Jz) + iyp * (nz_Jz) + (izd + 1)]).fetch_add(((coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp), std::memory_order_relaxed);
              reinterpret_cast<std::atomic<mini_float>&>(Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd]).fetch_add(((coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp), std::memory_order_relaxed);
              reinterpret_cast<std::atomic<mini_float>&>(Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)]).fetch_add(((coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp), std::memory_order_relaxed);
            }

        // Project on Jx
        /* {
         const mini_float coeffs[3] = {posxn - 0.5 - ixd, posyn - iyp, poszn - izp};

         std::atomic_ref(Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) + izp])+=((1 - coeffs[0]) * (1 -
         coeffs[1]) * (1 - coeffs[2]) * Jxp ) ; std::atomic_ref(Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) +
         (izp + 1)])+=((1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp) ;
         std::atomic_ref(Jx[ixd * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp])+=((1 - coeffs[0]) *
         (coeffs[1]) * (1 - coeffs[2]) * Jxp) ; std::atomic_ref(Jx[ixd * (nynz_Jx) + (iyp + 1) *
         (nz_Jx) + (izp + 1)])+=((1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp) ;
         std::atomic_ref(Jx[(ixd + 1) * (nynz_Jx) + iyp * (nz_Jx) + izp])+=((coeffs[0]) * (1 -
         coeffs[1]) * (1 - coeffs[2]) * Jxp) ; std::atomic_ref(Jx[(ixd + 1) * (nynz_Jx) + iyp *
         (nz_Jx) + (izp + 1)])+=((coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp) ;
         std::atomic_ref(Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp])+=((coeffs[0]) *
         (coeffs[1]) * (1 - coeffs[2]) * Jxp) ; std::atomic_ref(Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1)
         * (nz_Jx) + (izp + 1)])+=((coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp) ;
         }

         //Project on Jy
         {
           const mini_float coeffs[3] = {posxn - ixp, posyn - 0.5f - iyd, poszn - izp};

           std::atomic_ref(Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) + izp])+=((1 - coeffs[0]) * (1 -
         coeffs[1]) * (1 - coeffs[2]) * Jyp) ; std::atomic_ref(Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) +
         (izp + 1)])+=((1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp)  ;
           std::atomic_ref(Jy[ixp * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp])+=((1 - coeffs[0]) *
         (coeffs[1]) * (1 - coeffs[2]) * Jyp) ; std::atomic_ref(Jy[ixp * (nynz_Jy) + (iyd + 1) *
         (nz_Jy) + (izp + 1)])+=((1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp) ;
           std::atomic_ref(Jy[(ixp + 1) * (nynz_Jy) + iyd * (nz_Jy) + izp])+=((coeffs[0]) * (1 -
         coeffs[1]) * (1 - coeffs[2]) * Jyp) ; std::atomic_ref(Jy[(ixp + 1) * (nynz_Jy) + iyd *
         (nz_Jy) + (izp + 1)])+=((coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp) ;
           std::atomic_ref(Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp])+=((coeffs[0]) *
         (coeffs[1]) * (1 - coeffs[2]) * Jyp) ; std::atomic_ref(Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1)
         * (nz_Jy) + (izp + 1)])+=((coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp) ;

           }

           // Project on Jz
           {
             const mini_float coeffs[3] = {posxn - ixp, posyn - iyp, poszn - 0.5f - izd};

             std::atomic_ref(Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) + izd])+=((1 - coeffs[0]) * (1 -
         coeffs[1]) * (1 - coeffs[2]) * Jzp) ; std::atomic_ref(Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) +
         (izd + 1)])+= ((1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp) ;
             std::atomic_ref(Jz[ixp * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd])+=((1 - coeffs[0]) *
         (coeffs[1]) * (1 - coeffs[2]) * Jzp) ; std::atomic_ref(Jz[ixp * (nynz_Jz) + (iyp + 1) *
         (nz_Jz) + (izd + 1)])+=((1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp) ;
             std::atomic_ref(Jz[(ixp + 1) * (nynz_Jz) + iyp * (nz_Jz) + izd])+=((coeffs[0]) * (1 -
         coeffs[1]) * (1 - coeffs[2]) * Jzp) ; std::atomic_ref(Jz[(ixp + 1) * (nynz_Jz) + iyp *
         (nz_Jz) + (izd + 1)])+=((coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp) ;
             std::atomic_ref(Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd])+=((coeffs[0]) *
         (coeffs[1]) * (1 - coeffs[2]) * Jzp) ; std::atomic_ref(Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1)
         * (nz_Jz) + (izd + 1)])+=((coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp) ;

               //std::atomic<double>& atomicJ = std::atomic<double>(J);
           }*/
#endif
      }); // end for each particles

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
// NOT TESTED
auto project(Params &params, ElectroMagn &em, Patch &patch) -> void {

  const double dt = params.dt;

  const double inv_dx = params.inv_dx;
  const double inv_dy = params.inv_dy;
  const double inv_dz = params.inv_dz;

  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles                = patch.particles_m[is].size();
    const mini_float inv_cell_volume_x_q = params.inv_cell_volume * patch.particles_m[is].charge_m;

    mini_float *w = patch.particles_m[is].weight_.get_raw_pointer(minipic::device);

    mini_float *x = patch.particles_m[is].x_.get_raw_pointer(minipic::device);
    mini_float *y = patch.particles_m[is].y_.get_raw_pointer(minipic::device);
    mini_float *z = patch.particles_m[is].z_.get_raw_pointer(minipic::device);

    mini_float *Jx = patch.vec_Jx_m[is].get_raw_pointer(minipic::device);
    mini_float *Jy = patch.vec_Jy_m[is].get_raw_pointer(minipic::device);
    mini_float *Jz = patch.vec_Jz_m[is].get_raw_pointer(minipic::device);

    mini_float *mx = patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
    mini_float *my = patch.particles_m[is].my_.get_raw_pointer(minipic::device);
    mini_float *mz = patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

    const int nx_Jx = patch.vec_Jx_m[is].nx_m, ny_Jx = patch.vec_Jx_m[is].ny_m,
              nz_Jx = patch.vec_Jx_m[is].nz_m, nynz_Jx = ny_Jx * nz_Jx;
    const int nx_Jy = patch.vec_Jy_m[is].nx_m, ny_Jy = patch.vec_Jy_m[is].ny_m,
              nz_Jy = patch.vec_Jy_m[is].nz_m, nynz_Jy = ny_Jy * nz_Jy;
    const int nx_Jz = patch.vec_Jz_m[is].nx_m, ny_Jz = patch.vec_Jz_m[is].ny_m,
              nz_Jz = patch.vec_Jz_m[is].nz_m, nynz_Jz = ny_Jz * nz_Jz;

    std::for_each(policy, counting_iterator<size_t>(0), counting_iterator<size_t>(n_particles), [=](size_t ip) {
      // Delete if already compute by Pusher
      // mini_float usq = (moment[0]*moment[0] + moment[1]*moment[1] + moment[2]*moment[2]);
      // mini_float gamma = sqrt(1+usq);
      // gamma_inv = 1/gamma;

      const mini_float charge_weight = inv_cell_volume_x_q * w[ip];

      const mini_float gamma_inv =
        1 / sqrt(1 + mx[ip] * mx[ip] + my[ip] * my[ip] + mz[ip] * mz[ip]);

      const mini_float vx = mx[ip] * gamma_inv;
      const mini_float vy = my[ip] * gamma_inv;
      const mini_float vz = mz[ip] * gamma_inv;

      const mini_float Jxp = vx * charge_weight;
      const mini_float Jyp = vy * charge_weight;
      const mini_float Jzp = vz * charge_weight;

      // Calculate normalized positions
      // We come back 1/2 time step back in time for the position because of the leap frog scheme
      // As a consequence, we also have `+ 1` because the current grids have 2 additional ghost
      // cells (1 the min and 1 at the max border) when the direction is primal
      const mini_float posxn = (x[ip] - 0.5 * dt * vx) * inv_dx + 1;
      const mini_float posyn = (y[ip] - 0.5 * dt * vy) * inv_dy + 1;
      const mini_float poszn = (z[ip] - 0.5 * dt * vz) * inv_dz + 1;

      // Compute indexes in primal grid
      const int ixp = static_cast<int>(floor(posxn)); //- i_patch_topology_m * nx_cells_m;
      const int iyp = static_cast<int>(floor(posyn)); //- j_patch_topology_m * ny_cells_m;
      const int izp = static_cast<int>(floor(poszn)); //- k_patch_topology_m * nz_cells_m;

      // Compute indexes in dual grid
      const int ixd = static_cast<int>(floor(posxn - 0.5)); //- i_patch_topology_m * nx_cells_m;
      const int iyd = static_cast<int>(floor(posyn - 0.5)); //- j_patch_topology_m * ny_cells_m;
      const int izd = static_cast<int>(floor(poszn - 0.5)); //- k_patch_topology_m * nz_cells_m;

      // Projection particle on currant field
      // Compute interpolation coeff, p = primal, d = dual
      {
        mini_float coeffs[3] = {posxn - 0.5 - ixd, posyn - iyp, poszn - izp};

        Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) + izp] +=
          (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx[ixd * (nynz_Jx) + iyp * (nz_Jx) + (izp + 1)] +=
          (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
        Jx[ixd * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp] +=
          (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx[ixd * (nynz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)] +=
          (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;
        Jx[(ixd + 1) * (nynz_Jx) + iyp * (nz_Jx) + izp] +=
          (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx[(ixd + 1) * (nynz_Jx) + iyp * (nz_Jx) + (izp + 1)] +=
          (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
        Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1) * (nz_Jx) + izp] +=
          (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx[(ixd + 1) * (nynz_Jx) + (iyp + 1) * (nz_Jx) + (izp + 1)] +=
          (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;
      }

      {
        mini_float coeffs[3] = {posxn - ixp, posyn - 0.5 - iyd, poszn - izp};

        Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) + izp] +=
          (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy[ixp * (nynz_Jy) + iyd * (nz_Jy) + (izp + 1)] +=
          (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
        Jy[ixp * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp] +=
          (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy[ixp * (nynz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)] +=
          (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;
        Jy[(ixp + 1) * (nynz_Jy) + iyd * (nz_Jy) + izp] +=
          (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy[(ixp + 1) * (nynz_Jy) + iyd * (nz_Jy) + (izp + 1)] +=
          (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
        Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1) * (nz_Jy) + izp] +=
          (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy[(ixp + 1) * (nynz_Jy) + (iyd + 1) * (nz_Jy) + (izp + 1)] +=
          (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;
      }

      {
        mini_float coeffs[3] = {posxn - ixp, posyn - iyp, poszn - 0.5 - izd};

        Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) + izd] +=
          (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz[ixp * (nynz_Jz) + iyp * (nz_Jz) + (izd + 1)] +=
          (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
        Jz[ixp * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd] +=
          (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz[ixp * (nynz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)] +=
          (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
        Jz[(ixp + 1) * (nynz_Jz) + iyp * (nz_Jz) + izd] +=
          (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz[(ixp + 1) * (nynz_Jz) + iyp * (nz_Jz) + (izd + 1)] +=
          (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
        Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1) * (nz_Jz) + izd] +=
          (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz[(ixp + 1) * (nynz_Jz) + (iyp + 1) * (nz_Jz) + (izd + 1)] +=
          (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
      }
    }); // end for each particles
  }
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

  mini_float *Ex = em.Ex_m.get_raw_pointer(minipic::device);
  mini_float *Ey = em.Ey_m.get_raw_pointer(minipic::device);
  mini_float *Ez = em.Ez_m.get_raw_pointer(minipic::device);

  mini_float *Bx = em.Bx_m.get_raw_pointer(minipic::device);
  mini_float *By = em.By_m.get_raw_pointer(minipic::device);
  mini_float *Bz = em.Bz_m.get_raw_pointer(minipic::device);

  mini_float *Jx = em.Jx_m.get_raw_pointer(minipic::device);
  mini_float *Jy = em.Jy_m.get_raw_pointer(minipic::device);
  mini_float *Jz = em.Jz_m.get_raw_pointer(minipic::device);

  const auto nx_Ex = em.Ex_m.nx(), ny_Ex = em.Ex_m.ny(), nz_Ex = em.Ex_m.nz(),
             nynz_Ex = ny_Ex * nz_Ex;
  const auto nx_Ey = em.Ey_m.nx(), ny_Ey = em.Ey_m.ny(), nz_Ey = em.Ey_m.nz(),
             nynz_Ey = ny_Ey * nz_Ey;
  const auto nx_Ez = em.Ez_m.nx(), ny_Ez = em.Ez_m.ny(), nz_Ez = em.Ez_m.nz(),
             nynz_Ez = ny_Ez * nz_Ez;

  const auto nx_Bx = em.Bx_m.nx(), ny_Bx = em.Bx_m.ny(), nz_Bx = em.Bx_m.nz(),
             nynz_Bx = ny_Bx * nz_Bx;
  const auto nx_By = em.By_m.nx(), ny_By = em.By_m.ny(), nz_By = em.By_m.nz(),
             nynz_By = ny_By * nz_By;
  const auto nx_Bz = em.Bz_m.nx(), ny_Bz = em.Bz_m.ny(), nz_Bz = em.Bz_m.nz(),
             nynz_Bz = ny_Bz * nz_Bz;

  const auto nx_Jx = em.Jx_m.nx(), ny_Jx = em.Jx_m.ny(), nz_Jx = em.Jx_m.nz(),
             nynz_Jx = ny_Jx * nz_Jx;
  const auto nx_Jy = em.Jy_m.nx(), ny_Jy = em.Jy_m.ny(), nz_Jy = em.Jy_m.nz(),
             nynz_Jy = ny_Jy * nz_Jy;
  const auto nx_Jz = em.Jz_m.nx(), ny_Jz = em.Jz_m.ny(), nz_Jz = em.Jz_m.nz(),
             nynz_Jz = ny_Jz * nz_Jz;

  // Electric field Ex (d,p,p)
  std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Ex * nynz_Ex), [=](int idx) {
    const int ix = idx / (nynz_Ex);
    const int iy = (idx - ix * nynz_Ex) / nz_Ex;
    const int iz = idx - ix * nynz_Ex - iy * nz_Ex;

    Ex[idx] +=
      -dt * Jx[ix * (nynz_Jx) + (iy + 1) * nz_Jx + iz + 1] +
      dt_over_dy * (Bz[ix * nynz_Bz + (iy + 1) * nz_Bz + iz] - Bz[ix * nynz_Bz + iy * nz_Bz + iz]) -
      dt_over_dz * (By[ix * nynz_By + iy * nz_By + iz + 1] - By[ix * nynz_By + iy * nz_By + iz]);
  });

  // Electric field Ey (p,d,p)
  std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Ey * nynz_Ey), [=](int idx) {
    const int ix = idx / (nynz_Ey);
    const int iy = (idx - ix * nynz_Ey) / nz_Ey;
    const int iz = idx - ix * nynz_Ey - iy * nz_Ey;

    Ey[idx] +=
      -dt * Jy[(ix + 1) * nynz_Jy + iy * nz_Jy + iz + 1] -
      dt_over_dx * (Bz[(ix + 1) * nynz_Bz + iy * nz_Bz + iz] - Bz[ix * nynz_Bz + iy * nz_Bz + iz]) +
      dt_over_dz * (Bx[ix * nynz_Bx + iy * nz_Bx + iz + 1] - Bx[ix * nynz_Bx + iy * nz_Bx + iz]);
  });

  // Electric field Ez (p,p,d)
  std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Ez * nynz_Ez), [=](int idx) {
    const int ix = idx / (nynz_Ez);
    const int iy = (idx - ix * nynz_Ez) / nz_Ez;
    const int iz = idx - ix * nynz_Ez - iy * nz_Ez;

    Ez[idx] +=
      -dt * Jz[(ix + 1) * nynz_Jz + (iy + 1) * nz_Jz + iz] +
      dt_over_dx * (By[(ix + 1) * nynz_By + iy * nz_By + iz] - By[ix * nynz_By + iy * nz_By + iz]) -
      dt_over_dy * (Bx[ix * nynz_Bx + (iy + 1) * nz_Bx + iz] - Bx[ix * nynz_Bx + iy * nz_Bx + iz]);
  });

  /////     Solve Maxwell Faraday (B)

  // Magnetic field Bx (p,d,d)
  int nz   = nz_Bx - 2;
  int nynz = (ny_Bx - 2) * nz;
  std::for_each(policy,
                counting_iterator<int>(0),
                counting_iterator<int>(nx_Bx * (ny_Bx - 2) * (nz_Bx - 2)),
                [=](int idx) {
                  const int ix = idx / (nynz);
                  int iy       = (idx - ix * nynz) / nz;
                  int iz       = idx - ix * nynz - iy * nz;

                  iy += 1;
                  iz += 1;

                  Bx[ix * nynz_Bx + iy * nz_Bx + iz] +=
                    -dt_over_dy * (Ez[ix * nynz_Ez + iy * nz_Ez + iz] -
                                   Ez[ix * nynz_Ez + (iy - 1) * nz_Ez + iz]) +
                    dt_over_dz *
                      (Ey[ix * nynz_Ey + iy * nz_Ey + iz] - Ey[ix * nynz_Ey + iy * nz_Ey + iz - 1]);
                });

  // Magnetic field By (d,p,d)
  nz   = nz_By - 2;
  nynz = ny_By * nz;
  std::for_each(policy,
                counting_iterator<int>(0),
                counting_iterator<int>((nx_By - 2) * ny_By * (nz_By - 2)),
                [=](int idx) {
                  int ix       = idx / nynz;
                  const int iy = (idx - ix * nynz) / nz;
                  int iz       = idx - ix * nynz - iy * nz;

                  ix += 1;
                  iz += 1;

                  By[ix * nynz_By + iy * nz_By + iz] +=
                    -dt_over_dz * (Ex[ix * nynz_Ex + iy * nz_Ex + iz] -
                                   Ex[(ix)*nynz_Ex + iy * nz_Ex + iz - 1]) +
                    dt_over_dx * (Ez[ix * nynz_Ez + iy * nz_Ez + iz] -
                                  Ez[(ix - 1) * nynz_Ez + iy * nz_Ez + iz]);
                });

  // Magnetic field Bz (d,d,p)
  nz   = nz_Bz;
  nynz = (ny_Bz - 2) * nz;
  std::for_each(policy,
                counting_iterator<int>(0),
                counting_iterator<int>((nx_Bz - 2) * (ny_Bz - 2) * nz_Bz),
                [=](int idx) {
                  int ix       = idx / nynz;
                  int iy       = (idx - ix * nynz) / nz;
                  const int iz = idx - ix * nynz - iy * nz;

                  ix += 1;
                  iy += 1;

                  Bz[ix * nynz_Bz + iy * nz_Bz + iz] +=
                    -dt_over_dx * (Ey[ix * nynz_Ey + iy * nz_Ey + iz] -
                                   Ey[(ix - 1) * nynz_Ey + iy * nz_Ey + iz]) +
                    dt_over_dy * (Ex[ix * nynz_Ex + iy * nz_Ex + iz] -
                                  Ex[ix * nynz_Ex + (iy - 1) * nz_Ex + iz]);
                });

} // end solve

// _______________________________________________________________
//
//! \brief Boundaries condition on the global grid
//! \param[in] Params & params - global constant parameters
//! \param[in] ElectroMagn & em - global electromagnetic fields
// _______________________________________________________________
auto currentBC(Params &params, ElectroMagn &em) -> void {

  if (params.boundary_condition == "periodic") {

    mini_float *Jx = em.Jx_m.get_raw_pointer(minipic::device);
    mini_float *Jy = em.Jy_m.get_raw_pointer(minipic::device);
    mini_float *Jz = em.Jz_m.get_raw_pointer(minipic::device);

    const auto nx_Jx = em.Jx_m.nx(), ny_Jx = em.Jx_m.ny(), nz_Jx = em.Jx_m.nz(),
               nynz_Jx = ny_Jx * nz_Jx;
    const auto nx_Jy = em.Jy_m.nx(), ny_Jy = em.Jy_m.ny(), nz_Jy = em.Jy_m.nz(),
               nynz_Jy = ny_Jy * nz_Jy;
    const auto nx_Jz = em.Jz_m.nx(), ny_Jz = em.Jz_m.ny(), nz_Jz = em.Jz_m.nz(),
               nynz_Jz = ny_Jz * nz_Jz;

    // X
    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(ny_Jx * nz_Jx), [=](int idx) {
      // Method 1
      /*{
        const auto index_left  = idx;
        const auto index_right = (nx_Jx - 2) * ny_Jx * nz_Jx + idx;

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }

      {
        const auto index_left  = idx + ny_Jx * nz_Jx;
        const auto index_right = (nx_Jx - 1) * ny_Jx * nz_Jx + idx;

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }*/

  // Method 2
      const int iy = idx / nz_Jx;
      const int iz = idx - iy * nz_Jx;
      {
        auto index_left  = iy * nz_Jx + iz;
        auto index_right = (nx_Jx - 2) * nynz_Jx + iy * nz_Jx + iz;

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }
      {
        auto index_left  = nynz_Jx + iy * nz_Jx + iz;
        auto index_right = (nx_Jx - 1) * nynz_Jx + iy * nz_Jx + iz;

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }
    });

    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(ny_Jy * nz_Jy), [=](int idx) {
      // Method 1
      /*{
        const auto index_left  = idx;
        const auto index_right = (nx_Jy - 2) * ny_Jy * nz_Jy + idx;

        Jy[index_left] += Jy[index_right];
        Jy[index_right] = Jy[index_left];
      }

      {
        const auto index_left  = idx + ny_Jy * nz_Jy;
        const auto index_right = (nx_Jy - 1) * ny_Jy * nz_Jy + idx;

        Jy[index_left] += Jy[index_right];
        Jy[index_right] = Jy[index_left];
      }*/

  // Method 2
      const int iy = idx / nz_Jy;
      const int iz = idx - iy * nz_Jy;
      {
        auto index_left  = iy * nz_Jy + iz;
        auto index_right = (nx_Jy - 2) * nynz_Jy + iy * nz_Jy + iz;

        Jy[index_left] += Jy[index_right];
        Jy[index_right] = Jy[index_left];
      }
      {
        auto index_left  = nynz_Jy + iy * nz_Jy + iz;
        auto index_right = (nx_Jy - 1) * nynz_Jy + iy * nz_Jy + iz;

        Jy[index_left] += Jy[index_right];
        Jy[index_right] = Jy[index_left];
      }
    });

    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(ny_Jz * nz_Jz), [=](int idx) {
      // Method 1
      /* {
         const auto index_left  = idx;
         const auto index_right = (nx_Jz - 2) * ny_Jz * nz_Jz + idx;

         Jz[index_left] += Jz[index_right];
         Jz[index_right] = Jz[index_left];
       }

       {
         const auto index_left  = idx + ny_Jz * nz_Jz;
         const auto index_right = (nx_Jz - 1) * ny_Jz * nz_Jz + idx;

         Jz[index_left] += Jz[index_right];
         Jz[index_right] = Jz[index_left];
       }*/

  // Method 2
      const int iy = idx / nz_Jz;
      const int iz = idx - iy * nz_Jz;
      {
        auto index_left  = iy * nz_Jz + iz;
        auto index_right = (nx_Jz - 2) * nynz_Jz + iy * nz_Jz + iz;

        Jz[index_left] += Jz[index_right];
        Jz[index_right] = Jz[index_left];
      }
      {
        auto index_left  = nynz_Jz + iy * nz_Jz + iz;
        auto index_right = (nx_Jz - 1) * nynz_Jz + iy * nz_Jz + iz;

        Jz[index_left] += Jz[index_right];
        Jz[index_right] = Jz[index_left];
      }
    });

    // Y
    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Jx * nz_Jx), [=](int idx) {
      const int ix = idx / nz_Jx;
      const int iz = idx - ix * nz_Jx;
      {
        auto index_left  = ix * nynz_Jx + iz; // em.Jx_m(ix, 0, iz) += em.Jx_m(ix, ny_Jx - 2, iz);
        auto index_right = ix * nynz_Jx + (ny_Jx - 2) * nz_Jx + iz;

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }
      {
        auto index_left  = ix * nynz_Jx + nz_Jx + iz;
        auto index_right = ix * nynz_Jx + (ny_Jx - 1) * nz_Jx + iz;

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }
    });

    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Jy * nz_Jy), [=](int idx) {
      const int ix = idx / nz_Jy;
      const int iz = idx - ix * nz_Jy;
      {
        auto index_left  = ix * nynz_Jy + iz;
        auto index_right = ix * nynz_Jy + (ny_Jy - 2) * nz_Jy + iz;

        Jy[index_left] += Jy[index_right];
        Jy[index_right] = Jy[index_left];
      }
      {
        auto index_left  = ix * nynz_Jy + nz_Jy + iz;
        auto index_right = ix * nynz_Jy + (ny_Jy - 1) * nz_Jy + iz;

        Jy[index_left] += Jy[index_right];
        Jy[index_right] = Jy[index_left];
      }
    });

    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Jz * nz_Jz), [=](int idx) {
      const int ix = idx / nz_Jz;
      const int iz = idx - ix * nz_Jz;
      {
        auto index_left  = ix * nynz_Jz + iz;
        auto index_right = ix * nynz_Jz + (ny_Jz - 2) * nz_Jz + iz;

        Jz[index_left] += Jz[index_right];
        Jz[index_right] = Jz[index_left];
      }
      {
        auto index_left  = ix * nynz_Jz + nz_Jz + iz;
        auto index_right = ix * nynz_Jz + (ny_Jz - 1) * nz_Jz + iz;

        Jz[index_left] += Jz[index_right];
        Jz[index_right] = Jz[index_left];
      }
    });

    // Z
    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Jx * ny_Jx), [=](int idx) {
      const int ix = idx / ny_Jx;
      const int iy = idx - ix * ny_Jx;
      {
        auto index_left  = ix * nynz_Jx + iy * nz_Jx;
        auto index_right = ix * nynz_Jx + iy * nz_Jx + (nz_Jx - 2);

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }
      {
        auto index_left  = ix * nynz_Jx + iy * nz_Jx + 1;
        auto index_right = ix * nynz_Jx + iy * nz_Jx + (nz_Jx - 1);

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }
    });

    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Jy * ny_Jy), [=](int idx) {
      const int ix = idx / ny_Jy;
      const int iy = idx - ix * ny_Jy;
      {
        auto index_left  = ix * nynz_Jy + iy * nz_Jy;
        auto index_right = ix * nynz_Jy + iy * nz_Jy + (nz_Jy - 2);

        Jy[index_left] += Jy[index_right];
        Jy[index_right] = Jy[index_left];
      }
      {
        auto index_left  = ix * nynz_Jy + iy * nz_Jy + 1;
        auto index_right = ix * nynz_Jy + iy * nz_Jy + (nz_Jy - 1);

        Jy[index_left] += Jy[index_right];
        Jy[index_right] = Jy[index_left];
      }
    });

    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Jz * ny_Jz), [=](int idx) {
      const int ix = idx / ny_Jz;
      const int iy = idx - ix * ny_Jz;
      {
        auto index_left  = ix * nynz_Jz + iy * nz_Jz;
        auto index_right = ix * nynz_Jz + iy * nz_Jz + (nz_Jz - 2);

        Jz[index_left] += Jz[index_right];
        Jz[index_right] = Jz[index_left];
      }
      {
        auto index_left  = ix * nynz_Jz + iy * nz_Jz + 1;
        auto index_right = ix * nynz_Jz + iy * nz_Jz + (nz_Jz - 1);

        Jz[index_left] += Jz[index_right];
        Jz[index_right] = Jz[index_left];
      }
    });
  } // end if periodic
} // end currentBC

// _______________________________________________________________
//
//! \brief Boundaries condition on the global grid
//! \param[in] Params & params - global constant parameters
// _______________________________________________________________
auto solveBC(Params &params, ElectroMagn &em) -> void {

  mini_float *Bx = em.Bx_m.get_raw_pointer(minipic::device);
  mini_float *By = em.By_m.get_raw_pointer(minipic::device);
  mini_float *Bz = em.Bz_m.get_raw_pointer(minipic::device);

  const auto nx_Bx = em.Bx_m.nx(), ny_Bx = em.Bx_m.ny(), nz_Bx = em.Bx_m.nz(),
             nynz_Bx = ny_Bx * nz_Bx;
  const auto nx_By = em.By_m.nx(), ny_By = em.By_m.ny(), nz_By = em.By_m.nz(),
             nynz_By = ny_By * nz_By;
  const auto nx_Bz = em.Bz_m.nx(), ny_Bz = em.Bz_m.ny(), nz_Bz = em.Bz_m.nz(),
             nynz_Bz = ny_Bz * nz_Bz;

  if (params.boundary_condition == "periodic") {

    // X dim
    // By (d,p,d)
    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nynz_By), [=](int idx) {
      By[idx]                         = By[idx + (nx_By - 2) * nz_By * ny_By];
      By[idx + (nx_By - 1) * nynz_By] = By[idx + nynz_By];
    });

    // Bz (d,d,p)
    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nynz_Bz), [=](int idx) {
      Bz[idx]                         = Bz[idx + (nx_Bz - 2) * nz_Bz * ny_Bz];
      Bz[idx + (nx_Bz - 1) * nynz_Bz] = Bz[idx + nynz_Bz];
    });

    // Y dim
    // Bx (p,d,d)
    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Bx * nz_Bx), [=](int idx) {
      const int ix = idx / nz_Bx;
      const int iz = idx - ix * nz_Bx;

      Bx[ix * nynz_Bx + iz]                       = Bx[ix * nynz_Bx + (ny_Bx - 2) * nz_Bx + iz];
      Bx[ix * nynz_Bx + (ny_Bx - 1) * nz_Bx + iz] = Bx[ix * nynz_Bx + nz_Bx + iz];
    });

    // Bz (d,d,p)
    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Bz * nz_Bz), [=](int idx) {
      const int ix = idx / nz_Bz;
      const int iz = idx - ix * nz_Bz;

      Bz[ix * nynz_Bz + iz]                       = Bz[ix * nynz_Bz + (ny_Bz - 2) * nz_Bz + iz];
      Bz[ix * nynz_Bz + (ny_Bz - 1) * nz_Bz + iz] = Bz[ix * nynz_Bz + nz_Bz + iz];
    });

    // Z dim
    // Bx
    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_Bx * ny_Bx), [=](int idx) {
      const int ix = idx / ny_Bx;
      const int iy = idx - ix * ny_Bx;

      Bx[ix * nynz_Bx + iy * nz_Bx]               = Bx[ix * nynz_Bx + iy * nz_Bx + (nz_Bx - 2)];
      Bx[ix * nynz_Bx + iy * nz_Bx + (nz_Bx - 1)] = Bx[ix * nynz_Bx + iy * nz_Bx + 1];
    });

    // By
    std::for_each(policy, counting_iterator<int>(0), counting_iterator<int>(nx_By * ny_By), [=](int idx) {
      const int ix = idx / ny_By;
      const int iy = idx - ix * ny_By;

      By[ix * nynz_By + iy * nz_By]               = By[ix * nynz_By + iy * nz_By + (nz_By - 2)];
      By[ix * nynz_By + iy * nz_By + (nz_By - 1)] = By[ix * nynz_By + iy * nz_By + 1];
    });

  } else if (params.boundary_condition == "reflective") {
    // NOT TESTED car antenna, default et beam sont périodiques
    //  X dim
    //  By (d,p,d)
    for (unsigned int iy = 0; iy < ny_By; ++iy) {
      for (unsigned int iz = 0; iz < nz_By; ++iz) {
        // -X
        em.By_m(0, iy, iz) = em.By_m(1, iy, iz);
        // +X
        em.By_m(nx_By - 1, iy, iz) = em.By_m(nx_By - 2, iy, iz);
      }
    }

    // Bz (d,d,p)
    for (unsigned int iy = 0; iy < ny_Bz; iy++) {
      for (unsigned int iz = 0; iz < nz_Bz; iz++) {
        // -X
        em.Bz_m(0, iy, iz) = em.Bz_m(1, iy, iz);
        // +X
        em.Bz_m(nx_Bz - 1, iy, iz) = em.Bz_m(nx_Bz - 2, iy, iz);
      }
    }

    // Y dim
    // Bx (p,d,d)
    for (unsigned int ix = 0; ix < nx_Bx; ix++) {
      for (unsigned int iz = 0; iz < nz_Bx; iz++) {
        // -Y
        em.Bx_m(ix, 0, iz) = em.Bx_m(ix, 1, iz);
        // +Y
        em.Bx_m(ix, ny_Bx - 1, iz) = em.Bx_m(ix, ny_Bx - 2, iz);
      }
    }
    // Bz (-1 to avoid corner)
    for (unsigned int ix = 0; ix < nx_Bz; ix++) {
      for (unsigned int iz = 0; iz < nz_Bz; ++iz) {
        // -Y
        em.Bz_m(ix, 0, iz) = em.Bz_m(ix, 1, iz);
        // +Y
        em.Bz_m(ix, ny_Bz - 1, iz) = em.Bz_m(ix, ny_Bz - 2, iz);
      }
    }

    // Z dim
    // Bx
    for (unsigned int ix = 0; ix < nx_Bx; ix++) {
      for (unsigned int iy = 0; iy < ny_Bx; iy++) {
        // -Z
        em.Bx_m(ix, iy, 0) = em.Bx_m(ix, iy, 1);
        // +Z
        em.Bx_m(ix, iy, nz_Bx - 1) = em.Bx_m(ix, iy, nz_Bx - 2);
      }
    }
    // By
    for (unsigned int ix = 0; ix < nx_By; ix++) {
      for (unsigned int iy = 0; iy < ny_By; iy++) {
        // -Z
        em.By_m(ix, iy, 0) = em.By_m(ix, iy, 1);
        // +Z
        em.By_m(ix, iy, nz_By - 1) = em.By_m(ix, iy, nz_By - 2);
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

          const size_t buffer_size =
            vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].size();

          if (buffer_size > 0) {
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
  mini_float *Jx_0 = patch.vec_Jx_m[0].get_raw_pointer(minipic::device);
  mini_float *Jy_0 = patch.vec_Jy_m[0].get_raw_pointer(minipic::device);
  mini_float *Jz_0 = patch.vec_Jz_m[0].get_raw_pointer(minipic::device);

  for (int is = 1; is < patch.n_species_m; is++) {
    mini_float *Jx_is = patch.vec_Jx_m[is].get_raw_pointer(minipic::device);
    mini_float *Jy_is = patch.vec_Jy_m[is].get_raw_pointer(minipic::device);
    mini_float *Jz_is = patch.vec_Jz_m[is].get_raw_pointer(minipic::device);

    const auto nx_Jx = patch.vec_Jx_m[is].nx(), nx_Jy = patch.vec_Jy_m[is].nx(),
               nx_Jz = patch.vec_Jz_m[is].nx();
    const auto ny_Jx = patch.vec_Jx_m[is].ny(), ny_Jy = patch.vec_Jy_m[is].ny(),
               ny_Jz = patch.vec_Jz_m[is].ny();
    const auto nz_Jx = patch.vec_Jx_m[is].nz(), nz_Jy = patch.vec_Jy_m[is].nz(),
               nz_Jz = patch.vec_Jz_m[is].nz();

    // Only if particles projected
    if (patch.projected_[is]) {
      std::for_each(policy,
                    counting_iterator<int>(0),
                    counting_iterator<int>(nx_Jx * ny_Jx * nz_Jx),
                    [=](int idx) { Jx_0[idx] += Jx_is[idx]; });

      std::for_each(policy,
                    counting_iterator(0),
                    counting_iterator(nx_Jy * ny_Jy * nz_Jy),
                    [=](int idx) { Jy_0[idx] += Jy_is[idx]; });

      std::for_each(policy,
                    counting_iterator(0),
                    counting_iterator(nx_Jz * ny_Jz * nz_Jz),
                    [=](int idx) { Jz_0[idx] += Jz_is[idx]; });
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

    // std::cerr << i_global_p  << " " << j_global_p << " " << k_global_p << std::endl;
    // std::cerr << "nx: " << vec_Jx_m[0].nx()  << " " << em.Jx_m.nx() << " " << std::endl;

    const auto nx_Jx_loc = patch.vec_Jx_m[0].nx(), ny_Jx_loc = patch.vec_Jx_m[0].ny(),
               nz_Jx_loc = patch.vec_Jx_m[0].nz(), nynz_Jx_loc = ny_Jx_loc * nz_Jx_loc;
    const auto nx_Jy_loc = patch.vec_Jy_m[0].nx(), ny_Jy_loc = patch.vec_Jy_m[0].ny(),
               nz_Jy_loc = patch.vec_Jy_m[0].nz(), nynz_Jy_loc = ny_Jy_loc * nz_Jy_loc;
    const auto nx_Jz_loc = patch.vec_Jz_m[0].nx(), ny_Jz_loc = patch.vec_Jz_m[0].ny(),
               nz_Jz_loc = patch.vec_Jz_m[0].nz(), nynz_Jz_loc = ny_Jz_loc * nz_Jz_loc;

    const auto nx_Jx_glob = em.Jx_m.nx(), ny_Jx_glob = em.Jx_m.ny(), nz_Jx_glob = em.Jx_m.nz(),
               nynz_Jx_glob = ny_Jx_glob * nz_Jx_glob;
    const auto nx_Jy_glob = em.Jy_m.nx(), ny_Jy_glob = em.Jy_m.ny(), nz_Jy_glob = em.Jy_m.nz(),
               nynz_Jy_glob = ny_Jy_glob * nz_Jy_glob;
    const auto nx_Jz_glob = em.Jz_m.nx(), ny_Jz_glob = em.Jz_m.ny(), nz_Jz_glob = em.Jz_m.nz(),
               nynz_Jz_glob = ny_Jz_glob * nz_Jz_glob;

    mini_float *Jx   = em.Jx_m.get_raw_pointer(minipic::device);
    mini_float *Jy   = em.Jy_m.get_raw_pointer(minipic::device);
    mini_float *Jz   = em.Jz_m.get_raw_pointer(minipic::device);
    mini_float *Jx_0 = patch.vec_Jx_m[0].get_raw_pointer(minipic::device);
    mini_float *Jy_0 = patch.vec_Jy_m[0].get_raw_pointer(minipic::device);
    mini_float *Jz_0 = patch.vec_Jz_m[0].get_raw_pointer(minipic::device);

    std::for_each(
      policy,
      counting_iterator<int>(0),
      counting_iterator<int>(nx_Jx_loc * ny_Jx_loc * nz_Jx_loc),
      [=](int idx) {
        const int ix = idx / (nynz_Jx_loc);
        const int iy = (idx - ix * nynz_Jx_loc) / nz_Jx_loc;
        const int iz = idx - ix * nynz_Jx_loc - iy * nz_Jx_loc;

        Jx[(i_global_d + ix) * nynz_Jx_glob + (j_global_p + iy) * nz_Jx_glob + k_global_p + iz] +=
          Jx_0[idx];
      });

    std::for_each(
      policy,
      counting_iterator<int>(0),
      counting_iterator<int>(nx_Jy_loc * ny_Jy_loc * nz_Jy_loc),
      [=](int idx) {
        const int ix = idx / (nynz_Jy_loc);
        const int iy = (idx - ix * nynz_Jy_loc) / nz_Jy_loc;
        const int iz = idx - ix * nynz_Jy_loc - iy * nz_Jy_loc;

        Jy[(i_global_p + ix) * nynz_Jy_glob + (j_global_d + iy) * nz_Jy_glob + k_global_p + iz] +=
          Jy_0[idx];
      });

    std::for_each(
      policy,
      counting_iterator<int>(0),
      counting_iterator<int>(nx_Jz_loc * ny_Jz_loc * nz_Jz_loc),
      [=](int idx) {
        const int ix = idx / (nynz_Jz_loc);
        const int iy = (idx - ix * nynz_Jz_loc) / nz_Jz_loc;
        const int iz = idx - ix * nynz_Jz_loc - iy * nz_Jz_loc;

        Jz[(i_global_p + ix) * nynz_Jz_glob + (j_global_p + iy) * nz_Jz_glob + k_global_d + iz] +=
          Jz_0[idx];
      });
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
}

} // end namespace operators

#endif // OPERATORS_H
