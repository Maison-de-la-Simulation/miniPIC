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

  for (int is = 0; is < patch.n_species_m; is++) {

    const int n_particles = patch.particles_m[is].size();

#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
    for (unsigned int part = 0; part < n_particles; part++) {

      // // Calculate normalized positions
      const auto ixn = patch.particles_m[is].x_(part) * inv_dx_m;
      const auto iyn = patch.particles_m[is].y_(part) * inv_dy_m;
      const auto izn = patch.particles_m[is].z_(part) * inv_dz_m;

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

        const auto v00 =
          em.Ex_m(ixd, iyp, izp) * (1 - coeffs[0]) + em.Ex_m(ixd + 1, iyp, izp) * coeffs[0];
        const auto v01 =
          em.Ex_m(ixd, iyp, izp + 1) * (1 - coeffs[0]) + em.Ex_m(ixd + 1, iyp, izp + 1) * coeffs[0];
        const auto v10 =
          em.Ex_m(ixd, iyp + 1, izp) * (1 - coeffs[0]) + em.Ex_m(ixd + 1, iyp + 1, izp) * coeffs[0];
        const auto v11 = em.Ex_m(ixd, iyp + 1, izp + 1) * (1 - coeffs[0]) +
                         em.Ex_m(ixd + 1, iyp + 1, izp + 1) * coeffs[0];
        const auto v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const auto v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        patch.particles_m[is].Ex_(part) = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // Ey (p, d, p)
      {
        const mini_float coeffs[3] = {ixn, iyn + 0.5, izn};

        const mini_float v00 =
          em.Ey_m(ixp, iyd, izp) * (1 - coeffs[0]) + em.Ey_m(ixp + 1, iyd, izp) * coeffs[0];
        const mini_float v01 =
          em.Ey_m(ixp, iyd, izp + 1) * (1 - coeffs[0]) + em.Ey_m(ixp + 1, iyd, izp + 1) * coeffs[0];
        const mini_float v10 =
          em.Ey_m(ixp, iyd + 1, izp) * (1 - coeffs[0]) + em.Ey_m(ixp + 1, iyd + 1, izp) * coeffs[0];
        const mini_float v11 = em.Ey_m(ixp, iyd + 1, izp + 1) * (1 - coeffs[0]) +
                               em.Ey_m(ixp + 1, iyd + 1, izp + 1) * coeffs[0];
        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        patch.particles_m[is].Ey_(part) = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }

      // // //particles_m[is].Ey_.d_view(part) = compute_interpolation(ixp, b, izp, coeffs, Ey);
      // Ez (p, p, d)
      {
        const mini_float coeffs[3] = {ixn, iyn, izn + 0.5};

        const mini_float v00 =
          em.Ez_m(ixp, iyp, izd) * (1 - coeffs[0]) + em.Ez_m(ixp + 1, iyp, izd) * coeffs[0];
        const mini_float v01 =
          em.Ez_m(ixp, iyp, izd + 1) * (1 - coeffs[0]) + em.Ez_m(ixp + 1, iyp, izd + 1) * coeffs[0];
        const mini_float v10 =
          em.Ez_m(ixp, iyp + 1, izd) * (1 - coeffs[0]) + em.Ez_m(ixp + 1, iyp + 1, izd) * coeffs[0];
        const mini_float v11 = em.Ez_m(ixp, iyp + 1, izd + 1) * (1 - coeffs[0]) +
                               em.Ez_m(ixp + 1, iyp + 1, izd + 1) * coeffs[0];
        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        patch.particles_m[is].Ez_(part) = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }
      // particles_m[is].Ez_.d_view(part) = compute_interpolation(ixp, iyp, g, coeffs, Ez);

      // interpolation magnetic field
      // Bx (p, d, d)
      {
        const mini_float coeffs[3] = {ixn, iyn + 0.5, izn + 0.5};

        const mini_float v00 =
          em.Bx_m(ixp, iyd, izd) * (1 - coeffs[0]) + em.Bx_m(ixp + 1, iyd, izd) * coeffs[0];
        const mini_float v01 =
          em.Bx_m(ixp, iyd, izd + 1) * (1 - coeffs[0]) + em.Bx_m(ixp + 1, iyd, izd + 1) * coeffs[0];
        const mini_float v10 =
          em.Bx_m(ixp, iyd + 1, izd) * (1 - coeffs[0]) + em.Bx_m(ixp + 1, iyd + 1, izd) * coeffs[0];
        const mini_float v11 = em.Bx_m(ixp, iyd + 1, izd + 1) * (1 - coeffs[0]) +
                               em.Bx_m(ixp + 1, iyd + 1, izd + 1) * coeffs[0];
        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        patch.particles_m[is].Bx_(part) = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }
      // particles_m[is].Bx_.d_view(part) = compute_interpolation(ixp, b, g, coeffs, Bx);

      // By (d, p, d)
      {
        const mini_float coeffs[3] = {ixn + 0.5, iyn, izn + 0.5};

        const mini_float v00 =
          em.By_m(ixd, iyp, izd) * (1 - coeffs[0]) + em.By_m(ixd + 1, iyp, izd) * coeffs[0];
        const mini_float v01 =
          em.By_m(ixd, iyp, izd + 1) * (1 - coeffs[0]) + em.By_m(ixd + 1, iyp, izd + 1) * coeffs[0];
        const mini_float v10 =
          em.By_m(ixd, iyp + 1, izd) * (1 - coeffs[0]) + em.By_m(ixd + 1, iyp + 1, izd) * coeffs[0];
        const mini_float v11 = em.By_m(ixd, iyp + 1, izd + 1) * (1 - coeffs[0]) +
                               em.By_m(ixd + 1, iyp + 1, izd + 1) * coeffs[0];
        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        patch.particles_m[is].By_(part) = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }
      // particles_m[is].By_.d_view(part) = compute_interpolation(a, iyp, g, coeffs, By);

      // Bz (d, d, p)
      {
        const mini_float coeffs[3] = {ixn + 0.5, iyn + 0.5, izn};

        // Bzp(part)              = compute_interpolation(coeffs,
        //                                   Bz(ixd, iyd, izp),
        //                                   Bz(ixd, iyd, izp + 1),
        //                                   Bz(ixd, iyd + 1, izp),
        //                                   Bz(ixd, iyd + 1, izp + 1),
        //                                   Bz(ixd + 1, iyd, izp),
        //                                   Bz(ixd + 1, iyd, izp + 1),
        //                                   Bz(ixd + 1, iyd + 1, izp),
        //                                   Bz(ixd + 1, iyd + 1, izp + 1));

        const mini_float v00 =
          em.Bz_m(ixd, iyd, izp) * (1 - coeffs[0]) + em.Bz_m(ixd + 1, iyd, izp) * coeffs[0];
        const mini_float v01 =
          em.Bz_m(ixd, iyd, izp + 1) * (1 - coeffs[0]) + em.Bz_m(ixd + 1, iyd, izp + 1) * coeffs[0];
        const mini_float v10 =
          em.Bz_m(ixd, iyd + 1, izp) * (1 - coeffs[0]) + em.Bz_m(ixd + 1, iyd + 1, izp) * coeffs[0];
        const mini_float v11 = em.Bz_m(ixd, iyd + 1, izp + 1) * (1 - coeffs[0]) +
                               em.Bz_m(ixd + 1, iyd + 1, izp + 1) * coeffs[0];
        const mini_float v0 = v00 * (1 - coeffs[1]) + v10 * coeffs[1];
        const mini_float v1 = v01 * (1 - coeffs[1]) + v11 * coeffs[1];

        patch.particles_m[is].Bz_(part) = v0 * (1 - coeffs[2]) + v1 * coeffs[2];
      }
      // particles_m[is].Bz_.d_view(part) = compute_interpolation(a, b, izp, coeffs, Bz);
    } // End for each particle

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

    const int n_particles = patch.particles_m[is].size();

    // q' = dt * (q/2m)
    const mini_float qp = patch.particles_m[is].charge_m * dt * 0.5 / patch.particles_m[is].mass_m;

#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
    for (int ip = 0; ip < n_particles; ++ip) {

      // 1/2 E
      mini_float px = qp * patch.particles_m[is].Ex_(ip);
      mini_float py = qp * patch.particles_m[is].Ey_(ip);
      mini_float pz = qp * patch.particles_m[is].Ez_(ip);

      const mini_float ux = patch.particles_m[is].mx_(ip) + px;
      const mini_float uy = patch.particles_m[is].my_(ip) + py;
      const mini_float uz = patch.particles_m[is].mz_(ip) + pz;

      // gamma-factor
      mini_float gamma_inv = qp / sqrt(1 + (ux * ux + uy * uy + uz * uz));

      // B, T = Transform to rotate the particle
      const mini_float tx  = gamma_inv * patch.particles_m[is].Bx_(ip);
      const mini_float ty  = gamma_inv * patch.particles_m[is].By_(ip);
      const mini_float tz  = gamma_inv * patch.particles_m[is].Bz_(ip);
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
      patch.particles_m[is].mx_(ip) = px;
      patch.particles_m[is].my_(ip) = py;
      patch.particles_m[is].mz_(ip) = pz;

      // Update positions
      patch.particles_m[is].x_(ip) += px * dt * gamma_inv;
      patch.particles_m[is].y_(ip) += py * dt * gamma_inv;
      patch.particles_m[is].z_(ip) += pz * dt * gamma_inv;
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

    const int n_particles = patch.particles_m[is].size();

    // q' = dt * (q/2m)
    const mini_float qp = patch.particles_m[is].charge_m * dt * 0.5 / patch.particles_m[is].mass_m;

#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
    for (auto ip = 0; ip < n_particles; ++ip) {

      // 1/2 E
      mini_float px = qp * patch.particles_m[is].Ex_h(ip);
      mini_float py = qp * patch.particles_m[is].Ey_h(ip);
      mini_float pz = qp * patch.particles_m[is].Ez_h(ip);

      const mini_float ux = patch.particles_m[is].mx_h(ip) + px;
      const mini_float uy = patch.particles_m[is].my_h(ip) + py;
      const mini_float uz = patch.particles_m[is].mz_h(ip) + pz;

      // gamma-factor
      mini_float gamma_inv = qp / sqrt(1 + (ux * ux + uy * uy + uz * uz));

      // B, T = Transform to rotate the particle
      const mini_float tx  = gamma_inv * patch.particles_m[is].Bx_h(ip);
      const mini_float ty  = gamma_inv * patch.particles_m[is].By_h(ip);
      const mini_float tz  = gamma_inv * patch.particles_m[is].Bz_h(ip);
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
      patch.particles_m[is].mx_h(ip) = px;
      patch.particles_m[is].my_h(ip) = py;
      patch.particles_m[is].mz_h(ip) = pz;

    } // End for each particles
  }   // end for species
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

        unsigned int n_particles = patch.particles_m[is].size();

        for (unsigned int part = 0; part < n_particles; part++) {

          mini_float *pos[3] = {&patch.particles_m[is].x_(part),
                                &patch.particles_m[is].y_(part),
                                &patch.particles_m[is].z_(part)};

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
        } // End loop on particles

      } // End loop on species

      // Reflective conditions
    } else if (params.boundary_condition_code == 2) {
      for (int is = 0; is < patch.n_species_m; is++) {

        unsigned int n_particles = patch.particles_m[is].size();

        Vector<mini_float> &x = patch.particles_m[is].x_;
        Vector<mini_float> &y = patch.particles_m[is].y_;
        Vector<mini_float> &z = patch.particles_m[is].z_;

        Vector<mini_float> &mx = patch.particles_m[is].mx_;
        Vector<mini_float> &my = patch.particles_m[is].my_;
        Vector<mini_float> &mz = patch.particles_m[is].mz_;

        for (unsigned int part = 0; part < n_particles; part++) {

          mini_float *pos[3] = {&x(part), &y(part), &z(part)};

          mini_float *momentum[3] = {&mx(part), &my(part), &mz(part)};

          for (int d = 0; d < 3; d++) {

            if (*pos[d] >= sup_global[d]) {

              *pos[d]      = 2 * sup_global[d] - *pos[d];
              *momentum[d] = -*momentum[d];

            } else if (*pos[d] < inf_global[d]) {

              *pos[d]      = 2 * inf_global[d] - *pos[d];
              *momentum[d] = -*momentum[d];
            }
          }

        } // End loop on particles

      } // End loop on species
    }   // if type of conditions
  }     // if on border
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

    const int n_particles = patch.particles_m[is].size();

    for (int ip = 0; ip < n_particles; ++ip) {

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

    patch.vec_Jx_m[is].reset(minipic::host);
    patch.vec_Jy_m[is].reset(minipic::host);
    patch.vec_Jz_m[is].reset(minipic::host);

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

    const int n_particles = patch.particles_m[is].size();
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

      Field<mini_float> &Jx_loc = patch.vec_Jx_m[is];
      Field<mini_float> &Jy_loc = patch.vec_Jy_m[is];
      Field<mini_float> &Jz_loc = patch.vec_Jz_m[is];

      Vector<mini_float> &w = patch.particles_m[is].weight_;

      Vector<mini_float> &x = patch.particles_m[is].x_;
      Vector<mini_float> &y = patch.particles_m[is].y_;
      Vector<mini_float> &z = patch.particles_m[is].z_;

      for (int part = 0; part < n_particles; ++part) {

        const mini_float gamma_inv =
          1 / sqrt(1 + patch.particles_m[is].mx_h(part) * patch.particles_m[is].mx_h(part) +
                   patch.particles_m[is].my_h(part) * patch.particles_m[is].my_h(part) +
                   patch.particles_m[is].mz_h(part) * patch.particles_m[is].mz_h(part));

        const mini_float charge_weight = inv_cell_volume_x_q * w(part);

        const mini_float vx = patch.particles_m[is].mx_h(part) * gamma_inv;
        const mini_float vy = patch.particles_m[is].my_h(part) * gamma_inv;
        const mini_float vz = patch.particles_m[is].mz_h(part) * gamma_inv;

        // Current from the particle
        const mini_float Jxp = vx * charge_weight;
        const mini_float Jyp = vy * charge_weight;
        const mini_float Jzp = vz * charge_weight;

        // Calculate normalized position relative to the patch
        // ixn = (particles_m[is].x(part) ) * params.inv_dx;
        // iyn = (particles_m[is].y(part) ) * params.inv_dy;
        // izn = (particles_m[is].z(part) ) * params.inv_dz;
        const mini_float posxn = (x(part) - 0.5 * dt * vx - xmin) * inv_dx + 1;
        const mini_float posyn = (y(part) - 0.5 * dt * vy - ymin) * inv_dy + 1;
        const mini_float poszn = (z(part) - 0.5 * dt * vz - zmin) * inv_dz + 1;

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

        mini_float coeffs[3];

        coeffs[0] = posxn - 0.5 - ixd;
        coeffs[1] = posyn - iyp;
        coeffs[2] = poszn - izp;

        // Project on Jx
        Jx_loc(ixd, iyp, izp) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx_loc(ixd, iyp, izp + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
        Jx_loc(ixd, iyp + 1, izp) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx_loc(ixd, iyp + 1, izp + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;
        Jx_loc(ixd + 1, iyp, izp) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx_loc(ixd + 1, iyp, izp + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
        Jx_loc(ixd + 1, iyp + 1, izp) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx_loc(ixd + 1, iyp + 1, izp + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;

        coeffs[0] = posxn - ixp;
        coeffs[1] = posyn - 0.5 - iyd;
        coeffs[2] = poszn - izp;

        // compute_projection(ixp,
        //                    iyd,
        //                    izp,
        //                    vec_Jy_m[is].nx(),
        //                    vec_Jy_m[is].ny(),
        //                    vec_Jy_m[is].nz(),
        //                    coeffs,
        //                    vec_Jy_m[is].data(),
        //                    Jyp);

        // compute_projection(coeffs,
        //                    Jyp,
        //                    Jy_loc(ixp, iyd, izp),
        //                    Jy_loc(ixp, iyd, izp + 1),
        //                    Jy_loc(ixp, iyd + 1, izp),
        //                    Jy_loc(ixp, iyd + 1, izp + 1),
        //                    Jy_loc(ixp + 1, iyd, izp),
        //                    Jy_loc(ixp + 1, iyd, izp + 1),
        //                    Jy_loc(ixp + 1, iyd + 1, izp),
        //                    Jy_loc(ixp + 1, iyd + 1, izp + 1));

        Jy_loc(ixp, iyd, izp) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy_loc(ixp, iyd, izp + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
        Jy_loc(ixp, iyd + 1, izp) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy_loc(ixp, iyd + 1, izp + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;
        Jy_loc(ixp + 1, iyd, izp) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy_loc(ixp + 1, iyd, izp + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
        Jy_loc(ixp + 1, iyd + 1, izp) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy_loc(ixp + 1, iyd + 1, izp + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;

        coeffs[0] = posxn - ixp;
        coeffs[1] = posyn - iyp;
        coeffs[2] = poszn - 0.5 - izd;

        // compute_projection(ixp,
        //                    iyp,
        //                    izd,
        //                    vec_Jz_m[is].nx(),
        //                    vec_Jz_m[is].ny(),
        //                    vec_Jz_m[is].nz(),
        //                    coeffs,
        //                    vec_Jz_m[is].data(),
        //                    Jzp);

        // compute_projection(coeffs,
        //                    Jzp,
        //                    Jz_loc(ixp, iyp, izd),
        //                    Jz_loc(ixp, iyp, izd + 1),
        //                    Jz_loc(ixp, iyp + 1, izd),
        //                    Jz_loc(ixp, iyp + 1, izd + 1),
        //                    Jz_loc(ixp + 1, iyp, izd),
        //                    Jz_loc(ixp + 1, iyp, izd + 1),
        //                    Jz_loc(ixp + 1, iyp + 1, izd),
        //                    Jz_loc(ixp + 1, iyp + 1, izd + 1));

        Jz_loc(ixp, iyp, izd) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz_loc(ixp, iyp, izd + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
        Jz_loc(ixp, iyp + 1, izd) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz_loc(ixp, iyp + 1, izd + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
        Jz_loc(ixp + 1, iyp, izd) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz_loc(ixp + 1, iyp, izd + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
        Jz_loc(ixp + 1, iyp + 1, izd) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz_loc(ixp + 1, iyp + 1, izd + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
      } // end for each particles

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
auto project(Params &params, ElectroMagn &em, Patch &patch) -> void {

  Field<mini_float> &Jx = em.Jx_m;
  Field<mini_float> &Jy = em.Jy_m;
  Field<mini_float> &Jz = em.Jz_m;

  const double dt = params.dt;

  const double inv_dx = params.inv_dx;
  const double inv_dy = params.inv_dy;
  const double inv_dz = params.inv_dz;

  for (int is = 0; is < patch.n_species_m; is++) {

    const int n_particles                = patch.particles_m[is].size();
    const mini_float inv_cell_volume_x_q = params.inv_cell_volume * patch.particles_m[is].charge_m;
    // double m       = particles_m[is].mass_m;

    Vector<mini_float> &w = patch.particles_m[is].weight_;

    Vector<mini_float> &x = patch.particles_m[is].x_;
    Vector<mini_float> &y = patch.particles_m[is].y_;
    Vector<mini_float> &z = patch.particles_m[is].z_;

    Vector<mini_float> &mx = patch.particles_m[is].mx_;
    Vector<mini_float> &my = patch.particles_m[is].my_;
    Vector<mini_float> &mz = patch.particles_m[is].mz_;

    for (int part = 0; part < n_particles; ++part) {

      // Delete if already compute by Pusher
      // mini_float usq = (moment[0]*moment[0] + moment[1]*moment[1] + moment[2]*moment[2]);
      // mini_float gamma = sqrt(1+usq);
      // gamma_inv = 1/gamma;

      const mini_float charge_weight = inv_cell_volume_x_q * w(part);

      const mini_float gamma_inv =
        1 / sqrt(1 + mx(part) * mx(part) + my(part) * my(part) + mz(part) * mz(part));

      const mini_float vx = mx(part) * gamma_inv;
      const mini_float vy = my(part) * gamma_inv;
      const mini_float vz = mz(part) * gamma_inv;

      const mini_float Jxp = vx * charge_weight;
      const mini_float Jyp = vy * charge_weight;
      const mini_float Jzp = vz * charge_weight;

      // Calculate normalized positions
      // We come back 1/2 time step back in time for the position because of the leap frog scheme
      // As a consequence, we also have `+ 1` because the current grids have 2 additional ghost
      // cells (1 the min and 1 at the max border) when the direction is primal
      const mini_float posxn = (x(part) - 0.5 * dt * vx) * inv_dx + 1;
      const mini_float posyn = (y(part) - 0.5 * dt * vy) * inv_dy + 1;
      const mini_float poszn = (z(part) - 0.5 * dt * vz) * inv_dz + 1;

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

      mini_float coeffs[3];

      coeffs[0] = posxn - 0.5 - ixd;
      coeffs[1] = posyn - iyp;
      coeffs[2] = poszn - izp;

      Jx(ixd, iyp, izp) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
      Jx(ixd, iyp, izp + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
      Jx(ixd, iyp + 1, izp) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
      Jx(ixd, iyp + 1, izp + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;
      Jx(ixd + 1, iyp, izp) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
      Jx(ixd + 1, iyp, izp + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
      Jx(ixd + 1, iyp + 1, izp) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
      Jx(ixd + 1, iyp + 1, izp + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;

      coeffs[0] = posxn - ixp;
      coeffs[1] = posyn - 0.5 - iyd;
      coeffs[2] = poszn - izp;

      Jy(ixp, iyd, izp) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
      Jy(ixp, iyd, izp + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
      Jy(ixp, iyd + 1, izp) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
      Jy(ixp, iyd + 1, izp + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;
      Jy(ixp + 1, iyd, izp) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
      Jy(ixp + 1, iyd, izp + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
      Jy(ixp + 1, iyd + 1, izp) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
      Jy(ixp + 1, iyd + 1, izp + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;

      coeffs[0] = posxn - ixp;
      coeffs[1] = posyn - iyp;
      coeffs[2] = poszn - 0.5 - izd;

      Jz(ixp, iyp, izd) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
      Jz(ixp, iyp, izd + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
      Jz(ixp, iyp + 1, izd) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
      Jz(ixp, iyp + 1, izd + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
      Jz(ixp + 1, iyp, izd) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
      Jz(ixp + 1, iyp, izd + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
      Jz(ixp + 1, iyp + 1, izd) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
      Jz(ixp + 1, iyp + 1, izd + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
    } // end for each particles
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

  /////     Solve Maxwell Ampere (E)
  // Electric field Ex (d,p,p)

#if defined(__MINIPIC_OMP__)
#pragma omp for schedule(runtime)
#endif
  for (unsigned int ix = 0; ix < em.nx_d_m; ix++) {

    profiler.start(MAXWELL);

    for (unsigned int iy = 0; iy < em.ny_p_m; iy++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < em.nz_p_m; iz++) {

        em.Ex_m(ix, iy, iz) += -dt * em.Jx_m(ix, iy + 1, iz + 1) +
                               dt_over_dy * (em.Bz_m(ix, iy + 1, iz) - em.Bz_m(ix, iy, iz)) -
                               dt_over_dz * (em.By_m(ix, iy, iz + 1) - em.By_m(ix, iy, iz));
      }
    }

    profiler.stop();
  }

  // Electric field Ey (p,d,p)
#if defined(__MINIPIC_OMP__)
#pragma omp for schedule(runtime)
#endif
  for (unsigned int ix = 0; ix < em.nx_p_m; ix++) {

    profiler.start(MAXWELL);

    for (unsigned int iy = 0; iy < em.ny_d_m; iy++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < em.nz_p_m; iz++) {

        em.Ey_m(ix, iy, iz) += -dt * em.Jy_m(ix + 1, iy, iz + 1) -
                               dt_over_dx * (em.Bz_m(ix + 1, iy, iz) - em.Bz_m(ix, iy, iz)) +
                               dt_over_dz * (em.Bx_m(ix, iy, iz + 1) - em.Bx_m(ix, iy, iz));
      }
    }

    profiler.stop();
  }

  // Electric field Ez (p,p,d)

#if defined(__MINIPIC_OMP__)
#pragma omp for schedule(runtime)
#endif
  for (unsigned int ix = 0; ix < em.nx_p_m; ix++) {

    profiler.start(MAXWELL);

    for (unsigned int iy = 0; iy < em.ny_p_m; iy++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < em.nz_d_m; iz++) {

        em.Ez_m(ix, iy, iz) += -dt * em.Jz_m(ix + 1, iy + 1, iz) +
                               dt_over_dx * (em.By_m(ix + 1, iy, iz) - em.By_m(ix, iy, iz)) -
                               dt_over_dy * (em.Bx_m(ix, iy + 1, iz) - em.Bx_m(ix, iy, iz));
      }
    }

    profiler.stop();
  }

  /////     Solve Maxwell Faraday (B)

  // Magnetic field Bx (p,d,d)

#if defined(__MINIPIC_OMP__)
#pragma omp for schedule(runtime)
#endif
  for (unsigned int ix = 0; ix < em.nx_p_m; ix++) {

    profiler.start(MAXWELL);

    for (unsigned int iy = 1; iy < em.ny_d_m - 1; iy++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 1; iz < em.nz_d_m - 1; iz++) {

        em.Bx_m(ix, iy, iz) += -dt_over_dy * (em.Ez_m(ix, iy, iz) - em.Ez_m(ix, iy - 1, iz)) +
                               dt_over_dz * (em.Ey_m(ix, iy, iz) - em.Ey_m(ix, iy, iz - 1));
      }
    }

    profiler.stop();
  }

  // Magnetic field By (d,p,d)

#if defined(__MINIPIC_OMP__)
#pragma omp for schedule(runtime)
#endif
  for (unsigned int ix = 1; ix < em.nx_d_m - 1; ix++) {

    profiler.start(MAXWELL);

    for (unsigned int iy = 0; iy < em.ny_p_m; iy++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 1; iz < em.nz_d_m - 1; iz++) {
        em.By_m(ix, iy, iz) += -dt_over_dz * (em.Ex_m(ix, iy, iz) - em.Ex_m(ix, iy, iz - 1)) +
                               dt_over_dx * (em.Ez_m(ix, iy, iz) - em.Ez_m(ix - 1, iy, iz));
      }
    }

    profiler.stop();
  }

  // Magnetic field Bz (d,d,p)

#if defined(__MINIPIC_OMP__)
#pragma omp for schedule(runtime)
#endif
  for (unsigned int ix = 1; ix < em.nx_d_m - 1; ix++) {

    profiler.start(MAXWELL);

    for (unsigned int iy = 1; iy < em.ny_d_m - 1; iy++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < em.nz_p_m; iz++) {
        em.Bz_m(ix, iy, iz) += -dt_over_dx * (em.Ey_m(ix, iy, iz) - em.Ey_m(ix - 1, iy, iz)) +
                               dt_over_dy * (em.Ex_m(ix, iy, iz) - em.Ex_m(ix, iy - 1, iz));
      }
    }

    profiler.stop();
  }

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

  if (params.boundary_condition == "periodic") {

    // X dim
    // By (d,p,d)
    for (unsigned int iy = 0; iy < ny_By; ++iy) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_By; ++iz) {
        // -X
        em.By_m(0, iy, iz)         = em.By_m(nx_By - 2, iy, iz);
        em.By_m(nx_By - 1, iy, iz) = em.By_m(1, iy, iz);
      }
    }

    // Bz (d,d,p)
    for (unsigned int iy = 0; iy < ny_Bz; iy++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Bz; iz++) {
        // -X
        em.Bz_m(0, iy, iz)         = em.Bz_m(nx_Bz - 2, iy, iz);
        em.Bz_m(nx_Bz - 1, iy, iz) = em.Bz_m(1, iy, iz);
      }
    }

    // Y dim
    // Bx (p,d,d)
    for (unsigned int ix = 0; ix < nx_Bx; ix++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Bx; iz++) {
        // -Y
        em.Bx_m(ix, 0, iz)         = em.Bx_m(ix, ny_Bx - 2, iz);
        em.Bx_m(ix, ny_Bx - 1, iz) = em.Bx_m(ix, 1, iz);
      }
    }
    // Bz (d,d,p)
    for (unsigned int ix = 0; ix < nx_Bz; ix++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Bz; iz++) {
        // -Y
        em.Bz_m(ix, 0, iz)         = em.Bz_m(ix, ny_Bz - 2, iz);
        em.Bz_m(ix, ny_Bz - 1, iz) = em.Bz_m(ix, 1, iz);
      }
    }

    // Z dim
    // Bx
    for (unsigned int ix = 0; ix < nx_Bx; ix++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iy = 0; iy < ny_Bx; iy++) {
        // -Z
        em.Bx_m(ix, iy, 0) = em.Bx_m(ix, iy, nz_Bx - 2);
        // +Z
        em.Bx_m(ix, iy, nz_Bx - 1) = em.Bx_m(ix, iy, 1);
      }
    }
    // By
    for (unsigned int ix = 0; ix < nx_By; ++ix) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iy = 0; iy < ny_By; ++iy) {
        // -Z
        em.By_m(ix, iy, 0) = em.By_m(ix, iy, nz_By - 2);
        // +Z
        em.By_m(ix, iy, nz_By - 1) = em.By_m(ix, iy, 1);
      }
    }

  } else if (params.boundary_condition == "reflective") {

    // X dim
    // By (d,p,d)
    for (unsigned int iy = 0; iy < ny_By; ++iy) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_By; ++iz) {
        // -X
        em.By_m(0, iy, iz) = em.By_m(1, iy, iz);
        // +X
        em.By_m(nx_By - 1, iy, iz) = em.By_m(nx_By - 2, iy, iz);
      }
    }

    // Bz (d,d,p)
    for (unsigned int iy = 0; iy < ny_Bz; iy++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
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
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iz = 0; iz < nz_Bx; iz++) {
        // -Y
        em.Bx_m(ix, 0, iz) = em.Bx_m(ix, 1, iz);
        // +Y
        em.Bx_m(ix, ny_Bx - 1, iz) = em.Bx_m(ix, ny_Bx - 2, iz);
      }
    }
    // Bz (-1 to avoid corner)
    for (unsigned int ix = 0; ix < nx_Bz; ix++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
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
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
      for (unsigned int iy = 0; iy < ny_Bx; iy++) {
        // -Z
        em.Bx_m(ix, iy, 0) = em.Bx_m(ix, iy, 1);
        // +Z
        em.Bx_m(ix, iy, nz_Bx - 1) = em.Bx_m(ix, iy, nz_Bx - 2);
      }
    }
    // By
    for (unsigned int ix = 0; ix < nx_By; ix++) {
#if defined(__MINIPIC_SIMD__)
#pragma omp simd
#endif
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
    unsigned int n_particles = patch.particles_m[is].size();

    if (n_particles == 0)
      continue;

    // Count number of particles to move
    Vector<int> n_particle_to_move(26, 0, backend);

    // index for particle copy in the buffers
    Vector<unsigned int> ip_to_move(26, 0, backend);

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
    for (int ip = 0; ip < n_particles; ip++) {

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

    unsigned int total_particles_to_remove = 0;

    for (int ib = 0; ib < 26; ib++) {
      patch.particles_to_move_m[is][ib].resize(n_particle_to_move.h(ib), minipic::device);
      total_particles_to_remove += n_particle_to_move.h(ib);

      //      if (n_particle_to_move.h(ib) > 0) {
      //         std::cerr << ib << " " << n_particle_to_move.h(ib) << std::endl;
      //      }
    }

    // 3 -  Move tagged particles in the corresponding buffer

    if (total_particles_to_remove > 0) {

      for (int ip = 0; ip < n_particles; ++ip) {

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

          const int i = ip_to_move.h(ib);

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
      int ip = 0;

      // last particle index
      int last_ip = n_particles - 1;

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

    const int number_of_particles = patch.particles_m[is].size();

    // total number of particles coming from other patches
    int coming_number_of_particles = 0;

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
      }   // end for each neighbors
    }     // end for each neighbors

    if (coming_number_of_particles > 0) {
      patch.particles_m[is].resize(number_of_particles + coming_number_of_particles,
                                   minipic::device);
    }

    // Index where to start to copy the incoming particles in Particles
    int ip_buffer_start = number_of_particles;

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
            for (int ip = 0; ip < buffer_size; ++ip) {

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
      }   // end for each neighbors
    }     // end for each neighbors

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

      for (int ix = 0; ix < patch.vec_Jx_m[is].nx(); ix++) {
        for (int iy = 0; iy < patch.vec_Jx_m[is].ny(); iy++) {
          for (int iz = 0; iz < patch.vec_Jx_m[is].nz(); iz++) {
            patch.vec_Jx_m[0](ix, iy, iz) += patch.vec_Jx_m[is](ix, iy, iz);
          }
        }
      }

      for (int ix = 0; ix < patch.vec_Jy_m[is].nx(); ++ix) {
        for (int iy = 0; iy < patch.vec_Jy_m[is].ny(); ++iy) {
          for (int iz = 0; iz < patch.vec_Jy_m[is].nz(); ++iz) {
            patch.vec_Jy_m[0](ix, iy, iz) += patch.vec_Jy_m[is](ix, iy, iz);
          }
        }
      }

      for (int ix = 0; ix < patch.vec_Jz_m[0].nx(); ix++) {
        for (int iy = 0; iy < patch.vec_Jz_m[0].ny(); iy++) {
          for (int iz = 0; iz < patch.vec_Jz_m[0].nz(); iz++) {

            patch.vec_Jz_m[0](ix, iy, iz) += patch.vec_Jz_m[is](ix, iy, iz);
          }
        }
      }

    } // end check if particles
  }   // end for species
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

    for (int ix = 0; ix < patch.vec_Jx_m[0].nx(); ix++) {
      for (int iy = 0; iy < patch.vec_Jx_m[0].ny(); iy++) {
        for (int iz = 0; iz < patch.vec_Jx_m[0].nz(); iz++) {
          em.Jx_m(i_global_d + ix, j_global_p + iy, k_global_p + iz) +=
            patch.vec_Jx_m[0](ix, iy, iz);
        }
      }
    }

    for (int ix = 0; ix < patch.vec_Jy_m[0].nx(); ix++) {
      for (int iy = 0; iy < patch.vec_Jy_m[0].ny(); iy++) {
        for (int iz = 0; iz < patch.vec_Jy_m[0].nz(); iz++) {
          em.Jy_m(i_global_p + ix, j_global_d + iy, k_global_p + iz) +=
            patch.vec_Jy_m[0](ix, iy, iz);
        }
      }
    }

    for (int ix = 0; ix < patch.vec_Jz_m[0].nx(); ix++) {
      for (int iy = 0; iy < patch.vec_Jz_m[0].ny(); iy++) {
        for (int iz = 0; iz < patch.vec_Jz_m[0].nz(); iz++) {
          em.Jz_m(i_global_p + ix, j_global_p + iy, k_global_d + iz) +=
            patch.vec_Jz_m[0](ix, iy, iz);
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