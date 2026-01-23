/* _____________________________________________________________________ */
//! \file Operators.hpp

//! \brief contains generic kernels for the particle pusher

/* _____________________________________________________________________ */

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
auto interpolate(Params &params, ElectroMagn &em, Patch &patch) -> void {

  const auto inv_dx_m = em.inv_dx_m;
  const auto inv_dy_m = em.inv_dy_m;
  const auto inv_dz_m = em.inv_dz_m;

  const mini_float xmin = params.inf_x;
  const mini_float ymin = params.inf_y;
  const mini_float zmin = params.inf_z;

  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles = patch.particles_m[is].size();

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

    device_field_t Ex = em.Ex_m.data_m.d_view;
    device_field_t Ey = em.Ey_m.data_m.d_view;
    device_field_t Ez = em.Ez_m.data_m.d_view;

    device_field_t Bx = em.Bx_m.data_m.d_view;
    device_field_t By = em.By_m.data_m.d_view;
    device_field_t Bz = em.Bz_m.data_m.d_view;

    device_vector_t x = patch.particles_m[is].x_.data_.d_view;
    device_vector_t y = patch.particles_m[is].y_.data_.d_view;
    device_vector_t z = patch.particles_m[is].z_.data_.d_view;

    device_vector_t Exp = patch.particles_m[is].Ex_.data_.d_view;
    device_vector_t Eyp = patch.particles_m[is].Ey_.data_.d_view;
    device_vector_t Ezp = patch.particles_m[is].Ez_.data_.d_view;

    device_vector_t Bxp = patch.particles_m[is].Bx_.data_.d_view;
    device_vector_t Byp = patch.particles_m[is].By_.data_.d_view;
    device_vector_t Bzp = patch.particles_m[is].Bz_.data_.d_view;

#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

    device_field_t Ex = em.Ex_m.data_m;
    device_field_t Ey = em.Ey_m.data_m;
    device_field_t Ez = em.Ez_m.data_m;

    device_field_t Bx = em.Bx_m.data_m;
    device_field_t By = em.By_m.data_m;
    device_field_t Bz = em.Bz_m.data_m;

    device_vector_t x = patch.particles_m[is].x_.data_;
    device_vector_t y = patch.particles_m[is].y_.data_;
    device_vector_t z = patch.particles_m[is].z_.data_;

    device_vector_t Exp = patch.particles_m[is].Ex_.data_;
    device_vector_t Eyp = patch.particles_m[is].Ey_.data_;
    device_vector_t Ezp = patch.particles_m[is].Ez_.data_;

    device_vector_t Bxp = patch.particles_m[is].Bx_.data_;
    device_vector_t Byp = patch.particles_m[is].By_.data_;
    device_vector_t Bzp = patch.particles_m[is].Bz_.data_;

#endif

    Kokkos::parallel_for(
      n_particles,
      KOKKOS_LAMBDA(const size_t part) {

        // Calculate normalized positions
        const double ixn = (x(part) - xmin) * inv_dx_m;
        const double iyn = (y(part) - ymin) * inv_dy_m;
        const double izn = (z(part) - zmin) * inv_dz_m;

        // Compute indexes in global primal grid
        const unsigned int ixp = Kokkos::floor(ixn);
        const unsigned int iyp = Kokkos::floor(iyn);
        const unsigned int izp = Kokkos::floor(izn);

        // Compute indexes in global dual grid
        const unsigned int ixd = Kokkos::floor(ixn + 0.5);
        const unsigned int iyd = Kokkos::floor(iyn + 0.5);
        const unsigned int izd = Kokkos::floor(izn + 0.5);

        // Compute distances for interpolation
        const double dist_x_p = ixn - ixp;
        const double dist_y_p = iyn - iyp;
        const double dist_z_p = izn - izp;

        const double dist_x_d = (ixn + 0.5) - ixd;
        const double dist_y_d = (iyn + 0.5) - iyd;
        const double dist_z_d = (izn + 0.5) - izd;

        // Ex (d, p , p)
        {
          const double v00 =
            Ex(ixd, iyp, izp) * (1 - dist_x_d) + Ex(ixd + 1, iyp, izp) * dist_x_d;
          const double v01 =
            Ex(ixd, iyp, izp + 1) * (1 - dist_x_d) + Ex(ixd + 1, iyp, izp + 1) * dist_x_d;
          const double v10 =
            Ex(ixd, iyp + 1, izp) * (1 - dist_x_d) + Ex(ixd + 1, iyp + 1, izp) * dist_x_d;
          const double v11 =
            Ex(ixd, iyp + 1, izp + 1) * (1 - dist_x_d) + Ex(ixd + 1, iyp + 1, izp + 1) * dist_x_d;
          const double v0 = v00 * (1 - dist_y_p) + v10 * dist_y_p;
          const double v1 = v01 * (1 - dist_y_p) + v11 * dist_y_p;

          Exp(part) = v0 * (1 - dist_z_p) + v1 * dist_z_p;
        }

        // Ey (p, d, p)
        {
          const double v00 =
            Ey(ixp, iyd, izp) * (1 - dist_x_p) + Ey(ixp + 1, iyd, izp) * dist_x_p;
          const double v01 =
            Ey(ixp, iyd, izp + 1) * (1 - dist_x_p) + Ey(ixp + 1, iyd, izp + 1) * dist_x_p;
          const double v10 =
            Ey(ixp, iyd + 1, izp) * (1 - dist_x_p) + Ey(ixp + 1, iyd + 1, izp) * dist_x_p;
          const double v11 =
            Ey(ixp, iyd + 1, izp + 1) * (1 - dist_x_p) + Ey(ixp + 1, iyd + 1, izp + 1) * dist_x_p;
          const double v0 = v00 * (1 - dist_y_d) + v10 * dist_y_d;
          const double v1 = v01 * (1 - dist_y_d) + v11 * dist_y_d;

          Eyp(part) = v0 * (1 - dist_z_p) + v1 * dist_z_p;
        }

        // Ez (p, p, d)
        {
          const double v00 =
            Ez(ixp, iyp, izd) * (1 - dist_x_p) + Ez(ixp + 1, iyp, izd) * dist_x_p;
          const double v01 =
            Ez(ixp, iyp, izd + 1) * (1 - dist_x_p) + Ez(ixp + 1, iyp, izd + 1) * dist_x_p;
          const double v10 =
            Ez(ixp, iyp + 1, izd) * (1 - dist_x_p) + Ez(ixp + 1, iyp + 1, izd) * dist_x_p;
          const double v11 =
            Ez(ixp, iyp + 1, izd + 1) * (1 - dist_x_p) + Ez(ixp + 1, iyp + 1, izd + 1) * dist_x_p;
          const double v0 = v00 * (1 - dist_y_p) + v10 * dist_y_p;
          const double v1 = v01 * (1 - dist_y_p) + v11 * dist_y_p;

          Ezp(part) = v0 * (1 - dist_z_d) + v1 * dist_z_d;
        }

        // interpolation magnetic field
        // Bx (p, d, d)
        {
          const double v00 =
            Bx(ixp, iyd, izd) * (1 - dist_x_p) + Bx(ixp + 1, iyd, izd) * dist_x_p;
          const double v01 =
            Bx(ixp, iyd, izd + 1) * (1 - dist_x_p) + Bx(ixp + 1, iyd, izd + 1) * dist_x_p;
          const double v10 =
            Bx(ixp, iyd + 1, izd) * (1 - dist_x_p) + Bx(ixp + 1, iyd + 1, izd) * dist_x_p;
          const double v11 =
            Bx(ixp, iyd + 1, izd + 1) * (1 - dist_x_p) + Bx(ixp + 1, iyd + 1, izd + 1) * dist_x_p;
          const double v0 = v00 * (1 - dist_y_d) + v10 * dist_y_d;
          const double v1 = v01 * (1 - dist_y_d) + v11 * dist_y_d;

          Bxp(part) = v0 * (1 - dist_z_d) + v1 * dist_z_d;
        }

        // By (d, p, d)
        {
          const double v00 =
            By(ixd, iyp, izd) * (1 - dist_x_d) + By(ixd + 1, iyp, izd) * dist_x_d;
          const double v01 =
            By(ixd, iyp, izd + 1) * (1 - dist_x_d) + By(ixd + 1, iyp, izd + 1) * dist_x_d;
          const double v10 =
            By(ixd, iyp + 1, izd) * (1 - dist_x_d) + By(ixd + 1, iyp + 1, izd) * dist_x_d;
          const double v11 =
            By(ixd, iyp + 1, izd + 1) * (1 - dist_x_d) + By(ixd + 1, iyp + 1, izd + 1) * dist_x_d;
          const double v0 = v00 * (1 - dist_y_p) + v10 * dist_y_p;
          const double v1 = v01 * (1 - dist_y_p) + v11 * dist_y_p;

          Byp(part) = v0 * (1 - dist_z_d) + v1 * dist_z_d;
        }

        // Bz (d, d, p)
        {
          const double v00 =
            Bz(ixd, iyd, izp) * (1 - dist_x_d) + Bz(ixd + 1, iyd, izp) * dist_x_d;
          const double v01 =
            Bz(ixd, iyd, izp + 1) * (1 - dist_x_d) + Bz(ixd + 1, iyd, izp + 1) * dist_x_d;
          const double v10 =
            Bz(ixd, iyd + 1, izp) * (1 - dist_x_d) + Bz(ixd + 1, iyd + 1, izp) * dist_x_d;
          const double v11 =
            Bz(ixd, iyd + 1, izp + 1) * (1 - dist_x_d) + Bz(ixd + 1, iyd + 1, izp + 1) * dist_x_d;
          const double v0 = v00 * (1 - dist_y_d) + v10 * dist_y_d;
          const double v1 = v01 * (1 - dist_y_d) + v11 * dist_y_d;

          Bzp(part) = v0 * (1 - dist_z_p) + v1 * dist_z_p;
        }
      } // End for each particle

    ); // end KOKKOS PARALLEL

    Kokkos::fence();

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

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

    device_vector_t x = patch.particles_m[is].x_.data_.d_view;
    device_vector_t y = patch.particles_m[is].y_.data_.d_view;
    device_vector_t z = patch.particles_m[is].z_.data_.d_view;

    device_vector_t mx = patch.particles_m[is].mx_.data_.d_view;
    device_vector_t my = patch.particles_m[is].my_.data_.d_view;
    device_vector_t mz = patch.particles_m[is].mz_.data_.d_view;

    device_vector_t Exp = patch.particles_m[is].Ex_.data_.d_view;
    device_vector_t Eyp = patch.particles_m[is].Ey_.data_.d_view;
    device_vector_t Ezp = patch.particles_m[is].Ez_.data_.d_view;

    device_vector_t Bxp = patch.particles_m[is].Bx_.data_.d_view;
    device_vector_t Byp = patch.particles_m[is].By_.data_.d_view;
    device_vector_t Bzp = patch.particles_m[is].Bz_.data_.d_view;


#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

    device_vector_t x = patch.particles_m[is].x_.data_;
    device_vector_t y = patch.particles_m[is].y_.data_;
    device_vector_t z = patch.particles_m[is].z_.data_;

    device_vector_t mx = patch.particles_m[is].mx_.data_;
    device_vector_t my = patch.particles_m[is].my_.data_;
    device_vector_t mz = patch.particles_m[is].mz_.data_;

    device_vector_t Exp = patch.particles_m[is].Ex_.data_;
    device_vector_t Eyp = patch.particles_m[is].Ey_.data_;
    device_vector_t Ezp = patch.particles_m[is].Ez_.data_;

    device_vector_t Bxp = patch.particles_m[is].Bx_.data_;
    device_vector_t Byp = patch.particles_m[is].By_.data_;
    device_vector_t Bzp = patch.particles_m[is].Bz_.data_;

#endif

    Kokkos::parallel_for(
      n_particles,
      KOKKOS_LAMBDA(const size_t ip) {
        // 1/2 E
        mini_float px = qp * Exp(ip);
        mini_float py = qp * Eyp(ip);
        mini_float pz = qp * Ezp(ip);

        const mini_float ux = mx(ip) + px;
        const mini_float uy = my(ip) + py;
        const mini_float uz = mz(ip) + pz;

        // gamma-factor
        mini_float usq       = (ux * ux + uy * uy + uz * uz);
        mini_float gamma     = Kokkos::sqrt(1 + usq);
        mini_float gamma_inv = qp / gamma;

        // B, T = Transform to rotate the particle
        const mini_float tx  = gamma_inv * Bxp(ip);
        const mini_float ty  = gamma_inv * Byp(ip);
        const mini_float tz  = gamma_inv * Bzp(ip);
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
        usq   = (px * px + py * py + pz * pz);
        gamma = Kokkos::sqrt(1 + usq);

        // Update inverse gamma factor
        gamma_inv = 1 / gamma;

        // Update momentum
        mx(ip) = px;
        my(ip) = py;
        mz(ip) = pz;

        // Update positions
        x(ip) += mx(ip) * dt * gamma_inv;
        y(ip) += my(ip) * dt * gamma_inv;
        z(ip) += mz(ip) * dt * gamma_inv;
      });

    Kokkos::fence();

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

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

    device_vector_t mx = patch.particles_m[is].mx_.data_.d_view;
    device_vector_t my = patch.particles_m[is].my_.data_.d_view;
    device_vector_t mz = patch.particles_m[is].mz_.data_.d_view;

    device_vector_t Exp = patch.particles_m[is].Ex_.data_.d_view;
    device_vector_t Eyp = patch.particles_m[is].Ey_.data_.d_view;
    device_vector_t Ezp = patch.particles_m[is].Ez_.data_.d_view;

    device_vector_t Bxp = patch.particles_m[is].Bx_.data_.d_view;
    device_vector_t Byp = patch.particles_m[is].By_.data_.d_view;
    device_vector_t Bzp = patch.particles_m[is].Bz_.data_.d_view;

#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

    device_vector_t mx = patch.particles_m[is].mx_.data_;
    device_vector_t my = patch.particles_m[is].my_.data_;
    device_vector_t mz = patch.particles_m[is].mz_.data_;

    device_vector_t Exp = patch.particles_m[is].Ex_.data_;
    device_vector_t Eyp = patch.particles_m[is].Ey_.data_;
    device_vector_t Ezp = patch.particles_m[is].Ez_.data_;

    device_vector_t Bxp = patch.particles_m[is].Bx_.data_;
    device_vector_t Byp = patch.particles_m[is].By_.data_;
    device_vector_t Bzp = patch.particles_m[is].Bz_.data_;

#endif

    Kokkos::parallel_for(
      n_particles,
      KOKKOS_LAMBDA(const size_t ip) {
        // 1/2 E
        mini_float px = qp * Exp(ip);
        mini_float py = qp * Eyp(ip);
        mini_float pz = qp * Ezp(ip);

        const mini_float ux = mx(ip) + px;
        const mini_float uy = my(ip) + py;
        const mini_float uz = mz(ip) + pz;

        // gamma-factor
        mini_float usq       = (ux * ux + uy * uy + uz * uz);
        mini_float gamma     = Kokkos::sqrt(1 + usq);
        mini_float gamma_inv = qp / gamma;

        // B, T = Transform to rotate the particle
        const mini_float tx  = gamma_inv * Bxp(ip);
        const mini_float ty  = gamma_inv * Byp(ip);
        const mini_float tz  = gamma_inv * Bzp(ip);
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
        usq   = (px * px + py * py + pz * pz);
        gamma = Kokkos::sqrt(1 + usq);

        // Update inverse gamma factor
        gamma_inv = 1 / gamma;

        // Update momentum
        mx(ip) = px;
        my(ip) = py;
        mz(ip) = pz;
      });

    Kokkos::fence();
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

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

        device_vector_t x = patch.particles_m[is].x_.data_.d_view;
        device_vector_t y = patch.particles_m[is].y_.data_.d_view;
        device_vector_t z = patch.particles_m[is].z_.data_.d_view;

#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

        device_vector_t x = patch.particles_m[is].x_.data_;
        device_vector_t y = patch.particles_m[is].y_.data_;
        device_vector_t z = patch.particles_m[is].z_.data_;

#endif

        Kokkos::parallel_for(
          n_particles,
          KOKKOS_LAMBDA(const size_t part) {
            mini_float *pos[3] = {&x(part), &y(part), &z(part)};

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

        );

        Kokkos::fence();

      } // End loop on species

      // Reflective conditions
    } else if (params.boundary_condition_code == 2) {
      for (int is = 0; is < patch.n_species_m; is++) {

        size_t n_particles = patch.particles_m[is].size();

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

        device_vector_t x = patch.particles_m[is].x_.data_.d_view;
        device_vector_t y = patch.particles_m[is].y_.data_.d_view;
        device_vector_t z = patch.particles_m[is].z_.data_.d_view;

        device_vector_t mx = patch.particles_m[is].mx_.data_.d_view;
        device_vector_t my = patch.particles_m[is].my_.data_.d_view;
        device_vector_t mz = patch.particles_m[is].mz_.data_.d_view;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__) || defined(__MINIPIC_KOKKOS_VIEWS__)

        device_vector_t x = patch.particles_m[is].x_.data_;
        device_vector_t y = patch.particles_m[is].y_.data_;
        device_vector_t z = patch.particles_m[is].z_.data_;

        device_vector_t mx = patch.particles_m[is].mx_.data_;
        device_vector_t my = patch.particles_m[is].my_.data_;
        device_vector_t mz = patch.particles_m[is].mz_.data_;

#endif

        Kokkos::parallel_for(
          n_particles,
          KOKKOS_LAMBDA(const size_t part) {
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

        );

        Kokkos::fence();

      } // End loop on species
    } // if type of conditions
  } // if on border
}

// _______________________________________________________________________
//
//! \brief Current projection from global particles position to local grid
//! \param params global parameters
//! \param patch current patch
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

      const double inv_cell_volume_x_q = params.inv_cell_volume * patch.particles_m[is].charge_m;
      const double dt                  = params.dt;

      const double inv_dx = params.inv_dx;
      const double inv_dy = params.inv_dy;
      const double inv_dz = params.inv_dz;

      const double xmin = patch.inf_m[0];
      const double ymin = patch.inf_m[1];
      const double zmin = patch.inf_m[2];

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

      device_field_t Jx_species = patch.vec_Jx_m[is].data_m.d_view;
      device_field_t Jy_species = patch.vec_Jy_m[is].data_m.d_view;
      device_field_t Jz_species = patch.vec_Jz_m[is].data_m.d_view;

      device_vector_t w = patch.particles_m[is].weight_.data_.d_view;

      device_vector_t x = patch.particles_m[is].x_.data_.d_view;
      device_vector_t y = patch.particles_m[is].y_.data_.d_view;
      device_vector_t z = patch.particles_m[is].z_.data_.d_view;

      device_vector_t mx = patch.particles_m[is].mx_.data_.d_view;
      device_vector_t my = patch.particles_m[is].my_.data_.d_view;
      device_vector_t mz = patch.particles_m[is].mz_.data_.d_view;

#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

      device_field_t Jx_species = patch.vec_Jx_m[is].data_m;
      device_field_t Jy_species = patch.vec_Jy_m[is].data_m;
      device_field_t Jz_species = patch.vec_Jz_m[is].data_m;

      device_vector_t w = patch.particles_m[is].weight_.data_;

      device_vector_t x = patch.particles_m[is].x_.data_;
      device_vector_t y = patch.particles_m[is].y_.data_;
      device_vector_t z = patch.particles_m[is].z_.data_;

      device_vector_t mx = patch.particles_m[is].mx_.data_;
      device_vector_t my = patch.particles_m[is].my_.data_;
      device_vector_t mz = patch.particles_m[is].mz_.data_;

#endif

#if defined(__MINIPIC_KOKKOS_SCATTERVIEW__)
      Kokkos::Experimental::ScatterView<double ***> scatter_Jx_loc(Jx_species);
      Kokkos::Experimental::ScatterView<double ***> scatter_Jy_loc(Jy_species);
      Kokkos::Experimental::ScatterView<double ***> scatter_Jz_loc(Jz_species);
#else
      Kokkos::View<double ***, Kokkos::MemoryTraits<Kokkos::Atomic>> Jx_loc = Jx_species;
      Kokkos::View<double ***, Kokkos::MemoryTraits<Kokkos::Atomic>> Jy_loc = Jy_species;
      Kokkos::View<double ***, Kokkos::MemoryTraits<Kokkos::Atomic>> Jz_loc = Jz_species;
#endif

      Kokkos::parallel_for(
        n_particles,
        KOKKOS_LAMBDA(const size_t part) {

#if defined(__MINIPIC_KOKKOS_SCATTERVIEW__)
          auto Jx_loc = scatter_Jx_loc.access();
          auto Jy_loc = scatter_Jy_loc.access();
          auto Jz_loc = scatter_Jz_loc.access();
#endif

          const double charge_weight = inv_cell_volume_x_q * w(part);

          // gamma-factor
          const double gamma_inv =
            1 / Kokkos::sqrt(1 + (mx(part) * mx(part) + my(part) * my(part) + mz(part) * mz(part)));

          const double vx = mx(part) * gamma_inv;
          const double vy = my(part) * gamma_inv;
          const double vz = mz(part) * gamma_inv;

          // Current from the particle
          const double Jxp = vx * charge_weight;
          const double Jyp = vy * charge_weight;
          const double Jzp = vz * charge_weight;

          // Calculate normalized position relative to the patch
          // ixn = (particles_m[is].x(part) ) * params.inv_dx;
          // iyn = (particles_m[is].y(part) ) * params.inv_dy;
          // izn = (particles_m[is].z(part) ) * params.inv_dz;
          const double posxn = (x(part) - 0.5 * dt * vx - xmin) * inv_dx + 1;
          const double posyn = (y(part) - 0.5 * dt * vy - ymin) * inv_dy + 1;
          const double poszn = (z(part) - 0.5 * dt * vz - zmin) * inv_dz + 1;

          // Compute indexes in primal grid
          const int ixp = (int)(Kokkos::floor(posxn));
          const int iyp = (int)(Kokkos::floor(posyn));
          const int izp = (int)(Kokkos::floor(poszn));

          // Compute indexes in dual grid
          // For the current, the dual grid is 0.5 * dx shorter on each side of the grid (if dual
          // directions only)
          const int ixd = (int)(Kokkos::floor(posxn - 0.5));
          const int iyd = (int)(Kokkos::floor(posyn - 0.5));
          const int izd = (int)(Kokkos::floor(poszn - 0.5));

          // Projection particle on currant field
          // Compute interpolation coeff, p = primal, d = dual

          double coeffs[3];

          coeffs[0] = posxn - 0.5 - ixd;
          coeffs[1] = posyn - iyp;
          coeffs[2] = poszn - izp;

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

          Jz_loc(ixp, iyp, izd) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
          Jz_loc(ixp, iyp, izd + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
          Jz_loc(ixp, iyp + 1, izd) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
          Jz_loc(ixp, iyp + 1, izd + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
          Jz_loc(ixp + 1, iyp, izd) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
          Jz_loc(ixp + 1, iyp, izd + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
          Jz_loc(ixp + 1, iyp + 1, izd) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
          Jz_loc(ixp + 1, iyp + 1, izd + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
        } // end for each particles

      );

      Kokkos::fence();

#if defined(__MINIPIC_KOKKOS_SCATTERVIEW__)
      Kokkos::Experimental::contribute(Jx_species, scatter_Jx_loc);
      Kokkos::Experimental::contribute(Jy_species, scatter_Jy_loc);
      Kokkos::Experimental::contribute(Jz_species, scatter_Jz_loc);
#endif

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
void project(Params &params, ElectroMagn &em, Patch &patch) {

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

  device_field_t Jx_device = em.Jx_m.data_m.d_view;
  device_field_t Jy_device = em.Jy_m.data_m.d_view;
  device_field_t Jz_device = em.Jz_m.data_m.d_view;

#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

  device_field_t Jx_device = em.Jx_m.data_m;
  device_field_t Jy_device = em.Jy_m.data_m;
  device_field_t Jz_device = em.Jz_m.data_m;

#endif

#if defined(__MINIPIC_KOKKOS_SCATTERVIEW__)
  // Use ScatterView
  Kokkos::Experimental::ScatterView<double ***> scatter_Jx(Jx_device);
  Kokkos::Experimental::ScatterView<double ***> scatter_Jy(Jy_device);
  Kokkos::Experimental::ScatterView<double ***> scatter_Jz(Jz_device);
#else
  // Use atomic memory traits
  Kokkos::View<double ***, Kokkos::MemoryTraits<Kokkos::Atomic>> Jx(Jx_device);
  Kokkos::View<double ***, Kokkos::MemoryTraits<Kokkos::Atomic>> Jy(Jy_device);
  Kokkos::View<double ***, Kokkos::MemoryTraits<Kokkos::Atomic>> Jz(Jz_device);
#endif

  const double dt = params.dt;

  const double inv_dx = params.inv_dx;
  const double inv_dy = params.inv_dy;
  const double inv_dz = params.inv_dz;

#if (__MINIPIC_DEBUG__)
  int nx_Jx = em.Jx_m.nx_m;
  int ny_Jx = em.Jx_m.ny_m;
  int nz_Jx = em.Jx_m.nz_m;

  int nx_Jy = em.Jy_m.nx_m;
  int ny_Jy = em.Jy_m.ny_m;
  int nz_Jy = em.Jy_m.nz_m;

#endif

  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t n_particles            = patch.particles_m[is].size();
    const double inv_cell_volume_x_q = params.inv_cell_volume * patch.particles_m[is].charge_m;
    // double m       = particles_m[is].mass_m;

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

    device_vector_t w = patch.particles_m[is].weight_.data_.d_view;

    device_vector_t x = patch.particles_m[is].x_.data_.d_view;
    device_vector_t y = patch.particles_m[is].y_.data_.d_view;
    device_vector_t z = patch.particles_m[is].z_.data_.d_view;

    device_vector_t mx = patch.particles_m[is].mx_.data_.d_view;
    device_vector_t my = patch.particles_m[is].my_.data_.d_view;
    device_vector_t mz = patch.particles_m[is].mz_.data_.d_view;

#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

    device_vector_t w = patch.particles_m[is].weight_.data_;
    device_vector_t x = patch.particles_m[is].x_.data_;
    device_vector_t y = patch.particles_m[is].y_.data_;
    device_vector_t z = patch.particles_m[is].z_.data_;
    device_vector_t mx = patch.particles_m[is].mx_.data_;
    device_vector_t my = patch.particles_m[is].my_.data_;
    device_vector_t mz = patch.particles_m[is].mz_.data_;

#endif

    Kokkos::parallel_for(
      n_particles,
      KOKKOS_LAMBDA(const size_t part) {

#if defined(__MINIPIC_KOKKOS_SCATTERVIEW__)
        auto Jx = scatter_Jx.access();
        auto Jy = scatter_Jy.access();
        auto Jz = scatter_Jz.access();
#endif

        // Delete if already compute by Pusher
        // double usq = (moment[0]*moment[0] + moment[1]*moment[1] + moment[2]*moment[2]);
        // double gamma = sqrt(1+usq);
        // gamma_inv = 1/gamma;

        const double charge_weight = inv_cell_volume_x_q * w(part);

        const double gamma_inv =
          1 / Kokkos::sqrt(1 + (mx(part) * mx(part) + my(part) * my(part) + mz(part) * mz(part)));

        const double vx = mx(part) * gamma_inv;
        const double vy = my(part) * gamma_inv;
        const double vz = mz(part) * gamma_inv;

        const double Jxp = vx * charge_weight;
        const double Jyp = vy * charge_weight;
        const double Jzp = vz * charge_weight;

        // Calculate normalized positions
        // We come back 1/2 time step back in time for the position because of the leap frog scheme
        // As a consequence, we also have `+ 1` because the current grids have 2 additional ghost
        // cells (1 the min and 1 at the max border) when the direction is primal
        const double posxn = (x(part) - 0.5 * dt * vx) * inv_dx + 1;
        const double posyn = (y(part) - 0.5 * dt * vy) * inv_dy + 1;
        const double poszn = (z(part) - 0.5 * dt * vz) * inv_dz + 1;

        // Compute indexes in primal grid
        const int ixp = (int)(Kokkos::floor(posxn)); //- i_patch_topology_m * nx_cells_m;
        const int iyp = (int)(Kokkos::floor(posyn)); //- j_patch_topology_m * ny_cells_m;
        const int izp = (int)(Kokkos::floor(poszn)); //- k_patch_topology_m * nz_cells_m;

        // Compute indexes in dual grid
        const int ixd = (int)Kokkos::floor(posxn - 0.5); //- i_patch_topology_m * nx_cells_m;
        const int iyd = (int)Kokkos::floor(posyn - 0.5); //- j_patch_topology_m * ny_cells_m;
        const int izd = (int)Kokkos::floor(poszn - 0.5); //- k_patch_topology_m * nz_cells_m;

#if (__MINIPIC_DEBUG__)
        // Check if the indexes are in the correct range
        if (ixd < 0 || ixd + 1 >= nx_Jx || iyp < 0 || iyp + 1 >= ny_Jx || izp < 0 ||
            izp + 1 >= nz_Jx) {
          Kokkos::printf("Error: part = %d\n", part);
          Kokkos::printf("Error: ixp = %d, iyp = %d, izp = %d\n", ixp, iyp, izp);
          Kokkos::printf("Error: posxn = %f, posyn = %f, poszn = %f\n", posxn, posyn, poszn);
          Kokkos::printf("Error: x = %f, y = %f, z = %f\n", x(part), y(part), z(part));
          Kokkos::printf("Error: vx = %f, vy = %f, vz = %f\n", vx, vy, vz);
          Kokkos::printf("Error: gamma_inv = %f\n", gamma_inv);
        }
#endif

        // Projection particle on currant field
        // Compute interpolation coeff, p = primal, d = dual

        double coeffs[3];

        coeffs[0] = posxn - 0.5 - ixd;
        coeffs[1] = posyn - iyp;
        coeffs[2] = poszn - izp;

        // #if defined(__MINIPIC_KOKKOS_SCATTERVIEW__)
        Jx(ixd, iyp, izp) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx(ixd, iyp, izp + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
        Jx(ixd, iyp + 1, izp) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx(ixd, iyp + 1, izp + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;
        Jx(ixd + 1, iyp, izp) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx(ixd + 1, iyp, izp + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp;
        Jx(ixd + 1, iyp + 1, izp) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp;
        Jx(ixd + 1, iyp + 1, izp + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jxp;
        // #else
        //  // Use kokkos atomics
        //  Kokkos::atomic_add(&Jx_device(ixd, iyp, izp), (1 - coeffs[0]) * (1 - coeffs[1]) * (1 -
        //  coeffs[2]) * Jxp); Kokkos::atomic_add(&Jx_device(ixd, iyp, izp + 1), (1 - coeffs[0]) *
        //  (1 - coeffs[1]) * (coeffs[2]) * Jxp); Kokkos::atomic_add(&Jx_device(ixd, iyp + 1, izp),
        //  (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jxp);
        //  Kokkos::atomic_add(&Jx_device(ixd, iyp + 1, izp + 1), (1 - coeffs[0]) * (coeffs[1]) *
        //  (coeffs[2]) * Jxp); Kokkos::atomic_add(&Jx_device(ixd + 1, iyp, izp), (coeffs[0]) * (1 -
        //  coeffs[1]) * (1 - coeffs[2]) * Jxp); Kokkos::atomic_add(&Jx_device(ixd + 1, iyp, izp +
        //  1), (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jxp);
        //  Kokkos::atomic_add(&Jx_device(ixd + 1, iyp + 1, izp), (coeffs[0]) * (coeffs[1]) * (1 -
        //  coeffs[2]) * Jxp); Kokkos::atomic_add(&Jx_device(ixd + 1, iyp + 1, izp + 1), (coeffs[0])
        //  * (coeffs[1]) * (coeffs[2]) * Jxp);
        // #endif

        coeffs[0] = posxn - ixp;
        coeffs[1] = posyn - 0.5 - iyd;
        coeffs[2] = poszn - izp;

        // #if defined(__MINIPIC_KOKKOS_SCATTERVIEW__)
        Jy(ixp, iyd, izp) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy(ixp, iyd, izp + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
        Jy(ixp, iyd + 1, izp) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy(ixp, iyd + 1, izp + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;
        Jy(ixp + 1, iyd, izp) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy(ixp + 1, iyd, izp + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp;
        Jy(ixp + 1, iyd + 1, izp) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp;
        Jy(ixp + 1, iyd + 1, izp + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jyp;
        // #else
        //  Use kokkos atomics
        //  Kokkos::atomic_add(&Jy_device(ixp, iyd, izp), (1 - coeffs[0]) * (1 - coeffs[1]) * (1 -
        //  coeffs[2]) * Jyp); Kokkos::atomic_add(&Jy_device(ixp, iyd, izp + 1), (1 - coeffs[0]) *
        //  (1 - coeffs[1]) * (coeffs[2]) * Jyp); Kokkos::atomic_add(&Jy_device(ixp, iyd + 1, izp),
        //  (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jyp);
        //  Kokkos::atomic_add(&Jy_device(ixp, iyd + 1, izp + 1), (1 - coeffs[0]) * (coeffs[1]) *
        //  (coeffs[2]) * Jyp); Kokkos::atomic_add(&Jy_device(ixp + 1, iyd, izp), (coeffs[0]) * (1 -
        //  coeffs[1]) * (1 - coeffs[2]) * Jyp); Kokkos::atomic_add(&Jy_device(ixp + 1, iyd, izp +
        //  1), (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jyp);
        //  Kokkos::atomic_add(&Jy_device(ixp + 1, iyd + 1, izp), (coeffs[0]) * (coeffs[1]) * (1 -
        //  coeffs[2]) * Jyp); Kokkos::atomic_add(&Jy_device(ixp + 1, iyd + 1, izp + 1), (coeffs[0])
        //  * (coeffs[1]) * (coeffs[2]) * Jyp);
        // #endif

        coeffs[0] = posxn - ixp;
        coeffs[1] = posyn - iyp;
        coeffs[2] = poszn - 0.5 - izd;

        // #if defined(__MINIPIC_KOKKOS_SCATTERVIEW__)
        Jz(ixp, iyp, izd) += (1 - coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz(ixp, iyp, izd + 1) += (1 - coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
        Jz(ixp, iyp + 1, izd) += (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz(ixp, iyp + 1, izd + 1) += (1 - coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
        Jz(ixp + 1, iyp, izd) += (coeffs[0]) * (1 - coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz(ixp + 1, iyp, izd + 1) += (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp;
        Jz(ixp + 1, iyp + 1, izd) += (coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp;
        Jz(ixp + 1, iyp + 1, izd + 1) += (coeffs[0]) * (coeffs[1]) * (coeffs[2]) * Jzp;
        // #else
        //  // Use kokkos atomics
        //  Kokkos::atomic_add(&Jz_device(ixp, iyp, izd), (1 - coeffs[0]) * (1 - coeffs[1]) * (1 -
        //  coeffs[2]) * Jzp); Kokkos::atomic_add(&Jz_device(ixp, iyp, izd + 1), (1 - coeffs[0]) *
        //  (1 - coeffs[1]) * (coeffs[2]) * Jzp); Kokkos::atomic_add(&Jz_device(ixp, iyp + 1, izd),
        //  (1 - coeffs[0]) * (coeffs[1]) * (1 - coeffs[2]) * Jzp);
        //  Kokkos::atomic_add(&Jz_device(ixp, iyp + 1, izd + 1), (1 - coeffs[0]) * (coeffs[1]) *
        //  (coeffs[2]) * Jzp); Kokkos::atomic_add(&Jz_device(ixp + 1, iyp, izd), (coeffs[0]) * (1 -
        //  coeffs[1]) * (1 - coeffs[2]) * Jzp); Kokkos::atomic_add(&Jz_device(ixp + 1, iyp, izd +
        //  1), (coeffs[0]) * (1 - coeffs[1]) * (coeffs[2]) * Jzp);
        //  Kokkos::atomic_add(&Jz_device(ixp + 1, iyp + 1, izd), (coeffs[0]) * (coeffs[1]) * (1 -
        //  coeffs[2]) * Jzp); Kokkos::atomic_add(&Jz_device(ixp + 1, iyp + 1, izd + 1), (coeffs[0])
        //  * (coeffs[1]) * (coeffs[2]) * Jzp);
        // #endif
      } // end for each particles
    );

    Kokkos::fence();
#if defined(__MINIPIC_KOKKOS_SCATTERVIEW__)
    Kokkos::Experimental::contribute(Jx_device, scatter_Jx);
    Kokkos::Experimental::contribute(Jy_device, scatter_Jy);
    Kokkos::Experimental::contribute(Jz_device, scatter_Jz);
#endif
    // particles_m[is].sync(minipic::device, minipic::host);
  }
}

// _______________________________________________________
//
//! \brief Solve Maxwell equations to compute EM fields
//! \param params global parameters
// _______________________________________________________
auto solve_maxwell(const Params &params, ElectroMagn &em) -> void {

  const double dt         = params.dt;
  const double dt_over_dx = params.dt * params.inv_dx;
  const double dt_over_dy = params.dt * params.inv_dy;
  const double dt_over_dz = params.dt * params.inv_dz;

  /////     Solve Maxwell Ampere (E)
  // Electric field Ex (d,p,p)

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

  device_field_t Jx = em.Jx_m.data_m.d_view;
  device_field_t Jy = em.Jy_m.data_m.d_view;
  device_field_t Jz = em.Jz_m.data_m.d_view;

  device_field_t Ex = em.Ex_m.data_m.d_view;
  device_field_t Ey = em.Ey_m.data_m.d_view;
  device_field_t Ez = em.Ez_m.data_m.d_view;

  device_field_t Bx = em.Bx_m.data_m.d_view;
  device_field_t By = em.By_m.data_m.d_view;
  device_field_t Bz = em.Bz_m.data_m.d_view;

#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

  device_field_t Jx = em.Jx_m.data_m;
  device_field_t Jy = em.Jy_m.data_m;
  device_field_t Jz = em.Jz_m.data_m;

  device_field_t Ex = em.Ex_m.data_m;
  device_field_t Ey = em.Ey_m.data_m;
  device_field_t Ez = em.Ez_m.data_m;

  device_field_t Bx = em.Bx_m.data_m;
  device_field_t By = em.By_m.data_m;
  device_field_t Bz = em.Bz_m.data_m;

#endif

  typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;

  Kokkos::parallel_for(
    mdrange_policy({0, 0, 0}, {em.nx_d_m, em.ny_p_m, em.nz_p_m}),
    KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
      Ex(ix, iy, iz) += -dt * Jx(ix, iy + 1, iz + 1) +
                        dt_over_dy * (Bz(ix, iy + 1, iz) - Bz(ix, iy, iz)) -
                        dt_over_dz * (By(ix, iy, iz + 1) - By(ix, iy, iz));
    });

  // Electric field Ey (p,d,p)
  Kokkos::parallel_for(
    mdrange_policy({0, 0, 0}, {em.nx_p_m, em.ny_d_m, em.nz_p_m}),
    KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
      Ey(ix, iy, iz) += -dt * Jy(ix + 1, iy, iz + 1) -
                        dt_over_dx * (Bz(ix + 1, iy, iz) - Bz(ix, iy, iz)) +
                        dt_over_dz * (Bx(ix, iy, iz + 1) - Bx(ix, iy, iz));
    });

  // Electric field Ez (p,p,d)

  Kokkos::parallel_for(
    mdrange_policy({0, 0, 0}, {em.nx_p_m, em.ny_p_m, em.nz_d_m}),
    KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
      Ez(ix, iy, iz) += -dt * Jz(ix + 1, iy + 1, iz) +
                        dt_over_dx * (By(ix + 1, iy, iz) - By(ix, iy, iz)) -
                        dt_over_dy * (Bx(ix, iy + 1, iz) - Bx(ix, iy, iz));
    });

  Kokkos::fence();

  /////     Solve Maxwell Faraday (B)

  // Magnetic field Bx (p,d,d)

  Kokkos::parallel_for(
    mdrange_policy({0, 1, 1}, {em.nx_p_m, em.ny_d_m - 1, em.nz_d_m - 1}),
    KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
      Bx(ix, iy, iz) += -dt_over_dy * (Ez(ix, iy, iz) - Ez(ix, iy - 1, iz)) +
                        dt_over_dz * (Ey(ix, iy, iz) - Ey(ix, iy, iz - 1));
    });

  // Magnetic field By (d,p,d)

  Kokkos::parallel_for(
    mdrange_policy({1, 0, 1}, {em.nx_d_m - 1, em.ny_p_m, em.nz_d_m - 1}),
    KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
      By(ix, iy, iz) += -dt_over_dz * (Ex(ix, iy, iz) - Ex(ix, iy, iz - 1)) +
                        dt_over_dx * (Ez(ix, iy, iz) - Ez(ix - 1, iy, iz));
    });

  // Magnetic field Bz (d,d,p)

  Kokkos::parallel_for(
    mdrange_policy({1, 1, 0}, {em.nx_d_m - 1, em.ny_d_m - 1, em.nz_p_m}),
    KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
      Bz(ix, iy, iz) += -dt_over_dx * (Ey(ix, iy, iz) - Ey(ix - 1, iy, iz)) +
                        dt_over_dy * (Ex(ix, iy, iz) - Ex(ix, iy - 1, iz));
    });

  Kokkos::fence();

} // end solve

// _______________________________________________________________
//
//! \brief Boundaries condition on the global grid
//! \param[in] Params & params - global constant parameters
// _______________________________________________________________
void currentBC(Params &params, ElectroMagn &em) {

  if (params.boundary_condition == "periodic") {

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

    device_field_t Jx = em.Jx_m.data_m.d_view;
    device_field_t Jy = em.Jy_m.data_m.d_view;
    device_field_t Jz = em.Jz_m.data_m.d_view;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__) || defined(__MINIPIC_KOKKOS_VIEWS__)

    device_field_t Jx = em.Jx_m.data_m;
    device_field_t Jy = em.Jy_m.data_m;
    device_field_t Jz = em.Jz_m.data_m;

#endif

    const auto nx_Jx = em.Jx_m.nx();
    const auto ny_Jx = em.Jx_m.ny();
    const auto nz_Jx = em.Jx_m.nz();

    const auto nx_Jy = em.Jy_m.nx();
    const auto ny_Jy = em.Jy_m.ny();
    const auto nz_Jy = em.Jy_m.nz();

    const auto nx_Jz = em.Jz_m.nx();
    const auto ny_Jz = em.Jz_m.ny();
    const auto nz_Jz = em.Jz_m.nz();

    typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>> mdrange_policy;

    // X

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {ny_Jx, nz_Jx}),
      KOKKOS_LAMBDA(const int iy, const int iz) {
        Jx(0, iy, iz) += Jx(nx_Jx - 2, iy, iz);
        Jx(nx_Jx - 2, iy, iz) = Jx(0, iy, iz);

        Jx(1, iy, iz) += Jx(nx_Jx - 1, iy, iz);
        Jx(nx_Jx - 1, iy, iz) = Jx(1, iy, iz);
      });

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {ny_Jy, nz_Jy}),
      KOKKOS_LAMBDA(const int iy, const int iz) {
        Jy(0, iy, iz) += Jy(nx_Jy - 2, iy, iz);
        Jy(nx_Jy - 2, iy, iz) = Jy(0, iy, iz);

        Jy(1, iy, iz) += Jy(nx_Jy - 1, iy, iz);
        Jy(nx_Jy - 1, iy, iz) = Jy(1, iy, iz);
      });

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {ny_Jz, nz_Jz}),
      KOKKOS_LAMBDA(const int iy, const int iz) {
        Jz(0, iy, iz) += Jz(nx_Jz - 2, iy, iz);
        Jz(nx_Jz - 2, iy, iz) = Jz(0, iy, iz);

        Jz(1, iy, iz) += Jz(nx_Jz - 1, iy, iz);
        Jz(nx_Jz - 1, iy, iz) = Jz(1, iy, iz);
      });

    Kokkos::fence();

    // Y

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Jx, nz_Jx}),
      KOKKOS_LAMBDA(const int ix, const int iz) {
        Jx(ix, 0, iz) += Jx(ix, ny_Jx - 2, iz);
        Jx(ix, ny_Jx - 2, iz) = Jx(ix, 0, iz);

        Jx(ix, 1, iz) += Jx(ix, ny_Jx - 1, iz);
        Jx(ix, ny_Jx - 1, iz) = Jx(ix, 1, iz);
      });

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Jy, nz_Jy}),
      KOKKOS_LAMBDA(const int ix, const int iz) {
        Jy(ix, 0, iz) += Jy(ix, ny_Jy - 2, iz);
        Jy(ix, ny_Jy - 2, iz) = Jy(ix, 0, iz);

        Jy(ix, 1, iz) += Jy(ix, ny_Jy - 1, iz);
        Jy(ix, ny_Jy - 1, iz) = Jy(ix, 1, iz);
      });

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Jz, nz_Jz}),
      KOKKOS_LAMBDA(const int ix, const int iz) {
        Jz(ix, 0, iz) += Jz(ix, ny_Jz - 2, iz);
        Jz(ix, ny_Jz - 2, iz) = Jz(ix, 0, iz);

        Jz(ix, 1, iz) += Jz(ix, ny_Jz - 1, iz);
        Jz(ix, ny_Jz - 1, iz) = Jz(ix, 1, iz);
      });

    Kokkos::fence();

    // Z

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Jx, ny_Jx}),
      KOKKOS_LAMBDA(const int ix, const int iy) {
        Jx(ix, iy, 0) += Jx(ix, iy, nz_Jx - 2);
        Jx(ix, iy, nz_Jx - 2) = Jx(ix, iy, 0);

        Jx(ix, iy, 1) += Jx(ix, iy, nz_Jx - 1);
        Jx(ix, iy, nz_Jx - 1) = Jx(ix, iy, 1);
      });

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Jy, ny_Jy}),
      KOKKOS_LAMBDA(const int ix, const int iy) {
        Jy(ix, iy, 0) += Jy(ix, iy, nz_Jy - 2);
        Jy(ix, iy, nz_Jy - 2) = Jy(ix, iy, 0);

        Jy(ix, iy, 1) += Jy(ix, iy, nz_Jy - 1);
        Jy(ix, iy, nz_Jy - 1) = Jy(ix, iy, 1);
      });

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Jz, ny_Jz}),
      KOKKOS_LAMBDA(const int ix, const int iy) {
        Jz(ix, iy, 0) += Jz(ix, iy, nz_Jz - 2);
        Jz(ix, iy, nz_Jz - 2) = Jz(ix, iy, 0);

        Jz(ix, iy, 1) += Jz(ix, iy, nz_Jz - 1);
        Jz(ix, iy, nz_Jz - 1) = Jz(ix, iy, 1);
      });

    Kokkos::fence();

  } // end if periodic
} // end currentBC

// _______________________________________________________________
//
//! \brief Boundaries condition on the global grid
//! \param[in] Params & params - global constant parameters
// _______________________________________________________________
auto solveBC(Params &params, ElectroMagn &em) -> void {

  typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>> mdrange_policy;

  if (params.boundary_condition == "periodic") {

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

    device_field_t Bx = em.Bx_m.data_m.d_view;
    device_field_t By = em.By_m.data_m.d_view;
    device_field_t Bz = em.Bz_m.data_m.d_view;

#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

    device_field_t Bx = em.Bx_m.data_m;
    device_field_t By = em.By_m.data_m;
    device_field_t Bz = em.Bz_m.data_m;

#endif

    const auto nx_Bx = em.Bx_m.nx();
    const auto ny_Bx = em.Bx_m.ny();
    const auto nz_Bx = em.Bx_m.nz();

    const auto nx_By = em.By_m.nx();
    const auto ny_By = em.By_m.ny();
    const auto nz_By = em.By_m.nz();

    const auto nx_Bz = em.Bz_m.nx();
    const auto ny_Bz = em.Bz_m.ny();
    const auto nz_Bz = em.Bz_m.nz();

    // X dim
    // By (d,p,d)

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {ny_By, nz_By}),
      KOKKOS_LAMBDA(const int iy, const int iz) {
        // -X
        By(0, iy, iz)         = By(nx_By - 2, iy, iz);
        By(nx_By - 1, iy, iz) = By(1, iy, iz);
      });

    // Bz (d,d,p)
    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {ny_Bz, nz_Bz}),
      KOKKOS_LAMBDA(const int iy, const int iz) {
        // -X
        Bz(0, iy, iz)         = Bz(nx_Bz - 2, iy, iz);
        Bz(nx_Bz - 1, iy, iz) = Bz(1, iy, iz);
      });

    Kokkos::fence();

    // Y dim
    // Bx (p,d,d)

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Bx, nz_Bx}),
      KOKKOS_LAMBDA(const int ix, const int iz) {
        // -Y
        Bx(ix, 0, iz) = Bx(ix, ny_Bx - 2, iz);
        // +Y
        Bx(ix, ny_Bx - 1, iz) = Bx(ix, 1, iz);
      });

    // Bz (d,d,p)

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Bz, nz_Bz}),
      KOKKOS_LAMBDA(const int ix, const int iz) {
        // -Y
        Bz(ix, 0, iz) = Bz(ix, ny_Bz - 2, iz);
        // +Y
        Bz(ix, ny_Bz - 1, iz) = Bz(ix, 1, iz);
      });

    Kokkos::fence();

    // Z dim
    // Bx

    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Bx, ny_Bx}),
      KOKKOS_LAMBDA(const int ix, const int iy) {
        // -Z
        Bx(ix, iy, 0) = Bx(ix, iy, nz_Bx - 2);
        // +Z
        Bx(ix, iy, nz_Bx - 1) = Bx(ix, iy, 1);
      });

    // By
    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_By, ny_By}),
      KOKKOS_LAMBDA(const int ix, const int iy) {
        // -Z
        By(ix, iy, 0) = By(ix, iy, nz_By - 2);
        // +Z
        By(ix, iy, nz_By - 1) = By(ix, iy, 1);
      });

    Kokkos::fence();

  } else if (params.boundary_condition == "reflective") {

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

    device_field_t Bx = em.Bx_m.data_m.d_view;
    device_field_t By = em.By_m.data_m.d_view;
    device_field_t Bz = em.Bz_m.data_m.d_view;

#elif defined(__MINIPIC_KOKKOS_VIEWS__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

    device_field_t Bx = em.Bx_m.data_m;
    device_field_t By = em.By_m.data_m;
    device_field_t Bz = em.Bz_m.data_m;

#endif

    const auto nx_Bx = em.Bx_m.nx();
    const auto ny_Bx = em.Bx_m.ny();
    const auto nz_Bx = em.Bx_m.nz();

    const auto nx_By = em.By_m.nx();
    const auto ny_By = em.By_m.ny();
    const auto nz_By = em.By_m.nz();

    const auto nx_Bz = em.Bz_m.nx();
    const auto ny_Bz = em.Bz_m.ny();
    const auto nz_Bz = em.Bz_m.nz();

    // X dim
    // By (d,p,d)
    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {ny_By, nz_By}),
      KOKKOS_LAMBDA(const int iy, const int iz) {
        // -X
        By(0, iy, iz) = By(1, iy, iz);
        // +X
        By(nx_By - 1, iy, iz) = By(nx_By - 2, iy, iz);
      });

    // Bz (d,d,p)
    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {ny_Bz, nz_Bz}),
      KOKKOS_LAMBDA(const int iy, const int iz) {
        // -X
        Bz(0, iy, iz) = Bz(1, iy, iz);
        // +X
        Bz(nx_Bz - 1, iy, iz) = Bz(nx_Bz - 2, iy, iz);
      });

    // Y dim
    // Bx (p,d,d)
    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Bx, nz_Bx}),
      KOKKOS_LAMBDA(const int ix, const int iz) {
        // -Y
        Bx(ix, 0, iz) = Bx(ix, 1, iz);
        // +Y
        Bx(ix, ny_Bx - 1, iz) = Bx(ix, ny_Bx - 2, iz);
      });

    Kokkos::fence();

    // Bz (-1 to avoid corner)
    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Bz, nz_Bz}),
      KOKKOS_LAMBDA(const int ix, const int iz) {
        // -Y
        Bz(ix, 0, iz) = Bz(ix, 1, iz);
        // +Y
        Bz(ix, ny_Bz - 1, iz) = Bz(ix, ny_Bz - 2, iz);
      });

    // Z dim
    // Bx
    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_Bx, ny_Bx}),
      KOKKOS_LAMBDA(const int ix, const int iy) {
        // -Z
        Bx(ix, iy, 0) = Bx(ix, iy, 1);
        // +Z
        Bx(ix, iy, nz_Bx - 1) = Bx(ix, iy, nz_Bx - 2);
      });

    // By
    Kokkos::parallel_for(
      mdrange_policy({0, 0}, {nx_By, ny_By}),
      KOKKOS_LAMBDA(const int ix, const int iy) {
        // -Z
        By(ix, iy, 0) = By(ix, iy, 1);
        // +Z
        By(ix, iy, nz_By - 1) = By(ix, iy, nz_By - 2);
      });

    Kokkos::fence();
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

  const double patch_inf[3] = {patch.inf_m[0], patch.inf_m[1], patch.inf_m[2]};
  const double patch_sup[3] = {patch.sup_m[0], patch.sup_m[1], patch.sup_m[2]};

  const double inf[3]    = {params.inf_x, params.inf_y, params.inf_z};
  const double sup[3]    = {params.sup_x, params.sup_y, params.sup_z};
  const double length[3] = {params.Lx, params.Ly, params.Lz};

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

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

    device_vector_t w = patch.particles_m[is].weight_.data_.d_view;

    device_vector_t x = patch.particles_m[is].x_.data_.d_view;
    device_vector_t y = patch.particles_m[is].y_.data_.d_view;
    device_vector_t z = patch.particles_m[is].z_.data_.d_view;

    device_vector_t mx = patch.particles_m[is].mx_.data_.d_view;
    device_vector_t my = patch.particles_m[is].my_.data_.d_view;
    device_vector_t mz = patch.particles_m[is].mz_.data_.d_view;

    Kokkos::View<size_t *, Kokkos::MemoryTraits<Kokkos::Atomic>> n_particle_to_move_accessor =
      n_particle_to_move.data_.d_view;

    Kokkos::DualView<int *>::t_dev masks_accessor = masks.data_.d_view;

    Kokkos::DualView<size_t *>::t_dev ip_to_move_accessor = ip_to_move.data_.d_view;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

    device_vector_t w = patch.particles_m[is].weight_.data_;

    device_vector_t x = patch.particles_m[is].x_.data_;
    device_vector_t y = patch.particles_m[is].y_.data_;
    device_vector_t z = patch.particles_m[is].z_.data_;

    device_vector_t mx = patch.particles_m[is].mx_.data_;
    device_vector_t my = patch.particles_m[is].my_.data_;
    device_vector_t mz = patch.particles_m[is].mz_.data_;

    Kokkos::View<size_t *, Kokkos::SharedSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>
      n_particle_to_move_accessor = n_particle_to_move.data_;

    Kokkos::View<int *, Kokkos::SharedSpace> masks_accessor = masks.data_;

    Kokkos::View<size_t *, Kokkos::SharedSpace> ip_to_move_accessor = ip_to_move.data_;

#elif defined(__MINIPIC_KOKKOS_VIEWS__)

    device_vector_t w = patch.particles_m[is].weight_.data_;

    device_vector_t x = patch.particles_m[is].x_.data_;
    device_vector_t y = patch.particles_m[is].y_.data_;
    device_vector_t z = patch.particles_m[is].z_.data_;

    device_vector_t mx = patch.particles_m[is].mx_.data_;
    device_vector_t my = patch.particles_m[is].my_.data_;
    device_vector_t mz = patch.particles_m[is].mz_.data_;

    Kokkos::View<size_t *, Kokkos::MemoryTraits<Kokkos::Atomic>> n_particle_to_move_accessor = n_particle_to_move.data_;

    Kokkos::View<int *> masks_accessor = masks.data_;

    Kokkos::View<size_t *> ip_to_move_accessor = ip_to_move.data_;

#endif

    // 1 - Compute number of particles to move per buffer and tag them

    Kokkos::parallel_for(
      n_particles,
      KOKKOS_LAMBDA(const size_t ip) {
        double shift[3];

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
          int ib = (shift[0] + 1) * 9 + (shift[1] + 1) * 3 + (shift[2] + 1);
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

          n_particle_to_move_accessor(ib) += 1;

          // we store here the buffer id to use it later
          masks_accessor(ip) = ib;
        }
      } // end for particles
    );

    Kokkos::fence();

    // masks.sync(minipic::device, minipic::host);
    n_particle_to_move.sync(minipic::device, minipic::host);

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

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)
      Kokkos::View<double *> b_w[26];

      Kokkos::View<double *> b_x[26];
      Kokkos::View<double *> b_y[26];
      Kokkos::View<double *> b_z[26];

      Kokkos::View<double *> b_mx[26];
      Kokkos::View<double *> b_my[26];
      Kokkos::View<double *> b_mz[26];

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

      Kokkos::View<double *, Kokkos::SharedSpace> b_w[26];

      Kokkos::View<double *, Kokkos::SharedSpace> b_x[26];
      Kokkos::View<double *, Kokkos::SharedSpace> b_y[26];
      Kokkos::View<double *, Kokkos::SharedSpace> b_z[26];

      Kokkos::View<double *, Kokkos::SharedSpace> b_mx[26];
      Kokkos::View<double *, Kokkos::SharedSpace> b_my[26];
      Kokkos::View<double *, Kokkos::SharedSpace> b_mz[26];

#elif defined(__MINIPIC_KOKKOS_VIEWS__)

      Kokkos::View<double *> b_mz[26];
      Kokkos::View<double *> b_my[26];
      Kokkos::View<double *> b_mx[26];

      Kokkos::View<double *> b_z[26];
      Kokkos::View<double *> b_y[26];
      Kokkos::View<double *> b_x[26];

      Kokkos::View<double *> b_w[26];

#endif

      int buffer_size[26];

      for (int ib = 0; ib < 26; ib++) {

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

        b_x[ib] = patch.particles_to_move_m[is][ib].x_.data_.d_view;
        b_y[ib] = patch.particles_to_move_m[is][ib].y_.data_.d_view;
        b_z[ib] = patch.particles_to_move_m[is][ib].z_.data_.d_view;

        b_mx[ib] = patch.particles_to_move_m[is][ib].mx_.data_.d_view;
        b_my[ib] = patch.particles_to_move_m[is][ib].my_.data_.d_view;
        b_mz[ib] = patch.particles_to_move_m[is][ib].mz_.data_.d_view;

        b_w[ib] = patch.particles_to_move_m[is][ib].weight_.data_.d_view;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__) || defined(__MINIPIC_KOKKOS_VIEWS__)

        b_x[ib] = patch.particles_to_move_m[is][ib].x_.data_;
        b_y[ib] = patch.particles_to_move_m[is][ib].y_.data_;
        b_z[ib] = patch.particles_to_move_m[is][ib].z_.data_;

        b_mx[ib] = patch.particles_to_move_m[is][ib].mx_.data_;
        b_my[ib] = patch.particles_to_move_m[is][ib].my_.data_;
        b_mz[ib] = patch.particles_to_move_m[is][ib].mz_.data_;

        b_w[ib] = patch.particles_to_move_m[is][ib].weight_.data_;

#endif

        buffer_size[ib] = patch.particles_to_move_m[is][ib].size();
      }

      Kokkos::parallel_for(
        26,
        KOKKOS_LAMBDA(const int ib) {
          if (buffer_size[ib] > 0) {

            for (size_t ip = 0; ip < n_particles; ++ip) {

              if (masks_accessor(ip) == ib) {

                const size_t i = ip_to_move_accessor(ib);

                b_x[ib](i) = x(ip);
                b_y[ib](i) = y(ip);
                b_z[ib](i) = z(ip);

                b_mx[ib](i) = mx(ip);
                b_my[ib](i) = my(ip);
                b_mz[ib](i) = mz(ip);

                b_w[ib](i) = w(ip);

                ip_to_move_accessor(ib)++;
                //}
              }

            } // end particle loop
          }
        });
      Kokkos::fence();

    } // end if total_particles_to_remove

    // 4 - Move particles to remove at the end of the vector

    if (total_particles_to_remove > 0) {

      Kokkos::parallel_for(
        1,
        KOKKOS_LAMBDA(const int i) {
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
        });
      Kokkos::fence();

      // Delete tagged particles by resizing particles_m[is]
      patch.particles_m[is].resize(n_particles - total_particles_to_remove, minipic::device);

    } // end if total_particles_to_remove > 0

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

  // Get the current patch
  Patch &patch = vec_patch[id_patch];

  for (int is = 0; is < patch.n_species_m; is++) {

    const size_t number_of_particles = patch.particles_m[is].size();

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

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

            device_vector_t w = patch.particles_m[is].weight_.data_.d_view;

            device_vector_t x = patch.particles_m[is].x_.data_.d_view;
            device_vector_t y = patch.particles_m[is].y_.data_.d_view;
            device_vector_t z = patch.particles_m[is].z_.data_.d_view;

            device_vector_t mx = patch.particles_m[is].mx_.data_.d_view;
            device_vector_t my = patch.particles_m[is].my_.data_.d_view;
            device_vector_t mz = patch.particles_m[is].mz_.data_.d_view;

            device_vector_t b_w =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].weight_.data_.d_view;

            device_vector_t b_x =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].x_.data_.d_view;
            device_vector_t b_y =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].y_.data_.d_view;
            device_vector_t b_z =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].z_.data_.d_view;

            device_vector_t b_mx =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].mx_.data_.d_view;
            device_vector_t b_my =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].my_.data_.d_view;
            device_vector_t b_mz =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].mz_.data_.d_view;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__) || defined(__MINIPIC_KOKKOS_VIEWS__)

            device_vector_t w = patch.particles_m[is].weight_.data_;

            device_vector_t x = patch.particles_m[is].x_.data_;
            device_vector_t y = patch.particles_m[is].y_.data_;
            device_vector_t z = patch.particles_m[is].z_.data_;

            device_vector_t mx = patch.particles_m[is].mx_.data_;
            device_vector_t my = patch.particles_m[is].my_.data_;
            device_vector_t mz = patch.particles_m[is].mz_.data_;

            device_vector_t b_w =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].weight_.data_;

            device_vector_t b_x =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].x_.data_;
            device_vector_t b_y =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].y_.data_;
            device_vector_t b_z =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].z_.data_;

            device_vector_t b_mx =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].mx_.data_;
            device_vector_t b_my =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].my_.data_;
            device_vector_t b_mz =
              vec_patch[idx_neighbor].particles_to_move_m[is][idx_buffer].mz_.data_;

#endif

            Kokkos::parallel_for(
              buffer_size,
              KOKKOS_LAMBDA(const size_t ip) {
                x(ip_buffer_start + ip) = b_x(ip);
                y(ip_buffer_start + ip) = b_y(ip);
                z(ip_buffer_start + ip) = b_z(ip);

                mx(ip_buffer_start + ip) = b_mx(ip);
                my(ip_buffer_start + ip) = b_my(ip);
                mz(ip_buffer_start + ip) = b_mz(ip);

                w(ip_buffer_start + ip) = b_w(ip);

                // Ex.h(last_ip) = buffer.Ex_h(ip);
                // Ey.h(last_ip) = buffer.Ey_h(ip);
                // Ez.h(last_ip) = buffer.Ez_h(ip);

                // Bx.h(last_ip) = buffer.Bx_h(ip);
                // By.h(last_ip) = buffer.By_h(ip);
                // Bz/h(last_ip) = buffer.Bz_h(ip);

                // } // end ip loop
              });
            Kokkos::fence();

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

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)
      device_field_t Jx_0  = patch.vec_Jx_m[0].data_m.d_view;
      device_field_t Jx_is = patch.vec_Jx_m[is].data_m.d_view;

      device_field_t Jy_0  = patch.vec_Jy_m[0].data_m.d_view;
      device_field_t Jy_is = patch.vec_Jy_m[is].data_m.d_view;

      device_field_t Jz_0  = patch.vec_Jz_m[0].data_m.d_view;
      device_field_t Jz_is = patch.vec_Jz_m[is].data_m.d_view;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__) || defined(__MINIPIC_KOKKOS_VIEWS__)
      device_field_t Jx_0  = patch.vec_Jx_m[0].data_m;
      device_field_t Jx_is = patch.vec_Jx_m[is].data_m;

      device_field_t Jy_0  = patch.vec_Jy_m[0].data_m;
      device_field_t Jy_is = patch.vec_Jy_m[is].data_m;

      device_field_t Jz_0  = patch.vec_Jz_m[0].data_m;
      device_field_t Jz_is = patch.vec_Jz_m[is].data_m;

#endif

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
      Kokkos::parallel_for(
        mdrange_policy({0, 0, 0},
                       {patch.vec_Jx_m[is].nx(), patch.vec_Jx_m[is].ny(), patch.vec_Jx_m[is].nz()}),
        KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
          Jx_0(ix, iy, iz) += Jx_is(ix, iy, iz);
        }

      );

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
      Kokkos::parallel_for(
        mdrange_policy({0, 0, 0},
                       {patch.vec_Jy_m[is].nx(), patch.vec_Jy_m[is].ny(), patch.vec_Jy_m[is].nz()}),
        KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
          Jy_0(ix, iy, iz) += Jy_is(ix, iy, iz);
        });

      typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
      Kokkos::parallel_for(
        mdrange_policy({0, 0, 0},
                       {patch.vec_Jz_m[is].nx(), patch.vec_Jz_m[is].ny(), patch.vec_Jz_m[is].nz()}),
        KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
          Jz_0(ix, iy, iz) += Jz_is(ix, iy, iz);
        });
      Kokkos::fence();

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

#if defined(__MINIPIC_KOKKOS_DUALVIEW_COMMON__)

    device_field_t Jx_0 = patch.vec_Jx_m[0].data_m.d_view;
    device_field_t Jx   = em.Jx_m.data_m.d_view;

    device_field_t Jy_0 = patch.vec_Jy_m[0].data_m.d_view;
    device_field_t Jy   = em.Jy_m.data_m.d_view;

    device_field_t Jz_0 = patch.vec_Jz_m[0].data_m.d_view;
    device_field_t Jz   = em.Jz_m.data_m.d_view;

#elif defined(__MINIPIC_KOKKOS_UNIFIED__) || defined(__MINIPIC_KOKKOS_VIEWS__)

    device_field_t Jx_0 = patch.vec_Jx_m[0].data_m;
    device_field_t Jx   = em.Jx_m.data_m;

    device_field_t Jy_0 = patch.vec_Jy_m[0].data_m;
    device_field_t Jy   = em.Jy_m.data_m;

    device_field_t Jz_0 = patch.vec_Jz_m[0].data_m;
    device_field_t Jz   = em.Jz_m.data_m;

#endif

    typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
    Kokkos::parallel_for(
      mdrange_policy({0, 0, 0},
                     {patch.vec_Jx_m[0].nx(), patch.vec_Jx_m[0].ny(), patch.vec_Jx_m[0].nz()}),
      KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
        Jx(i_global_d + ix, j_global_p + iy, k_global_p + iz) += Jx_0(ix, iy, iz);
      });

    typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
    Kokkos::parallel_for(
      mdrange_policy({0, 0, 0},
                     {patch.vec_Jy_m[0].nx(), patch.vec_Jy_m[0].ny(), patch.vec_Jy_m[0].nz()}),
      KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
        Jy(i_global_p + ix, j_global_d + iy, k_global_p + iz) += Jy_0(ix, iy, iz);
      });

    typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
    Kokkos::parallel_for(
      mdrange_policy({0, 0, 0},
                     {patch.vec_Jz_m[0].nx(), patch.vec_Jz_m[0].ny(), patch.vec_Jz_m[0].nz()}),
      KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
        Jz(i_global_p + ix, j_global_p + iy, k_global_d + iz) += Jz_0(ix, iy, iz);
      });
    Kokkos::fence();

  } // end if total_particles
}

// ____________________________________________________________________________
//! \brief Emit a laser field in the x direction using an antenna
//! \param[in] Params & params - global constant parameters
//! \param[in] profile - (std::function<double(double y, double z, double t)>) profile of the
//! antenna \param[in] x - (double) position of the antenna \param[in] double t - (double) current
//! time
// ____________________________________________________________________________
auto antenna(Params &params,
             ElectroMagn &em,
             std::function<double(double, double, double)> profile,
             double x,
             double t) -> void {

  em.Jz_m.sync(minipic::device, minipic::host);

  Field<mini_float> *J = &em.Jz_m;

  const int ix = Kokkos::floor((x - params.inf_x - J->dual_x_m * 0.5 * params.dx) / params.dx);

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

} // end antenna

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
                        std::function<double(double, double, double, double)> func_weight) -> void {

} // end function

} // end namespace operators

#endif // OPERATORS_H