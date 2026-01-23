/* _____________________________________________________________________ */
//! \file Operators.hpp

//! \brief contains generic kernels for the particle pusher

/* _____________________________________________________________________ */

#ifndef OPERATORS_H
#define OPERATORS_H

#include "ElectroMagn.hpp"
#include "Patch.hpp"
#include "Profiler.hpp"

namespace operators {

// ______________________________________________________________________________
//
//! \brief Interpolation operator at the patch level :
//! interpolate EM fields from global grid for each particle
//! \param[in] em  global electromagnetic fields
//! \param[in] patch  patch data structure
// ______________________________________________________________________________
auto interpolate(Params &params, ElectroMagn &em, Patch &patch) -> void {

  const mini_float inv_dx_m = em.inv_dx_m;
  const mini_float inv_dy_m = em.inv_dy_m;
  const mini_float inv_dz_m = em.inv_dz_m;

  const mini_float xmin = params.inf_x;
  const mini_float ymin = params.inf_y;
  const mini_float zmin = params.inf_z;

  const auto ny_Ex = em.Ex_m.ny_m, nz_Ex = em.Ex_m.nz_m;
  const auto ny_Ey = em.Ey_m.ny_m, nz_Ey = em.Ey_m.nz_m;
  const auto ny_Ez = em.Ez_m.ny_m, nz_Ez = em.Ez_m.nz_m;

  const auto ny_Bx = em.Bx_m.ny_m, nz_Bx = em.Bx_m.nz_m;
  const auto ny_By = em.By_m.ny_m, nz_By = em.By_m.nz_m;
  const auto ny_Bz = em.Bz_m.ny_m, nz_Bz = em.Bz_m.nz_m;

  const auto nynz_Ex = nz_Ex * ny_Ex;
  const auto nynz_Ey = nz_Ey * ny_Ey;
  const auto nynz_Ez = nz_Ez * ny_Ez;

  const auto nynz_Bx = nz_Bx * ny_Bx;
  const auto nynz_By = nz_By * ny_By;
  const auto nynz_Bz = nz_Bz * ny_Bz;

  const mini_float *const __restrict__ Ex = em.Ex_m.get_raw_pointer(minipic::device);
  const mini_float *const __restrict__ Ey = em.Ey_m.get_raw_pointer(minipic::device);
  const mini_float *const __restrict__ Ez = em.Ez_m.get_raw_pointer(minipic::device);

  const mini_float *const __restrict__ Bx = em.Bx_m.get_raw_pointer(minipic::device);
  const mini_float *const __restrict__ By = em.By_m.get_raw_pointer(minipic::device);
  const mini_float *const __restrict__ Bz = em.Bz_m.get_raw_pointer(minipic::device);

  for (int is = 0; is < patch.n_species_m; is++) {

#if defined(__MINIPIC_THRUST_COUNTING__)

    const mini_float *const __restrict__ part_x =
      patch.particles_m[is].x_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ part_y =
      patch.particles_m[is].y_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ part_z =
      patch.particles_m[is].z_.get_raw_pointer(minipic::device);

    mini_float *const __restrict__ part_Ex =
      patch.particles_m[is].Ex_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ part_Ey =
      patch.particles_m[is].Ey_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ part_Ez =
      patch.particles_m[is].Ez_.get_raw_pointer(minipic::device);

    mini_float *const __restrict__ part_Bx =
      patch.particles_m[is].Bx_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ part_By =
      patch.particles_m[is].By_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ part_Bz =
      patch.particles_m[is].Bz_.get_raw_pointer(minipic::device);

#endif

    const size_t n_particles = patch.particles_m[is].size();

#if defined(__MINIPIC_THRUST_ZIP__)

    // Use thrust for each with zip iterator and lambda function
    thrust::for_each(
      thrust::make_zip_iterator(thrust::make_tuple(patch.particles_m[is].x_.device_data_.begin(),
                                                   patch.particles_m[is].y_.device_data_.begin(),
                                                   patch.particles_m[is].z_.device_data_.begin(),
                                                   patch.particles_m[is].Ex_.device_data_.begin(),
                                                   patch.particles_m[is].Ey_.device_data_.begin(),
                                                   patch.particles_m[is].Ez_.device_data_.begin(),
                                                   patch.particles_m[is].Bx_.device_data_.begin(),
                                                   patch.particles_m[is].By_.device_data_.begin(),
                                                   patch.particles_m[is].Bz_.device_data_.begin())),
      thrust::make_zip_iterator(
        thrust::make_tuple(patch.particles_m[is].x_.device_data_.begin() + n_particles,
                           patch.particles_m[is].y_.device_data_.begin() + n_particles,
                           patch.particles_m[is].z_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Ex_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Ey_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Ez_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Bx_.device_data_.begin() + n_particles,
                           patch.particles_m[is].By_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Bz_.device_data_.begin() + n_particles)),
      [=] __device__(thrust::tuple<mini_float &,
                                   mini_float &,
                                   mini_float &, // x, y, z
                                   mini_float &,
                                   mini_float &,
                                   mini_float &, // Ex, Ey, Ez
                                   mini_float &,
                                   mini_float &,
                                   mini_float & // Bx, By, Bz
                                   > data) {
        mini_float &x = thrust::get<0>(data);
        mini_float &y = thrust::get<1>(data);
        mini_float &z = thrust::get<2>(data);

        mini_float &Exp = thrust::get<3>(data);
        mini_float &Eyp = thrust::get<4>(data);
        mini_float &Ezp = thrust::get<5>(data);

        mini_float &Bxp = thrust::get<6>(data);
        mini_float &Byp = thrust::get<7>(data);
        mini_float &Bzp = thrust::get<8>(data);

        // Calculate normalized positions
        const mini_float ixn = (x - xmin) * inv_dx_m;
        const mini_float iyn = (y - ymin) * inv_dy_m;
        const mini_float izn = (z - zmin) * inv_dz_m;

#elif defined(__MINIPIC_THRUST_COUNTING__)

    // Use thrust for each with lambda function
    thrust::for_each(
      thrust::counting_iterator<size_t>(0),
      thrust::counting_iterator<size_t>(n_particles),
      [=] __device__(size_t ip) {
        // Calculate normalized positions
        const mini_float ixn = (part_x[ip] - xmin) * inv_dx_m;
        const mini_float iyn = (part_y[ip] - ymin) * inv_dy_m;
        const mini_float izn = (part_z[ip] - zmin) * inv_dz_m;

#endif

        // Compute indexes in global primal grid
        const unsigned int ixp = static_cast<unsigned int>(floor(ixn));
        const unsigned int iyp = static_cast<unsigned int>(floor(iyn));
        const unsigned int izp = static_cast<unsigned int>(floor(izn));

        // Compute indexes in global dual grid
        const unsigned int ixd = static_cast<unsigned int>(floor(ixn + 0.5));
        const unsigned int iyd = static_cast<unsigned int>(floor(iyn + 0.5));
        const unsigned int izd = static_cast<unsigned int>(floor(izn + 0.5));

        // Compute distances

        const mini_float dist_x_p = ixn - static_cast<mini_float>(ixp);
        const mini_float dist_y_p = iyn - static_cast<mini_float>(iyp);
        const mini_float dist_z_p = izn - static_cast<mini_float>(izp);

        const mini_float dist_x_d = (ixn + 0.5) - static_cast<mini_float>(ixd);
        const mini_float dist_y_d = (iyn + 0.5) - static_cast<mini_float>(iyd);
        const mini_float dist_z_d = (izn + 0.5) - static_cast<mini_float>(izd);

        // interpolation electric field
        // Ex (d, p , p)
        const auto v00 = Ex[ixd * (nynz_Ex) + iyp * (nz_Ex) + izp] * (1 - dist_x_d) +
                         Ex[(ixd + 1) * (nynz_Ex) + iyp * (nz_Ex) + izp] * dist_x_d;
        const auto v01 = Ex[ixd * (nynz_Ex) + iyp * (nz_Ex) + (izp + 1)] * (1 - dist_x_d) +
                         Ex[(ixd + 1) * (nynz_Ex) + iyp * (nz_Ex) + (izp + 1)] * dist_x_d;
        const auto v10 = Ex[ixd * (nynz_Ex) + (iyp + 1) * (nz_Ex) + izp] * (1 - dist_x_d) +
                         Ex[(ixd + 1) * (nynz_Ex) + (iyp + 1) * (nz_Ex) + izp] * dist_x_d;
        const auto v11 = Ex[ixd * (nynz_Ex) + (iyp + 1) * (nz_Ex) + (izp + 1)] * (1 - dist_x_d) +
                         Ex[(ixd + 1) * (nynz_Ex) + (iyp + 1) * (nz_Ex) + (izp + 1)] * dist_x_d;

#if defined(__MINIPIC_THRUST_ZIP__)
        Exp =
#elif defined(__MINIPIC_THRUST_COUNTING__)
        part_Ex[ip] =
#endif
          (v00 * (1 - dist_y_p) + v10 * dist_y_p) * (1 - dist_z_p) +
          (v01 * (1 - dist_y_p) + v11 * dist_y_p) * dist_z_p;

        // Ey (p, d, p)
        {
          const auto v00 = Ey[ixp * (nynz_Ey) + iyd * (nz_Ey) + izp] * (1 - dist_x_p) +
                           Ey[(ixp + 1) * (nynz_Ey) + iyd * (nz_Ey) + izp] * dist_x_p;
          const auto v01 = Ey[ixp * (nynz_Ey) + iyd * (nz_Ey) + (izp + 1)] * (1 - dist_x_p) +
                           Ey[(ixp + 1) * (nynz_Ey) + iyd * (nz_Ey) + (izp + 1)] * dist_x_p;
          const auto v10 = Ey[ixp * (nynz_Ey) + (iyd + 1) * (nz_Ey) + izp] * (1 - dist_x_p) +
                           Ey[(ixp + 1) * (nynz_Ey) + (iyd + 1) * (nz_Ey) + izp] * dist_x_p;
          const auto v11 = Ey[ixp * (nynz_Ey) + (iyd + 1) * (nz_Ey) + (izp + 1)] * (1 - dist_x_p) +
                           Ey[(ixp + 1) * (nynz_Ey) + (iyd + 1) * (nz_Ey) + (izp + 1)] * dist_x_p;
          const auto v0 = v00 * (1 - dist_y_d) + v10 * dist_y_d;
          const auto v1 = v01 * (1 - dist_y_d) + v11 * dist_y_d;

#if defined(__MINIPIC_THRUST_ZIP__)
          Eyp =
#elif defined(__MINIPIC_THRUST_COUNTING__)
          part_Ey[ip] =
#endif
            v0 * (1 - dist_z_p) + v1 * dist_z_p;
        }

        // Ez (p, p, d)
        {
          const auto v00 = Ez[ixp * (nynz_Ez) + iyp * (nz_Ez) + izd] * (1 - dist_x_p) +
                           Ez[(ixp + 1) * (nynz_Ez) + iyp * (nz_Ez) + izd] * dist_x_p;
          const auto v01 = Ez[ixp * (nynz_Ez) + iyp * (nz_Ez) + (izd + 1)] * (1 - dist_x_p) +
                           Ez[(ixp + 1) * (nynz_Ez) + iyp * (nz_Ez) + (izd + 1)] * dist_x_p;
          const auto v10 = Ez[ixp * (nynz_Ez) + (iyp + 1) * (nz_Ez) + izd] * (1 - dist_x_p) +
                           Ez[(ixp + 1) * (nynz_Ez) + (iyp + 1) * (nz_Ez) + izd] * dist_x_p;
          const auto v11 = Ez[ixp * (nynz_Ez) + (iyp + 1) * (nz_Ez) + (izd + 1)] * (1 - dist_x_p) +
                           Ez[(ixp + 1) * (nynz_Ez) + (iyp + 1) * (nz_Ez) + (izd + 1)] * dist_x_p;
          const auto v0 = v00 * (1 - dist_y_p) + v10 * dist_y_p;
          const auto v1 = v01 * (1 - dist_y_p) + v11 * dist_y_p;
#if defined(__MINIPIC_THRUST_ZIP__)
          Ezp =
#elif defined(__MINIPIC_THRUST_COUNTING__)
          part_Ez[ip] =
#endif
            v0 * (1 - dist_z_d) + v1 * dist_z_d;
        }

        // interpolation magnetic field
        // Bx (p, d, d)
        {

          const auto v00 = Bx[ixp * (nynz_Bx) + iyd * (nz_Bx) + izd] * (1 - dist_x_p) +
                           Bx[(ixp + 1) * (nynz_Bx) + iyd * (nz_Bx) + izd] * dist_x_p;
          const auto v01 = Bx[ixp * (nynz_Bx) + iyd * (nz_Bx) + (izd + 1)] * (1 - dist_x_p) +
                           Bx[(ixp + 1) * (nynz_Bx) + iyd * (nz_Bx) + (izd + 1)] * dist_x_p;
          const auto v10 = Bx[ixp * (nynz_Bx) + (iyd + 1) * (nz_Bx) + izd] * (1 - dist_x_p) +
                           Bx[(ixp + 1) * (nynz_Bx) + (iyd + 1) * (nz_Bx) + izd] * dist_x_p;
          const auto v11 = Bx[ixp * (nynz_Bx) + (iyd + 1) * (nz_Bx) + (izd + 1)] * (1 - dist_x_p) +
                           Bx[(ixp + 1) * (nynz_Bx) + (iyd + 1) * (nz_Bx) + (izd + 1)] * dist_x_p;
          const auto v0 = v00 * (1 - dist_y_d) + v10 * dist_y_d;
          const auto v1 = v01 * (1 - dist_y_d) + v11 * dist_y_d;
#if defined(__MINIPIC_THRUST_ZIP__)
          Bxp =
#elif defined(__MINIPIC_THRUST_COUNTING__)
          part_Bx[ip] =
#endif
            v0 * (1 - dist_z_d) + v1 * dist_z_d;
        }

        // By (d, p, d)
        {
          const auto v00 = By[ixd * (nynz_By) + iyp * (nz_By) + izd] * (1 - dist_x_d) +
                           By[(ixd + 1) * (nynz_By) + iyp * (nz_By) + izd] * dist_x_d;
          const auto v01 = By[ixd * (nynz_By) + iyp * (nz_By) + (izd + 1)] * (1 - dist_x_d) +
                           By[(ixd + 1) * (nynz_By) + iyp * (nz_By) + (izd + 1)] * dist_x_d;
          const auto v10 = By[ixd * (nynz_By) + (iyp + 1) * (nz_By) + izd] * (1 - dist_x_d) +
                           By[(ixd + 1) * (nynz_By) + (iyp + 1) * (nz_By) + izd] * dist_x_d;
          const auto v11 = By[ixd * (nynz_By) + (iyp + 1) * (nz_By) + (izd + 1)] * (1 - dist_x_d) +
                           By[(ixd + 1) * (nynz_By) + (iyp + 1) * (nz_By) + (izd + 1)] * dist_x_d;
          const auto v0 = v00 * (1 - dist_y_p) + v10 * dist_y_p;
          const auto v1 = v01 * (1 - dist_y_p) + v11 * dist_y_p;
#if defined(__MINIPIC_THRUST_ZIP__)
          Byp =
#elif defined(__MINIPIC_THRUST_COUNTING__)
          part_By[ip] =
#endif
            v0 * (1 - dist_z_d) + v1 * dist_z_d;
        }

        // Bz (d, d, p)
        {

          const auto v00 = Bz[ixd * (nynz_Bz) + iyd * (nz_Bz) + izp] * (1 - dist_x_d) +
                           Bz[(ixd + 1) * (nynz_Bz) + iyd * (nz_Bz) + izp] * dist_x_d;
          const auto v01 = Bz[ixd * (nynz_Bz) + iyd * (nz_Bz) + (izp + 1)] * (1 - dist_x_d) +
                           Bz[(ixd + 1) * (nynz_Bz) + iyd * (nz_Bz) + (izp + 1)] * dist_x_d;
          const auto v10 = Bz[ixd * (nynz_Bz) + (iyd + 1) * (nz_Bz) + izp] * (1 - dist_x_d) +
                           Bz[(ixd + 1) * (nynz_Bz) + (iyd + 1) * (nz_Bz) + izp] * dist_x_d;
          const auto v11 = Bz[ixd * (nynz_Bz) + (iyd + 1) * (nz_Bz) + (izp + 1)] * (1 - dist_x_d) +
                           Bz[(ixd + 1) * (nynz_Bz) + (iyd + 1) * (nz_Bz) + (izp + 1)] * dist_x_d;
          const auto v0 = v00 * (1 - dist_y_d) + v10 * dist_y_d;
          const auto v1 = v01 * (1 - dist_y_d) + v11 * dist_y_d;
#if defined(__MINIPIC_THRUST_ZIP__)
          Bzp =
#elif defined(__MINIPIC_THRUST_COUNTING__)
          part_Bz[ip] =
#endif
            v0 * (1 - dist_z_p) + v1 * dist_z_p;
        }
      });

  } // Species loop
}

// ______________________________________________________________________________
//
//! \brief Lambda function to push the particle momentum
// ______________________________________________________________________________

struct momentum_pusher_lambda {
  const mini_float dt;
  const mini_float qp;

#if defined(__MINIPIC_THRUST_ZIP__)

  momentum_pusher_lambda(mini_float dt, mini_float qp) : dt(dt), qp(qp) {}

  __host__ __device__ void operator()(thrust::tuple<mini_float &,
                                                    mini_float &,
                                                    mini_float &, // mx, my, mz
                                                    mini_float &,
                                                    mini_float &,
                                                    mini_float &, // Ex, Ey, Ez
                                                    mini_float &,
                                                    mini_float &,
                                                    mini_float & // Bx, By, Bz
                                                    > data) const {

    mini_float &mx = thrust::get<0>(data);
    mini_float &my = thrust::get<1>(data);
    mini_float &mz = thrust::get<2>(data);

    mini_float &Ex = thrust::get<3>(data);
    mini_float &Ey = thrust::get<4>(data);
    mini_float &Ez = thrust::get<5>(data);

    mini_float &Bx = thrust::get<6>(data);
    mini_float &By = thrust::get<7>(data);
    mini_float &Bz = thrust::get<8>(data);

    // 1/2 E
    mini_float px = qp * Ex;
    mini_float py = qp * Ey;
    mini_float pz = qp * Ez;

    const mini_float ux = mx + px;
    const mini_float uy = my + py;
    const mini_float uz = mz + pz;

    // gamma-factor
    mini_float gamma_inv = qp / sqrt(1 + (ux * ux + uy * uy + uz * uz));

    // B, T = Transform to rotate the particle
    const mini_float tx  = gamma_inv * Bx;
    const mini_float ty  = gamma_inv * By;
    const mini_float tz  = gamma_inv * Bz;
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

    // Update inverse gamma factor
    gamma_inv = 1 / sqrt(1 + (px * px + py * py + pz * pz));

    // Update momentum
    mx = px;
    my = py;
    mz = pz;
  }

#elif defined(__MINIPIC_THRUST_COUNTING__)

  const mini_float *part_Ex;
  const mini_float *part_Ey;
  const mini_float *part_Ez;
  const mini_float *part_Bx;
  const mini_float *part_By;
  const mini_float *part_Bz;
  mini_float *part_mx;
  mini_float *part_my;
  mini_float *part_mz;

  momentum_pusher_lambda(mini_float dt,
                         mini_float qp,
                         const mini_float *part_Ex,
                         const mini_float *part_Ey,
                         const mini_float *part_Ez,
                         const mini_float *part_Bx,
                         const mini_float *part_By,
                         const mini_float *part_Bz,
                         mini_float *part_mx,
                         mini_float *part_my,
                         mini_float *part_mz)
    : dt(dt), qp(qp), part_Ex(part_Ex), part_Ey(part_Ey), part_Ez(part_Ez), part_Bx(part_Bx),
      part_By(part_By), part_Bz(part_Bz), part_mx(part_mx), part_my(part_my), part_mz(part_mz) {}

  __device__ void operator()(size_t ip) const {

    // 1/2 E
    mini_float px = qp * part_Ex[ip];
    mini_float py = qp * part_Ey[ip];
    mini_float pz = qp * part_Ez[ip];

    const mini_float ux = part_mx[ip] + px;
    const mini_float uy = part_my[ip] + py;
    const mini_float uz = part_mz[ip] + pz;

    // gamma-factor
    mini_float gamma_inv = qp / sqrt(1 + (ux * ux + uy * uy + uz * uz));

    // B, T = Transform to rotate the particle
    const mini_float tx  = gamma_inv * part_Bx[ip];
    const mini_float ty  = gamma_inv * part_By[ip];
    const mini_float tz  = gamma_inv * part_Bz[ip];
    const mini_float tsq = 1.0 + (tx * tx + ty * ty + tz * tz);
    mini_float tsq_inv   = 1.0 / tsq;

    px += ((1.0 + tx * tx - ty * ty - tz * tz) * ux + 2.0 * (tx * ty + tz) * uy +
           2.0 * (tz * tx - ty) * uz) *
          tsq_inv;

    py += (2.0 * (tx * ty - tz) * ux + (1.0 - tx * tx + ty * ty - tz * tz) * uy +
           2.0 * (ty * tz + tx) * uz) *
          tsq_inv;

    pz += (2.0 * (tz * tx + ty) * ux + 2.0 * (ty * tz - tx) * uy +
           (1.0 - tx * tx - ty * ty + tz * tz) * uz) *
          tsq_inv;

    // Update inverse gamma factor
    gamma_inv = 1.0 / sqrt(1 + (px * px + py * py + pz * pz));

    // Update momentum
    part_mx[ip] = px;
    part_my[ip] = py;
    part_mz[ip] = pz;
  }
#endif
};

// //! \brief Lambda function to push the particle momentum
// struct position_pusher_lambda
// {
//   const mini_float dt;

//   position_pusher_lambda(mini_float dt) : dt(dt) {}

//   __host__ __device__
//   void operator()(thrust::tuple<
//           mini_float&, mini_float&, mini_float&,  // x, y, z
//           mini_float&, mini_float&, mini_float&  // mx, my, mz
//           > data) const
//   {

//     mini_float& x = thrust::get<0>(data);
//     mini_float& y = thrust::get<1>(data);
//     mini_float& z = thrust::get<2>(data);

//     mini_float& mx = thrust::get<3>(data);
//     mini_float& my = thrust::get<4>(data);
//     mini_float& mz = thrust::get<5>(data);

//     const mini_float gamma_inv = 1 / sqrt(1 + (mx * mx + my * my + mz * mz));

//     // Update position
//     x += mx * dt * gamma_inv;
//     y += my * dt * gamma_inv;
//     z += mz * dt * gamma_inv;

//   }
// };

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

#if defined(__MINIPIC_THRUST_COUNTING__)

    const mini_float *const __restrict__ part_Ex =
      patch.particles_m[is].Ex_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ part_Ey =
      patch.particles_m[is].Ey_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ part_Ez =
      patch.particles_m[is].Ez_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ part_Bx =
      patch.particles_m[is].Bx_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ part_By =
      patch.particles_m[is].By_.get_raw_pointer(minipic::device);
    const mini_float *const __restrict__ part_Bz =
      patch.particles_m[is].Bz_.get_raw_pointer(minipic::device);

    mini_float *const __restrict__ part_x =
      patch.particles_m[is].x_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ part_y =
      patch.particles_m[is].y_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ part_z =
      patch.particles_m[is].z_.get_raw_pointer(minipic::device);

    mini_float *const __restrict__ part_mx =
      patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ part_my =
      patch.particles_m[is].my_.get_raw_pointer(minipic::device);
    mini_float *const __restrict__ part_mz =
      patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

#endif

#if defined(__MINIPIC_THRUST_ZIP__)

    // Push the momentum
    thrust::for_each(
      thrust::make_zip_iterator(thrust::make_tuple(patch.particles_m[is].mx_.device_data_.begin(),
                                                   patch.particles_m[is].my_.device_data_.begin(),
                                                   patch.particles_m[is].mz_.device_data_.begin(),
                                                   patch.particles_m[is].Ex_.device_data_.begin(),
                                                   patch.particles_m[is].Ey_.device_data_.begin(),
                                                   patch.particles_m[is].Ez_.device_data_.begin(),
                                                   patch.particles_m[is].Bx_.device_data_.begin(),
                                                   patch.particles_m[is].By_.device_data_.begin(),
                                                   patch.particles_m[is].Bz_.device_data_.begin())),
      thrust::make_zip_iterator(
        thrust::make_tuple(patch.particles_m[is].mx_.device_data_.begin() + n_particles,
                           patch.particles_m[is].my_.device_data_.begin() + n_particles,
                           patch.particles_m[is].mz_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Ex_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Ey_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Ez_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Bx_.device_data_.begin() + n_particles,
                           patch.particles_m[is].By_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Bz_.device_data_.begin() + n_particles)),
      momentum_pusher_lambda(dt, qp));

    // Update the position
    thrust::for_each(
      thrust::make_zip_iterator(thrust::make_tuple(patch.particles_m[is].x_.device_data_.begin(),
                                                   patch.particles_m[is].y_.device_data_.begin(),
                                                   patch.particles_m[is].z_.device_data_.begin(),
                                                   patch.particles_m[is].mx_.device_data_.begin(),
                                                   patch.particles_m[is].my_.device_data_.begin(),
                                                   patch.particles_m[is].mz_.device_data_.begin())),
      thrust::make_zip_iterator(
        thrust::make_tuple(patch.particles_m[is].x_.device_data_.begin() + n_particles,
                           patch.particles_m[is].y_.device_data_.begin() + n_particles,
                           patch.particles_m[is].z_.device_data_.begin() + n_particles,
                           patch.particles_m[is].mx_.device_data_.begin() + n_particles,
                           patch.particles_m[is].my_.device_data_.begin() + n_particles,
                           patch.particles_m[is].mz_.device_data_.begin() + n_particles)),
      [=] __device__(
        thrust::
          tuple<mini_float &, mini_float &, mini_float &, mini_float &, mini_float &, mini_float &>
            data) {
        mini_float &x = thrust::get<0>(data);
        mini_float &y = thrust::get<1>(data);
        mini_float &z = thrust::get<2>(data);

        mini_float &mx = thrust::get<3>(data);
        mini_float &my = thrust::get<4>(data);
        mini_float &mz = thrust::get<5>(data);

        const mini_float gamma_inv = 1 / sqrt(1 + (mx * mx + my * my + mz * mz));

        // Update position
        x += mx * dt * gamma_inv;
        y += my * dt * gamma_inv;
        z += mz * dt * gamma_inv;
      });

#elif defined(__MINIPIC_THRUST_COUNTING__)

    // Push the momentum
    thrust::for_each(thrust::make_counting_iterator<size_t>(0),
                     thrust::make_counting_iterator<size_t>(n_particles),
                     momentum_pusher_lambda(dt,
                                            qp,
                                            part_Ex,
                                            part_Ey,
                                            part_Ez,
                                            part_Bx,
                                            part_By,
                                            part_Bz,
                                            part_mx,
                                            part_my,
                                            part_mz));

    // Update the position
    thrust::for_each(thrust::make_counting_iterator<size_t>(0),
                     thrust::make_counting_iterator<size_t>(n_particles),
                     [=] __device__(size_t ip) {
                       const mini_float gamma_inv =
                         1 / sqrt(1 + (part_mx[ip] * part_mx[ip] + part_my[ip] * part_my[ip] +
                                       part_mz[ip] * part_mz[ip]));

                       // Update position
                       part_x[ip] += part_mx[ip] * dt * gamma_inv;
                       part_y[ip] += part_my[ip] * dt * gamma_inv;
                       part_z[ip] += part_mz[ip] * dt * gamma_inv;
                     });
#endif
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

#if defined(__MINIPIC_THRUST_ZIP__)

    // Push the momentum
    thrust::for_each(
      thrust::make_zip_iterator(thrust::make_tuple(patch.particles_m[is].mx_.device_data_.begin(),
                                                   patch.particles_m[is].my_.device_data_.begin(),
                                                   patch.particles_m[is].mz_.device_data_.begin(),
                                                   patch.particles_m[is].Ex_.device_data_.begin(),
                                                   patch.particles_m[is].Ey_.device_data_.begin(),
                                                   patch.particles_m[is].Ez_.device_data_.begin(),
                                                   patch.particles_m[is].Bx_.device_data_.begin(),
                                                   patch.particles_m[is].By_.device_data_.begin(),
                                                   patch.particles_m[is].Bz_.device_data_.begin())),
      thrust::make_zip_iterator(
        thrust::make_tuple(patch.particles_m[is].mx_.device_data_.begin() + n_particles,
                           patch.particles_m[is].my_.device_data_.begin() + n_particles,
                           patch.particles_m[is].mz_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Ex_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Ey_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Ez_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Bx_.device_data_.begin() + n_particles,
                           patch.particles_m[is].By_.device_data_.begin() + n_particles,
                           patch.particles_m[is].Bz_.device_data_.begin() + n_particles)),
      momentum_pusher_lambda(dt, qp));

#elif defined(__MINIPIC_THRUST_COUNTING__)

    // Push the momentum
    thrust::for_each(
      thrust::make_counting_iterator<size_t>(0),
      thrust::make_counting_iterator<size_t>(n_particles),
      momentum_pusher_lambda(dt,
                             qp,
                             patch.particles_m[is].Ex_.get_raw_pointer(minipic::device),
                             patch.particles_m[is].Ey_.get_raw_pointer(minipic::device),
                             patch.particles_m[is].Ez_.get_raw_pointer(minipic::device),
                             patch.particles_m[is].Bx_.get_raw_pointer(minipic::device),
                             patch.particles_m[is].By_.get_raw_pointer(minipic::device),
                             patch.particles_m[is].Bz_.get_raw_pointer(minipic::device),
                             patch.particles_m[is].mx_.get_raw_pointer(minipic::device),
                             patch.particles_m[is].my_.get_raw_pointer(minipic::device),
                             patch.particles_m[is].mz_.get_raw_pointer(minipic::device)));
#endif

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

    const mini_float domain_x_min = params.inf_x;
    const mini_float domain_y_min = params.inf_y;
    const mini_float domain_z_min = params.inf_z;
    const mini_float domain_x_max = params.sup_x, domain_y_max = params.sup_y,
                     domain_z_max = params.sup_z;

    // Periodic conditions
    if (params.boundary_condition_code == 1) {

      const int nx_patch = patch.nx_patchs_m, ny_patch = patch.ny_patchs_m,
                nz_patch = patch.nz_patchs_m;

      const mini_float Lx = params.Lx;
      const mini_float Ly = params.Ly;
      const mini_float Lz = params.Lz;

      for (int is = 0; is < patch.n_species_m; is++) {

        const size_t n_particles = patch.particles_m[is].size();

#if defined(__MINIPIC_THRUST_ZIP__)

        thrust::for_each(
          thrust::make_zip_iterator(
            thrust::make_tuple(patch.particles_m[is].x_.device_data_.begin(),
                               patch.particles_m[is].y_.device_data_.begin(),
                               patch.particles_m[is].z_.device_data_.begin())),
          thrust::make_zip_iterator(
            thrust::make_tuple(patch.particles_m[is].x_.device_data_.begin() + n_particles,
                               patch.particles_m[is].y_.device_data_.begin() + n_particles,
                               patch.particles_m[is].z_.device_data_.begin() + n_particles)),
          [=] __device__(thrust::tuple<mini_float &, mini_float &, mini_float &> data) {
            auto &x = thrust::get<0>(data);
            auto &y = thrust::get<1>(data);
            auto &z = thrust::get<2>(data);

            // Only relevant if there is just 1 patch in this direction
            // Else the patch exchange with periodicity is managed in the dedicated function
            if (nx_patch == 1) {
              if (x >= domain_x_max) {
                x -= Lx;
              } else if (x < domain_x_min) {
                x += Lx;
              }
            }
            // y direction
            if (ny_patch == 1) {
              if (y >= domain_y_max) {
                y -= Ly;
              } else if (y < domain_y_min) {
                y += Ly;
              }
            }
            // z direction
            if (nz_patch == 1) {
              if (z >= domain_z_max) {
                z -= Lz;
              } else if (z < domain_z_min) {
                z += Lz;
              }
            }
          });

#elif defined(__MINIPIC_THRUST_COUNTING__)

        mini_float *const __restrict__ part_x =
          patch.particles_m[is].x_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ part_y =
          patch.particles_m[is].y_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ part_z =
          patch.particles_m[is].z_.get_raw_pointer(minipic::device);

        thrust::for_each(
          // thrust::device,
          thrust::make_counting_iterator<size_t>(0),
          thrust::make_counting_iterator(n_particles),
          [=] __device__(size_t ip) {
            // Only relevant if there is just 1 patch in this direction
            // Else the patch exchange with periodicity is managed in the dedicated function
            if (nx_patch == 1) {
              if (part_x[ip] >= domain_x_max) {
                part_x[ip] -= Lx;
              } else if (part_x[ip] < domain_x_min) {
                part_x[ip] += Lx;
              }
            }
            // y direction
            if (ny_patch == 1) {
              if (part_y[ip] >= domain_y_max) {
                part_y[ip] -= Ly;
              } else if (part_y[ip] < domain_y_min) {
                part_y[ip] += Ly;
              }
            }
            // z direction
            if (nz_patch == 1) {
              if (part_z[ip] >= domain_z_max) {
                part_z[ip] -= Lz;
              } else if (part_z[ip] < domain_z_min) {
                part_z[ip] += Lz;
              }
            }
          });
#endif

      } // End loop on species

      // Reflective conditions
    } else if (params.boundary_condition_code == 2) {

      for (int is = 0; is < patch.n_species_m; is++) {

        const size_t n_particles = patch.particles_m[is].size();

#if defined(__MINIPIC_THRUST_ZIP__)

        thrust::for_each(thrust::make_zip_iterator(
                           thrust::make_tuple(patch.particles_m[is].x_.device_data_.begin(),
                                              patch.particles_m[is].y_.device_data_.begin(),
                                              patch.particles_m[is].z_.device_data_.begin(),
                                              patch.particles_m[is].mx_.device_data_.begin(),
                                              patch.particles_m[is].my_.device_data_.begin(),
                                              patch.particles_m[is].mz_.device_data_.begin())),
                         thrust::make_zip_iterator(thrust::make_tuple(
                           patch.particles_m[is].x_.device_data_.begin() + n_particles,
                           patch.particles_m[is].y_.device_data_.begin() + n_particles,
                           patch.particles_m[is].z_.device_data_.begin() + n_particles,
                           patch.particles_m[is].mx_.device_data_.begin() + n_particles,
                           patch.particles_m[is].my_.device_data_.begin() + n_particles,
                           patch.particles_m[is].mz_.device_data_.begin() + n_particles)),
                         [=] __device__(thrust::tuple<mini_float &,
                                                      mini_float &,
                                                      mini_float &,
                                                      mini_float &,
                                                      mini_float &,
                                                      mini_float &> data) {
                           auto &x  = thrust::get<0>(data);
                           auto &y  = thrust::get<1>(data);
                           auto &z  = thrust::get<2>(data);
                           auto &mx = thrust::get<3>(data);
                           auto &my = thrust::get<4>(data);
                           auto &mz = thrust::get<5>(data);

                           if (x >= domain_x_max) {
                             x  = 2 * domain_x_max - x;
                             mx = -mx;
                           } else if (x < domain_x_min) {
                             x  = 2 * domain_x_min - x;
                             mx = -mx;
                           }

                           if (y >= domain_y_max) {
                             y  = 2 * domain_y_max - y;
                             my = -my;
                           } else if (y < domain_y_min) {
                             y  = 2 * domain_y_min - y;
                             my = -my;
                           }

                           if (z >= domain_z_max) {
                             z  = 2 * domain_z_max - z;
                             mz = -mz;
                           } else if (z < domain_z_min) {
                             z  = 2 * domain_z_min - z;
                             mz = -mz;
                           }
                         });

#elif defined(__MINIPIC_THRUST_COUNTING__)

        mini_float *const __restrict__ part_x =
          patch.particles_m[is].x_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ part_y =
          patch.particles_m[is].y_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ part_z =
          patch.particles_m[is].z_.get_raw_pointer(minipic::device);

        mini_float *const __restrict__ part_mx =
          patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ part_my =
          patch.particles_m[is].my_.get_raw_pointer(minipic::device);
        mini_float *const __restrict__ part_mz =
          patch.particles_m[is].mz_.get_raw_pointer(minipic::device);

        thrust::for_each(thrust::make_counting_iterator<size_t>(0),
                         thrust::make_counting_iterator<size_t>(n_particles),
                         [=] __device__(size_t ip) {
                           if (part_x[ip] >= domain_x_max) {
                             part_x[ip]  = 2 * domain_x_max - part_x[ip];
                             part_mx[ip] = -part_mx[ip];
                           } else if (part_x[ip] < domain_x_min) {
                             part_x[ip]  = 2 * domain_x_min - part_x[ip];
                             part_mx[ip] = -part_mx[ip];
                           }

                           if (part_y[ip] >= domain_y_max) {
                             part_y[ip]  = 2 * domain_y_max - part_y[ip];
                             part_my[ip] = -part_my[ip];
                           } else if (part_y[ip] < domain_y_min) {
                             part_y[ip]  = 2 * domain_y_min - part_y[ip];
                             part_my[ip] = -part_my[ip];
                           }

                           if (part_z[ip] >= domain_z_max) {
                             part_z[ip]  = 2 * domain_z_max - part_z[ip];
                             part_mz[ip] = -part_mz[ip];
                           } else if (part_z[ip] < domain_z_min) {
                             part_z[ip]  = 2 * domain_z_min - part_z[ip];
                             part_mz[ip] = -part_mz[ip];
                           }
                         });

#endif

      } // End loop on species
    } // if type of conditions
  } // if on border
}

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
                nz_Jx = patch.vec_Jx_m[is].nz_m;
      const int nx_Jy = patch.vec_Jy_m[is].nx_m, ny_Jy = patch.vec_Jy_m[is].ny_m,
                nz_Jy = patch.vec_Jy_m[is].nz_m;
      const int nx_Jz = patch.vec_Jz_m[is].nx_m, ny_Jz = patch.vec_Jz_m[is].ny_m,
                nz_Jz = patch.vec_Jz_m[is].nz_m;

      const auto nynz_Jx = patch.vec_Jx_m[is].nynz_;
      const auto nynz_Jy = patch.vec_Jy_m[is].nynz_;
      const auto nynz_Jz = patch.vec_Jz_m[is].nynz_;

      const mini_float inv_cell_volume_x_q =
        params.inv_cell_volume * patch.particles_m[is].charge_m;
      const mini_float dt = params.dt;

      const mini_float inv_dx = params.inv_dx;
      const mini_float inv_dy = params.inv_dy;
      const mini_float inv_dz = params.inv_dz;

      const mini_float xmin = patch.inf_m[0];
      const mini_float ymin = patch.inf_m[1];
      const mini_float zmin = patch.inf_m[2];

      mini_float *Jx = patch.vec_Jx_m[is].get_raw_pointer(minipic::device);
      mini_float *Jy = patch.vec_Jy_m[is].get_raw_pointer(minipic::device);
      mini_float *Jz = patch.vec_Jz_m[is].get_raw_pointer(minipic::device);

#if defined(__MINIPIC_THRUST_COUNTING__)

      mini_float *const __restrict__ part_w =
        patch.particles_m[is].weight_.get_raw_pointer(minipic::device);

      mini_float *const __restrict__ part_x =
        patch.particles_m[is].x_.get_raw_pointer(minipic::device);
      mini_float *const __restrict__ part_y =
        patch.particles_m[is].y_.get_raw_pointer(minipic::device);
      mini_float *const __restrict__ part_z =
        patch.particles_m[is].z_.get_raw_pointer(minipic::device);

      mini_float *const __restrict__ part_mx =
        patch.particles_m[is].mx_.get_raw_pointer(minipic::device);
      mini_float *const __restrict__ part_my =
        patch.particles_m[is].my_.get_raw_pointer(minipic::device);
      mini_float *const __restrict__ part_mz =
        patch.particles_m[is].mz_.get_raw_pointer(minipic::device);
#endif

#if defined(__MINIPIC_THRUST_ZIP__)

      thrust::for_each(
        thrust::make_zip_iterator(
          thrust::make_tuple(patch.particles_m[is].weight_.device_data_.begin(),
                             patch.particles_m[is].x_.device_data_.begin(),
                             patch.particles_m[is].y_.device_data_.begin(),
                             patch.particles_m[is].z_.device_data_.begin(),
                             patch.particles_m[is].mx_.device_data_.begin(),
                             patch.particles_m[is].my_.device_data_.begin(),
                             patch.particles_m[is].mz_.device_data_.begin())),
        thrust::make_zip_iterator(
          thrust::make_tuple(patch.particles_m[is].weight_.device_data_.begin() + n_particles,
                             patch.particles_m[is].x_.device_data_.begin() + n_particles,
                             patch.particles_m[is].y_.device_data_.begin() + n_particles,
                             patch.particles_m[is].z_.device_data_.begin() + n_particles,
                             patch.particles_m[is].mx_.device_data_.begin() + n_particles,
                             patch.particles_m[is].my_.device_data_.begin() + n_particles,
                             patch.particles_m[is].mz_.device_data_.begin() + n_particles)),
        [=] __device__(thrust::tuple<mini_float &,
                                     mini_float &,
                                     mini_float &,
                                     mini_float &,
                                     mini_float &,
                                     mini_float &,
                                     mini_float &> data) {
          mini_float &w  = thrust::get<0>(data);
          mini_float &x  = thrust::get<1>(data);
          mini_float &y  = thrust::get<2>(data);
          mini_float &z  = thrust::get<3>(data);
          mini_float &mx = thrust::get<4>(data);
          mini_float &my = thrust::get<5>(data);
          mini_float &mz = thrust::get<6>(data);

          const mini_float gamma_inv = 1 / sqrt(1 + mx * mx + my * my + mz * mz);

          const mini_float charge_weight = inv_cell_volume_x_q * w;

#elif defined(__MINIPIC_THRUST_COUNTING__)

      thrust::for_each(
        thrust::make_counting_iterator<size_t>(0),
        thrust::make_counting_iterator(n_particles),
        [=] __device__(size_t ip) {
          const mini_float gamma_inv =
            1 / sqrt(1 + part_mx[ip] * part_mx[ip] + part_my[ip] * part_my[ip] +
                     part_mz[ip] * part_mz[ip]);

          const mini_float charge_weight = inv_cell_volume_x_q * part_w[ip];
#endif

#if defined(__MINIPIC_THRUST_ZIP__)
          const mini_float vx = mx * gamma_inv;
          const mini_float vy = my * gamma_inv;
          const mini_float vz = mz * gamma_inv;
#elif defined(__MINIPIC_THRUST_COUNTING__)
          const mini_float vx = part_mx[ip] * gamma_inv;
          const mini_float vy = part_my[ip] * gamma_inv;
          const mini_float vz = part_mz[ip] * gamma_inv;
#endif

          // Current from the particle
          const mini_float Jxp = vx * charge_weight;
          const mini_float Jyp = vy * charge_weight;
          const mini_float Jzp = vz * charge_weight;

        // Calculate normalized position relative to the patch
#if defined(__MINIPIC_THRUST_ZIP__)
          const mini_float posxn = (x - 0.5 * dt * vx - xmin) * inv_dx + 1;
          const mini_float posyn = (y - 0.5 * dt * vy - ymin) * inv_dy + 1;
          const mini_float poszn = (z - 0.5 * dt * vz - zmin) * inv_dz + 1;
#elif defined(__MINIPIC_THRUST_COUNTING__)
          const mini_float posxn = (part_x[ip] - 0.5 * dt * vx - xmin) * inv_dx + 1;
          const mini_float posyn = (part_y[ip] - 0.5 * dt * vy - ymin) * inv_dy + 1;
          const mini_float poszn = (part_z[ip] - 0.5 * dt * vz - zmin) * inv_dz + 1;
#endif

          // Compute indexes in primal grid
          const int ixp = static_cast<int>(floor(posxn));
          const int iyp = static_cast<int>(floor(posyn));
          const int izp = static_cast<int>(floor(poszn));

          // Compute indexes in dual grid
          const int ixd = static_cast<int>(floor(posxn - 0.5));
          const int iyd = static_cast<int>(floor(posyn - 0.5));
          const int izd = static_cast<int>(floor(poszn - 0.5));

          // Projection particle on currant field
          // Compute interpolation coeff, p = primal, d = dual
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
        });

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

    const size_t n_particles                = patch.particles_m[is].size();
    const mini_float inv_cell_volume_x_q = params.inv_cell_volume * patch.particles_m[is].charge_m;
    // double m       = particles_m[is].mass_m;

    Vector<mini_float> &w = patch.particles_m[is].weight_;

    Vector<mini_float> &x = patch.particles_m[is].x_;
    Vector<mini_float> &y = patch.particles_m[is].y_;
    Vector<mini_float> &z = patch.particles_m[is].z_;

    Vector<mini_float> &mx = patch.particles_m[is].mx_;
    Vector<mini_float> &my = patch.particles_m[is].my_;
    Vector<mini_float> &mz = patch.particles_m[is].mz_;

    for (size_t part = 0; part < n_particles; ++part) {

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
//! \param em electromagnetic fields
//! \param profiler profiler to handle time profiler
// _______________________________________________________
auto solve_maxwell(const Params &params, ElectroMagn &em, Profiler &profiler) -> void {

  const auto dt         = params.dt;
  const auto dt_over_dx = params.dt * params.inv_dx;
  const auto dt_over_dy = params.dt * params.inv_dy;
  const auto dt_over_dz = params.dt * params.inv_dz;

  mini_float *Ex   = em.Ex_m.get_raw_pointer(minipic::device);
  const auto nx_Ex = em.Ex_m.nx(), ny_Ex = em.Ex_m.ny(), nz_Ex = em.Ex_m.nz();

  mini_float *Ey   = em.Ey_m.get_raw_pointer(minipic::device);
  const auto nx_Ey = em.Ey_m.nx(), ny_Ey = em.Ey_m.ny(), nz_Ey = em.Ey_m.nz();

  mini_float *Ez   = em.Ez_m.get_raw_pointer(minipic::device);
  const auto nx_Ez = em.Ez_m.nx(), ny_Ez = em.Ez_m.ny(), nz_Ez = em.Ez_m.nz();

  mini_float *Bx   = em.Bx_m.get_raw_pointer(minipic::device);
  const auto nx_Bx = em.Bx_m.nx(), ny_Bx = em.Bx_m.ny(), nz_Bx = em.Bx_m.nz();

  mini_float *By   = em.By_m.get_raw_pointer(minipic::device);
  const auto nx_By = em.By_m.nx(), ny_By = em.By_m.ny(), nz_By = em.By_m.nz();

  mini_float *Bz   = em.Bz_m.get_raw_pointer(minipic::device);
  const auto nx_Bz = em.Bz_m.nx(), ny_Bz = em.Bz_m.ny(), nz_Bz = em.Bz_m.nz();

  mini_float *Jx   = em.Jx_m.get_raw_pointer(minipic::device);
  const auto nx_Jx = em.Jx_m.nx(), ny_Jx = em.Jx_m.ny(), nz_Jx = em.Jx_m.nz();

  mini_float *Jy   = em.Jy_m.get_raw_pointer(minipic::device);
  const auto nx_Jy = em.Jy_m.nx(), ny_Jy = em.Jy_m.ny(), nz_Jy = em.Jy_m.nz();

  mini_float *Jz   = em.Jz_m.get_raw_pointer(minipic::device);
  const auto nx_Jz = em.Jz_m.nx(), ny_Jz = em.Jz_m.ny(), nz_Jz = em.Jz_m.nz();

  const auto nynz_Jx = ny_Jx * nz_Jx;
  const auto nynz_Jy = ny_Jy * nz_Jy;
  const auto nynz_Jz = ny_Jz * nz_Jz;

  const auto nynz_Ex = ny_Ex * nz_Ex;
  const auto nynz_Ey = ny_Ey * nz_Ey;
  const auto nynz_Ez = ny_Ez * nz_Ez;

  const auto nynz_Bx = ny_Bx * nz_Bx;
  const auto nynz_By = ny_By * nz_By;
  const auto nynz_Bz = ny_Bz * nz_Bz;

  thrust::counting_iterator<int> index(0);

  /////     Solve Maxwell Ampere (E)
  // Electric field Ex (d,p,p)

  thrust::for_each(thrust::device, index, index + nx_Ex * nynz_Ex, [=] __device__(int idx) {
    const int ix = idx / (nynz_Ex);
    const int iy = (idx - ix * nynz_Ex) / nz_Ex;
    const int iz = idx - ix * nynz_Ex - iy * nz_Ex;

    Ex[idx] +=
      -dt * Jx[ix * (nynz_Jx) + (iy + 1) * (nz_Jx) + iz + 1] +
      dt_over_dy *
        (Bz[ix * (nynz_Bz) + (iy + 1) * (nz_Bz) + iz] - Bz[ix * (nynz_Bz) + iy * (nz_Bz) + iz]) -
      dt_over_dz *
        (By[ix * (nynz_By) + iy * (nz_By) + iz + 1] - By[ix * (nynz_By) + iy * (nz_By) + iz]);
  });

  // Electric field Ey (p,d,p)

  thrust::for_each(thrust::device, index, index + nx_Ey * nynz_Ey, [=] __device__(int idx) {
    const int ix = idx / (nynz_Ey);
    const int iy = (idx - ix * nynz_Ey) / nz_Ey;
    const int iz = idx - ix * nynz_Ey - iy * nz_Ey;

    Ey[idx] +=
      -dt * Jy[(ix + 1) * (nynz_Jy) + iy * (nz_Jy) + iz + 1] -
      dt_over_dx *
        (Bz[(ix + 1) * (nynz_Bz) + iy * (nz_Bz) + iz] - Bz[ix * (nynz_Bz) + iy * (nz_Bz) + iz]) +
      dt_over_dz *
        (Bx[ix * (nynz_Bx) + iy * (nz_Bx) + iz + 1] - Bx[ix * (nynz_Bx) + iy * (nz_Bx) + iz]);
  });

  // Electric field Ez (p,p,d)

  thrust::for_each(thrust::device, index, index + nx_Ez * nynz_Ez, [=] __device__(int idx) {
    const int ix = idx / (nynz_Ez);
    const int iy = (idx - ix * nynz_Ez) / nz_Ez;
    const int iz = idx - ix * nynz_Ez - iy * nz_Ez;

    Ez[idx] +=
      -dt * Jz[(ix + 1) * (nynz_Jz) + (iy + 1) * (nz_Jz) + iz] +
      dt_over_dx *
        (By[(ix + 1) * (nynz_By) + iy * (nz_By) + iz] - By[ix * (nynz_By) + iy * (nz_By) + iz]) -
      dt_over_dy *
        (Bx[ix * (nynz_Bx) + (iy + 1) * (nz_Bx) + iz] - Bx[ix * (nynz_Bx) + iy * (nz_Bx) + iz]);
  });

  /////     Solve Maxwell Faraday (B)

  // Magnetic field Bx (p,d,d)

  int nz   = nz_Bx - 2;
  int nynz = (ny_Bx - 2) * nz;
  thrust::for_each(thrust::device,
                   index,
                   index + nx_Bx * (ny_Bx - 2) * (nz_Bx - 2),
                   [=] __device__(int idx) {
                     const int ix = idx / nynz;
                     int iy       = (idx - ix * nynz) / nz;
                     int iz       = idx - ix * nynz - iy * nz;

                     iy += 1;
                     iz += 1;

                     Bx[ix * nynz_Bx + iy * nz_Bx + iz] +=
                       -dt_over_dy * (Ez[ix * (nynz_Ez) + iy * (nz_Ez) + iz] -
                                      Ez[ix * (nynz_Ez) + (iy - 1) * (nz_Ez) + iz]) +
                       dt_over_dz * (Ey[ix * (nynz_Ey) + iy * (nz_Ey) + iz] -
                                     Ey[ix * (nynz_Ey) + iy * (nz_Ey) + iz - 1]);
                   });

  // Magnetic field By (d,p,d)

  nz   = nz_By - 2;
  nynz = ny_By * nz;
  thrust::for_each(thrust::device,
                   index,
                   index + (nx_By - 2) * ny_By * (nz_By - 2),
                   [=] __device__(int idx) {
                     int ix       = idx / nynz;
                     const int iy = (idx - ix * nynz) / nz;
                     int iz       = idx - ix * nynz - iy * nz;

                     ix += 1;
                     iz += 1;

                     By[ix * nynz_By + iy * nz_By + iz] +=
                       -dt_over_dz * (Ex[ix * (nynz_Ex) + iy * (nz_Ex) + iz] -
                                      Ex[(ix) * (nynz_Ex) + iy * (nz_Ex) + iz - 1]) +
                       dt_over_dx * (Ez[ix * (nynz_Ez) + iy * (nz_Ez) + iz] -
                                     Ez[(ix - 1) * (nynz_Ez) + iy * (nz_Ez) + iz]);
                   });

  // Magnetic field Bz (d,d,p)
  nynz = (ny_Bz - 2) * nz_Bz;
  thrust::for_each(thrust::device,
                   index,
                   index + (nx_Bz - 2) * (ny_Bz - 2) * nz_Bz,
                   [=] __device__(int idx) {
                     int ix       = idx / nynz;
                     int iy       = (idx - ix * nynz) / nz_Bz;
                     const int iz = idx - ix * nynz - iy * nz_Bz;

                     ix += 1;
                     iy += 1;

                     Bz[ix * nynz_Bz + iy * nz_Bz + iz] +=
                       -dt_over_dx * (Ey[ix * (nynz_Ey) + iy * (nz_Ey) + iz] -
                                      Ey[(ix - 1) * (nynz_Ey) + iy * (nz_Ey) + iz]) +
                       dt_over_dy * (Ex[ix * (nynz_Ex) + iy * (nz_Ex) + iz] -
                                     Ex[ix * (nynz_Ex) + (iy - 1) * (nz_Ex) + iz]);
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

    const auto nx_Jx = em.Jx_m.nx();
    const auto ny_Jx = em.Jx_m.ny();
    const auto nz_Jx = em.Jx_m.nz();

    const auto nx_Jy = em.Jy_m.nx();
    const auto ny_Jy = em.Jy_m.ny();
    const auto nz_Jy = em.Jy_m.nz();

    const auto nx_Jz = em.Jz_m.nx();
    const auto ny_Jz = em.Jz_m.ny();
    const auto nz_Jz = em.Jz_m.nz();

    const auto nynz_Jx = ny_Jx * nz_Jx;
    const auto nynz_Jy = ny_Jy * nz_Jy;
    const auto nynz_Jz = em.Jz_m.nynz_;

    mini_float *Jx = em.Jx_m.get_raw_pointer(minipic::device);
    mini_float *Jy = em.Jy_m.get_raw_pointer(minipic::device);
    mini_float *Jz = em.Jz_m.get_raw_pointer(minipic::device);

    // _______________________________________________________________
    // X faces

    thrust::counting_iterator<int> index(0);

    thrust::for_each(thrust::device, index, index + nynz_Jx, [=] __device__(int i) {
      // const int iy = i / nz_loc;
      // const int iz = i - iy * nz_loc;

      {
        const auto index_left  = i;
        const auto index_right = (nx_Jx - 2) * nynz_Jx + i;

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }

      {
        const auto index_left  = i + nynz_Jx;
        const auto index_right = (nx_Jx - 1) * nynz_Jx + i;

        Jx[index_left] += Jx[index_right];
        Jx[index_right] = Jx[index_left];
      }
    });

    thrust::for_each(thrust::device, index, index + nynz_Jy, [=] __device__(int i) {
      // const int ix = i / nz_loc;
      // const int iz = i - ix * nz_loc;

      auto index_left  = i;
      auto index_right = i + (nx_Jy - 2) * nynz_Jy;

      Jy[index_left] += Jy[index_right];
      Jy[index_right] = Jy[index_left];

      index_left  = i + nynz_Jy;
      index_right = i + (nx_Jy - 1) * nynz_Jy;

      Jy[index_left] += Jy[index_right];
      Jy[index_right] = Jy[index_left];
    });

    thrust::for_each(thrust::device, index, index + nynz_Jz, [=] __device__(int i) {
      auto index_left  = i;
      auto index_right = i + (nx_Jz - 2) * nynz_Jz;

      Jz[index_left] += Jz[index_right];
      Jz[index_right] = Jz[index_left];

      index_left  = i + nynz_Jz;
      index_right = i + (nx_Jz - 1) * nynz_Jz;

      Jz[index_left] += Jz[index_right];
      Jz[index_right] = Jz[index_left];
    });

    // _______________________________________________________________
    // Y faces

    thrust::for_each(thrust::device, index, index + nx_Jx * nz_Jx, [=] __device__(int i) {
      const int ix = i / nz_Jx;
      const int iz = i - ix * nz_Jx;

      auto index_left  = ix * nynz_Jx + iz;                       // iy = 0
      auto index_right = ix * nynz_Jx + (ny_Jx - 2) * nz_Jx + iz; // iy = ny_Jx-2

      Jx[index_left] += Jx[index_right];
      Jx[index_right] = Jx[index_left];

      index_left  = ix * nynz_Jx + nz_Jx + iz;               // iy = 1
      index_right = ix * nynz_Jx + (ny_Jx - 1) * nz_Jx + iz; // iy = ny_Jx-1

      Jx[index_left] += Jx[index_right];
      Jx[index_right] = Jx[index_left];
    });

    thrust::for_each(thrust::device, index, index + nx_Jy * nz_Jy, [=] __device__(int i) {
      const int ix = i / nz_Jy;
      const int iz = i - ix * nz_Jy;

      auto index_left  = ix * nynz_Jy + iz;                       // iy = 0
      auto index_right = ix * nynz_Jy + (ny_Jy - 2) * nz_Jy + iz; // iy = ny_Jy-2

      Jy[index_left] += Jy[index_right];
      Jy[index_right] = Jy[index_left];

      index_left  = ix * nynz_Jy + nz_Jy + iz;               // iy = 1
      index_right = ix * nynz_Jy + (ny_Jy - 1) * nz_Jy + iz; // iy = ny_Jy-1

      Jy[index_left] += Jy[index_right];
      Jy[index_right] = Jy[index_left];
    });

    thrust::for_each(thrust::device, index, index + nx_Jz * nz_Jz, [=] __device__(int i) {
      const int ix = i / nz_Jz;
      const int iz = i - ix * nz_Jz;

      auto index_left  = ix * nynz_Jz + iz;                       // iy = 0
      auto index_right = ix * nynz_Jz + (ny_Jz - 2) * nz_Jz + iz; // iy = ny_Jz-2

      Jz[index_left] += Jz[index_right];
      Jz[index_right] = Jz[index_left];

      index_left  = ix * nynz_Jz + nz_Jz + iz;               // iy = 1
      index_right = ix * nynz_Jz + (ny_Jz - 1) * nz_Jz + iz; // iy = ny_Jz-1

      Jz[index_left] += Jz[index_right];
      Jz[index_right] = Jz[index_left];
    });

    // _______________________________________________________________
    // Z faces

    thrust::for_each(thrust::device, index, index + nx_Jx * ny_Jx, [=] __device__(int i) {
      const int ix = i / ny_Jx;
      const int iy = i - ix * ny_Jx;

      auto index_left  = ix * nynz_Jx + iy * nz_Jx;               // iz = 0
      auto index_right = ix * nynz_Jx + iy * nz_Jx + (nz_Jx - 2); // iz = nz_Jx-2

      Jx[index_left] += Jx[index_right];
      Jx[index_right] = Jx[index_left];

      index_left  = ix * nynz_Jx + iy * nz_Jx + 1;           // iz = 1
      index_right = ix * nynz_Jx + iy * nz_Jx + (nz_Jx - 1); // iz = nz_Jx-1

      Jx[index_left] += Jx[index_right];
      Jx[index_right] = Jx[index_left];
    });

    thrust::for_each(thrust::device, index, index + nx_Jy * ny_Jy, [=] __device__(int i) {
      const int ix = i / ny_Jy;
      const int iy = i - ix * ny_Jy;

      auto index_left  = ix * nynz_Jy + iy * nz_Jy;               // iz = 0
      auto index_right = ix * nynz_Jy + iy * nz_Jy + (nz_Jy - 2); // iz = nz_Jy-2

      Jy[index_left] += Jy[index_right];
      Jy[index_right] = Jy[index_left];

      index_left  = ix * nynz_Jy + iy * nz_Jy + 1;           // iz = 1
      index_right = ix * nynz_Jy + iy * nz_Jy + (nz_Jy - 1); // iz = nz_Jy-1

      Jy[index_left] += Jy[index_right];
      Jy[index_right] = Jy[index_left];
    });

    thrust::for_each(thrust::device, index, index + nx_Jz * ny_Jz, [=] __device__(int i) {
      const int ix = i / ny_Jz;
      const int iy = i - ix * ny_Jz;

      auto index_left  = ix * nynz_Jz + iy * nz_Jz;               // iz = 0
      auto index_right = ix * nynz_Jz + iy * nz_Jz + (nz_Jz - 2); // iz = nz_Jz-2

      Jz[index_left] += Jz[index_right];
      Jz[index_right] = Jz[index_left];

      index_left  = ix * nynz_Jz + iy * nz_Jz + 1;           // iz = 1
      index_right = ix * nynz_Jz + iy * nz_Jz + (nz_Jz - 1); // iz = nz_Jz-1

      Jz[index_left] += Jz[index_right];
      Jz[index_right] = Jz[index_left];
    });

  } // end if periodic
} // end currentBC

// _______________________________________________________________
//
//! \brief Boundaries condition on the global grid
//! \param[in] Params & params - global constant parameters
// _______________________________________________________________
auto solveBC(Params &params, ElectroMagn &em) -> void {

  mini_float *Bx   = em.Bx_m.get_raw_pointer(minipic::device);
  const auto nx_Bx = em.Bx_m.nx();
  const auto ny_Bx = em.Bx_m.ny();
  const auto nz_Bx = em.Bx_m.nz();

  mini_float *By   = em.By_m.get_raw_pointer(minipic::device);
  const auto nx_By = em.By_m.nx();
  const auto ny_By = em.By_m.ny();
  const auto nz_By = em.By_m.nz();

  mini_float *Bz   = em.Bz_m.get_raw_pointer(minipic::device);
  const auto nx_Bz = em.Bz_m.nx();
  const auto ny_Bz = em.Bz_m.ny();
  const auto nz_Bz = em.Bz_m.nz();

  const auto nynz_Bx = ny_Bx * nz_Bx;
  const auto nynz_By = ny_By * nz_By;
  const auto nynz_Bz = ny_Bz * nz_Bz;

  thrust::counting_iterator<int> index(0);

  if (params.boundary_condition == "periodic") {

    // X dim
    // By (d,p,d)
    thrust::for_each(thrust::device, index, index + nynz_By, [=] __device__(int i) {
      // -X
      By[i] = By[(nx_By - 2) * nynz_By + i];
      // +X
      By[(nx_By - 1) * nynz_By + i] = By[nynz_By + i];
    });

    // Bz (d,d,p)
    thrust::for_each(thrust::device, index, index + nynz_Bz, [=] __device__(int i) {
      // -X
      Bz[i] = Bz[(nx_Bz - 2) * nynz_Bz + i];
      // +X
      Bz[(nx_Bz - 1) * nynz_Bz + i] = Bz[nynz_Bz + i];
    });

    // Y dim
    // Bx (p,d,d)
    thrust::for_each(thrust::device, index, index + nx_Bx * nz_Bx, [=] __device__(int i) {
      const int ix = i / nz_Bx;
      const int iz = i - ix * nz_Bx;

      // -Y
      Bx[ix * nynz_Bx + iz] = Bx[ix * nynz_Bx + (ny_Bx - 2) * nz_Bx + iz];
      // +Y
      Bx[ix * nynz_Bx + (ny_Bx - 1) * nz_Bx + iz] = Bx[ix * nynz_Bx + nz_Bx + iz];
    });

    // Bz (d,d,p)
    thrust::for_each(thrust::device, index, index + nx_Bz * nz_Bz, [=] __device__(int i) {
      const int ix = i / nz_Bz;
      const int iz = i - ix * nz_Bz;

      // -Y
      Bz[ix * nynz_Bz + iz] = Bz[ix * nynz_Bz + (ny_Bz - 2) * nz_Bz + iz];
      // +Y
      Bz[ix * nynz_Bz + (ny_Bz - 1) * nz_Bz + iz] = Bz[ix * nynz_Bz + nz_Bz + iz];
    });

    // Z dim
    // Bx
    thrust::for_each(thrust::device, index, index + nx_Bx * ny_Bx, [=] __device__(int i) {
      const int ix = i / ny_Bx;
      const int iy = i - ix * ny_Bx;

      // -Z
      Bx[ix * nynz_Bx + iy * nz_Bx] = Bx[ix * nynz_Bx + iy * nz_Bx + (nz_Bx - 2)];
      // +Z
      Bx[ix * nynz_Bx + iy * nz_Bx + (nz_Bx - 1)] = Bx[ix * nynz_Bx + iy * nz_Bx + 1];
    });

    // By
    thrust::for_each(thrust::device, index, index + nx_By * ny_By, [=] __device__(int i) {
      const int ix = i / ny_By;
      const int iy = i - ix * ny_By;

      // -Z
      By[ix * nynz_By + iy * nz_By] = By[ix * nynz_By + iy * nz_By + (nz_By - 2)];
      // +Z
      By[ix * nynz_By + iy * nz_By + (nz_By - 1)] = By[ix * nynz_By + iy * nz_By + 1];
    });

  } else if (params.boundary_condition == "reflective") {

    // X dim
    // By (d,p,d)
    thrust::for_each(thrust::device, index, index + nynz_By, [=] __device__(int i) {
      // -X
      By[i] = By[nynz_By + i];
      // +X
      By[(nx_By - 1) * nynz_By + i] = By[(nx_By - 2) * nynz_By + i];
    });

    // Bz (d,d,p)
    thrust::for_each(thrust::device, index, index + nynz_Bz, [=] __device__(int i) {
      // -X
      Bz[i] = Bz[nynz_Bz + i];
      // +X
      Bz[(nx_Bz - 1) * nynz_Bz + i] = Bz[(nx_Bz - 2) * nynz_Bz + i];
    });

    // Y dim
    // Bx (p,d,d)
    thrust::for_each(thrust::device, index, index + nx_Bx * nz_Bx, [=] __device__(int i) {
      const int ix = i / nz_Bx;
      const int iz = i - ix * nz_Bx;

      // -Y
      Bx[ix * nynz_Bx + iz] = Bx[ix * nynz_Bx + nz_Bx + iz];
      // +Y
      Bx[ix * nynz_Bx + (ny_Bx - 1) * nz_Bx + iz] = Bx[ix * nynz_Bx + (ny_Bx - 2) * nz_Bx + iz];
    });

    // Bz (-1 to avoid corner)
    thrust::for_each(thrust::device, index, index + nx_Bz * nz_Bz, [=] __device__(int i) {
      const int ix = i / nz_Bz;
      const int iz = i - ix * nz_Bz;

      // -Y
      Bz[ix * nynz_Bz + iz] = Bz[ix * nynz_Bz + nz_Bz + iz];
      // +Y
      Bz[ix * nynz_Bz + (ny_Bz - 1) * nz_Bz + iz] = Bz[ix * nynz_Bz + (ny_Bz - 2) * nz_Bz + iz];
    });

    // Z dim
    // Bx
    thrust::for_each(thrust::device, index, index + nx_Bx * ny_Bx, [=] __device__(int i) {
      const int ix = i / ny_Bx;
      const int iy = i - ix * ny_Bx;

      // -Z
      Bx[ix * nynz_Bx + iy * nz_Bx] = Bx[ix * nynz_Bx + iy * nz_Bx + 1];
      // +Z
      Bx[ix * nynz_Bx + iy * nz_Bx + (nz_Bx - 1)] = Bx[ix * nynz_Bx + iy * nz_Bx + (nz_Bx - 2)];
    });

    // By
    thrust::for_each(thrust::device, index, index + nx_By * ny_By, [=] __device__(int i) {
      const int ix = i / ny_By;
      const int iy = i - ix * ny_By;

      // -Z
      By[ix * nynz_By + iy * nz_By] = By[ix * nynz_By + iy * nz_By + 1];
      // +Z
      By[ix * nynz_By + iy * nz_By + (nz_By - 1)] = By[ix * nynz_By + iy * nz_By + (nz_By - 2)];
    });

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
auto reduc_current(Patch &patch) -> void {

  for (int is = 1; is < patch.n_species_m; is++) {

    // Only if particles projected
    if (patch.projected_[is]) {

#if defined(__MINIPIC_THRUST__)

      thrust::transform(thrust::device,
                        patch.vec_Jx_m[is].device_data_.begin(),
                        patch.vec_Jx_m[is].device_data_.end(),
                        patch.vec_Jx_m[0].device_data_.begin(),
                        patch.vec_Jx_m[0].device_data_.begin(),
                        thrust::plus<mini_float>());

      thrust::transform(thrust::device,
                        patch.vec_Jy_m[is].device_data_.begin(),
                        patch.vec_Jy_m[is].device_data_.end(),
                        patch.vec_Jy_m[0].device_data_.begin(),
                        patch.vec_Jy_m[0].device_data_.begin(),
                        thrust::plus<mini_float>());

      thrust::transform(thrust::device,
                        patch.vec_Jz_m[is].device_data_.begin(),
                        patch.vec_Jz_m[is].device_data_.end(),
                        patch.vec_Jz_m[0].device_data_.begin(),
                        patch.vec_Jz_m[0].device_data_.begin(),
                        thrust::plus<mini_float>());

#elif defined(__MINIPIC_THRUST_UNIFIED__)

      thrust::transform(thrust::device,
                        patch.vec_Jx_m[is].data_.begin(),
                        patch.vec_Jx_m[is].data_.end(),
                        patch.vec_Jx_m[0].data_.begin(),
                        patch.vec_Jx_m[0].data_.begin(),
                        thrust::plus<mini_float>());

      thrust::transform(thrust::device,
                        patch.vec_Jy_m[is].data_.begin(),
                        patch.vec_Jy_m[is].data_.end(),
                        patch.vec_Jy_m[0].data_.begin(),
                        patch.vec_Jy_m[0].data_.begin(),
                        thrust::plus<mini_float>());

      thrust::transform(thrust::device,
                        patch.vec_Jz_m[is].data_.begin(),
                        patch.vec_Jz_m[is].data_.end(),
                        patch.vec_Jz_m[0].data_.begin(),
                        patch.vec_Jz_m[0].data_.begin(),
                        thrust::plus<mini_float>());
#endif

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
    const auto ix_origin_p = patch.ix_origin_m;
    const auto iy_origin_p = patch.iy_origin_m;
    const auto iz_origin_p = patch.iz_origin_m;
    const auto ix_origin_d = patch.ix_origin_m;
    const auto iy_origin_d = patch.iy_origin_m;
    const auto iz_origin_d = patch.iz_origin_m;

    // Use thrust transform and iterators from 0 to nx_Jx * ny_Jx * nz_Jx

    thrust::counting_iterator<int> index(0);

    {

      const auto nx = em.Jx_m.nx();
      const auto ny = em.Jx_m.ny();
      const auto nz = em.Jx_m.nz();

      const auto nx_loc     = patch.vec_Jx_m[0].nx();
      const auto ny_loc     = patch.vec_Jx_m[0].ny();
      const auto nz_loc     = patch.vec_Jx_m[0].nz();
      const auto nynz_loc   = patch.vec_Jx_m[0].nynz();
      const auto nz_loc_inv = 1. / patch.vec_Jx_m[0].nz();

      mini_float *J_loc = patch.vec_Jx_m[0].get_raw_pointer(minipic::device);
      mini_float *J     = em.Jx_m.get_raw_pointer(minipic::device);

      thrust::for_each(thrust::device,
                       index,
                       index + nx_loc * nynz_loc,
                       [=] __device__(int i) {
                         const int ix = i / (nynz_loc);
                         const int iy = (i - ix * nynz_loc) * nz_loc_inv;
                         const int iz = i - ix * nynz_loc - iy * nz_loc;

                         const auto global_index = (ix + ix_origin_d) * ny * nz +
                                                   (iy + iy_origin_p) * nz + (iz + iz_origin_p);

                         J[global_index] = J_loc[i];
                       });
    }

    {
      const auto nx = em.Jy_m.nx();
      const auto ny = em.Jy_m.ny();
      const auto nz = em.Jy_m.nz();

      const auto nx_loc = patch.vec_Jy_m[0].nx();
      const auto ny_loc = patch.vec_Jy_m[0].ny();
      const auto nz_loc = patch.vec_Jy_m[0].nz();

      mini_float *J_loc = patch.vec_Jy_m[0].get_raw_pointer(minipic::device);
      mini_float *J     = em.Jy_m.get_raw_pointer(minipic::device);

      thrust::for_each(thrust::device,
                       index,
                       index + nx_loc * ny_loc * nz_loc,
                       [=] __device__(int i) {
                         const int ix = i / (ny_loc * nz_loc);
                         const int iy = (i - ix * ny_loc * nz_loc) / nz_loc;
                         const int iz = i - ix * ny_loc * nz_loc - iy * nz_loc;

                         const auto global_index = (ix + ix_origin_p) * ny * nz +
                                                   (iy + iy_origin_d) * nz + (iz + iz_origin_p);

                         J[global_index] = J_loc[i];
                       });
    }

    {
      const auto nx = em.Jz_m.nx();
      const auto ny = em.Jz_m.ny();
      const auto nz = em.Jz_m.nz();

      const auto nx_loc = patch.vec_Jz_m[0].nx();
      const auto ny_loc = patch.vec_Jz_m[0].ny();
      const auto nz_loc = patch.vec_Jz_m[0].nz();

      mini_float *J_loc = patch.vec_Jz_m[0].get_raw_pointer(minipic::device);
      mini_float *J     = em.Jz_m.get_raw_pointer(minipic::device);

      thrust::for_each(thrust::device,
                       index,
                       index + nx_loc * ny_loc * nz_loc,
                       [=] __device__(int i) {
                         const int ix = i / (ny_loc * nz_loc);
                         const int iy = (i - ix * ny_loc * nz_loc) / nz_loc;
                         const int iz = i - ix * ny_loc * nz_loc - iy * nz_loc;

                         const auto global_index = (ix + ix_origin_p) * ny * nz +
                                                   (iy + iy_origin_p) * nz + (iz + iz_origin_d);

                         J[global_index] = J_loc[i];
                       });
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
