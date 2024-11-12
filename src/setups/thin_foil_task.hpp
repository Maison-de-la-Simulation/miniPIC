//! \file thin_foil.hpp

//! \brief This setup simulate the interaction of a thin foil with a laser beam

/* _____________________________________________________________________ */

#include "Params.hpp"

//! \brief Functiun to setup all input parameters
void setup(Params &params) {

  // Simulation name
  params.name = "Laser thin foil interaction";

  // Physics parameters
  const double temperature  = 0 / 511.;
  const double n0           = 10;
  const double debye_length = std::sqrt(temperature / n0);

  // Space
  params.inf_x = 0.;
  params.inf_y = 0.;
  params.inf_z = 0.;
  params.sup_x = 3.;
  params.sup_y = 1.;
  params.sup_z = 1.;

  // Decomp
  params.n_subdomains = 1;

  //  Bin number
#if defined(__MINIPIC_OMP_TASK__) || defined(__MINIPIC_EVENTIFY__)
  params.bin_size = 16000;
#endif

  // Number of patches
  params.nx_patch = 32;
  params.ny_patch = 16;
  params.nz_patch = 16;

  // Cells per patch per direction
  params.nx_cells_by_patch = 8;
  params.ny_cells_by_patch = 8;
  params.nz_cells_by_patch = 8;

  // // Number of patches
  // params.nx_patch = 32;
  // params.ny_patch = 8;
  // params.nz_patch = 8;

  // // Cells per patch per direction
  // params.nx_cells_by_patch = 12;
  // params.ny_cells_by_patch = 8;
  // params.nz_cells_by_patch = 64;

  // Time

  // const double dx = (params.sup_x - params.inf_x) / (params.nx_cells_by_patch * params.nx_patch);
  // const double dy = (params.sup_y - params.inf_y) / (params.ny_cells_by_patch * params.ny_patch);
  // const double dz = (params.sup_z - params.inf_z) / (params.nz_cells_by_patch * params.nz_patch);

  params.dt = 0.95; // Fraction of the CFL

  params.simulation_time = 1000 * params.dt;

  // Antenna profile to generate a gaussian laser beam
  auto laser_profile = [](double y, double z, double t) -> double {
    // Intensity
    const double E0 = 1;

    // Full width at half maximum of the focal spot
    const double fwhm_focal_spot = 0.1;

    // Period of a laser oscillation
    const double period = 0.025;

    // Full width at half maximum of the time envelope
    const double fwhm_time = 0.5;

    // Time at maximum intensity
    const double t0 = 0.75;

    // Compute the waist of the focal spot
    const double waist_focal_spot_square =
      fwhm_focal_spot * fwhm_focal_spot / (2.0 * std::log(2.0));

    // Compute the waist of the time envelope
    const double waist_time_square = fwhm_time * fwhm_time / (2.0 * std::log(2.0));

    const double alpha = 1. / waist_time_square;
    const double beta  = 2 * M_PI / period;
    const double trel  = t - t0;

    const double focal_spot = std::exp(-(y * y + z * z) / (waist_focal_spot_square));

    // Derivative of a E field of the form E0 * exp(-alpha * t^2) * cos(beta * t)

    return E0 * std::exp(-alpha * trel * trel) *
           (beta * std::cos(beta * trel) - 2 * alpha * trel * sin(beta * trel)) * focal_spot;
  };

  params.add_antenna(laser_profile, params.inf_x + 2 * dx);

  // foil
  // auto profile = [n0](double x, double y, double z) -> double {
  //   if (x > 0.4 && x < 0.6) {
  //     return n0;
  //   } else {
  //     return 0;
  //   }
  // };

  // sphere
  auto profile = [n0](double x, double y, double z) -> double {
    const double zc = (z * 3 - 1.5);
    const double R  = zc * zc + (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
    if (R < 0.25 * 0.25) {
      return n0;
    } else {
      return 0;
    }
  };

  // Hydrogen Plasma
  params.add_species("electron", 1, -1, temperature, profile, {0, 0, 0}, 128, "random", "cell");
  params.add_species("proton", 1836.125, 1, temperature, profile, {0, 0, 0}, 128, "electron", "cell");

  // Momentum correction at init
  params.momentum_correction = false;

  // Bourndary conditions
  params.boundary_condition = "periodic";
  params.no_diagnostics_at_init = true;

  // Display
  params.print_period = 5000;

  // Random seed
  params.seed = 0;

  // Field Diagnostics
  params.field_diagnostics_period = 50000;
  params.field_diagnostics_format = "vtk";

  // Cloud Diagnostics
  params.particle_cloud_period = 50000;
  params.particle_cloud_format = "binary";

  // Binning Diagnostics
  /*params.add_particle_binning("diag_x_y_z_d",
                              "density",
                              {"x", "y", "z"},
                              {params.nx_patch * params.nx_cells_by_patch,
                               params.ny_patch * params.ny_cells_by_patch,
                               params.nz_patch * params.nz_cells_by_patch},
                              {0., 0., 0.},
                              {params.sup_x, params.sup_y, params.sup_z},
                              {0, 1},
                              50,
                              "vtk");*/
}