//! \file antenna.hpp

//! \brief This setup generates a Gaussian laser field using an antenna

/* _____________________________________________________________________ */

#include "Params.hpp"

//! \brief Functiun to setup all input parameters
void setup(Params &params) {

  // Simulation name
  params.name = "antenna";

  // Space
  params.inf_x = 0.;
  params.inf_y = 0.;
  params.inf_z = 0.;
  params.sup_x = 3.;
  params.sup_y = 1.;
  params.sup_z = 1.;

  // Decomp
  params.n_subdomains = 1;

  // Number of patches
  params.nx_patch = 1;
  params.ny_patch = 1;
  params.nz_patch = 1;

  // Cells per patch per direction
  params.nx_cells_by_patch = 384;
  params.ny_cells_by_patch = 64;
  params.nz_cells_by_patch = 64;

  // Time

  // const double dx = (params.sup_x - params.inf_x) / (params.nx_cells_by_patch * params.nx_patch);
  // const double dy = (params.sup_y - params.inf_y) / (params.ny_cells_by_patch * params.ny_patch);
  // const double dz = (params.sup_z - params.inf_z) / (params.nz_cells_by_patch * params.nz_patch);

  params.dt = 0.95;

  params.simulation_time = 600 * params.dt;

  // Antenna profile to generate a gaussian laser beam
  auto profile = [](double y, double z, double t) -> double {
    // Intensity
    const double E0 = 1;

    // Full width at half maximum of the focal spot
    const double fwhm_focal_spot = 0.1;

    // Period of a laser oscillation
    const double period = 0.05;

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

  params.add_antenna(profile, params.inf_x + 0.5 * (params.sup_x - params.inf_x));

  // Momentum correction at init
  params.momentum_correction = false;

  // Current projection
  params.current_projection = false;

  // Bourndary conditions
  params.boundary_condition = "periodic";

  // Display
  params.print_period = 50;

  // Random seed
  params.seed = 0;

  // Scalar Diagnostics
  params.scalar_diagnostics_period = 10;

  // Field Diagnostics
  params.field_diagnostics_period = 100;
  params.field_diagnostics_format = "bin";

  // Timers
  params.save_timers_period      = 50;
  params.save_timers_start       = 50;
  params.bufferize_timers_output = false;
}
