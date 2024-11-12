/* _____________________________________________________________________ */
//! \file ppc_study.hpp

//! \brief Used for the parametric study focusing on the number of particles per cell
//! CASE: homogeneous uniform thermalized plasma with periodic boundary conditions

/* _____________________________________________________________________ */

#ifndef SETUP
#define SETUP

#include "Params.hpp"

//! \brief Functiun to setup all input parameters
void setup(Params &params) {

  // Space
  params.inf_x = 0.;
  params.inf_y = 0.;
  params.inf_z = 0.;
  params.sup_x = 1.;
  params.sup_y = 1.;
  params.sup_z = 1.;

  // Decomp
  params.n_subdomains = 1;

  // Number of patches
  params.nx_patch = 1;
  params.ny_patch = 1;
  params.nz_patch = 1;

  // Cells per patch per direction
  params.nx_cells_by_patch = 64;
  params.ny_cells_by_patch = 64;
  params.nz_cells_by_patch = 64;

  // Time

  const double dx = (params.sup_x - params.inf_x) / (params.nx_cells_by_patch * params.nx_patch);
  const double dy = (params.sup_y - params.inf_y) / (params.ny_cells_by_patch * params.ny_patch);
  const double dz = (params.sup_z - params.inf_z) / (params.nz_cells_by_patch * params.nz_patch);

  params.dt = 0.9 ; // Fraction of the CFL condition

  params.simulation_time = 100 * params.dt;

  // Species

  // custom density profile
  auto profile = [](double x, double y, double z) -> double {
    // if ((x > 0.2) && (x < 0.8)) {
    return 1e-5;
    // } else {
    //   return 0;
    // }
  };

  // name, mass, charge, temperature, density profile, drift velocity, particles per cell,
  // position initialization
  params.add_species("electron", 1, -1, 1e-2, profile, {0, 0, 0}, 8, "random", "cell");
  params.add_species("proton", 1836.125, 1, 1e-2, profile, {0, 0, 0}, 8, "electron", "cell");

  // Momentum correction at init
  params.momentum_correction = false;

  // Bourndary conditions
  params.boundary_condition = "periodic";

  // Display
  params.print_period = 50;

  // Random seed
  params.seed = 0;

  // No diags at init
  params.no_diagnostics_at_init = true;

  // Scalar Diagnostics
  params.scalar_diagnostics_period = 0;

  // Field Diagnostics
  params.field_diagnostics_period = 0;
  params.field_diagnostics_format = "vtk";

  // Cloud Diagnostics
  params.particle_cloud_period = 0;
  params.particle_cloud_format = "binary";
}

#endif // SETUP
