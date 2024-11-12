/* _____________________________________________________________________ */
//! \file default.hpp

//! \brief This is the default setup used to initialize MiniPic
//! CASE: homogeneous uniform thermalized plasma with periodic boundary conditions

/* _____________________________________________________________________ */

#include "Params.hpp"

//! \brief Functiun to setup all input parameters
void setup(Params &params) {

  // Simulation name
  params.name = "Constant magnetic field";

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
  params.nx_patch = 4;
  params.ny_patch = 4;
  params.nz_patch = 4;

  // Cells per patch per direction
  params.nx_cells_by_patch = 8;
  params.ny_cells_by_patch = 8;
  params.nz_cells_by_patch = 8;

  // Time

  // const double dx = (params.sup_x - params.inf_x) / (params.nx_cells_by_patch * params.nx_patch);
  // const double dy = (params.sup_y - params.inf_y) / (params.ny_cells_by_patch * params.ny_patch);
  // const double dz = (params.sup_z - params.inf_z) / (params.nz_cells_by_patch * params.nz_patch);

  params.dt = 0.9;

  params.simulation_time = 50000 * params.dt;

  // Species

  // custom density profile
  auto profile = [](double x, double y, double z) -> double {
    // if ((x > 0.2) && (x < 0.8)) {
    return 1e-5;
    // } else {
    //   return 0;
    // }
  };

  // name, mass, charge, density, temperature, density profile, drift velocity, particles per cell,
  // position initialization
  params.add_species("electron", 1, -1, 1e-2, profile, {0, 0, 0}, 0, "random", "cell");

  // Add a single particle
  params.add_particle(0, 1e-3, 0.99, 0.51, 0.51, 0., 4., 0.);

  // Initiale electric field values
  params.initialize_magnetic_field(0, 0, 9);

  // Bourndary conditions
  params.boundary_condition = "periodic";

  // Display
  params.print_period = 10000;

  // Random seed
  params.seed = 0;

  // Field Diagnostics
  params.particle_cloud_period = 1000;

  // Scalar Diagnostics
  params.scalar_diagnostics_period = 100;

  // Field Diagnostics
  params.field_diagnostics_period = 0;

  // Momentum correction at init
  params.momentum_correction = true;

  // Current projection
  params.current_projection = false;

  // Maxwell solver
  params.maxwell_solver = false;
}
