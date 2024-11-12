/* _____________________________________________________________________ */
//! \file default.hpp

//! \brief This is the default setup used to initialize MiniPic
//! CASE: homogeneous uniform thermalized plasma with periodic boundary conditions

/* _____________________________________________________________________ */

#include "Params.hpp"

//! \brief Functiun to setup all input parameters
void setup(Params &params) {

  // Simulation name
  params.name = "Thermalized plasma";

  // Space
  params.inf_x = 0.;
  params.inf_y = 0.;
  params.inf_z = 0.;
  params.sup_x = 1.;
  params.sup_y = 1.;
  params.sup_z = 1.;

  // Decomp
  params.n_subdomains = 1;

  int nx_cells = 32;
  int ny_cells = 32;
  int nz_cells = 32;

  // const double dx = (params.sup_x - params.inf_x) / nx_cells;
  // const double dy = (params.sup_y - params.inf_y) / ny_cells;
  // const double dz = (params.sup_z - params.inf_z) / nz_cells;

  // Number of patches
  params.nx_patch = 2;
  params.ny_patch = 2;
  params.nz_patch = 2;

  // Cells per patch per direction
  params.nx_cells_by_patch = nx_cells / params.nx_patch;
  params.ny_cells_by_patch = ny_cells / params.ny_patch;
  params.nz_cells_by_patch = nz_cells / params.nz_patch;

  // Time

  params.dt = 0.9;

  params.simulation_time = 300 * params.dt;

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

  // Timers
  params.save_timers_period = 50;

  // Scalar Diagnostics
  params.scalar_diagnostics_period = 10;

  // Field Diagnostics
  params.field_diagnostics_period = 50;
  params.field_diagnostics_format = "vtk";

  // Cloud Diagnostics
  params.particle_cloud_period = 50;
  params.particle_cloud_format = "binary";

  // Binning Diagnostics
  params.add_particle_binning("diag_w_gamma",
                              "weight",
                              {"gamma"},
                              {32},
                              {0},
                              {0},
                              {0, 1},
                              50,
                              "binary");

  params.add_particle_binning("diag_x_y_z_d",
                              "density",
                              {"x", "y", "z"},
                              {params.nx_patch * params.nx_cells_by_patch,
                               params.ny_patch * params.ny_cells_by_patch,
                               params.nz_patch * params.nz_cells_by_patch},
                              {0., 0., 0.},
                              {params.sup_x, params.sup_y, params.sup_z},
                              {0, 1},
                              50,
                              "vtk");

  params.add_particle_binning("diag_px_py_pz_d",
                              "density",
                              {"px", "py", "pz"},
                              {params.nx_patch * params.nx_cells_by_patch,
                               params.ny_patch * params.ny_cells_by_patch,
                               params.nz_patch * params.nz_cells_by_patch},
                              {0., 0., 0.},
                              {0, 0, 0},
                              {0, 1},
                              50,
                              "vtk");
}
