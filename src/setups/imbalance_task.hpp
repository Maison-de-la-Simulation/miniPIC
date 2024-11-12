/* _____________________________________________________________________ */
//! \file imbalance_task.hpp

//! \brief This is the imbalance setup used to initialize MiniPic
//! CASE: homogeneous uniform thermalized plasma with periodic boundary conditions with a artifical
//! physical operator

/* _____________________________________________________________________ */

#include "Params.hpp"

//! \brief Functiun to setup all input parameters
void setup(Params &params) {

  // Simulation name
  params.name = "Imbalance";

  // Space
  params.inf_x = 0.;
  params.inf_y = 0.;
  params.inf_z = 0.;
  params.sup_x = 1.;
  params.sup_y = 1.;
  params.sup_z = 1.;

  // Decomp
  params.n_subdomains = 1;

  //  Bin number
#if defined(__MINIPIC_OMP_TASK__) || defined(__MINIPIC_EVENTIFY__)
  params.bin_size = 16000;
#endif

  // int nx_cells = 32;
  // int ny_cells = 32;
  // int nz_cells = 32;

  int nx_cells = 128;
  int ny_cells = 128;
  int nz_cells = 128;

  const double dx = (params.sup_x - params.inf_x) / nx_cells;
  const double dy = (params.sup_y - params.inf_y) / ny_cells;
  const double dz = (params.sup_z - params.inf_z) / nz_cells;

  // Number of patches
  //params.nx_patch = 16;
  //params.ny_patch = 16;
  //params.nz_patch = 16;

  params.nx_patch = 16;
  params.ny_patch = 16;
  params.nz_patch = 16;

  // Cells per patch per direction
  params.nx_cells_by_patch = nx_cells / params.nx_patch;
  params.ny_cells_by_patch = ny_cells / params.ny_patch;
  params.nz_cells_by_patch = nz_cells / params.nz_patch;

  // Time

  params.dt = 0.9; // Use the CFL to always adjust the dt to the space delta

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
  // params.add_species("electron", 1, -1, 1e-2, profile, {0, 0, 0}, 2, "random");
  // params.add_species("proton", 1836.125, 1, 1e-2, profile, {0, 0, 0}, 2, "electron");

  params.add_species("electron", 1, -1, 1e-2, profile, {0, 0, 0}, 16, "random", "cell");
  params.add_species("proton", 1836.125, 1, 1e-2, profile, {0, 0, 0}, 16, "electron", "cell");

  // Momentum correction at init
  params.momentum_correction = false;

  // Bourndary conditions
  params.boundary_condition = "periodic";
  params.no_diagnostics_at_init = true;

  // Display
  params.print_period = 50000;

  // Random seed
  params.seed = 0;

  // Scalar Diagnostics
  params.scalar_diagnostics_period = 10000;

  // Field Diagnostics
  params.field_diagnostics_period = 50000;
  params.field_diagnostics_format = "vtk";

  // Cloud Diagnostics
  params.particle_cloud_period = 50000;
  params.particle_cloud_format = "binary";

  // Binning Diagnostics
  /*params.add_particle_binning("diag_w_gamma",
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
                              "vtk");*/

  // Add imbalance operator
  auto func_weight = [](double x, double y, double z, double t) -> double {
    // Intensity
    double E0  = 10;
    double res = 0;

    double r = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5);
    res = E0 * std::exp(-r / 0.1);
    // Derivative of a E field of the form E0 * exp(-alpha * t^2) * cos(beta * t)

    return res;
  };

  params.add_imbalance(func_weight);
}