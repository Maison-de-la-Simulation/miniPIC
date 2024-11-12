/* _____________________________________________________________________ */
//! \file beam.hpp

//! \brief Simulation of a drifting neutral beam

/* _____________________________________________________________________ */

#include "Params.hpp"

//! \brief Functiun to setup all input parameters
void setup(Params &params) {

  // Simulation name
  params.name = "beam";

  // Physics parameters
  const double temperature  = 0.1 / 511.;
  const double n0           = 1;
  const double debye_length = sqrt(temperature / n0);
  const double v_drift      = 0.9;

  const double dx = debye_length / 8;
  const double dy = debye_length / 8;
  const double dz = debye_length / 8;

  // Space: compute the domain size from the number of cells and number of patches
  params.inf_x = 0.;
  params.inf_y = 0.;
  params.inf_z = 0.;
  params.sup_x = 32 * dx;
  params.sup_y = 32 * dy;
  params.sup_z = 32 * dz;

  // Decomp
  params.n_subdomains = 1;

  int nx_cells = static_cast<int>((params.sup_x - params.inf_x) / dx);
  int ny_cells = static_cast<int>((params.sup_y - params.inf_y) / dy);
  int nz_cells = static_cast<int>((params.sup_z - params.inf_z) / dz);

  // Number of patches
  params.nx_patch = 4;
  params.ny_patch = 4;
  params.nz_patch = 4;

  // #else

  //   // Number of patches
  //   params.nx_patch = 4;
  //   params.ny_patch = 4;
  //   params.nz_patch = 4;

  // #endif

  // Cells per patch per direction
  params.nx_cells_by_patch = nx_cells / params.nx_patch;
  params.ny_cells_by_patch = ny_cells / params.ny_patch;
  params.nz_cells_by_patch = nz_cells / params.nz_patch;

  // Time
  params.dt = 0.9;

  params.simulation_time = 500 * params.dt;

  // Species

  // custom density profile
  auto profile = [n0](double x, double y, double z) -> double {
    const double R = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5);
    if (R < 0.25 * 0.25) {
      return n0;
    } else {
      return 0;
    }
  };

  // Plasma
  params.add_species("electron", 1, -1, temperature, profile, {v_drift, 0, 0}, 8, "random", "cell");
  params.add_species("proton",
                     1836.125,
                     1,
                     temperature,
                     profile,
                     {v_drift, 0, 0},
                     8,
                     "electron",
                     "cell");

  // Bourndary conditions
  params.boundary_condition = "periodic";

  // Display
  params.print_period = 50;

  // Random seed
  params.seed = 0;

  // Scalar Diagnostics
  params.scalar_diagnostics_period = 10;

  // Field Diagnostics
  params.field_diagnostics_period = 50;

  // Particle binning
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

  // Timers
  params.save_timers_start       = 10;
  params.save_timers_period      = 100;
  params.bufferize_timers_output = true;
}

