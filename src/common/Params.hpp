/* _____________________________________________________________________ */
//! \file Params.hpp

//! \brief description of the class Params that contains all
//! global simulation parameters

/* _____________________________________________________________________ */

#pragma once

#include "Headers.hpp"
#include "Random.hpp"
#include "Tools.hpp"
#include <chrono>
#include <csignal>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <math.h>
#include <vector>

// ___________________________________________________________________
//
//! structure to store the properties of a particle binning diagnostics
// ___________________________________________________________________
struct ParticleBinningProperties {
  //! Name of the diagnostic (used for output)
  std::string name_;
  //! Name of the parameter to project
  std::string projected_parameter_;
  //! vector of axis names
  std::vector<std::string> axis_;
  //! Number of cells for each axis
  std::vector<int> n_cells_;
  //! min values for each axis
  std::vector<double> min_;
  //! max values for each axis
  std::vector<double> max_;
  //! Species number
  std::vector<int> species_indexes_;
  //! Output period
  int period_;
  //! format of the output file
  std::string format_;

  //! Constructor
  ParticleBinningProperties(std::string name,
                            std::string projected_parameter,
                            std::vector<std::string> axis,
                            std::vector<int> n_cells,
                            std::vector<double> min,
                            std::vector<double> max,
                            std::vector<int> species_indexes,
                            int period,
                            std::string format)
    : name_(name), projected_parameter_(projected_parameter), axis_(axis), n_cells_(n_cells),
      min_(min), max_(max), species_indexes_(species_indexes), period_(period), format_(format){};
};

// ___________________________________________________________________
//
//! structure to store a single particle
// ___________________________________________________________________
struct Particle {

  unsigned int is_;
  double w_;
  double x_;
  double y_;
  double z_;
  double mx_;
  double my_;
  double mz_;

  //! Constructor
  Particle(unsigned int is, double w, double x, double y, double z, double mx, double my, double mz)
    : is_(is), w_(w), x_(x), y_(y), z_(z), mx_(mx), my_(my), mz_(mz){};
};

// ___________________________________________________________________
//
//! Structure to store all input parameters and those that can be computed
// ___________________________________________________________________
class Params {
public:
  Params() {}
  ~Params() {}

  // Simulation name
  std::string name;

  //      Discretisation param
  // Space
  //! Domain boundaries
  double inf_x, inf_y, inf_z;
  double sup_x, sup_y, sup_z;

  //! Domain decomposition
  int n_subdomains;
  int nx_patch, ny_patch, nz_patch;
  int nx_cells_by_patch, ny_cells_by_patch, nz_cells_by_patch;

#if defined(__MINIPIC_OMP_TASK__) || defined(__MINIPIC_EVENTIFY__)
  //! Particles bin numbers per patch
  int bin_size = 32000;
#endif

  //! Input computed
  int N_patches;
  int nx_cells, ny_cells, nz_cells;
  double Lx, Ly, Lz;
  double dx, dy, dz;
  double inv_dx, inv_dy, inv_dz;
  double dx_sq, dy_sq, dz_sq;
  double cell_volume, inv_cell_volume;

  int nx_p, ny_p, nz_p;
  int nx_d, ny_d, nz_d;

  //! Number of primal cells per patch
  int nx_p_by_patch, ny_p_by_patch, nz_p_by_patch;

  // Number of dual cells per patch
  int nx_d_by_patch, ny_d_by_patch, nz_d_by_patch;

  //! String to store the name of the boundary condition
  std::string boundary_condition;

  //! Code to shortcut the type of boundary conditions
  //! 0 - free
  //! 1 - periodic
  //! 2 - reflective
  int boundary_condition_code;

  //! Iteration period for information display in the terminal
  int print_period;

  // Time
  //! Start, stop and step for the simulation time
  double dt;
  //! Time total, computed
  double simulation_time;
  //! Number of iteration, computed
  int n_it;
  //! CFL
  double dt_cfl;

  // ____________________________________________________________
  // Species parameters

  //! Name of the species
  std::vector<std::string> species_names_;
  //! Particle per cell for each species at init
  std::vector<int> ppc_;
  //! Normalized density, temperature, mass, charge for each species at init
  std::vector<double> temp_, mass_, charge_;
  //! Density profile for each species at init
  //! Lambda of 3 doubles: x, y, z that representes a normalized position (between 0 and 1)
  std::vector<std::function<double(double, double, double)>> density_profiles_;
  //! Mean velocity for each species at init (vector of 3 doubles)
  std::vector<std::vector<double>> drift_velocity_;
  //! Method to use for position init
  std::vector<std::string> position_initialization_method_;
  //! Position init level (patch or cell)
  std::vector<std::string> position_initialization_level_;

  //! Number of particles in one patch for each species at init, computed
  std::vector<int> n_particles_by_species_in_patch;
  //! Number of particles total for each species at init, computed
  std::vector<int> n_particles_by_species;
  //! Number of particles at init, computed
  int n_particles;

  //! list of particles to add at init
  std::vector<Particle> particles_to_add_;

  // ____________________________________________________________
  // Antenna parameters

  //! Antenna profile (x, y, t)
  std::vector<std::function<double(double, double, double)>> antenna_profiles_;
  //! Antenna position
  std::vector<double> antenna_positions_;

  // ____________________________________________________________
  // Imbalance parameters

  //! imbalance profile (x, y, z)
  std::vector<std::function<double(double, double, double, double)>> imbalance_function_;

  //! Random number
  int seed = 0;

  //! Number of threads
  int number_of_threads;

  // ____________________________________________________________
  // Initial Electric and magnetic field values

  double E0_[3] = {0, 0, 0};
  double B0_[3] = {0, 0, 0};

  // ____________________________________________________________
  // Operators

  //! Momentum correction
  bool momentum_correction = false;

  //! Projector
  bool current_projection = true;

  //! Projector
  bool maxwell_solver = true;

  // ____________________________________________________________
  // Parallelism

  bool on_gpu_ = false;

  // ____________________________________________________________
  // Diagnostics parameters

  //! If true, no diagnostics at iteration 0
  bool no_diagnostics_at_init = false;

  //! Output directory relative path
  std::string output_directory = "diags";

  //! vector of ParticleBinningProperties
  std::vector<ParticleBinningProperties> particle_binning_properties_;

  //! Period of the particle cloud diagnostic
  unsigned int particle_cloud_period = 0;

  //! Format of the particle cloud diagnostic
  std::string particle_cloud_format = "binary";

  //! Period of the field diagnostic
  unsigned int field_diagnostics_period = 0;

  //! Format of the field diagnostic
  std::string field_diagnostics_format = "vtk";

  //! Period for scalar diagnostics
  unsigned int scalar_diagnostics_period = 0;

  //! Digits for iteration in output file names
  unsigned int max_it_digits = 0;

  //! Period for timers saving
  unsigned int save_timers_period = 0;

  //! start for timers saving
  unsigned int save_timers_start = 0;

  //! bufferize the output
  bool bufferize_timers_output = true;

  // _________________________________________
  // Methods

  // _______________________________________________________________
  //
  //! \brief Compute all global parameters using the user input
  //! Check that the input parameters are coherent
  // _______________________________________________________________
  void compute() {

    N_patches = nx_patch * ny_patch * nz_patch;

    nx_cells = nx_cells_by_patch * nx_patch;
    ny_cells = ny_cells_by_patch * ny_patch;
    nz_cells = nz_cells_by_patch * nz_patch;

    Lx              = sup_x - inf_x;
    Ly              = sup_y - inf_y;
    Lz              = sup_z - inf_z;

    dx              = Lx / nx_cells;
    dy              = Ly / ny_cells;
    dz              = Lz / nz_cells;

    inv_dx          = 1. / dx;
    inv_dy          = 1. / dy;
    inv_dz          = 1. / dz;

    dx_sq           = dx * dx;
    dy_sq           = dy * dy;
    dz_sq           = dz * dz;

    cell_volume     = dx * dy * dz;
    inv_cell_volume = inv_dx * inv_dy * inv_dz;

    nx_p = nx_cells + 1;
    ny_p = ny_cells + 1;
    nz_p = nz_cells + 1;
    nx_d = nx_cells + 2;
    ny_d = ny_cells + 2;
    nz_d = nz_cells + 2;

    nx_p_by_patch = nx_cells_by_patch + 1;
    ny_p_by_patch = ny_cells_by_patch + 1;
    nz_p_by_patch = nz_cells_by_patch + 1;
    nx_d_by_patch = nx_cells_by_patch + 2;
    ny_d_by_patch = ny_cells_by_patch + 2;
    nz_d_by_patch = nz_cells_by_patch + 2;

    // boundary conditions

    if (boundary_condition == "periodic") {
      boundary_condition_code = 1;
    } else if (boundary_condition == "reflective") {
      boundary_condition_code = 2;
    } else {
      ERROR(" Boundary condition " << boundary_condition << " is not supported");
      std::raise(SIGABRT);
    }

    // Time

    // Compute the exacte CFL condition
    dt_cfl = std::sqrt(1 / (1 / dx_sq + 1 / dy_sq + 1 / dz_sq));

    if (dt > 1 || dt <= 0) {
      ERROR("ERROR in setup: dt (fraction of the CFL) must be between 0 and 1 to comply with the CFL condition")
      std::raise(SIGABRT);
    }

    // Get the number of time steps
    n_it = static_cast<int>(std::round(simulation_time / dt));

    // Convert dt into a fraction of the CFL
    dt = dt * dt_cfl;
    simulation_time = n_it * dt;

    // const double cfl = (dt * dt / (1 / dx_sq + 1 / dy_sq + 1 / dz_sq));
    // if (cfl > 1) {
    //   std::cerr << " CFL condition is not respected, you must have : 1/dt**2 <= 1/dx**2 + 1/dy**2 "
    //                "+ 1/dz**2 "
    //             << std::endl;
    //   std::raise(SIGABRT);
    // }

    // Physics
    n_particles = 0;
    n_particles_by_species_in_patch.resize(species_names_.size());
    n_particles_by_species.resize(species_names_.size());
    for (int is = 0; is < species_names_.size(); is++) {
      n_particles_by_species_in_patch[is] =
        nx_cells_by_patch * ny_cells_by_patch * nz_cells_by_patch * ppc_[is];
      n_particles_by_species[is] =
        nx_patch * ny_patch * nz_patch * n_particles_by_species_in_patch[is];
      n_particles += n_particles_by_species[is];
    }

    // Check species initialization
    for (auto is = 0; is < species_names_.size(); ++is) {

      bool passed = false;

      if (position_initialization_method_[is] == "random") {

        passed = true;

      } else {

        // We check that the position init is one of the previous species
        for (auto is2 = 0; is2 < is; ++is2) {
          if (position_initialization_method_[is] == species_names_[is2]) {
            passed = true;
          }
        }
      }
      // if not passed, return an error
      if (!passed) {
        ERROR(" Position initialization " << position_initialization_method_[is]
                  << " is not supported");
        std::raise(SIGABRT);
      }
    }

    // Check species init level (should be "cell" or "patch")
    for (auto is = 0; is < species_names_.size(); ++is) {
      if (position_initialization_level_[is] != "cell" &&
          position_initialization_level_[is] != "patch") {
        ERROR(" Position initialization level " << position_initialization_level_[is]
                  << " is not supported");
        std::raise(SIGABRT);
      }
    }

    // _________________________________________
    // Diagnostics

    // get number of digit for the number of iterations
    // Used for diagnostics names
    max_it_digits = 0;
    int n_it_tmp  = n_it;
    while (n_it_tmp > 0) {
      n_it_tmp /= 10;
      max_it_digits++;
    }

    // if the period is 0, then no outputs
    if (particle_cloud_period == 0) {
      particle_cloud_period = n_it + 1;
    }

    if (field_diagnostics_period == 0) {
      field_diagnostics_period = n_it + 1;
    }

    if (scalar_diagnostics_period == 0) {
      scalar_diagnostics_period = n_it + 1;
    }

    if (save_timers_period == 0) {
      save_timers_period = n_it + 1;
    }

    // _________________________________________
    // Parallelism
#if defined(__MINIPIC_OMP__) || defined(__MINIPIC_OMP_TASK__) || defined(__MINIPIC_EVENTIFY__)
    number_of_threads = omp_get_max_threads();
#else
    number_of_threads = 1;
#endif
  }

  // _________________________________________________________________________________________________
  //! \brief Add an antenna
  //! \param antenna_profile Lambda std::function<double(double, double, double)> representing the
  //! antenna profile
  //! \param x Position of the antenna
  // _________________________________________________________________________________________________
  void add_antenna(std::function<double(double, double, double)> antenna_profile, double x);

  // _________________________________________________________________________________________________
  //! \brief Add imbalance function
  //! \param function_profile function
  // _________________________________________________________________________________________________
  void add_imbalance(std::function<double(double, double, double, double)> function_profile);

  // _________________________________________________________________________________________________
  //! \brief Add a species
  //! \param name Name of the species
  //! \param mass Mass of the species
  //! \param charge Charge of the species
  //! \param temp Temperature of the species
  //! \param density_profile Lambda std::function<double(double, double, double)> representing the
  //! density profile of the species \param drift_velocity Drift velocity of the species \param ppc
  //! Number of particles per cell of the species \param position_initiatization Method to use for
  //! position init
  // _________________________________________________________________________________________________
  void add_species(std::string name,
                   double mass,
                   double charge,
                   double temp,
                   std::function<double(double, double, double)> density_profile,
                   std::vector<double> drift_velocity,
                   double ppc,
                   std::string position_initiatization,
                   std::string position_initialization_level);

  // _____________________________________________________
  //
  //! \brief Add a particle binning object (diagnostic) to the vector
  //! of particle binning properties
  //! \param[in] diag_name - string used to initiate the diag name
  //! \param[in] projected_parameter - property projected on the grid, can be
  //! `weight`
  //! \param[in] axis - axis to use for the grid. The axis vector size
  //! determines the dimension of the diag. Axis can be `gamma`, `weight`, `x`,
  //! `y`, `z`, `px`, `py`, `pz`
  //! \param[in] n_cells - number of cells for each axis
  //! \param[in] min - min value for each axis
  //! \param[in] max - max value for each axis
  //! \param[in] is - species
  //! \param[in] period - output period
  //! \param[in] format - (optional argument) - determine the output format, can
  //! be `binary` (default) or `vtk`
  // _____________________________________________________
  void add_particle_binning(std::string diag_name,
                            std::string projected_parameter,
                            std::vector<std::string> axis,
                            std::vector<int> n_cells,
                            std::vector<double> min,
                            std::vector<double> max,
                            std::vector<int> species_indexes,
                            int period,
                            std::string format = "binary");

  // _____________________________________________________
  //
  //! \brief Add a single particle to the corresponding species
  //! \param[in] is species index
  //! \param[in] w  particle weight
  //! \param[in] x  particle weight
  //! \param[in] y  particle weight
  //! \param[in] z  particle weight
  //! \param[in] mx  particle weight
  //! \param[in] my  particle weight
  //! \param[in] mz  particle weight
  // _____________________________________________________
  void add_particle(unsigned int is,
                    double w,
                    double x,
                    double y,
                    double z,
                    double mx,
                    double my,
                    double mz);

  // _____________________________________________________
  //
  //! \brief Set the intial values of the electric field
  //! \param[in] Ex initial value for Ex
  //! \param[in] Ez initial value for Ey
  //! \param[in] Ez initial value for Ez
  // _____________________________________________________
  void initialize_electric_field(double Ex, double Ey, double Ez);

  // _____________________________________________________
  //
  //! \brief Set the intial values of the magnetic field
  //! \param[in] Bx initial value for Bx
  //! \param[in] By initial value for By
  //! \param[in] Bz initial value for Bz
  // _____________________________________________________
  void initialize_magnetic_field(double Bx, double By, double Bz);

  // _____________________________________________________
  //
  //! \brief return the number of species
  // _____________________________________________________
  inline auto get_species_number() const { return species_names_.size(); }

  // _____________________________________________________
  //
  //! \brief return the number of particle binning objects
  // _____________________________________________________
  inline auto get_particle_binning_number() const { return particle_binning_properties_.size(); }

  // _____________________________________________________
  //
  //! \brief print the help for command line options
  // _____________________________________________________
  void help() const {
    std::cout << " \n"
              << " Help for command line options \n"
              << " Note: command line options overwrite program parameters.\n\n"
              << " -h   (--help): print the help page for command line options\n"
              << " -it  (--iterations) int: change the number of iterations\n"
              << " -dmin (--domain_min) double double double: change the domain minimum boundaries\n"
              << " -dmax (--domain_max) double double double: change the domain maximum boundaries\n"
              << " -p   (--patches) int int int: number of patches per direction\n"
              << " -cpp (--cells_per_patch) int int int: number of cells per patch per direction\n"
              << " -rs  (--random_seed) int: seed for random generator\n"
              << " -pp  (--print_period) int: iteration period for terminal printing\n"
              << " -stp (--save_timers_period) int: iteration period for timers saving\n"
              << " -sts (--save_timers_start) int: iteration start for timers saving\n"
              << std::endl;
    std::_Exit(EXIT_SUCCESS);
  }

  // _____________________________________________________
  //
  //! \brief Simple parser to read some input parameters
  //! from command line arguments
  // _____________________________________________________
  void read_from_command_line_arguments(int argc, char *argv[]);

  // _____________________________________________________
  //
  //! \brief Create a line composed of `size` characters
  //! \param[in] size - number of characters for a single line
  // _____________________________________________________
  std::string seperator(int size) const {
    std::string line = " ";
    for (int i = 0; i < size; ++i) {
      line += "_";
    }
    return line;
  }

  // _____________________________________________________
  //
  //! \brief print the title
  // _____________________________________________________
  void title();

  // _____________________________________________________
  //
  //! \brief Print input parameters summary
  // _____________________________________________________
  void info();

  // _____________________________________________________________
  //
  //! \brief Give the index patch topology from 3D indexes,
  //!       -1 if out of domain
  // _____________________________________________________________
  INLINE int get_patch_index(int i, int j, int k) {

    int ixp = i;
    int iyp = j;
    int izp = k;

    // Periodic management of the topology
    if (ixp < 0) {
      ixp = nx_patch - 1;
    } else if (ixp >= nx_patch) {
      ixp = 0;
    }

    if (iyp < 0) {
      iyp = ny_patch - 1;
    } else if (iyp >= ny_patch) {
      iyp = 0;
    }

    if (izp < 0) {
      izp = nz_patch - 1;
    } else if (izp >= nz_patch) {
      izp = 0;
    }

    // Reflective management of the topology
    // if(i<0 || i>=nx_patchs_m)
    //   return -1;
    // if(j<0 || j>=ny_patchs_m)
    //   return -1;
    // if(k<0 || k>=nz_patchs_m)
    //   return -1;

    return ixp * nz_patch * ny_patch + iyp * nz_patch + izp;
  }
};
