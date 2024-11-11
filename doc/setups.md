# Setups

Input parameters to run a miniPIC simulation must be provided before compiling via a specific header containing the `setup` function. the setup function then is used to initiate global parameters in the `Params` class.

A list of benchmarks can be found in the directory `src/setups/`.

## List of parameters

- `params.name` (`std::string`): name of the setup
- `params.inf_x`, `params.inf_y`, `params.inf_z` (`double`): lower bounds of the domain
- `params.sup_x`, `params.sup_y`, `params.sup_z` (`double`): upper bounds of the domain
- `params.nx_patch`, `params.ny_patch`, `params.nz_patch` (`int`): number of patches in each direction
- `params.nx_cells_by_patch`, `params.ny_cells_by_patch`, `params.nz_cells_by_patch` (`int`): number of cells per patch in each direction
- `params.dt` (`double`): time step
- `params.simulation_time` (`double`): simulation time (the total number of time steps is automatically determined `simulation_time / dt`)

## Add a plasma

Species can be added using the `add_species` method:

- `params.add_species`(`std::string name`, `double mass`, `int charge`, `double temperature`, `std::function<double(double xn, double yn, double zn)> profile`, `std::array<double, 3> drift_velocity`, `int number_of_particles_per_cell`, `std::string position_initialization`)

    - `name` (`std::string`): name of the species
    - `mass` (`double`): mass of the species
    - `charge` (`int`): charge of the species
    - `temperature` (`double`): temperature of the species (for the Maxwellian distribution initialization)
    - `profile` (`std::function<double(double xn, double yn, double zn)>`): density profile of the species, `xn`, `yn`, `zn` are the normalized coordinates (between 0 and 1)
    - `drift_velocity` (`std::array<double, 3>`): drift velocity of the species
    - `number_of_particles` (`int`): number of particles per cell for initialization
    - `position_initialization` (`std::string`): initialization of the position of the particles (`random` or a species name). If `random`, the particles are initialized randomly in the cell. If a species name, the particles are initialized using the position of the specified species.
    - `params.boundary_condition` (`std::string`): boundary conditions
    - `params.print_period` (`int`): period for printing the state of the simulation
    - `params.seed` (`int`): random seed for the initialization of the particles

```C++
// Example of added species

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
  params.add_species("electron", 1, -1, temperature, profile, {v_drift, 0, 0}, 8, "random");
  params.add_species("proton", 1836.125, 1, temperature, profile, {v_drift, 0, 0}, 8, "electron");
```


## Add diags

List of available diags:


| Diags                 | Description                         | Parameters  | Format      |
|-----------------------|-------------------------------------| ----------- | ----------- |
| Scalars               | Scalar diags are list of reduced values for a given timestep, for instance global kinetic energy for each species, global field energy per component and more | iteration period | ascii |
| Fields                | Output of the global field grids (for each component Ex, Ey, Ez, Bx, By, Bz) for a given timestep |  iteration period | custom binary or vtk |
| Particle Clouds       | Output of all particles (`x`, `y`, `z`, `px`, `py`, `pz`, `gamma`)                     | iteration period | custom binary or vtk |
| Particle binnings     | Projection of the particle properties on a 1D, 2D or 3D grid | grid dimension, list of parameters, iteration period, starting iteration | custom binary or vtk |


For scalars:

- `params.scalar_diagnostics_period` (`int`): period for scalar diags (if 0, no scalar diags)

For fields:

- `params.field_diagnostics_period` (`int`): period for field diags (if 0, no field diags)
- `params.field_diagnostics_format` (`std::string`): format for field diags (binary or vtk)

For particle clouds:

- `params.particle_cloud_period` (`int`): period for particle cloud diags (if 0, no particle cloud diags)
- `params.particle_cloud_format` (`std::string`): format for particle cloud diags (binary or vtk)

For particle binning:

Binning diagnostics are defined using the `add_particle_binning` method:

- `params.add_particle_binning`(`std::string name`, `std::string projected_parameter`, `std::vector<std::string> axis_parameters`, `std::vector<int> grid_dimensions`, `std::vector<double> inf`, `std::vector<double> sup`, `std::vector<int> periodic`, `int iteration_period`, `std::string format`)

    - `name` (`std::string`): name of the binning
    - `projected_parameter` (`std::string`): projected parameter of the binning (`density` or `weight`)
    - `axis_parameters` (`std::vector<std::string>`): list of parameters for axis
    - `grid_dimensions` (`std::vector<int>`): dimensions of the grid
    - `inf` (`std::vector<double>`): lower bounds of the grid
    - `sup` (`std::vector<double>`): upper bounds of the grid
    - `periodic` (`std::vector<int>`): periodicity of the grid
    - `iteration_period` (`int`): period for the binning diags
    - `format` (`std::string`): format for the binning diags (binary or vtk)


```C++
//1D diag binning example
params.add_particle_binning("diag_w_gamma",
                            "weight",
                            {"gamma"},
                            {32},
                            {0},
                            {0},
                            {0, 1},
                            50,
                            "binary");
```

```C++
//3D diag binning example
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
```