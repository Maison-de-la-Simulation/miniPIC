# Plotting Diags

Some python scripts are provided to plot the diagnostics files.
The scripts are located in the `scripts` folder.

## Available scripts

- `plot_energy_balance.py`: plot the energy balance of the simulation using the scalars

| Option | Long Option | Description |
| --- | --- | --- |
| `-f` | `--folder` | Path to the diags folder |

```bash
# Example
python plot_energy_balance.py -f diags/
```

<img title="field diag example" alt="field diag example" src="./images/energy_balance.svg" height="300">

- `plot_particle_binning.py`: plot particle binning diags (1D, 2D or 3D) only for the custom binary format

| Option | Long Option | Description |
| --- | --- | --- |
| `-f` | `--file` | Path to the file to plot |
| `-c` | `--colormap` | Colormap to use |


```bash
# Example
python plot_particle_binning.py -f diags/diag_w_gamma_s00_0300.bin
```

- `plot_particle_cloud.py`: plot particle cloud using Matplotlib and the binary custom format. We recommend to use Paraview with the vtk format for better results.

| Option | Long Option | Description |
| --- | --- | --- |
| `-f` | `--file` | Path to the file to plot |

- `plot_field.py`: plot a field diag using Matplotlib (2D slice)

| Option | Long Option | Description |
| --- | --- | --- |
| `-f` | `--file` | Path to the file to plot |
| `-c` | `--colormap` | Colormap to use |

```bash
# Example
python plot_field.py -f diags/Ez_200.bin
```

<img title="field diag example" alt="field diag example" src="./images/Ez_200.png" height="300">

- `plot_profiler.py`: plot the profiler diagnostics