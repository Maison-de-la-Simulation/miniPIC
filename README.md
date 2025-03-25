

<img title="Title" alt="title" src="./doc/images/title.png" height="100">

## Presentation

miniPIC is a playground for computer science and HPC experiments applied to the Particle-In-Cell method.

> [!WARNING]
> miniPIC is not a code intended to be used for numerical simulation of physical cases.

### CPU backends

<table>
    <tr>
        <th></th>
        <th style="background-color: #CDD8E0; color: black">Internal Backend</th>
        <th style="background-color: #CDD8E0; color: black">CPU Intel</th>
        <th style="background-color: #CDD8E0; color: black">CPU AMD</th>
        <th style="background-color: #CDD8E0; color: black">ARM based CPU</th>
    </tr>
    <tr>
        <td>OpenMP loop</td>
        <td></td>
        <td style="background-color: #B1E0CB">LLVM, GCC, OneAPI</td>
        <td style="background-color: #B1E0CB">LLVM, GCC, OneAPI</td>
        <td style="background-color: #B1E0CB">LLVM 14, GCC 10</td>
    </tr>
    <tr>
        <td>SYCL</td>
        <td></td>
        <td style="background-color: #B1E0CB">LLVM, GCC, OneAPI</td>
        <td style="background-color: #B1E0CB">OneAPI</td>
        <td style="background-color: #E0B4B2">Not tested</td>
    </tr>
    <tr>
        <td>Kokkos</td>
        <td>OpenMP</td>
        <td style="background-color: #B1E0CB">LLVM, GCC, OneAPI</td>
        <td style="background-color: #E0B4B2">Not tested</td>
        <td style="background-color: #B1E0CB">LLVM, GCC</td>
    </tr>
    <tr>
        <td>Eventify</td>
        <td></td>
        <td style="background-color: #B1E0CB"></td>
        <td style="background-color: #B1E0CB"></td>
        <td style="background-color: #E0B4B2">Not tested</td>
    </tr>
    <tr>
        <td>OpenMP task</td>
        <td></td>
        <td style="background-color: #B1E0CB"></td>
        <td style="background-color: #B1E0CB"></td>
        <td style="background-color: #E0B4B2">Not tested</td>
    </tr>
</table>

### GPU backends

<table>
    <tr>
        <th></th>
        <th style="background-color: #CDD8E0; color: black">Internal Backend</th>
        <th style="background-color: #CDD8E0; color: black">GPU NVIDIA</th>
        <th style="background-color: #CDD8E0; color: black">GPU AMD</th>
        <th style="background-color: #CDD8E0; color: black">GPU Intel</th>
    </tr>
    <tr>
        <td>Kokkos</td>
        <td>CUDA, HIP, SYCL</td>
        <td style="background-color: #B1E0CB">CUDA 11 (tested on V100, A100)</td>
        <td style="background-color: #B1E0CB">HIP (tested on MI250)</td>
        <td style="background-color: #B1E0CB">Tested on Intel MAX 1550</td>
    </tr>
    <tr>
        <td>SYCL</td>
        <td></td>
        <td style="background-color: #B1E0CB">tested on V100, A100</td>
        <td style="background-color: #B1E0CB">Tested on MI250</td>
        <td style="background-color: #B1E0CB">Tested on Intel MAX 1550</td>
    </tr>
    <tr>
        <td>THRUST</td>
        <td></td>
        <td style="background-color: #B1E0CB">tested on V100, A100 (nvcc 11)</td>
        <td style="background-color: #B1E0CB">Tested on MI250</td>
        <td style="background-color: #E0B4B2">Not tested</td>
    </tr>
    <tr>
        <td>OpenACC</td>
        <td></td>
        <td style="background-color: #B1E0CB">tested on V100, A100 (nvhpc 23.7)</td>
        <td style="background-color: #E0B4B2">Not tested</td>
        <td style="background-color: #E0B4B2">Not tested</td>
    </tr>
    <tr>
        <td>OpenMP target</td>
        <td></td>
        <td style="background-color: #E0B4B2">Tested (nvhpc 23.7)</td>
        <td style="background-color: #E0B4B2">Not tested</td>
        <td style="background-color: #E0B4B2">Not tested</td>
    </tr>  
</table>

## Repository structure

- `doc`: documentation pages
- `src`: C++ source
  - `setups`: headers used to initialize the physical parameters
  - `common`: source files common to all backends
  - backend specific folder (`kokkos`, `thrust`, etc): backend specific operators
- `lib`: Python libraries for miniPIC python tools
- `script`: Python scripts to read and plot diags
- `tests`: simulation tests used to validate the code
- `validation`: Python validation scripts
- `unit_tests`: unit tests


## How to use miniPIC

- [Compilation](./doc/compilation.md)
- [Understand and create your own setups](./doc/setups.md)
- [plot diags](./doc/diags.md)
- [Timers](./doc/timers.md)
- [Python tools](./doc/python_tools.md)

## How to contribute

- [Run the validation tests](./doc/validation.md)
- [Continuous Integration](./doc/ci.md)
- [Code structure](./doc/code_structure.md)

