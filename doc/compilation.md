# Compilation

## General

miniPIC uses CMake as a build system.

```bash
mkdir build
cd build
cmake ../ 
make
```

<img title="Warning" alt="Warning" src="./doc/images/warning.png" height="20"> Building in the root directory is not supported.

<img title="Warning" alt="Warning" src="./doc/images/warning.png" height="20"> By default, the code is compiled in sequential mode.


## Options

CMake useful options:

- `-DCMAKE_CXX_COMPILER=<compiler choice>`: specify the compiler to use

Backends:

- `-DBACKEND`: enable to choose the backend

| CPU backends      | Description                         |
|-------------------|-------------------------------------|
| sequential        | Sequential CPU version              |
| openmp            | OpenMP CPU version                  |
| stdpar_cpu        | CPU-oriented Stdpar version         |

| CPU TASK backends    | Description                    |
|----------------------|--------------------------------|
| eventify             | Eventify version               |
| openmp_task          | OpenMP task version            |

| GPU backends            | Description                                          |
|-------------------------|------------------------------------------------------|
| kokkos                  | GPU-oriented Kokkos version using dual views         |
| kokkos_dualview_unified | Kokkos dualview using unified memory                 |
| kokkos_unified          | Kokkos using normal views and unified memory         |
| thrust                  | GPU-oriented Thrust version         |
| thrust_unified          | GPU-oriented Thrust unified version |
| sycl                    | SYCL version                        |
| openacc                 | OpenACC version                     |
| openmp_target           | OpenMP target version               |
| stdpar                  | GPU-oriented Stdpar version         |

Tools:

- `-DSHAMAN=ON/OFF`: enable/disable Shaman tool (only available in OpenMP mode or sequential mode)

Others:

- `-DDEBUG=ON/OFF`: enable/disable debug mode (`OFF` by default)
- `-DTEST=ON/OFF`: enable/disable tests mode (for CI, `OFF` by default)
- `-DWARNING=ON/OFF`: enable/disable warnings (`OFF` by default)

- `-DDEVICE`: enable to tune the code for a specific device (required for some backends)

| CPU devices   | Description                         |
|---------------|-------------------------------------|
| nvidia_grace  | Nvidia Grace CPU                    |
| amd_genoa     | AMD Genoa CPU                       |

| GPU devices   | Description                         |
|---------------|-------------------------------------|
| nvidia_v100   | Nvidia V100 GPU                     |
| nvidia_a100   | Nvidia A100 GPU                     |
| nvidia_h100   | Nvidia H100 GPU                     |
| nvidia_gh200  | Nvidia GH200 GPU                    |
| amd_mi250     | AMD MI250 GPU                       |
| amd_mi300     | AMD MI300 GPU                       |
| intel_pvc     | Intel Ponte Vecchio GPU             |

- `-DMINIPIC=ON/OFF`: enable/disable minipic compilation (`ON` by default)
- `-DUNIT_TESTS=ON/OFF`: enable/disable unit tests compilation (`OFF` by default)

## Examples

- Sequential compilation

```bash
cmake ../ 
make
```

- OpenMP compilation using g++

```bash
cmake ../ -DCMAKE_CXX_COMPILER=g++ -DBACKEND=openmp
make
```

- Kokkos compilation using clang++

```bash
cmake ../ -DCMAKE_CXX_COMPILER=clang++ -DBACKEND=kokkos 
make
```

- THRUST compilation using nvcc for Nvidia V100

```bash
cmake ../ -DCMAKE_CXX_COMPILER=nvcc -DBACKEND=thrust -DDEVICE=nvidia_v100
make
```

- SYCL compilation using icpx for Intel Ponte Vecchio

```bash
cmake ../ -DCMAKE_CXX_COMPILER=icpx -DBACKEND=sycl -DDEVICE=intel_pvc -D CMAKE_BUILD_TYPE=Release
make
```

- STDPAR compilation using nvc++ for Nvidia V100

```bash
cmake ../ -DCMAKE_CXX_COMPILER=nvc++ -DBACKEND=stdpar -DDEVICE=nvidia_v100
make
```

- STDPAR compilation on CPU using g++ or nvc++

```bash
cmake ../ -DCMAKE_CXX_COMPILER=g++ -DBACKEND=stdpar_cpu
make
OR
cmake ../ -DCMAKE_CXX_COMPILER=nvc++ -DBACKEND=stdpar_cpu -DCMAKE_CXX_FLAGS="-stdpar=multicore"
make
```

