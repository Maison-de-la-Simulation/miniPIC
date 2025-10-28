# Developer zone

## Domain decomposition

For the moment, miniPIC does not support distributed memory parallelism.

The domain is decomposed into patches using a 3D cartesian decomposition. Each patch is designed to be independent and can be computed in parallel. It represents a piece of the domain that contains the particles and the local current field grid.

Maxwell's equations are solved at the domain scale.
Therefore, electromagnetic fields grids are global.
A reduction operation is performed to compute the global current grid from the local current grids.

## PIC loop steps

<img title="pic loop" alt="pic loop" src="./images/pic_loop.png" height="500">

## Code design

The figure below illustrates schematically the code design. It shows how the different classes are organized and how they interact with each other.

<img title="code design" alt="code design" src="./images/code_design.png" height="700">

Each file provides either a set of functions, a namespace or a data container (class).

| File                   | Where  |Description                                                                                 |
|------------------------|--------|---------------------------------------------------------------------------------------------|
| Headers                | common | Determine the best headers to use depending on the selected backend                         |
| Backend                | common | Data container that contains backend specific parameters (often global) for parallelism     |
| Vector                 | common | Vector class that mimics the std::vector with backend abstraction for both CPU and GPU      |
| Field                  | common | Class that provides 3D arrays with backend abstraction for both CPU and GPU                 |
| Particle               | common | Class that provides a particle container with backend abstraction for both CPU and GPU      |
| ElectroMagn            | common | Class that provide a data container for electromagnetic and current grids                   |
| Patch                  | common | Data container representing a patch entity (see patch decomposition)                        |
| SubDomain              | model specific folders | SubDomain is a data container representing a domain piece                   |
| Diagnostics            | common | Function to perform diagnostic output                                                       |
| Operators              | model specific folders | Functions to perform the Particle-In-Cell loop (such as interpolator, pusher, projection, etc), this header is duplicated for each programming models |
| Timers                 | common | Class that provide timer functionality to monitor the time and make statistics              |
| Profilers              | common | Class that provide a home-made profiler                                                     |
| Main                   | src    | Main source file for the global code structure                                              |

## Macros

| Macros                 | Description                                            |
|------------------------|--------------------------------------------------------|
| `__MINIPIC_SIMD__`     | Activate specific SIMD pragmas                         |
| `__MINIPIC_OMP__`      | Activate specific OPENMP operators                     |
| `__MINIPIC_THRUST__`   | Activate specific THRUST operators                     |
| `__MINIPIC_THRUST_COMMON__` | Activate THRUST common backend (used for Fields and vectors to manage memory)                    |
| `__MINIPIC_KOKKOS__`   | Activate specific KOKKOS operators                     |
| `__MINIPIC_OMP_TASK__` | Activate specific OPENMP task operators                |
| `__MINIPIC_EVENTIFY__` | Activate specific eventify operators                   |
| `__MINIPIC_OMP_TARGET__` | Activate specific OPENMP target operators (use the )             |
| `__MINIPIC_OPENACC__`   | Activate specific SHAMAN operators                     |
| `__MINIPIC_STDPAR__`   | Activate specific STDPAR operators                     |


## Known issues

| Backend                | Issue                                                         |
|------------------------|---------------------------------------------------------------|
| OpenMP target          | Reduction on device does not work properly (gives 0)          |
| OpenMP target          | Compilation issue on GH200                                    |