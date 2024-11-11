# Validation

## Benchmark validation

A validation mechanism can be used locally and on super-computer to validate the code modification.
A python script `run.py` located in the `tests` folder is designed to test the whole code against reference data for a list of benchmarks.
This script needs Python libraries located in `lib`.

You can access the help page by doing:
```bash
python run.py -h
```

### Available options:

Here is a list of available options:

| Option | Long Option | Description |
| --- | --- | --- |
| `-h` | `--help` | Show this help message and exit |
| `-g CONFIG` | `--config CONFIG` | Configuration choice: sequential, openmp (default), kokkos, kokkos_gpu, thrust |
| `-c COMPILER` | `--compiler COMPILER` | Custom compiler choice |
| `-b BENCHMARKS` | `--benchmarks BENCHMARKS` | Specific benchmark, you can specify several benchmarks with a coma. For instance "default,beam" |
| `-t THREADS` | `--threads THREADS` | Default number of threads |
| `-a ARGUMENTS` | `--arguments ARGUMENTS` | Default arguments |
| | `--clean` | Whether to delete or not the generated files |
| | `--no-evaluate` | If used, do not evaluate against the reference |
| | `--compile-only` | If used, only compile the tests |
| | `--threshold THRESHOLD` | Threshold for the validation |
| | `--save-timers` | Save the timers for each benchmark |

### Configurations

Here is a list of possible configurations:

| Configuration | Description | Compiler | CMake Options | Run Prefix |
| --- | --- | --- | --- | --- |
| sequential | | | `-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="sequential"` | |
| openmp | OpenMP for version (for CPU) |  | `-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="openmp"` | `OMP_PROC_BIND=spread` |
| openmp_task || | `-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="openmp_task"` | `OMP_PROC_BIND=spread OMP_MAX_ACTIVE_LEVELS=10` |
| eventify || | `-DCMAKE_VERBOSE_MAKEFILE=ON -DBACKEND="eventify"` | `KMP_AFFINITY=` |


### Usage examples:

- default run

```bash
python run.py
```

- specific configuration

```bash
python run.py -g eventify
```

- custom compiler

```bash
python run.py -c clang++
```