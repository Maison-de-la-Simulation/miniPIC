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

| CPU TASK backends    | Description                    |
|----------------------|--------------------------------|
| eventify             | Eventify version               |
| openmp_task          | OpenMP task version            |


Tools:

- `-DSHAMAN=ON/OFF`: enable/disable Shaman tool (only available in OpenMP mode or sequential mode)

Others:

- `-DDEBUG=ON/OFF`: enable/disable debug mode (`OFF` by default)
- `-DTEST=ON/OFF`: enable/disable tests mode (for CI, `OFF` by default)

- `-DMINIPIC=ON/OFF`: enable/disable minipic compilation (`ON` by default)
- `-UNIT_TESTS=ON/OFF`: enable/disable unit tests compilation (`OFF` by default)

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


