# Timers

## Introduction

MiniPIC has an internal timing system that can be used to measure the time taken by different parts of the code.
The timers are updated all along the simulation.
The results are printed at the end of the simulation and can be saved in a file `timers.json`.

## Parameters

### Setup parameters

In the setup file :

- `params.save_timers_start` (unsigned int): First time step at which the timers are saved in the file `timers.json`. Default value is 0.
- `params.save_timers_period` (unsigned int): Period at which the timers are saved in the file `timers.json` from the start parameter `params.save_timers_start`. Default value is 0 (no output).
- `params.bufferize_timers_output` (bool) : If true, the timers are buffered and only output at the end of the simulation else, the file is appended at each time step. Default value is `true`.


### Command line parameters

Some parameters can be passed to MiniPIC as command line arguments:

- `-sts` or `--save-timers-start` (unsigned int): First time step at which the timers are saved in the file `timers.json`. Default value is 0.
- `-stp` or `--save-timers-period` (unsigned int): Period at which the timers are saved in the file `timers.json` from the start parameter `params.save_timers_start`. Default value is 0 (no output).


## Undestand the File `timers.json`

### Global file structure

The file `timers.json` contains the timing information in JSON format.
Depending on the input parameters, the timing at intermediate steps can be saved in the file as well to follow the time evolution all along the simulation.

The file has a first key `parameters` that contains some input parameters useful for analysing the timers:

```json
  "parameters" : {
    "number_of_patches" : 1,
    "number_of_threads" : 1,
    "iterations" : 25,
    "save_timers_period" : 30,
    "save_timers_start" : 5
  },
```

Then, the next key is `initialization` and represents the time to initialize the simulation:

```json
  "initialization" : 4.38563E1,
```

Then, if a timer output is requested at intermediate steps (using parameter `save_timers_period`), the file contains a key representing the corresponding time step as show below:


```json
{
  "50" : {
    "pic iteration" : 1.884484e-01,
    "diags" : 2.000652e-01,
    "all diags" : [2.000652e-01, 0.000000e+00],
    "interpolate" : [0.000000e+00, 1.533500e-05],
    "push" : [0.000000e+00, 3.668000e-06],
    "pushBC" : [0.000000e+00, 7.496000e-06],
    "id_parts_to_move" : [0.000000e+00, 1.334000e-06],
    "exchange" : [0.000000e+00, 0.000000e+00],
    "reset_current" : [0.000000e+00, 0.000000e+00],
    "projection" : [0.000000e+00, 0.000000e+00],
    "current_local_reduc" : [0.000000e+00, 0.000000e+00],
    "current_global_reduc" : [0.000000e+00, 0.000000e+00],
    "currentBC" : [0.000000e+00, 0.000000e+00],
    "maxwell_solver" : [1.587832e-01, 0.000000e+00],
    "maxwellBC" : [2.397584e-02, 0.000000e+00],
    "diags_sync" : [0.000000e+00, 0.000000e+00],
    "diags_binning" : [3.204291e-03, 0.000000e+00],
    "diags_cloud" : [0.000000e+00, 0.000000e+00],
    "diags_scalar" : [8.773675e-02, 0.000000e+00],
    "diags_field" : [1.052382e-01, 0.000000e+00],
    "imbalance" : [0.000000e+00, 0.000000e+00]
  },
```

At the end of the simulation, a key `final` contains the timing information for the last time step:

```json
  "final" : {
    "pic iteration" : 2.203330e+00,
    "diags" : 1.790520e+00,
    "all diags" : [1.790520e+00, 0.000000e+00],
    "interpolate" : [0.000000e+00, 1.932510e-04],
    "push" : [0.000000e+00, 4.426300e-05],
    "pushBC" : [0.000000e+00, 8.508800e-05],
    "id_parts_to_move" : [0.000000e+00, 2.546700e-05],
    "exchange" : [0.000000e+00, 0.000000e+00],
    "reset_current" : [0.000000e+00, 0.000000e+00],
    "projection" : [0.000000e+00, 0.000000e+00],
    "current_local_reduc" : [0.000000e+00, 0.000000e+00],
    "current_global_reduc" : [0.000000e+00, 0.000000e+00],
    "currentBC" : [0.000000e+00, 0.000000e+00],
    "maxwell_solver" : [1.852589e+00, 0.000000e+00],
    "maxwellBC" : [2.834123e-01, 0.000000e+00],
    "diags_sync" : [0.000000e+00, 0.000000e+00],
    "diags_binning" : [3.636822e-02, 0.000000e+00],
    "diags_cloud" : [0.000000e+00, 0.000000e+00],
    "diags_scalar" : [8.885825e-01, 0.000000e+00],
    "diags_field" : [8.198371e-01, 0.000000e+00],
    "imbalance" : [0.000000e+00, 0.000000e+00]
  }
}
```

And the file is terminated by the last key `main_loop` that represents the total time in the main loop:

```json
  "main_loop" : 5.99304E1
}
```

### How to understand the timers for a specific time step

The timers for a specific time setp are decomposed into different parts of the code.

Depending on the code section, timers can be placed inside parallel loop (in OpenMP for or OpenMP task for instance ) or outside the parallel loop.
Sometime, a timer for a specific code section is the combination of different subparts with one being outside a parallel loop and the other inside a parallel loop.
Therefore, timers are a list of values that can be interpreted as follows:

- The first element of the list is the time spent outside the parallel loop.
- Then the rest of the list is the time spent in each patch inside the parallel loop.

In the example below, we use 64 patches in OpenMP for mode. For the particle pusher, the timer is only inside the parallel loop so that the first element of the list is 0 and the rest of the list contains the time spent in each patch.

```json
    "push" : [0.000000e+00, 5.168000e-06, 8.956000e-06, 9.364000e-06, 4.551000e-06, 1.095900e-05, 2.726333e-03, 2.543290e-03, 1.075100e-05, 9.871000e-06, 2.663749e-03, 2.549078e-03, 1.417200e-05, 4.456000e-06, 1.200600e-05, 1.121300e-05, 4.623000e-06, 4.755000e-06, 9.091000e-06, 1.270500e-05, 3.989000e-06, 1.141800e-05, 2.740460e-03, 2.652875e-03, 1.412800e-05, 1.308300e-05, 2.773991e-03, 2.485336e-03, 1.658100e-05, 4.251000e-06, 1.303900e-05, 1.512400e-05, 4.291000e-06, 5.628000e-06, 8.633000e-06, 9.588000e-06, 3.876000e-06, 7.911000e-06, 2.761343e-03, 2.708544e-03, 8.757000e-06, 1.625400e-05, 2.833870e-03, 2.697086e-03, 1.203300e-05, 4.002000e-06, 1.132900e-05, 1.029300e-05, 4.618000e-06, 5.500000e-06, 1.016400e-05, 9.788000e-06, 4.327000e-06, 1.421000e-05, 2.924826e-03, 2.616838e-03, 8.212000e-06, 1.068000e-05, 2.687716e-03, 2.610037e-03, 1.066900e-05, 3.791000e-06, 1.029900e-05, 9.622000e-06, 4.168000e-06],
```

Keys `pic iteration` and `diags` are never in a parallel region.

## How to interprete the timers

### Computing the average time spent in a specific code section

The average time per thread spent in a specific code section can be computed by summing the time spent in each patch and dividing by the number of threads and summing with the time outside the parallel region.

For instance, the average time spent in the particle pusher can be computed as follows:

```python
import json

with open('timers.json') as json_file:
    timers = json.load(json_file)

mean_push_time = timers["final"]["push"][0] + sum(timers["final"]["push"][1:]) / timers["parameters"]["number_of_threads"]

print(f"Mean time spent in the particle pusher: {mean_push_time}")
```