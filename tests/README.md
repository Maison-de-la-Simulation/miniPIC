# Tests

## Description

This directory contains scripts for testing and for continuous integration.

The script `run.py` is used by `.gitlab-ci.yml` to compile and run the following list of becnhmarks:

- `default.hpp`
- `antenna.hpp`
- `Bcst.hpp`
- `Ecst.hpp`
- `beam.hpp`

For a given benchmark, the output results can be analyzed to validate the run.
This only happens if a script of the same name is present in the directory `validate` on the root of the repository.
For instance, `default.py` being implemented in `validate`, the run of `default.hpp` will be analyzed at the end of the simulation.

## How to use

Some configurations are alreday implemented at the top of `run.py`.

You can get some help with:

```bash
python run.py -h
```