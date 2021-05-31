# cluster-editing-2021

This repository stores the source code of my solver for the Exact track of the [PACE 2021](https://pacechallenge.org/2021/) challenge.

### DOI of Version 1

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4877899.svg)](https://doi.org/10.5281/zenodo.4877899)

### Solver Description

TBD

### Requirements

- [GNU Make](https://www.gnu.org/software/make/)
- [CMake](https://cmake.org/) version 3.5 or higher
- C++ compiler that supports the C++14 standard.

### How to build

In this directory, run the following command.

```
make build
```

And the executable file `Exact` will be generated under the `dist` directory.

Notes: The submission file (`dist/Exact.tgz`) will be created by `make publish`.
