# N-Body

N-body simulation written in C that currently models the planets of the solar system. The program fetcges initial vectors from JPL Horizons and integrates them using selectable integration methods (Euler, Verlet, RK4).

## Requirements

The current version is tested on CachyOS with the latest versions of the compilers (Clang 20.1.8 and GCC 15.1.1). It should work without issues using WSL on Windows. It should also function natively on Windows if compiled using MSVC, although the makefile with not work and it has not been tested.

## Quick Start

#### 1. Build:
```bash
make
```

#### 2. Run:
```bash
make run
```

#### 3. Clean build and program outputs
```bash
make clean
```

## Configuration

The simulation parameters can be changed in src/main.h

- SYSTEM_METHOD: EULER / VERLET / RK4
- SYSTEM_YEAR: epoch year used for initial vectors
- SYSTEM_SIZE: number of planets (excludes Sun)
- SYSTEM_STEP: integration step in days
- SYSTEM_TIME: total simulated time
- NAME: output filename

To change the used compiler between GCC and Clang, change the value for the `CC` variables in the makefile.

## TODO

- Parallelize API calls (multithread curl section)
- Use one thread to write to file (speed up writing for large simulations)
- Add Barnes-Hut and other force approximations
- Add algorithms to initialize other simulation types (galaxies)
- Add multithreading and/or GPU acceleration (OpenMP, HIP, or OpenGL/OpenCL)
- Add UI to make it easier to use and visualize data
