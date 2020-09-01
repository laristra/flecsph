![logo](doc/flecsph_logo_bg.png)

[![Build Status](https://travis-ci.com/laristra/flecsph.svg?branch=master)](https://travis-ci.com/laristra/flecsph)
<!---
[![codecov.io](https://codecov.io/github/laristra/flecsph/coverage.svg?branch=master)](https://codecov.io/github/laristra/flecsph?branch=master)
[![Quality Gate](https://sonarqube.com/api/badges/gate?key=flecsph%3A%2Fmaster)](https://sonarqube.com/dashboard?id=flecsph%3A%2Fmaster)
--->

# Introduction 

FleCSPH is a multi-physics compact application that exercises FleCSI parallel
data structures for tree-based particle methods. In particular, FleCSPH
implements a smoothed-particle hydrodynamics (SPH) solver for the solution of
Lagrangian problems in astrophysics and cosmology. FleCSPH includes support
for gravitational forces using the fast multipole method (FMM).
Currently, particle affinity and gravitation is handled using the parallel
implementation of the octree data structure provided by FleCSI.

We provide several examples of physics problems in 1D, 2D and 3D:

- Sod shock tubes in 1D/2D/3D;
- Noh shock test in 2D/3D;
- Sedov blast waves 2D and 3D;
- airfoil flow in a wind tunnel (2D/3D);
- pressure-induced spherical implosion (2D/3D);
- single and binary stars with Newtonian gravity in 3D.

# Building FleCSPH with Spack

FleCSPH can now be installed as a Spack package:
- Download spack at: github.com/spack/spack 
- Follow installation instructions 
- Use the following command to install _module_ support for spack and load the module. The second line can be added in your bash_profile.sh
```{engine=sh}
spack install lmod
. $(spack location -i lmod)/lmod/lmod/init/bash
```
- Run:
```{engine=sh}
spack install flecsph 
```
This will build all the dependencies, compile and install FleCSPH. 
In order to use FleCSPH executables simply run: 
```{engine=sh}
spack load flecsph 
```

You will then have access to the generators and the drivers: 
- sodtube\_{1-2-3}d\_generator, sedov\_{2-3}d\_generator...
- hydro\_{1-2-3}d, newtonian\_3d...

You can access to pre-configured parameter files and examples by downloading this repository: 
```{engine=sh}
git clone --recursive git@github.com:laristra/flecsph.git
cd flecsph
```
Sample parameter files and the intial data can be found in the `data` subdirectory.

For the developper guideline, please refer to this page: 
[Development Guidelines](https://github.com/laristra/flecsph/blob/master/doc/development.md)

# Building FleCSPH manually

Below we assume that FleCSPH is installed in FLECSPH root directory 
`${FLECSPH_ROOT}`.

## Suggested directory structure

We recommend to use an isolated installation of FleCSPH and FleCSI, such that
the software and all the dependencies in a separate directory, with the 
following directory structure:

```{engine=sh}
  ${FLECSPH_ROOT}
  ├── flecsi
  │   └── build
  ├── flecsph
  │   ├── build
  │   └── third-party-libraries
  └── local
      ├── bin
      ├── include
      ├── lib
      ├── lib64
      └── share
```

All the build happens in `build` subdirectories, and compiled dependencies are
installed in `local` subdirectory. 
Make sure to set your CMAKE prefix to this location:

    % export CMAKE_PREFIX_PATH=${FLECSPH_ROOT}/local

## Prerequisites

You will need the following tools:

- C++17 - capable compiler, such as gcc version >= 7;
- git version > 2.14;
- MPI libraries;
- cmake version >= 3.15;
- boost library version > 1.59;
- Python version >= 3.6.
- HDF5 compiled with parallel flag version > 1.8
- GSL library 

## FleCSI

Clone FleCSI repo at the `stable/flecsph` branch.
Checkout submodules recursively, then configure as below:

```{engine=sh}    
   export CMAKE_PREFIX_PATH=${FLECSPH_ROOT}/local
   cd ${FLECSPH_ROOT}
   git clone --recursive git@github.com:laristra/flecsi.git
   cd flecsi
   git checkout stable/flecsph
   git submodule update --recursive
   mkdir build ; cd build
   cmake .. \
       -DENABLE_OPENMP=OFF \
       -DCXX_CONFORMANCE_STANDARD=c++17 \
       -DENABLE_METIS=ON \
       -DENABLE_PARMETIS=ON \
       -DENABLE_COLORING=ON \
       -DENABLE_DEVEL_TARGETS=ON \
       -DENABLE_LOG=ON           \
       -DFLECSI_RUNTIME_MODEL=mpi
```    

In this configuration, FleCSI is installed with the MPI backend.

Build as a final step:

    % make -j

## FleCSPH

Clone FleCSPH git repo:
```{engine=sh}
   cd ${FLECSPH_ROOT}
   git clone --recursive git@github.com:laristra/flecsph.git
```    

### Building FleCSPH

Configure and build commands:

```{engine=sh}
   # in ${FLECSPH_ROOT}/build:
   export CMAKE_PREFIX_PATH=${FLECSPH_ROOT}/local
   cmake .. \
       -DCMAKE_BUILD_TYPE=debug \
       -DENABLE_UNIT_TESTS=ON   \
       -DENABLE_DEBUG=OFF       \
       -DLOG_STRIP_LEVEL=1
   make -j
   make install
```


### Building FleCSPH on various architectures
Architecture-/machine-specific notes for building FleCSPH are collected in
[doc/machines](https://github.com/laristra/flecsph/tree/master/doc/machines).
If you succeeded in compiling and running FleCSPH on new architectures,
please do not hesitate to share your recipe.
We appreciate user contributions.

# Running FleCSPH applications
Current FleCSPH contains several initial data generators and two evolution
drivers: `hydro` and `newtonian`. Initial data generators are located in
`app/id_generators/`:

- `sodtube`: 1D/2D/3D sodtube shock test;
- `sedov`: 2D and 3D Sedov blast wave;
- `noh`: 2D and 3D Noh implosion test.
- etc.

Evolution drivers are located in `app/drivers`:

- `hydro`: 1D/2D/3D hydro evolution without gravity;
- `newtonian`: 3D hydro evolution with self-gravity.

To run a test, you also need an input parameter file, specifying parameters of the
problem. Parameter files are located in `data/` subdirectory. Running an
application consists of two steps:

- generating initial data;
- running evolution code.

For instance, to run a `sodtube` 1D shock test, do the following (assuming
you are in your build directory after having successfully built FleCSPH):
```{engine=sh}
  cp ../data/sodtube_t1_n1000.par sodtube.par
  # edit the file sodtube.par to adjust the number of particles etc.
  app/id_generators/sodtube_1d_generator sodtube.par
  app/driver/hydro_1d sodtube.par
```
Our [wiki](https://github.com/laristra/flecsph/wiki#several-examples-you-may-want-to-try) 
page contains more examples that you may want to try.

# Creating your own initial data or drivers
You can add your own initial data generator or a new evolution module under
`app/id_generators` or `app/drivers` directories. Create a directory with
unique name for your project and modify `CMakeLists.txt` to inform the cmake
system that your project needs to be built.

A new initial data generator usually has a single `main.cc` file and an optional
include file. You can use existing interfaces for lattice generators or equations
of state in the `include/` directory.
The file `app/drivers/include/user.h` defines the dimensions of your problem, both
for initial data generators and for the evolution drivers.
This is done via a compile-time macro `EXT_GDIMENSION`, which allows users to have
the same source code for different problem dimensions. Actual dimension is set at
compile time via the `target_compile_definitions` directive of cmake, e.g.:
```
   target_compile_definitions(sodtube_1d_generator PUBLIC -DEXT_GDIMENSION=1)
   target_compile_definitions(sodtube_2d_generator PUBLIC -DEXT_GDIMENSION=2)
```

A new evolution driver must have a `main.cc` and `main_driver.cc` files. Do not edit
`main.cc`, because FleCSI expects certain format of this file. It is easier to start
by copying existing files to your folder under `app/drivers`. Include cmake
targets with different dimensions using examples in `app/drivers/CMakeLists.txt`.

Make sure to document your subproject in a corresponding `README.md` file
that describes the problem you want to run. In order to get all files easily and
correctly, you can copy them from other subprojects such as `sodtube` or `hydro`.

# Logs

In the code, you can set the level of output from trace(0) and info(1) to warn(2), error(3) and fatal(4).
You can then control the level of output at compile time by setting the flag `LOG_STRIP_LEVEL`:
by default it is set to 0 (trace), but for simulations it is perhaps preferrable to set it to 1 (info).
```cpp
  log_one(trace) << "This is verbose output  (level 0)" << std::endl;
  log_one(info) << "This is essential output (level 1)" << std::endl;
  log_one(warn) << "This is a warning output (level 2)" << std::endl;
  log_one(fatal) << "Farewell!" << std::endl;
```

For further details, refer to the documentation at:
[Cinch: Logging](https://github.com/laristra/cinch/blob/master/logging/README.md)

# Contacts

If you have any questions or concerns regarding FleCSPH, please contact us 
via the mailing list flecsph-support@lanl.gov
