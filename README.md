![logo](doc/flecsph_logo_bg.png)

[![Build Status](https://travis-ci.org/laristra/flecsph.svg?branch=master)](https://travis-ci.org/laristra/flecsph)
[![codecov.io](https://codecov.io/github/laristra/flecsph/coverage.svg?branch=master)](https://codecov.io/github/laristra/flecsph?branch=master)
[![Quality Gate](https://sonarqube.com/api/badges/gate?key=flecsph%3A%2Fmaster)](https://sonarqube.com/dashboard?id=flecsph%3A%2Fmaster)


<aside class="warning">
The distributed version of gravitation with FMM is not working yet
Working on it
</aside>

# SPH on FleCSI

This project is an implementation of SPH problem using FleCSI framework.
This code intent to provide distributed and parallel implementation of the octree data structure provide by FleCSI.
The Binary, Quad and Oct tree respectively for 1, 2 and 3 dimensions is developed here for Smoothed Particle Hydrodynamics problems.

For the first version of the code we intent to provide several basic physics problems:

- Sod Shock Tube Tests in 1D
- Sedov Blast Wave 2D
- Binary Compact Object Merger 3D

# Getting the Code

FleCSPH can be installed anywhere in your system; to be particular, below we
assume that all repositories are downloaded in FleCSPH root directory `${HOME}/FleCSPH`.
All dependencies are installed in `${HOME}/FleCSPH/local`.

    % mkdir -p $HOME/FleCSPH/local
    % cd $HOME/FleCSPH
    % git clone --recursive git@github.com:laristra/flecsph.git

# Requirements

The code requires:

- FleCSI third party library
- FleCSI

## Install the dependencies

On DARWIN supercomputer load the modules:

    % module load gcc/6.2.0 openmpi/2.0.1-gcc_6.2.0 git/2.8.0 cinch/cinch-utils cmake/3.7.1 boost/1.59.0_gcc-6.2.0

### FleCSI third part libraries

Clone the FleCSI repo with third party libraries and check out FleCSPH-compatible branch `flecsph`:

    % cd $HOME/FleCSPH
    % git clone --recursive https://github.com/laristra/flecsi-third-party.git
    % git checkout flecsph
    % git submodule update
    % mkdir build ; cd build
    % ccmake ..

Let all the flags ON, make sure the conduit of GASNET is MPI.
Set `CMAKE_INSTALL_PREFIX -> ~/FleCSPH/local`.
Run configure and make using several cores:

    % cmake ..
    % make -j8

### FleCSI

Clone FleCSI repo and change to FlecSPH branch:

    % git clone --recursive https://github.com/laristra/flecsi.git
    % git checkout FleCSPH
    % git submodule update
    % mkdir build ; cd build
    % ccmake ..

Press `c` to do initial configuring.
- Set `CMAKE_INSTALL_PREFIX -> ~/FleCSPH/local`.

Press `c` to reconfigure. Press `t` for advanced options; scroll down to select:
Here add:
- `ENABLE_MPI`: ON
- `ENABLE_MPI_CXX_BINDINGS`: ON
- `ENABLE_OPENMP`: ON
- `ENABLE_LEGION`: ON
- `FLECSI_RUNTIME_MODEL`: legion

Press `c` to reconfigure and `g` to generate configurations scripts.
If no errors have emerged, configure, make and install:

      % ccmake ..
      % make -j8
      % make install

In case of errors: if you are rebuilding everything from scratch, 
make sure that your installation directory (`$HOME/FleCSPH/local` 
in our example) is empty.

# Build

## Dependencies

In order to build flecsph some other dependencies can be found in the third-party-libraries/ directory.
Use the two script to install:
- hdf5 parallel
- h5hut
- ScalingFramework is available in LANL property right now, soon open-source

## Build FleCSPH

    % mkdir build
    % cd build
    % ccmake ../

- ENABLE_MPI: ON
- ENABLE_OPENMPI: ON
- ENABLE_LEGION: ON

Then make:

    % make

 # Running test cases

 You can find each test case in the corresponding app/ directory. Currently, we have below cases:

 - sodtube : 1D Sod shock tube problem
 - sedov : 2D circular Sedov blast wave expansion
 - fluid : 3D simple fluid drop simulation
 - bwd : 3D binary white dwarf merger with analytic zero temperature EOS

 Also, app/miscell directory contains python scripts that perform different aspects for initial data generators and converter for h5part format.

 # Contacts

 If you have any questions or find any worries about FleCSPH, please contact Julien Loiseau (julien.loiseau@univ-reims.fr) and/or Hyun Lim (hylim1988@gmail.com)

