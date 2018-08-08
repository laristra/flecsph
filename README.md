![logo](doc/flecsph_logo_bg.png)

[![Build Status](https://travis-ci.org/laristra/flecsph.svg?branch=master)](https://travis-ci.org/laristra/flecsph)
[![codecov.io](https://codecov.io/github/laristra/flecsph/coverage.svg?branch=master)](https://codecov.io/github/laristra/flecsph?branch=master)
<!---
[![Quality Gate](https://sonarqube.com/api/badges/gate?key=flecsph%3A%2Fmaster)](https://sonarqube.com/dashboard?id=flecsph%3A%2Fmaster)
--->

# SPH on FleCSI

This project is an implementation of SPH problem using FleCSI framework.
This code intent to provide distributed and parallel implementation of the octree data structure provide by FleCSI.
The Binary, Quad and Oct tree respectively for 1, 2 and 3 dimensions is developed here for Smoothed Particle Hydrodynamics problems.

For the first version of the code we intent to provide several basic physics problems:

- Sod Shock Tube Tests in 1D
- Sedov Blast Wave 2D
- Binary Compact Object Merger 3D

You can find detail ingredients of FleCSPH such as formulations and algorithm under `doc/notes.pdf`. The document is constantly updated to contain latest development of FleCSPH

# Getting the Code

FleCSPH can be installed anywhere in your system; to be particular, below we
assume that all repositories are downloaded in FLECSPH root directory `${HOME}/FLECSPH`.
The code requires:

- FleCSI third party library
- FleCSI
- shared local directory

### Suggested directory structure

```{engine=sh}
   mkdir -p $HOME/FLECSPH/local
   cd $HOME/FLECSPH
   git clone --recursive git@github.com:laristra/flecsph.git
```    

```{engine=sh}
  ${HOME}/FLECSPH
  ├── flecsi
  │   └── build
  ├── flecsi-third-party
  │   └── build
  ├── flecsph
  │   ├── build
  │   └── third-party-libraries
  │       ├── install_h5hut.sh
  │       └── install_hdf5_parallel.sh
  └── local
      ├── bin
      ├── include
      ├── lib
      ├── lib64
      └── share
```

In this configuration, all dependencies are installed in `${HOME}/FLECSPH/local`.

## Install the dependencies

On DARWIN supercomputer load the modules:

    % module load gcc/6.2.0 openmpi/2.0.1-gcc_6.2.0 git/2.8.0 cinch/cinch-utils cmake/3.7.1 boost/1.59.0_gcc-6.2.0

### FleCSI third part libraries

Clone the FleCSI repo with third party libraries and check out FleCSPH-compatible branch `flecsph`:

```{engine=sh}    
   cd $HOME/FLECSPH
   git clone --recursive git@github.com:laristra/flecsi-third-party.git
   cd flecsi-third-party
   git checkout FleCSPH
   git submodule update
   mkdir build ; cd build
   ccmake ..
```    

Let all the flags ON, make sure the conduit of GASNET is MPI.
Set `CMAKE_INSTALL_PREFIX -> ~/FLECSPH/local`.
Build the libraries using several cores (note that no install step is required):

    % make -j8

### FleCSI

Clone FleCSI repo and change to FlecSPH branch:

```{engine=sh}    
   cd $HOME/FLECSPH
   git clone --recursive git@github.com:laristra/flecsi.git
   cd flecsi
   git checkout stable/flecsph-compatible
   git submodule update
   mkdir build ; cd build
   ccmake ..
```    

Press `c` to do initial configuring.
- Set `CMAKE_INSTALL_PREFIX -> ~/FLECSPH/local`.

Press `c` to reconfigure. Press `t` for advanced options; scroll down to select:
Here add:
- `ENABLE_MPI`: ON
- `ENABLE_MPI_CXX_BINDINGS`: ON
- `ENABLE_OPENMP`: ON
- `ENABLE_LEGION`: ON
- `FLECSI_RUNTIME_MODEL`: legion

Press `c` to reconfigure and `g` to generate configurations scripts.

You can also supply CMakeCache.txt to avoid multiple reconfigures:

```{engine=sh}
cat > CMakeCache.txt << EOF
  CMAKE_INSTALL_PREFIX:PATH=$HOME/FLECSPH/local
  ENABLE_MPI:BOOL=ON
  ENABLE_MPI_CXX_BINDINGS:BOOL=ON
  ENABLE_OPENMP:BOOL=ON
  ENABLE_LEGION:BOOL=ON
  FLECSI_RUNTIME_MODEL:STRING=legion
EOF
ccmake ..
```

If no errors appeared, build and install:

    % make -j8 install

In case of errors: if you are rebuilding everything from scratch, 
make sure that your installation directory (`$HOME/FLECSPH/local` 
in our example) is empty.

# Build

## Dependencies

In order to build flecsph some other dependencies can be found in the third-party-libraries/ directory.
- Use the scripts to install HDF5 and H5Hut from within build/ directory:

```{engine=sh}
   cd ~/FLECSPH/flecsph
   mkdir build; cd build
   ../third-party-libraries/install_hdf5_parallel.sh
   ../third-party-libraries/install_h5hut.sh
```    

- ScalingFramework is available in LANL property right now, soon open-source

## Build FleCSPH

Continue with the build:

    % ccmake ..

Set the following options:
- `CMAKE_INSTALL_PREFIX`: ~/FLECSPH/local
- `ENABLE_MPI`: ON
- `ENABLE_OPENMPI`: ON
- `ENABLE_LEGION`: ON
- `ENABLE_UNIT_TESTS`: ON
- `HDF5_C_LIBRARY_hdf5`: `~/FLECSPH/local/lib/libhdf5.so`

You can also use the following command to setup cmake cache:

```{engine=sh}
cat > CMakeCache.txt << EOF
  CMAKE_INSTALL_PREFIX:PATH=$HOME/FLECSPH/local
  ENABLE_LEGION:BOOL=ON
  ENABLE_MPI:BOOL=ON
  ENABLE_MPI_CXX_BINDINGS:BOOL=ON
  ENABLE_OPENMP:BOOL=ON
  ENABLE_UNIT_TESTS:BOOL=ON
  HDF5_C_LIBRARY_hdf5:FILEPATH=$HOME/FLECSPH/local/lib/libhdf5.so
  VERSION_CREATION:STRING=
EOF
ccmake ..
```

Configure, build and install:

    % make -j8 install

 # Running test cases

 You can find each test case in the corresponding app/ directory. Currently, we have below cases:

 - sodtube : 1D Sod shock tube problem
 - sedov : 2D circular Sedov blast wave expansion
 - fluid : 3D simple fluid drop simulation
 - bwd : 3D binary white dwarf merger with analytic zero temperature EOS

 Also, app/miscell directory contains python scripts that perform different aspects for initial data generators and converter for h5part format.

# Adding your own projects
FleCSPH can handle different projects. To add your own project, you first need to create a corresponding directory in the `app` folder, for example `myproject`. In `myproject`, you must have the following files: `CMakeLists.txt`, `generator/main.cc`, `/include/user.h`, `main.cc`, and `main_driver.cc`. It is also advisable to have a README.md file that describes the problem you want to run. In order to get all files easily and correctly, you can copy them from other application directories such as `sodtube`, paste them into `myproject` and modify according to your needs.

You need to modify the file `CMakeLists.txt` to have the correct name of your executable. The file `generator/main.cc` should contain the initialization of your problem. All values and functions that you want to put into your project and/or change from their default values should go into `main_driver.cc`. While `generator/main.cc` initializes your problem, `main_driver.cc` describes how it is run. The file `/include/user.h` defines the dimensions of your problem. Do not edit `main.cc`.

To compile your project, you need to add it in `flecsph/config/project.cmake`: At the end of this file you have a section called `Add application targets`. Here, you need to add:

```cinch add application directory("app/myproject")```

Then, `cmake` links your project and you can compile it. 

# Style guide

FleCSPH follows the FleCSI coding style, which in turn follows (in general) the Google coding conventions.
FleCSI coding style is documented here:
https://github.com/laristra/flecsi/blob/master/flecsi/style.md

# Development Logistic
If you would are new person for development of the FleCSPH, please check development logistic under doc/DEVELOPMENT.md

# Logs 

Cinch Log is the logging tool for this project. 
In order to display log set the environement variable as: 
```bash 
export CLOG_ENABLE_STDLOG=1
```

In the code, you can set the level of output from trace(0) and info(1) to warn(2), error(3) and fatal(4).
You can then control the level of output at compile time by setting the flag `CLOG_STRIP_LEVEL`:
by default it is set to 0 (trace), but for simulations it is perhaps preferrable to set it to 1 (info).
In the code, if you only want a single rank (e.g. rank 0) to output, you need to use `clog_one` command,
and indicate the output rank using the `clog_set_output_rank(rank)` macro right after MPI initialization:
```cpp
clog_set_output_rank(0);
clog_one(info) << "This is essential output (level 1)" << std::endl;
clog_one(trace) << "This is verbose output  (level 0)" << std::endl;
clog_one(fatal) << "Farewell!" << std::endl; 
```

For further details, refer to the documentation at: 
https://github.com/laristra/cinch/blob/master/logging/README.md

 # Contacts

 If you have any questions or find any worries about FleCSPH, please contact Julien Loiseau (julien.loiseau@univ-reims.fr) and/or Hyun Lim (hylim1988@gmail.com)

