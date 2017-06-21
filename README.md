![logo](doc/flecsph_logo_bg.png)


<aside class="warning">
The distributed version of gravitation with FMM is not working yet
Working on it
</aside>

# SPH on FleCSI 

This project is an implementation of SPH problem using FleCSI framework.
This code intent to provide distributed and parallel implementation of the octree data structure provide by FleCSI. 
The Binary, Quad and Oct tree respectively for 1, 2 and 3 dimensions is developed here for Smoothed Particle Hydrodynamics problems.

For the first version of the code we intent to provide several basic physics problems: 

- Sod Shock Tube in 1D
- Sedov Blast Wave 2D
- Binary Neutron Star Merger 3D  

# Getting the Code 

    % git clone --recursive git@github.com:laristra/flecsph.git

# Requirements

The code requires:

- FleCSI third part library 
- FleCSI 

## Install the dependencies 

On DARWIN supercomputer load the modules: 

    % module load gcc/6.2.0
    % module load openmpi/2.0.1-gcc_6.2.0
    % module load git/2.8.0
    % module load cinch/cinch-utils
    % module load cmake/3.7.1
    % module load boost/1.59.0_gcc-6.2.0

### FleCSI third part libraries

    % git clone --recursive https://github.com/laristra/flecsi-third-party.git
    % mkdir build ; cd build
    % ccmake ../

Let all the flags ON, make sure the conduit of GASNET is MPI. 
If not administrator: set path for the CMAKE_INSTALL_PREFIX like /home/XXX/local/

    % make ; make install 

### FleCSI 

    % git clone --recursive https://github.com/laristra/flecsi.git

Here we need to change to the FleCSPH branch 

    % git checkout flecsph 
    % mkdir build ; cd build 
    % ccmake ../

Here add:
- ENABLE_MPI 
- ENABLE_OPENMP 
- FLECSI_RUNTIME_MODEL legion
If not administrator:  
- CMAKE_INSTALL_PREFIX like /home/XXX/local/

    % make ; make install 

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

    % make 
