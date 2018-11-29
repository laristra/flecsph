# Install FLECSPH on MacOSX

## Install gcc 8.2.0 and OpenMPI 3.1.3:
- Install the xcode requirements
```
$ xcode-select --install
$ cd /Library/Developer/CommandLineTools/Packages/
$ open .
```
- Install the pkg present in the windows that pops up
- Use the installation script specifying the installation directory.
- Prepare some coffee, sit comfortably, this can be long.
```
$ ./install_gcc_openpmi.sh /home/XX/local
```
- The script will create and gcc/ and openmpi/ suddirectory.
- Add them to your PATH in bach_profile. Directly or as function:
```
# File ~/.bach_profile
load_gcc_openmpi(){
  export PATH=...:$PATH
  export LD_LIBRARY_PATH=...:$LD_LIBRARY_PATH      
}
# Terminal
$ load_gcc_openpmi
```
- Be sure gcc and g++ are in the path. You should see:  
```
$ g++ --version
g++ (GCC) 8.2.0
Copyright (C) 2018 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

## Install boost
- Download latest boost version: https://sourceforge.net/projects/boost/files/latest/download
- untar and compile/install
```
$ CC=gcc CXX=g++ ./bootstrap.sh
$ CC=gcc CXX=g++ ./b2 install --prefix=/home/XX/local/boost
```

## Install HDF5 from the script
- Use the script to install the latest version of HDF5:
```
./install_hdf5.sh ~/local/
```

## Install FleCSI, branch feature/flecsph
- Install FleCSI with the flag DBUILD_SHARED_LIBS=OFF
- Be sure to use the gcc/g++: CMAKE_CXX_COMPILER=g++ and CMAKE_C_COMPILER=gcc
```
$ git clone --recursive git@github.com:laristra/flecsi.git
$ cd flecsi
$  git checkout feature/flecsph
$  git submodule update --recursive
$  mkdir build ; cd build
$  export CMAKE_PREFIX_PATH=~/local/flecsph
$  cmake .. \
       -DCMAKE_INSTALL_PREFIX=$CMAKE_PREFIX_PATH  \
       -DENABLE_MPI=ON                            \
       -DENABLE_OPENMP=ON                         \
       -DCXX_CONFORMANCE_STANDARD=c++17           \
       -DENABLE_CLOG=ON                           \
       -DFLECSI_RUNTIME_MODEL=mpi                 \
       -DENABLE_FLECSIT=OFF                       \
       -DENABLE_FLECSI_TUTORIAL=OFF               \
       -DBUILD_SHARED_LIBS=OFF                    \
       -DCMAKE_CXX_COMPILER=g++                   \
       -DCMAKE_C_COMPILER=gcc
$ make install -j4
```

# Install FleCSPH
- Be sure to use the correct gcc/g++ like previously
