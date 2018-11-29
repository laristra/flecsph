#!/bin/bash

# Stop script on error
set -e

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "1 argument required: INSTALLATION PATH, $# provided"

echo "INSTALLING GCC IN:     $1gcc"
echo "INSTALLING OPENMPI IN: $1openmpi"

export INSTALL_PATH=$1

# Clean previous
rm -rf gmp*
rm -rf mpfr*
rm -rf mpc*
rm -rf isl*
rm -rf $INSTALL_PATH/gcc
rm -rf gcc*/build
echo "Cleaned"

#Download latest
curl -L ftp://gcc.gnu.org/pub/gcc/releases/gcc-8.2.0/gcc-8.2.0.tar.gz | tar xjf -

# Dependencises

curl -L ftp://gcc.gnu.org/pub/gcc/infrastructure/gmp-6.1.0.tar.bz2 | tar xf -
curl -L ftp://gcc.gnu.org/pub/gcc/infrastructure/mpfr-3.1.4.tar.bz2 | tar xf -
curl -L ftp://gcc.gnu.org/pub/gcc/infrastructure/mpc-1.0.3.tar.gz | tar xf -
curl -L ftp://gcc.gnu.org/pub/gcc/infrastructure/isl-0.18.tar.bz2 | tar xf -

# GMP
cd gmp*
mkdir build && cd build
../configure --prefix=$INSTALL_PATH/gcc --enable-cxx
make -j4
make install-strip

# MPFR
cd ../..
cd mpfr*
mkdir build && cd build
../configure --prefix=$INSTALL_PATH/gcc --with-gmp=$INSTALL_PATH/gcc
make -j4
make install-strip

# MPC
cd ../..
cd mpc*
mkdir build && cd build
../configure --prefix=$INSTALL_PATH/gcc --with-gmp=$INSTALL_PATH/gcc \
  --with-mpfr=$INSTALL_PATH/gcc
make -j4
make install-strip

# ISL
cd ../..
cd isl*
mkdir build && cd build
../configure --prefix=$INSTALL_PATH/gcc --with-gmp-prefix=$INSTALL_PATH/gcc
make -j4
make install

# GCC
cd ../..
cd gcc*
mkdir build && cd build
../configure --prefix=$INSTALL_PATH/gcc --enable-checking=release \
  --with-gmp=$INSTALL_PATH/gcc --with-mpfr=$INSTALL_PATH/gcc \
  --with-mpc=$INSTALL_PATH/gcc --enable-languages=c,c++,fortran \
	--with-isl=$INSTALL_PATH/gcc
make -j4
make install-strip

# Export GCC
export PATH=$INSTALL_PATH/gcc/bin:$PATH
export LD_LIBRARY_PATH=$INSTALL_PATH/gcc/lib:$LD_LIBRARY_PATH

which gcc
which g++

# OPENMPI
cd ../..
rm -rf openmpi*

wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.3.tar.gz
tar xf openmpi-3.1.3.tar.gz
cd openmpi*

./configure CC=gcc CXX=g++ --prefix=$INSTALL_PATH/openmpi --enable-shared=yes --enable-static=yes
make -j4
make install

# Export OPENMPI
export PATH=$INSTALL_PATH/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$INSTALL_PATH/openmpi/lib:$LD_LIBRARY_PATH
