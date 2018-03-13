#!/bin/bash

VERSION="hdf5-1.8.19"
VERS=${VERSION%.[0-9]*}
TARBALL="${VERSION}.tar.gz"
CDIR=`pwd`

# TODO: use cmake to specify installation directory
INSTALLDIR=${CDIR%/*/*}/local

case $CDIR in
*/flecsph/build)
  echo "Building HDF5 in flecsph/build/hdf5_parallel directory.."
  ;;
*)
  echo "ERROR: this script must be run in build/ subdir of flecsph"
  exit 1
esac

echo "Downloading HDF5 code..."
rm -f ${TARBALL}
wget https://support.hdfgroup.org/ftp/HDF5/releases/${VERS}/${VERSION}/src/${TARBALL}
echo "done."

printf "Wiping hdf5_parallel..."
rm -rf hdf5_parallel
echo "done."

echo "Unpacking to hdf5_parallel..."
mkdir  hdf5_parallel
tar -zxvf ${TARBALL} -C hdf5_parallel --strip 1
echo "done."

cd hdf5_parallel/

echo "Configure..."
CC=mpicc CXX=mpicxx ./configure --enable-parallel --prefix="$INSTALLDIR"

make -j8 
make -j8 check 

echo "Installing to $INSTALLDIR..."
make install
echo "... finished."
