#!/bin/bash

VERSION="hdf5-1.8.19.tar"

echo "Downloading HDF5 code..."
wget https://support.hdfgroup.org/ftp/HDF5/current18/src/${VERSION}
echo "done."

echo "Untar to hdf5_parallel..."
mkdir hdf5_parallel
tar -xvf ${VERSION} -C hdf5_parallel --strip 1
echo "done."

cd hdf5_parallel/

echo "Configure..."
CC=mpicc ./configure --enable-parallel --prefix="$PWD/../local/"
make 
make check 
make install
