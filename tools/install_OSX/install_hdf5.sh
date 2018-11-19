#!/bin/bash

# Stop script on error
set -e

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "1 argument required: INSTALLATION PATH, $# provided"

echo "Installing HDF5 in $1/hdf5"

wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz
tar xf hdf5*
cd hdf5*

CC=mpicc CXX=mpicxx ./configure --prefix=$1/hdf5 --enable-parallel
make -j4
make install
