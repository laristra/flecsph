#!/bin/bash 

export PATH=$PWD/local/bin:$PATH
export LD_LIBRARY_PATH=$PWD/local/lib:$LD_LIBRARY_PATH
echo $PATH

VERSION="H5hut-1.99.12.tar.gz"

echo "Downloading H5hut..."
wget --no-check-certificate https://amas.psi.ch/H5hut/raw-attachment/wiki/DownloadSources/${VERSION}
echo "done."

echo "Untar..."
mkdir h5hut
tar -xvzf ${VERSION} -C h5hut --strip=1
echo "done."

cd h5hut/

CC=mpicc CXX=mpicxx ./configure --enable-parallel --with-hdf5="$PWD/../local" --prefix="$PWD/../local/" --enable-shared=yes

echo "Replace the non valid function to be able to compile"
sed -e 's/if ( H5Pset_fapl_mpiposix (/if ( H5Pset_fapl_mpio (/g' ./src/h5core/h5_hdf5_private.h -i


make -j
make install 
