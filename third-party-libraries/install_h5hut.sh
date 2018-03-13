#!/bin/bash 

VERSION="H5hut-1.99.12"
TARBALL="${VERSION}.tar.gz"
CDIR=`pwd`

# TODO: use cmake to specify installation directory
INSTALLDIR=${CDIR%/*/*}/local

case $CDIR in
*/flecsph/build)
  echo "Building H5hut in flecsph/build/h5hut directory.."
  ;;
*)
  echo "ERROR: this script must be run in build/ subdir of flecsph"
  exit 1
esac

export PATH=$INSTALLDIR/bin:$PATH
export LD_LIBRARY_PATH=$INSTALLDIR/lib:$LD_LIBRARY_PATH
echo $PATH

echo "Downloading H5hut..."
rm -f ${TARBALL}
wget --no-check-certificate https://amas.psi.ch/H5hut/raw-attachment/wiki/DownloadSources/${TARBALL}
echo "done."

printf "Wiping build/h5hut directory..."
rm -rf h5hut
echo "done."

echo "Unpacking..."
mkdir h5hut
tar -xvzf ${TARBALL} -C h5hut --strip=1
echo "done."

cd h5hut/

CC=mpicc CXX=mpicxx ./configure --enable-parallel --with-hdf5="$INSTALLDIR" --prefix="$INSTALLDIR" --enable-shared=yes

echo "Replace the non valid function to be able to compile"
sed -e 's/if ( H5Pset_fapl_mpiposix (/if ( H5Pset_fapl_mpio (/g' ./src/h5core/h5_hdf5_private.h -i

echo "Building... "
make -j8
echo "done."

echo "Installing to $INSTALLDIR..."
make install 
echo "... finished."
