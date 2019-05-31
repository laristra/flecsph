#!/bin/bash
#
# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved
#
#  This script extracts 1D cuts from a single HDF5 for all timesteps
#  in a gnuplot-readable format, and writes them in two files:
#  - xdp.x.dat: coordinate (x), density (d) and pressure(p);
#  - xuv.x.dat: coordinate (x), internal energy (u) and velocity in the x-direction
#


# - Options:
#   path to h5part slicer: modify this for your system
h5slicer=../../flecsph/tools/h5part_ascii_slicer.py  

if [ $# -ne 1 ]; then
  echo "Creates 1D slices of an HDF5 evolution file. "
  echo "Usage: ./$0 <evolution_file.h5part>"
  exit -1
fi

evfile=$1
if [ ! -f $evfile ]; then
  echo "ERROR: file $evfile does not exist!"
  exit -2
fi
tmpfile=tmp.h5part # make a copy just in case
echo "making a copy of the evolution h5part file... "
cp $evfile $tmpfile

# remove previous versions of these files
xdpfile="xdp.x.dat"
xuvfile="xuv.x.dat"
rm -f $xdpfile $xuvfile

# get the number of timesteps in $evfile
Nsteps=`h5ls $tmpfile | wc -l`
maxstep=$[ $Nsteps - 1 ]

# main loop
for step in `seq 0 $maxstep`; do
  echo "processing iteration $step .. "
  $h5slicer -x -d 0.1 $tmpfile -s $step | sort -gk1 >> $xdpfile
  echo >> $xdpfile
  echo >> $xdpfile

  $h5slicer -x -d 0.1 $tmpfile -s $step -v1 u -v2 vx | sort -gk1 >> $xuvfile
  echo >> $xuvfile
  echo >> $xuvfile
done

# get rid of the evidence
rm $tmpfile

