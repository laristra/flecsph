#!/usr/bin/env python

# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved

import sys, h5py, argparse
import numpy as np

########################
my_description = """
Extracts single particle trajectory from a FleCSPH H5part file."""
my_usage="""
   %(prog)s ifile -i|--id <particle-id>
"""

parser = argparse.ArgumentParser ( description=my_description,
    usage=my_usage,
    epilog="example: ./%(prog)s sodtube_evolution.h5part -i 42")

parser.add_argument('ifile', type=str,
    help='input file the FleCSPH hdf5 format')
parser.add_argument('-i', '--id', action='store', type=int, default=0,
    help='particle id', dest='pid')
args = parser.parse_args()

# read the input file
try:
  h5file = h5py.File(args.ifile)
except:
  sys.exit ("ERROR: cannot read input file %s" % args.ifile)

# get dimensionality of the file
try:
  ndim = h5file.attrs['ndim'][0]
except:
  try:
    ndim = h5file.attrs['dimension'][0]
  except:
    sys.exit ("ERROR: cannot determine file dimensionality")

# access the data
key_step = h5file.keys()[0]
Nsteps = len(h5file.keys())
dset = h5file[key_step]
Npart = dset['x'].len()
if args.pid < 1 or args.pid > Npart:
  sys.exit ("ERROR: wrong particle ID '%d'" % args.pid)

# allocate arrays
ids = np.zeros(Npart,dtype=np.int64)
xyz = np.zeros(Npart)

# write header
print ("# Particle %d" % args.pid)
print ("# 1:iteration 2:time 3:x 4:y 5:z 6:rho 7:P 8:u 9:vx 10:vy 11:vz")
print ("# 12:ax 13:ay 14:az 15:h 16:m 17:dt 18:pnum 19:rank 20:type")

# main output loop
for step in range(Nsteps):
  key_step = ("Step#%d" % step)
  dset = h5file[key_step]
  tm = dset.attrs['time'][0]
  it = dset.attrs['iteration'][0]
  dset['id'].read_direct(ids)
  pn = np.where(ids == args.pid)[0][0]
  
  dset['x'].read_direct(xyz); x = xyz[pn]
  dset['y'].read_direct(xyz); y = xyz[pn]
  dset['z'].read_direct(xyz); z = xyz[pn]
  
  dset['rho'].read_direct(xyz); rho = xyz[pn]
  dset['P'  ].read_direct(xyz); P = xyz[pn]
  dset['u'  ].read_direct(xyz); u = xyz[pn]
  
  dset['vx'].read_direct(xyz); vx = xyz[pn]
  dset['vy'].read_direct(xyz); vy = xyz[pn]
  dset['vz'].read_direct(xyz); vz = xyz[pn]
  
  dset['ax'].read_direct(xyz); ax = xyz[pn]
  dset['ay'].read_direct(xyz); ay = xyz[pn]
  dset['az'].read_direct(xyz); az = xyz[pn]
  
  dset['h'].read_direct(xyz); h = xyz[pn]
  dset['m'].read_direct(xyz); m = xyz[pn]
  dset['dt'].read_direct(xyz); dt = xyz[pn]

  dset['rank'].read_direct(ids); rank = ids[pn]
  dset['type'].read_direct(ids); ptype= ids[pn]

  print ((" % 9d % 14.7e    % 14.7e % 14.7e % 14.7e % 14.7e % 14.7e % 14.7e"+ 
          " % 14.7e % 14.7e % 14.7e % 14.7e % 14.7e % 14.7e % 14.7e % 14.7e"+
          " % 14.7e % 9d % 5d % 3d")
  % (it, tm, x,y,z, rho,P,u, vx,vy,vz, ax,ay,az, 
     h,m,dt, pn,rank,ptype))
  
