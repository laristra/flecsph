#!/usr/bin/env python

# Copyright (c) 2017 Los Alamos National Security, LLC
# All rights reserved

# Converter from FleCSPH hdf5-output to 1D data. Type:
#   ./conv1d.py    # or ./conv1d.py -h
# for usage.

import sys, h5py, argparse
import numpy as np

########################
my_description = """
Reads in a FleCSPH hdf5 file and produces its 1D ASCII slice."""
my_usage="""
   %(prog)s ifile [-h] [-x|-y|-z|-xy|-yz|-xz]
              [-s|--step <val>]
              [-v1|--var1 <val>] [-v2|--var2 <val>]
              [-d <val>] [-x0 <val>] [-y0 <val>] [-z0 <val>]"""

parser = argparse.ArgumentParser ( description=my_description,
    usage=my_usage,
    epilog="example: ./%(prog)s -x sodtube_evolution.h5part -s 10 -v2 u")

parser.add_argument('ifile', type=str,
    help='input file the FleCSPH hdf5 format')
parser.add_argument('-x', action='store_const', const=True, default=False,
    help='output along the x-axis (default)')
parser.add_argument('-y', action='store_const', const=True, default=False,
    help='output along the y-axis')
parser.add_argument('-z', action='store_const', const=True, default=False,
    help='output along the z-axis')
parser.add_argument('-xy', action='store_const', const=True, default=False,
    help='output xy-parallel plane')
parser.add_argument('-yz', action='store_const', const=True, default=False,
    help='output yz-parallel plane')
parser.add_argument('-xz', action='store_const', const=True, default=False,
    help='output xz-parallel plane')
parser.add_argument('-s', '--step', action='store', type=int, default=0,
    help='timestep number (default: 0)', dest='step')
parser.add_argument('-v1', '--var1', action='store', type=str, default='rho',
    help='which variable to output (default: rho)', dest='var1')
parser.add_argument('-v2', '--var2', action='store', type=str, default='P',
    help='which variable to output (default: P)', dest='var2')
parser.add_argument('-d', action='store', type=float, default=0.0,
    help="""only output particles within stripe (in 2D) or a bar (in 3D)
          within +-D from the line/plane (default: output all particles)""")
parser.add_argument('-x0', action='store', type=float, default=0.0,
    help='origin on the x-axis (default: 0)')
parser.add_argument('-y0', action='store', type=float, default=0.0,
    help='origin on the y-axis (default: 0)')
parser.add_argument('-z0', action='store', type=float, default=0.0,
    help='origin on the z-axis (default: 0)')
args = parser.parse_args()

# arguments compatibility check
if not (args.x or args.y or args.z or args.xy or args.yz or args.xz):
  args.x = True

if (args.x + args.y + args.z + args.xy + args.yz + args.xz != 1):
  sys.exit ("ERROR: only one of the [-x|-y|-z|-xy|-yz|-xz] is allowed")

# flag to indicate whether this will be a 1D or 2D output
out1d = args.x or args.y or args.z

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

# sanity checks
if ndim == 1 and (args.y or args.z):
  sys.exit ("ERROR: -y or -z cannot be used with 1D files")

if ndim == 2 and args.z:
  sys.exit ("ERROR: -z cannot be used with 2D files")

if ndim == 1 and not out1d:
  sys.exit ("ERROR: cannot take a 2D slice of a 1D file")

if ndim == 2 and (args.yz or args.xz):
  sys.exit ("ERROR: 2D files can only be sliced with -xy option")

# by default, read the first timestep
key_step = ("Step#%d" % args.step)
if not key_step in h5file:
  sys.exit ("ERROR: cannot find step '%s' in your file" % key_step)

# access the data
dset = h5file[key_step]
try:
  tm = dset.attrs['time'][0]
except:
  tm = 0.0
Npart = dset['x'].len()

x = np.zeros(Npart)
dset['x'].read_direct(x)

if ndim > 1:
  y = np.zeros(Npart)
  dset['y'].read_direct(y)

if ndim > 2:
  z = np.zeros(Npart)
  dset['z'].read_direct(z)

var1 = np.zeros(Npart)
dset[args.var1].read_direct(var1)

var2 = np.zeros(Npart)
dset[args.var2].read_direct(var2)

# main output
print ("# step=%d  time=%12.5e" % (args.step, tm))
if out1d: # 1D slice --------------------------------------------------------
  if args.d == 0.0 or ndim == 1:
    if args.x:
      print ("# 1:x 2:%s 3:%s" % (args.var1, args.var2))
      for i in range(Npart):
        print ("%12.5e %12.5e %12.5e" % (x[i], var1[i], var2[i]))

    if args.y:
      print ("# 1:y 2:%s 3:%s" % (args.var1, args.var2))
      for i in range(Npart):
        print ("%12.5e %12.5e %12.5e" % (y[i], var1[i], var2[i]))

    if args.z:
      print ("# 1:z 2:%s 3:%s" % (args.var1, args.var2))
      for i in range(Npart):
        print ("%12.5e %12.5e %12.5e" % (z[i], var1[i], var2[i]))

  else: # thick slice = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    x0 = args.x0
    y0 = args.y0
    z0 = args.z0
    d  = args.d

    if args.x:
      print ("# 1:x 2:%s 3:%s" % (args.var1, args.var2))
      if ndim == 2:
        for i in range(Npart):
          if abs(y[i] - y0) < d:
            print ("%12.5e %12.5e %12.5e" % (x[i], var1[i], var2[i]))
      else:
        for i in range(Npart):
          if abs(y[i]-y0) < d and abs(z[i]-z0) < d:
            print ("%12.5e %12.5e %12.5e" % (x[i], var1[i], var2[i]))

    if args.y:
      print ("# 1:y 2:%s 3:%s" % (args.var1, args.var2))
      if ndim == 2:
        for i in range(Npart):
          if abs(x[i] - x0) < d:
            print ("%12.5e %12.5e %12.5e" % (y[i], var1[i], var2[i]))
      else:
        for i in range(Npart):
          if abs(x[i]-x0) < d and abs(z[i]-z0) < d:
            print ("%12.5e %12.5e %12.5e" % (y[i], var1[i], var2[i]))

    if args.z:
      print ("# 1:z 2:%s 3:%s" % (args.var1, args.var2))
      for i in range(Npart):
        if abs(x[i]-x0) < d and abs(y[i]-y0) < d:
          print ("%12.5e %12.5e %12.5e" % (z[i], var1[i], var2[i]))

else: # 2D slice -------------------------------------------------------------
  if args.d == 0.0 or ndim == 2:
    if args.xy:
      print ("# 1:x 2:y 3:%s 4:%s" % (args.var1, args.var2))
      for i in range(Npart):
        print ("%12.5e %12.5e %12.5e %12.5e" % (x[i],y[i], var1[i],var2[i]))

    if args.yz:
      print ("# 1:y 2:z 3:%s 4:%s" % (args.var1, args.var2))
      for i in range(Npart):
        print ("%12.5e %12.5e %12.5e %12.5e" % (y[i],z[i], var1[i],var2[i]))

    if args.xz:
      print ("# 1:x 2:z 3:%s 4:%s" % (args.var1, args.var2))
      for i in range(Npart):
        print ("%12.5e %12.5e %12.5e %12.5e" % (x[i],z[i], var1[i],var2[i]))

  else: # thick slice = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    x0 = args.x0
    y0 = args.y0
    z0 = args.z0
    d  = args.d

    if args.xy:
      print ("# 1:x 2:y 3:%s 4:%s" % (args.var1, args.var2))
      for i in range(Npart):
        if abs(z[i] - z0) < d:
          print ("%12.5e %12.5e %12.5e %12.5e" % (x[i],y[i],var1[i],var2[i]))

    if args.yz:
      print ("# 1:y 2:z 3:%s 4:%s" % (args.var1, args.var2))
      for i in range(Npart):
        if abs(x[i] - x0) < d:
          print ("%12.5e %12.5e %12.5e %12.5e" % (y[i],z[i],var1[i],var2[i]))

    if args.xz:
      print ("# 1:x 2:z 3:%s 4:%s" % (args.var1, args.var2))
      for i in range(Npart):
        if abs(y[i] - y0) < d:
          print ("%12.5e %12.5e %12.5e %12.5e" % (x[i],z[i],var1[i],var2[i]))
