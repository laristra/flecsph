#!/usr/bin/env python3

# Author : Hyun Lim
# Date : Jul.13.2017
# Script that converts ASCII file to h5part

import numpy as np
import h5py
import sys

filename = 'dwd_id.h5part'

print("Reading in files:\n\t{}".format(sys.argv[1]))
f1 = np.loadtxt(sys.argv[1])

# Read the data
x = np.array(f1['x'])
y = np.array(f1['y'])
z = np.array(f1['z'])
m = np.array(f1['m'])
vx = np.array(f1['vx'])
vy = np.array(f1['vy'])
vz = np.array(f1['vz'])
u = np.array(f1['u'])
rho = np.array(f1['rho'])

# Right the output
print("Generating data...")
with h5py.File(filename,'w') as g:
      f = g.create_group("Step#0")
      xdset = f.create_dataset('x', data = x)
      ydset = f.create_dataset('y', data = x)
      zdset = f.create_dataset('z', data = x)
      mdset = f.create_dataset('m', data = x)
      vxdset = f.create_dataset('vx', data = x)
      vydset = f.create_dataset('vy', data = x)
      vzdset = f.create_dataset('vz', data = x)
      udset = f.create_dataset('u', data = x)
      rhodset = f.create_dataset('rho', data = x)
print("Thank you for using!!")
