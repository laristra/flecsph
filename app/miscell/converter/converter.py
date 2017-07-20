# Copyright (c) 2017 Los Alamos National Security, LLC
# All rights reserved. 

#!/usr/bin/env python3

# Author : Hyun Lim
# Date : Jul.13.2017
# Script that converts ASCII file to h5part file

import numpy as np
import h5py
import sys

filename = 'dwd_id.h5part'

print("Reading in file:\n\t{}".format(sys.argv[1]))
f1 = np.loadtxt(sys.argv[1])

# Global variables
global dim
dim = 3 # Depend on the problem

# Extract data by cols
def column(matrix, i):
  return [row[i] for row in matrix]

x = column(f1,0)
y = column(f1,1)
z = column(f1,2)
m = column(f1,3)
vx = column(f1,4)
vy = column(f1,5)
vz = column(f1,6)
u = column(f1,7)
rho = column(f1,8)

# Write the output
print("Generating data...")
with h5py.File(filename,'w') as g:
      f = g.create_group("Step#0")
      f.attrs['dimension'] = 'dim'
      xdset = f.create_dataset('x', data = x)
      ydset = f.create_dataset('y', data = y)
      zdset = f.create_dataset('z', data = z)
      mdset = f.create_dataset('m', data = m)
      vxdset = f.create_dataset('vx', data = vx)
      vydset = f.create_dataset('vy', data = vy)
      vzdset = f.create_dataset('vz', data = vz)
      udset = f.create_dataset('u', data = u)
      rhodset = f.create_dataset('rho', data = rho)
print("Thank you for using!!")
