# Copyright (c) 2017 Triad National Security, LLC
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

# Code units
p_cgs_to_cu = 1.0e-22
p_cu_to_cgs = 1.0e-22
rho_cgs_to_cu = 1.0e-12
rho_cu_to_cgs = 1.0e12
u_cgs_to_cu = 1.0e-10
u_cu_to_cgs = 1.0e10
distance_cu_to_cgs = 1.0e5
distance_cgs_to_cu = 1.0e-5
v_cu_to_cgs = 1.0e5
v_cgs_to_cu = 1.0e-5
mass_cu_to_cgs = 1.0e27
mass_cgs_to_cu = 1.0e-27

# Extract data by cols
def column(matrix, i):
  return [row[i] for row in matrix]

x_col = column(f1,0)
y_col = column(f1,1)
z_col = column(f1,2)
m_col = column(f1,3)
vx_col = column(f1,4)
vy_col = column(f1,5)
vz_col = column(f1,6)
u_col = column(f1,7)
rho_col = column(f1,8)

x = np.array(x_col)*distance_cgs_to_cu
y = np.array(y_col)*distance_cgs_to_cu
z = np.array(z_col)*distance_cgs_to_cu
m = np.array(m_col)*mass_cgs_to_cu
vx = np.array(vx_col)*v_cgs_to_cu
vy = np.array(vy_col)*v_cgs_to_cu
vz = np.array(vz_col)*v_cgs_to_cu
u = np.array(u_col)*u_cgs_to_cu
rho = np.array(rho_col)*rho_cgs_to_cu


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
