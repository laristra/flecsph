#!/usr/bin/env python

# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved


#
# ./1D_reduction.ph -ifile <input.h5part> -step Step#<x>
#
# where <input.h5part> is the input file in hdf5 format and 
# x is step number


import numpy as np
import h5py
import sys
import argparse

my_description = """
Short script to get radial averages of thermodynamic quantities 
from FleCSPH h5part output files. 
Columns in the output: 
 - column1: particle x coord.; 
 - column2: particle y coord.; 
 - column3: particle radial distance; 
 - column4: Density; 
 - column5: Pressure; 
 - column6: Specific internal energy."""
my_usage = "./1D_reduction.py -ifile input.h5part -step Step#10"

def create_output(dataset):
    print (dataset)

    # Get timestep number
    s = dataset[5:]
    i = int(s)

    # Get data from datasets 
    x_data   = file[k]['x']
    y_data   = file[k]['y']
    z_data   = file[k]['z']
    rho_data = file[k]['rho']
    P_data   = file[k]['P']
    u_data   = file[k]['u']
    
    # Coordinates of center 
    x_c= (file['Step#0']['x'].value.min() + file['Step#0']['x'].value.max())/2
    y_c= (file['Step#0']['y'].value.min() + file['Step#0']['y'].value.max())/2
    z_c= (file['Step#0']['z'].value.min() + file['Step#0']['z'].value.max())/2
    
    # Calculate radial distance
    x_rc = x_data - x_c
    y_rc = y_data - y_c
    z_rc = z_data - z_c
    rad_distance = np.sqrt(np.multiply(x_rc, x_rc)
                 + np.multiply(y_rc, y_rc) 
                 + np.multiply(z_rc, z_rc))

    print (i, np.max(rad_distance))

    # Create table
    head = 'column1:x, column2:y, column3:z, column4:r, column5:rho, column6:P, column7:u'
    data = np.column_stack((x_data, y_data, z_data, rad_distance, rho_data, P_data, u_data))

    # Write table
    with open("1d_output_{0:05d}.dat".format(i), 'wb+') as datafile_id:
        np.savetxt(datafile_id, data, fmt='%04.03e', header = head)
    return;


########################

parser = argparse.ArgumentParser(description=my_description, 
         usage=my_usage, 
         formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-i', '--ifile', dest='ifile', type=str, 
                    help='Input FlecSPH hdf5 file', required=True)
parser.add_argument('-s', '--step', dest='step', type=str, 
                    help='Timestep x in the format Step#x')
args = parser.parse_args()


# Read file from input
file = h5py.File(args.ifile,'r')

# Read time step from input
if args.step is not None:
    k = args.step
    create_output(k)
# Otherwise, iterate over time steps
else:
    for k in iter(file.keys()):
        create_output(k)
