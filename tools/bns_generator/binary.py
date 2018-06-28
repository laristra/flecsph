# Copyright (c) 2017 Los Alamos National Security, LLC
# All rights reserved.

#!/usr/bin/env python3

# Author : Hyun Lim
# Date : Jun3.27.2017
# This is a python script that generates binary system in Keperlian orbit

# Here is the packages that I use normally for this
from __future__ import print_function
import numpy as np
import h5py
import sys
import multiprocessing as multi
from matplotlib import pylab
from scipy.interpolate import interp1d
from random import uniform
from scipy import optimize

USE_PN=True
USE_TIDALLY_LOCKED=True
FLIP_Z_1=True
FLIP_Z_2=True

# code units
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

# Important constants
G_newt = 6.67259e-8
c_light = 2.99792458e10

def get_reduced_mass(m1,m2):
    # Reduced mass
    return (m1/(m1+m2))*m2

def get_kepplerian_velocity(m_mine,m_other,r_mine,separation):
    # The tangential velocity for stable circular orbit
    angular_velocity = np.sqrt(G_newt*(m_mine+m_other)/separation**3)
    return angular_velocity*r_mine

def get_center_of_mass(m,x,y,z):
    # Find the center of mass of a star
    mtot = np.sum(m)
    x_com = np.sum(x*m)/mtot
    y_com = np.sum(y*m)/mtot
    z_com = np.sum(z*m)/mtot
    return x_com,y_com,z_com

def center_star(m,x,y,z):
    # Re-centers star so it sits at the origin
    x_com,y_com,z_com = get_center_of_mass(m,x,y,z)
    xnew = x - x_com
    ynew = y - y_com
    znew = z - z_com
    return xnew,ynew,znew

def get_pn_velocity_angular(m_mine,m_other,r_mine,separation):
    # The tangential velocity for initial data for PN2.5
    mtot = m_mine+m_other
    eta = (m_mine/mtot)*(m_other/mtot)
    factor=G_newt*mtot/(separation*c_light**2)
    omega_kepp = np.sqrt(G_newt*mtot/separation**3)
    omega = omega_kepp*(1 - 0.5*factor*(3-eta)
                        + (factor**2)*(15.+47.*eta+3.*(eta**2))/8.)
    vtan = omega*r_mine
    #print("\tvtan = ",vtan[0])
    return vtan

def get_pn_rdot(m_mine,m_other,separation):
    mtot=m_mine+m_other
    eta = (m_mine/mtot)*(m_other/mtot)
    rdot = -(64./5.)*(G_newt**3)*(mtot**2)*m_other*eta/((separation**3)*(c_light**5))
    #print("\trdot = ",rdot)
    return rdot

# DO NOT USE. Should use polar coordinates, not spherical
def get_radius(x,y,z):
    # Radius of a particle at coordinates x,y,z
    return np.sqrt(x**2+y**2+z**2)

def get_offsets(m1,m2,separation):
    # Returns radii of m1 and m2 in binary
    # r1,r2 are DISPLACEMENTS, not magnitudes
    r2 = separation*m1/(m1+m2)
    r1 = r2 - separation 
    return r1,r2

def separate(r1,r2,x1,y1,z1,x2,y2,z2):
    #Separate stars by radius and return new position arrays.
    #By default we separate along x-axis. 
    
    x1new = x1 + r1
    y1new = y1
    z1new = z1
    x2new = x2 + r2
    y2new = y2
    z2new = z2
    return x1new,y1new,z1new,x2new,y2new,z2new

def get_velocities(m_mine,m_other,separation,x,y,z,rstar):
    # Get the velocities of a particle.
    # Tangential velocity is along y-axis.
    # For counterclockwise rotation, vy > 0 iff x > 0
    if USE_TIDALLY_LOCKED:
        r = np.sqrt(x**2+y**2)
        theta = np.arctan2(y,x)
    else:
        r = np.abs(rstar)
        theta = 0 if np.sign(rstar) > 0 else np.pi
    if USE_PN:
        vtan = get_pn_velocity_angular(m_mine,m_other,r,separation)
        vrad = get_pn_rdot(m_mine,m_other,separation)
    else:
        vtan = get_kepplerian_velocity(m_mine,m_other,r,separation)
        vrad = 0
    vx = vrad*np.cos(theta) - vtan*np.sin(theta)
    vy = vrad*np.sin(theta) + vtan*np.cos(theta)
    vz = 0
    return vx,vy,vz

def get_initial_data(separation,file1,file2,output_name):
    # when reading in data, assume particles roughly at COM
    # this may need to change
    print("Reading in files:\n\t{}\n\t{}".format(file1,file2))
    f1 = h5py.File(file1)
    f2 = h5py.File(file2)

    # NOTE: Assumes center of mass is origin for each file. 
    # TODO: Fix this if need be.

    # read in data and convert to cgs
    # currently everything goes into memory
    # TODO: change this if need be
    x1 = np.array(f1['x'])*distance_cu_to_cgs
    x2 = np.array(f2['x'])*distance_cu_to_cgs
    y1 = np.array(f1['y'])*distance_cu_to_cgs
    y2 = np.array(f2['y'])*distance_cu_to_cgs
    z1 = np.array(f1['z'])*distance_cu_to_cgs
    z2 = np.array(f2['z'])*distance_cu_to_cgs
    h1 = np.array(f1['h'])*distance_cu_to_cgs
    h2 = np.array(f2['h'])*distance_cu_to_cgs
    u1 = np.array(f1['u'])*u_cu_to_cgs
    u2 = np.array(f2['u'])*u_cu_to_cgs
    m1 = np.array(f1['mass'])*mass_cu_to_cgs
    m2 = np.array(f2['mass'])*mass_cu_to_cgs

    # Read in particle type and electron fraction if they exist
    try:
        ye1 = np.array(f1['ye'])
    except:
        ye1 = 0.5*np.ones_like(x1)
        print("Warning. No electron fraction for star 1. "
              +"Using a value of ye = 0.5")
    try:
        ye2 = np.array(f2['ye'])
    except:
        ye2 = 0.5*np.ones_like(x2)
        print("Warning. No electron fraction for star 2. "
              +"Using a value of ye = 0.5")
    try:
        parttype1 = np.array(f1['parttype'])
    except:
        parttype1 = np.ones(len(x1),dtype=int)
        print("Warning. No particle type for star 1. "
              +"Using a value of 1.")
    try:
        parttype2 = np.array(f2['parttype'])
    except:
        parttype2 = np.ones(len(x2),dtype=int)
        print("Warning. No particle type for star 2. "
              +"Using a value of 1.")

    # recenter both stars so the are centered exactly at the origin
    # of their respective coordinate systems
    x1,y1,z1 = center_star(m1,x1,y1,z1)
    x2,y2,z2 = center_star(m2,x2,y2,z2)

    # Reflect star 1 about its origin so the stars 
    # are perfectly symmetric about their 
    # center of mass
    x1 *= -1
    y1 *= -1
    if FLIP_Z_1:
        z1 *= -1
    if FLIP_Z_2:
        z2 *= -1

    # total masses and reduced mass
    mtot1,mtot2 = np.sum(m1),np.sum(m2)
    print("Total masses are {} and {} g".format(mtot1,mtot2))
    reduced_mass = get_reduced_mass(mtot1,mtot2)
    print("Reduced mass is {} g".format(reduced_mass))

    # separate the stars
    r1,r2 = get_offsets(m1,m2,separation)
    x1,y1,z1,x2,y2,z2 = separate(r1,r2,x1,y1,z1,x2,y2,z2)

    # get velocities
    print ("Getting Kepplerian velocities.")
    vx1,vy1,vz1 = get_velocities(mtot1,mtot2,separation,x1,y1,z1,r1)
    vx2,vy2,vz2 = get_velocities(mtot2,mtot1,separation,x2,y2,z2,r2)
    vnorm1 = np.sqrt(vx1**2+vy1**2+vz1**2)
    vnorm2 = np.sqrt(vx2**2+vy2**2+vz2**2)
    print("Average velocity per star is {}".format(np.abs(np.mean(vnorm1))))
    print("min velocity is {} and max is {}".format(np.min(np.hstack((vnorm1,vnorm2))),
                                                    np.max(np.hstack((vnorm1,vnorm2)))))

    # get new variables for combined particles
    x = np.hstack((x1,x2))
    y = np.hstack((y1,y2))
    z = np.hstack((z1,z2))
    vx = np.hstack((vx1,vx2))
    vy = np.hstack((vy1,vy2))
    vz = np.hstack((vz1,vz2))
    h = np.hstack((h1,h2))
    m = np.hstack((m1,m2))
    u = np.hstack((u1,u2))
    ye = np.hstack((ye1,ye2))
    parttype = np.hstack((parttype1,parttype2))

    # output to file
    print("Writing data")
    with h5py.File(output_name,'w') as g:
	f = g.create_group("Step#0")
	xdset = f.create_dataset('x', data = x)
	ydset = f.create_dataset('y', data = y)
	zdset = f.create_dataset('z', data = z)
	vxdset = f.create_dataset('vx', data = vx)
	vydset = f.create_dataset('vy', data = vy)
	vzdset = f.create_dataset('vz', data = vz)
	udset = f.create_dataset('u', data = u)
	mdset = f.create_dataset('m', data = m)
	hdset = f.create_dataset('h', data = h)
	yedset = f.create_dataset('ye', data = ye)
	partdset = f.create_dataset('particel-type', data = parttype)

if __name__ == "__main__":
    separation = float(sys.argv[1])*distance_cu_to_cgs
    file1 = sys.argv[2]
    file2 = sys.argv[3]
    output_name = sys.argv[4]
    print("Making binary initial data..")
    get_initial_data(separation,file1,file2,output_name)
    print("Finish! Have fun with your simulation with this data!")

