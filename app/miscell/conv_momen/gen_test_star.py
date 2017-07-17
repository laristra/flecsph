# Copyright (c) 2017 Los Alamos National Security, LLc
# All rights reserved.

#!/usr/bin/env python

# Author : Hyun Lim
# Date : Jul.14.2017
# This is a python script for generating initial data that is a unit
# sphere with unit density and pressure. This is a diagnostic tool for
# checking conservation of momentum.

import sys
import numpy as np
import h5py
from numpy.random import uniform

# The polytrope terms 
POLY_GAMMA= 1.5
POLY_K = 0.5
TOTAL_MASS = 1.0
SPHERE_RADIUS=1.0
NUM_NEIGHBORS=200

nparticles=int(sys.argv[1])
try:
    v_init = float(sys.argv[2])
except:
    v_init = 0

filename='test_star_{}_v_{}_particles.h5part'.format(v_init,nparticles)

# Define properties. We use simple polytrope for pressure etc

def get_p(rho):
    return POLY_K*rho**POLY_GAMMA

def get_v_sound(rho):
    return np.sqrt(get_p(rho)/rho)

def get_u(rho):
    return POLY_K*(rho**(POLY_GAMMA-1))/(POLY_GAMMA-1)

def get_particle_mass(num_particles):
    return TOTAL_MASS/num_particles

def get_volume(r):
    return (4./3.)*np.pi*(r**3)

def get_density(num_particles):
    volume = get_volume(SPHERE_RADIUS)
    density = TOTAL_MASS/volume
    return density

def get_position(x):
    while True:
        rsample = uniform(rmin,rmax)
        rhosample = uniform(rho_1d_min,rho_1d_max)
        if rhosample <= rho1d_interp(rsample):
            return rsample

def get_particle_positions_polar(num_particles):
    rpart = ((3./4.)*(1./np.pi)*uniform(0,SPHERE_RADIUS,num_particles))**(1./3.)
    thetapart=np.arccos(2*uniform(0,1,num_particles)-1.0)
    phipart=uniform(0,2*np.pi,num_particles)
    return rpart,thetapart,phipart

def get_particle_positions_cartesian(rpart,thetapart,phipart):
    xpart=rpart*np.cos(phipart)*np.sin(thetapart)
    ypart=rpart*np.sin(phipart)*np.sin(thetapart)
    zpart=rpart*np.cos(thetapart)
    return xpart,ypart,zpart

def approximate_h(num_particles):
    V = get_volume(SPHERE_RADIUS)
    particle_density = num_particles/V
    h = ((NUM_NEIGHBORS/particle_density)*(3./4.)/np.pi)**(1./3.)
    return h

rpart,thetapart,phipart=get_particle_positions_polar(nparticles)
xpart,ypart,zpart=get_particle_positions_cartesian(rpart,thetapart,phipart)
mpart=get_particle_mass(nparticles)*np.ones_like(xpart)
rhopart=get_density(nparticles)*np.ones_like(xpart)
upart=get_u(rhopart)
ppart=get_p(rhopart)
vrpart = v_init*np.ones_like(rhopart)
vxpart = vrpart*np.sin(thetapart)*np.cos(phipart)
vypart = vrpart*np.sin(thetapart)*np.sin(phipart)
vzpart = vrpart*np.cos(thetapart)
hpart = approximate_h(nparticles)*np.ones_like(xpart)

print("Generating initial data for testing conservation of momentum")
with h5py.File(filename,'w') as g:
    f = g.create_group("Step#0")
    xdset = f.create_dataset('x',data=xpart)
    ydset = f.create_dataset('y',data=xpart)
    zdset = f.create_dataset('z',data=xpart)
    vxdset = f.create_dataset('vx',data=vxpart)
    vydset = f.create_dataset('vy',data=vxpart)
    vzdset = f.create_dataset('vz',data=vxpart)
    mdset = f.create_dataset('m',data=mpart)
    hdset = f.create_dataset('h',data=hpart)
    udset = f.create_dataset('u',data=upart)

print("Finish")
