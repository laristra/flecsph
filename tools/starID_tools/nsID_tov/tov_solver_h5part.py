#!/usr/bin/env python

# Author : Hyun Lim
# Date : June.26.2017
# This is a python script that solves the TOV equation (in Newtonian sense)
# with piecewise polytrope EOS and generates SPH data in H5part format.

import numpy as np
import h5py
import sys
import multiprocessing as multi
from matplotlib import pylab
from scipy.constants import pi, G, c, hbar, m_n #Physical constants (Not needed maybe but keep anyway)
from scipy.interpolate import interp1d
from random import uniform
#from scipy import optimize # UNUSED

# User input parameters
num_particles = int(sys.argv[1]) # Number of particles that you want to obtain in the data
target_num_neighbors = 50 # Number of neighbors. We choose fifty
rho_center = 10**(18) # Center of pressure in MKS
dr = 5.0 # Radial step
filename = 'tov_pp_{}_particles.h5part'.format(num_particles) #Generate h5part file
isentropic = False # Do not much care about isentropic EOS
use_code_units = False
Msun = 1.98892e30 # Mass of Sun in MKS

# Piecewise polytrope equations of state
Gamma0 = 5.0/3.0 # For low density 
Gamma1 = 2.5 # For high density (May use with 3.0)
rho1 = 5e17

K0 = (3.0*pi**2)**(2.0/3.0)*hbar**2/(5.0*m_n**(8.0/3.0))
KR = (3.0*pi**2)**(1.0/3.0)*hbar*c/(4.0*m_n**(4.0/3.0))

P1 = K0*rho1**Gamma0
K1 = P1/rho1**Gamma1

# Some unit conversions to avoid over/under flow problems 
# (HL, TODO : I just made assumptions. We may need to change 
#  the value of this based on our code or not need this part 
#  if we have normalized everything in the code)

# Here is the conversion factors (memory from 2HOT)
# mass_cu = 1.0e20g
# length_cu = 1.0e7cm

length_cu_to_cgs = 1.0e7
length_cgs_to_cu = 1.0e-7
mass_cu_to_cgs = 1.0e20
mass_cgs_to_cu = 1.0e-20
density_cu_to_cgs = 10.0
density_cgs_to_cu = 0.1
specific_energy_cu_to_cgs = 1.0e14
specific_energy_cgs_to_cu = 1.0e-14
pressure_cu_to_cgs = 1.0e13
pressure_cgs_to_cu = 1.0e-13

############ START ##############

print("Here is the beginning")
print("We solve newtonian-TOV equation (simple stellar structure equation a.k.a Poisson form)")
print("with dr = {} and number of particles = {}".format(dr,num_particles))

print("Piecewise polytrope EOS is used. The polytropic values are:")
print(K0, KR, K1)

# 4th order Runge-Kutta routines for integration
def rk4(f,y,x,h):
	k1 = f(y,x)*h
	k2 = f(y+0.5*k1, x+0.5*h)*h
	k3 = f(y+0.5*k2, x+0.5*h)*h
	k4 = f(y+k3, x+h)*h
	return y+k1/6.0+k2/3.0+k3/3.0+k4/6.0

# EOS and its inverse part
def eos(rho):
	if rho < rho1:
	 return K0*rho**Gamma0
	else:
	 return K1*rho**Gamma1

def inv_eos(P):
	if P < P1:
	 return (P/K0)**(1.0/Gamma0)
	else:
	 return (P/K1)**(1.0/Gamma1)

# Define Tolman-Oppenheimer-Volkoff equation in Newtonian gravity
def tov(y,r):
	P, m = y[0], y[1] # Pressure and mass
	rho = inv_eos(P) # Density can be obtained by inverse of EOS
	dPdr = -G*rho*m/r**2
	dmdr = 4.0*pi*rho*r**2
	return pylab.array([dPdr, dmdr])

# For calculation purpose, define derivative of quantities with respect to r
def get_dfdr(f,dr):
	dfdr = np.empty_like(f)
	dfdr[1:-1] = (f[2:]-f[:-2])/(2*dr)
	dfdr[0] = (f[2]-f[0])/dr
	dfdr[-1] = (f[-1]-f[-2])/dr
	return dfdr
# Solve TOV equation within a given central density
def sol_tov(rho_center):
	rmin = dr
	rmax = 20000.0
	r = pylab.arange(rmin,rmax+dr,dr)
	m = pylab.zeros_like(r)
	P = pylab.zeros_like(r)
	m[0] = (4.0/3.0)*pi*rho_center*r[0]**3
	P[0] = eos(rho_center)
	y = pylab.array([P[0],m[0]])
	i = 0
	while P[i]>0.0 and i<len(r)-1:
	 y = rk4(tov,y,r[i],dr)
	 P[i+1] = y[0]
	 m[i+1] = y[1]
	 i = i+1
	rho = pylab.array(list(map(lambda p: inv_eos(p),P)))

	m,r,rho,P = m[:i],r[:i],rho[:i],P[:i] # Give restriction for region

	if isentropic: # Consider the isentropic case
	 # We use conventional finite differencing to handle this case
	 drhodr = get_dfdr(rho,dr)
	 rprime = (r[-1]-r)[::-1] # Inward integration variable
	 # Specific internal energy from thermo identities
	 u_integrand = (-P*drhodr/(rho**2))[::-1]
	 # Interpolation
	 integrand_interpol = iterp1d(rprime,u_integrand,kind='cubic')
	 urhs = lambda u, rprime: integrand_interpol(rprime)
	 u = np.zeros_like(r)
	 for i in range(0,len(rprime)-1):
	  u[i+1] = rk4(urhs,u[i],rprime[i],dr) # Using RK4 to integrate
	  u = u[::-1] # Reverse u
	else:
	 u = np.zeros_like(r)
	return m, m[-1]/Msun, r, rho, P, u # Return the mass and radius of star

# Solving the equations
print("Solving the equations...Please wait...")
mball, M, r, rho, P, u = sol_tov(rho_center)
print("Problem solved!")

# Converting this into cgs unit
print("Converting into cgs")
mball *= 1000.0 #kg to g
r *= 100.0 # m to cm
rho *= 1.0e-3 # kg/m^3 to g/cm^3
P *= 10.0 # pascal to barye
u *= 1.0e10 # J/kg to erg/g

print("Checking... M = {} solar masses and {} grams".format(M, M*Msun*1000))

# Converting with code unit system (HL : we may not need this)
if use_code_units:
 print("Converting into code units")
 mball *= mass_cgs_to_cu
 r *= length_cgs_to_cu
 rho *= density_cgs_to_cu
 P *= pressure_cgs_to_cu
 u *= specific_energy_cgs_to_cu
M = mball[-1]


# Print the results
print("Here is maximum values in chosen units")
print("Mass = {}".format(M))
print("Length = {}".format(r[-1]))
print("Density = {}".format(np.max(rho)))
print("Pressure = {}".format(np.max(P)))
print("Specific energy = {}".format(np.max(u)))

# Rejection sampling part
rmin = r[0]
rmax = r[-1]
rhomin = rho[-1]
rhomax = rho_center
rho_interpol = interp1d(r,rho,kind='cubic')
u_interpol = interp1d(r,u,kind='cubic')
mpart = M/num_particles
rpart = np.empty(num_particles)

rho_1d = 4*np.pi*rho*r**2
rho_1d_interpol = interp1d(r,rho_1d,kind='cubic')
rho_1d_max = np.max(rho_1d)
rho_1d_min = np.min(rho_1d)

def get_position(x):
	while True:
	 rsample = uniform(rmin,rmax)
	 rhosample = uniform(rho_1d_min,rho_1d_max)
	 if rhosample <= rho_1d_interpol(rsample):
	  return rsample

# Make it as particle data to use SPH calculation
p = multi.Pool()
rpart = p.map(get_position, range(num_particles))

# Properties of particles
mparticles = mpart*np.ones_like(rpart) # Mass part
upart = u_interpol(rpart) # Internal energy
vxpart = np.zeros_like(rpart)
vypart = vxpart.copy()
vzpart = vypart.copy()
thetapart = np.array([np.arccos(2*uniform(0,1)-1) for i in range(num_particles)])
phipart = np.array([uniform(0,2*np.pi) for i in range(num_particles)])
xpart = rpart*np.cos(phipart)*np.sin(thetapart)
ypart = rpart*np.sin(phipart)*np.sin(thetapart)
zpart = rpart*np.cos(thetapart)

def solve_for_h(p): # Smoothing length
	x = xpart[p]
	y = ypart[p]
	z = zpart[p]
	rel_distance = np.sqrt((x-xpart)**2 + (y-ypart)**2 + (z-zpart)**2)
	rel_distance.sort(kind = 'quicksort') 
	h = rel_distance[target_num_neighbors-1]
	return h

p = multi.Pool()
hparticles = p.map(solve_for_h, range(num_particles))
norm_hparticles = np.divide(hparticles,rpart)

# Calculating mass in a certain ball
def get_mass_ball(r):
	m = 0 # Initialization
	for p in range(num_particles):
	 if rpart[p] <= r:
	  m += mparticles[p]
	return m
p = multi.Pool()
particle_mball = p.map(get_mass_ball,r)

# Make output
print("Making to output file Please wait...")
with h5py.File(filename,'w') as g:
	f = g.create_group("Step#0")
	xdset = f.create_dataset('x',data = xpart)
	ydset = f.create_dataset('y',data = ypart)
	zdset = f.create_dataset('z',data = zpart)
	udset = f.create_dataset('u',data = upart)
	vxdset = f.create_dataset('vx',data = vxpart)
	vydset = f.create_dataset('vy',data = vypart)
	vydset = f.create_dataset('vz',data = vzpart)
	mdset = f.create_dataset('m',data = mparticles)
	hdset = f.create_dataset('h',data = norm_hparticles)
	mbdset = f.create_dataset('mb',data = particle_mball)
	#r_a_dset = f.create_dataset('r',data = r)
	#rho_a_dset = f.create_dataset('rho',data = rho)
	#m_a_dset = f.create_dataset('gmb',data = mball)
	#P_a_dset = f.create_dataset('P',data = P)
	#u_a_dset = f.create_dataset('gu',data = u)
	
print("Finish! Have fun with your simulations with this data!")
