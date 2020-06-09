#!/usr/bin/env python3

import numpy as np
import h5py
import sys
import multiprocessing as multi
import scipy.optimize as opt

from matplotlib import pylab
from scipy.constants import pi, G, c, hbar, m_n, m_e     # Physical constants
from scipy.interpolate import interp1d
from random import uniform

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

# User input parameters in CGS
rho_center     = float(sys.argv[1]) #1e7                            # Central density [g/cm^3], max. stable value 832742
dr             = float(sys.argv[2]) #1.e6                           # Radial step [cm]
TOV            = str2bool(sys.argv[3])  #False                      # TOV or Newtonian

# Convert physical constants from MKS into CGS
Gr   = G*1.e3                                 # Grav. constant from m^3/(kg*s^2) to cm^3/(g*s^2)
hb   = hbar*1.e7                              # from kg*m^2/s to g*cm^2/s
me   = m_e*1.e3                               # from kg to g
mn   = m_n*1.e3                               # from kg to g            
cs   = c*1.e2                                 # from m/s to cm/s
Msun = 1.98892e33                             # Solar mass [g]

# Proton fraction
Ye = 0.5                                      # Proton fraction

# Polytrope equations of state
#Gamma = 5.0/3.0                             # For non-rel. polytropes
#Gamma = 4.0/3.0                             # For rel. polytropes

# Non-rel K
#K = (hb**2/(15*me*pi**2)) * (3*Ye*pi**2/(mn*cs**2))**(5.0/3.0)
# Rel K
#K = (hb*cs/(12*pi**2)) * (3*Ye*pi**2/(mn*cs**2))**(4.0/3.0)

############ START ##############
print("")
print("")
print("")
print("---- White Dwarf Stellar Structure Solver ----")
print("")
print("Using:")
print("Gravitational constant: {} cm^3/(g*s^2)".format(Gr))
print("Neutron mass:           {} g".format(mn))
print("Electron mass:          {} g".format(me))
print("Speed of light:         {} cm/s".format(cs))
print("Planck constant:        {} g*cm^2/s".format(hb))
print("TOV:                    {}".format(TOV))
#print "Prefactor K in cgs:     {} ".format(K)
print("")

# pressure (difference) for given mass density 
def eos(rho, P):
	 #return K*(rho*cs**2)**Gamma - P
	 A_wd = 6.00288e22
	 B_wd = 9.81011e5/Ye
	 x_wd = (rho/B_wd)**(1.0/3.0)
	 pressure = A_wd*(x_wd*(2.0*x_wd*x_wd-3.0)*
              np.sqrt(x_wd*x_wd+1.0)+3.0*np.arcsinh(x_wd))
	 return (pressure - P)

# energy density for given mass density 
def eos_e(rho):
	A_wd = 6.00288e22
	B_wd = 9.81011e5/Ye
	x_wd = (rho/B_wd)**(1.0/3.0)
	energy_nucleons  = rho*cs*cs                                # nucleons with rho = n*mn*Ye
	energy_electrons = 3.0*A_wd*(x_wd*(2.0*x_wd*x_wd+1.0)*      # electrons
	       np.sqrt(x_wd*x_wd+1.0) - np.arcsinh(x_wd))
	return (energy_nucleons+energy_electrons)

# bisection routine
def bisection(xl, xr, tol, P):
	n = 1
	NMAX = 1e4
	while abs((xl-xr)/2.0) > tol and n < NMAX: 
	 xmid = (xl+xr)/2.0
	 n = n+1
	 a = eos(xmid,P)
	 b = eos(xl,P)
	 if np.sign(a) == np.sign(b):
	  xl = xmid
	 else:
	  xr = xmid
	return xmid 

# 4th order Runge-Kutta routines for integration
def rk4(f,y,x,h,rho):
	k1 = f(y,x,rho)*h
	k2 = f(y + 0.5*k1, x + 0.5*h, rho)*h
	k3 = f(y + 0.5*k2, x + 0.5*h, rho)*h
	k4 = f(y + k3, x + h, rho)*h
	return y + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0

# Define Tolman-Oppenheimer-Volkoff equation in Newtonian gravity
def tov(y,r,rho):
	P, m = y[0], y[1]                                   # pressure and mass
	rho  = bisection(2.0*rho, 0.1*rho, 1.e-8, P)        # find density for given pressure
	en   = eos_e(rho)                                   # determine total energy density 
	if TOV == True:
	  GR1  = (1.0 + P/en)
	  GR2  = (1.0 + 4.0*pi*r**3*P/(m*cs*cs))
	  GR3  = 1.0/(1.0 + 2.0*Gr*m/(cs*cs*r))
	else: 
	  GR1 = 1
	  GR2 = 1
	  GR3 = 1

	dPdr = -Gr*en*m/(cs**2*r**2)*GR1*GR2*GR3            # change in pressure, TOV
	dmdr = 4.0*pi*en*r**2/cs**2                         # change in mass
	return pylab.array([dPdr, dmdr])

# For calculation purpose, define derivative of quantities with respect to r
def get_dfdr(f,dr):
	dfdr       = np.empty_like(f)
	dfdr[1:-1] = (f[2:] - f[:-2])/(2*dr)
	dfdr[0]    = (f[2] - f[0])/dr
	dfdr[-1]   = (f[-1] - f[-2])/dr
	return dfdr

# Solve TOV equation within a given central density
def sol_tov(rho_center):
	rmin = dr
	rmax = 2.e10

	r   = pylab.arange(rmin,rmax+dr,dr)
	m   = pylab.zeros_like(r)
	P   = pylab.zeros_like(r)
	rho = pylab.zeros_like(r)

	i = 0
	rho[i] = rho_center
	m[i]   = (4.0/3.0)*pi*dr**3
	P[i]   = eos(rho_center, 0.0)

	y = pylab.array([P[i],m[i]])    

	while P[i]>0.0 and i<len(r)-1:
	 y = rk4(tov, y, r[i], dr, rho[i])     
	 m[i+1] = y[1]
	 P[i+1] = y[0]
	 rho[i+1] = bisection(2*rho[i], 0.1*rho[i], 1.e-8, P[i+1])
	 i = i+1
     
	if P[i]<0.0:
	 P[i]   = 0.0
	 rho[i] = 0.0

	m,r,rho,P = m[:i],r[:i],rho[:i],P[:i]                   # Give restriction for region

	return m, m[-1]/Msun, r, rho, P                         # Return the mass and radius of star

# Solving the equations
mball, M, r, rho, P = sol_tov(rho_center)

j=0
drhodr = get_dfdr(rho,dr)
i = len(r)-1
print_every = int(round(i/40))
Rtot = r[i]
Mtot = mball[i]

print("Resulting Star:")
print("Mtot: {} g; Rtot {} cm".format(Mtot, Rtot))
print("rho_central: {} g/cm^3".format(rho[0]))
print("P_central: {} g/cm^3".format(P[0]))
print(" ")
print("c1: r/Rtot, c2: rho*Rtot^3/Mtot, c3: m/Mtot, c4: (drho/dr)*Rtot^4/Mtot")
print('{:<18.11E}{:<18.11E}{:<18.11E}{:<18.11E}'.format(0.0, rho[0]*Rtot**3/Mtot, 0.0, 0.0))
for j in range(0, i, print_every):
 print('{:<18.11E}{:<18.11E}{:<18.11E}{:<18.11E}'.format(r[j]/Rtot, rho[j]*Rtot**3/Mtot, mball[j]/Mtot, drhodr[j]*Rtot**4/Mtot))
print('{:<18.11E}{:<18.11E}{:<18.11E}{:<18.11E}'.format(r[i]/Rtot, rho[i]*Rtot**3/Mtot, mball[i]/Mtot, drhodr[i]*Rtot**4/Mtot))
print("Done")

