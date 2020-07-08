/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * EOS_CONSTS.H                                                               *
 *                                                                            *
 * GLOBAL CONSTANTS                                                           *
 *                                                                            *
 ******************************************************************************/

#pragma once 

#include <math.h>
#include <stdlib.h>

// Fundamental constants in CGS
const double M_SUN_CGS   = 1.98847e33;        // Solar mass [g]
const double C_LIGHT_CGS = 2.99792458e10;     // Speed of light [cm/s]
const double EE          = 4.80320680e-10;    // Electron charge [CGS]
const double ME          = 9.1093826e-28;     // Electron mass [g]
const double MP          = 1.67262171e-24;    // Proton mass [g]
const double MN          = 1.67492728e-24;    // Neutron mass [g]
const double AMU         = 1.66053878283e-24; // Atomic Mass Unit [g/baryon]
const double HPL         = 6.62607015e-27;    // Planck constant [erg*s]
const double HBAR        = HPL/(2*M_PI);      // Reduced Planck constant [erg*s]
const double KBOL        = 1.3806505e-16;     // Boltzmann constant [erg/K]
const double GNEWT       = 6.6742e-8;         // Gravitational constant [cm^3 g^-1 s^-2]
const double SIG         = 5.670400e-5;       // Stefan-Boltzmann constant [erg cm^-2 s^-1 K^-4]
const double AR          = 4*SIG/C_LIGHT_CGS; // Radiation constant [erg cm^-3 K^-4]
const double THOMSON     = 0.665245873e-24;   // Thomson cross section [cm^2]
const double COULOMB_LOG = 20.;               // Coulomb logarithm [.]
const double ALPHAFS     = 0.007299270073;    // Fine structure constant ~ 1./137. [.]
const double GFERM       = 1.435850814e-49;   // Fermi constant [??? TODO: check]
const double GA          = -1.272323;         // Axial-vector coupling [.]
const double GA2         = GA*GA;             // Axial-vector coupling squared [.]
const double S2THW       = 0.222321;          // sin^2(Theta_W), Theta_W = Weinberg angle [.]
const double S4THW       = S2THW*S2THW;       // sin^4(Theta_W), Theta_W = Weinberg angle [.]
const double NUSIGMA0    = 1.7611737037e-44;  // Fundamental neutrino cross section [cm^2]
const double AVO         = 6.0221417930e23;   // Avogadro's number [mol^-1]

// Unit Conversion factors
const double EV   = 1.60217653e-12;   // Electron-volt [erg]
const double MEV  = 1.0e6 * EV;       // Mega-Electron-Volt [erg]
const double GEV  = 1.0e9 * EV;       // Giga-Electron-Volt [erg]
const double JY   = 1.e-23;           // Jansky [erg cm^-2 s^-1 Hz^-1]
const double PC   = 3.085678e18;      // Parsec [cm]
const double AU   = 1.49597870691e13; // Astronomical unit [cm]
const double RSUN = 6.957e+10;        // Solar radius [cm]
const double HOUR = 3600.;            // hour [s]
const double DAY  = 86400.;           // day [s]
const double YEAR = 3.15576e+7;       // Julian year = 365.25 d [s] 

