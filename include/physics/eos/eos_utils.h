/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * EOS_UTILS.H                                                                     *
 *                                                                            *
 * GLOBAL MACROS, FUNCTION DEFINITIONS, INCLUDES, AND DECLARATIONS            *
 *                                                                            *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#if 1
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#endif

// Fundamental constants in CGS
#define EE          (4.80320680e-10  ) // Electron charge 
#define CL          (2.99792458e10   ) // Speed of light
#define ME          (9.1093826e-28   ) // Electron mass
#define MP          (1.67262171e-24  ) // Proton mass
#define MN          (1.67492728e-24  ) // Neutron mass
#define HPL         (6.6260693e-27   ) // Planck constant
#define HBAR        (HPL/(2.*M_PI)   ) // Reduced Planck constant
#define KBOL        (1.3806505e-16   ) // Boltzmann constant
#define GNEWT       (6.6742e-8       ) // Gravitational constant
#define SIG         (5.670400e-5     ) // Stefan-Boltzmann constant
#define AR          (4*SIG/CL        ) // Radiation constant
#define THOMSON     (0.665245873e-24 ) // Thomson cross section
#define COULOMB_LOG (20.             ) // Coulomb logarithm
#define ALPHAFS     (0.007299270073  ) // Fine structure constant ~ 1./137.
#define GFERM       (1.435850814e-49 ) // Fermi constant
#define GA          (-1.272323       ) // Axial-vector coupling
#define GA2         (GA*GA)
#define S2THW       (0.222321) // sin^2(Theta_W), Theta_W = Weinberg angle
#define S4THW       (S2THW*S2THW)
#define NUSIGMA0    (1.7611737037e-44) // Fundamental neutrino cross section

// Unit Conversion factors
#define EV   (1.60217653e-12  ) // Electron-volt
#define MEV  (1.0e6*EV        ) // Mega-Electron-Volt
#define GEV  (1.0e9*EV        ) // Giga-Electron-Volt
#define JY   (1.e-23          ) // Jansky
#define PC   (3.085678e18     ) // Parsec
#define AU   (1.49597870691e13) // Astronomical unit
#define YEAR (31536000.       )
#define DAY  (86400.          )
#define HOUR (3600.           )
#define MSUN (1.989e33        ) // Solar mass


// Macros
// ----------------------------------------------------------------------

// EOS settings
#define EOS (EOS_TYPE_TABLE)
#define GAMMA_FALLBACK (0)
#define POLYTROPE_FALLBACK (1)
#define NVAR_PASSIVE (1) // Ye
#define RADIATION (0)

// Primitive and conserved variables
#define RHO (0)
#define UU  (1)
#define U1  (2)
#define U2  (3)
#define U3  (4)
#define B1  (5)
#define B2  (6)
#define B3  (7)
#define NVAR_BASE (B3 + 1)

// Passive variables (if present)
#define PASSIVE_START (NVAR_BASE)
#define PASSIVE_STOP (NVAR_BASE + NVAR_PASSIVE)
#define PASSTYPE_INTRINSIC (0)
#define PASSTYPE_NUMBER    (1)
#if EOS == EOS_TYPE_TABLE
#define YE (PASSIVE_START)
#endif // EOS_TYPE_TABLE

// EOS
#define EOS_TYPE_GAMMA     (0)
#define EOS_TYPE_POLYTROPE (1)
#define EOS_TYPE_TABLE     (2)
#if EOS == EOS_TYPE_GAMMA
#define EOS_NUM_EXTRA (0)
#define POLYTROPE_FALLBACK (0)
#elif EPS == EOS_TYPE_POLYTROPE
#define EOS_NUM_EXTRA (0)
#define POLYTROPE_FALLBACK (0)
#elif EOS == EOS_TYPE_TABLE
#if GAMMA_FALLBACK
#define POLYTROPE_FALLBACK (0)
#else
#define POLYTROPE_FALLBACK (1)
#endif // GAMMA_FALLBACK
#define EOS_NUM_EXTRA (3)
#define EOS_LRHO (0)
#define EOS_LT   (1)
#define EOS_YE   (2)
// mass fractions
#define NUM_MASS_FRACTIONS (4)
#define MF_XA    (0)
#define MF_XH    (1)
#define MF_XN    (2)
#define MF_XP    (3)
#endif // EOS

// Fixup parameters
// may only apply for EOS GAMMA
#define RHOMINLIMIT (1.e-17)
#define UUMINLIMIT  (1.e-20)
#define RHOMIN      (1.e-5)
#define UUMIN       (1.e-8)
#define BSQORHOMAX  (50.)
#define BSQOUMAX    (2500.)
#define RHOEPS      (2.0)
#define UORHOMAX    (50.)

// Root finding
#define ROOT_SUCCESS (1)
#define ROOT_FAIL    (0)
#define FCOUNT_NBINS (6)
#define FCOUNT_MORE  (FCOUNT_NBINS-1)

// Numerical convenience to represent a small (<< 1) non-zero quantity
#define SMALL (1.e-20)

// Loop over primitive variables
#define PLOOP for(int ip = 0; ip < NVAR; ip++)
#define BASELOOP for (int ip = 0; ip < NVAR_BASE; ip++)

// Loop over extra variables
// TODO: Figure out how to make this conditionally defined.
#define EOS_ELOOP for (int e = 0; e < EOS_NUM_EXTRA; e++)

// Simple macro
#define MY_MAX(a, b) ( ( a > b) ? a : b ) 


// ----------------------------------------------------------------------


// Structs
// ----------------------------------------------------------------------
#if EOS == EOS_TYPE_TABLE
struct of_adiabat {
  double s, ye;
  double lrho_min, lrho_max;
  int imin, imax;
  double hm1_min, hm1_max;
  double* lT;
};
#endif

#if EOS == EOS_TYPE_TABLE
struct of_tablebounds {
  int Nrho, NT, NYe;
  double lrho_min, lrho_max, dlrho;
  double lT_min, lT_max, dlT;
  double Ye_min, Ye_max, dYe;
};
#endif

// ----------------------------------------------------------------------


// Important global variables
// ----------------------------------------------------------------------
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
extern double gam;
#endif
extern double M_unit;
extern double Reh;
extern double Risco;
#if NEED_UNITS
extern double mbh, Mbh, L_unit, T_unit, M_unit, RHO_unit, U_unit, B_unit;
#endif
#if EOS == EOS_TYPE_TABLE
extern double TEMP_unit;
#endif


//This will merge with other stuffs
#if 0
// ----------------------------------------------------------------------


// Public function APIs
// ----------------------------------------------------------------------

// eos.c
void init_EOS();
double EOS_bad_eos_error();
double EOS_get_gamma(const double* extra); 
double EOS_pressure_rho0_u(double rho, double u, const double* extra);

double EOS_enthalpy_rho0_u(double rho, double u, const double* extra);
double EOS_sound_speed_rho0_u(double rho, double u, const double* extra);
void EOS_set_floors(double scale, double rho, double u, double bsq,
  double* rhoflr, double* uflr, const double* extra);
double EOS_entropy_rho0_u(double rho, double u, const double* extra);
double EOS_adiabatic_constant(double rho, double u, const double* extra);
double EOS_temperature(double rho, double u, const double* extra);
double EOS_u_press(double press, double rho, double* extra);

// eos_stellar_collapse.c
#if EOS == EOS_TYPE_TABLE
void EOS_SC_init(char *name);
void EOS_SC_fill(double* rhoIn, double* uIn, double* yeIn, double* restrict eos);
double EOS_SC_pressure_rho0_u(double lrho, double lT, double ye);
double EOS_SC_pressure_rho0_w(double rho, double w, double ye, double *lTold);
double EOS_SC_specific_enthalpy_rho0_u(double lrho, double lT, double ye);
double EOS_SC_sound_speed(double lrho, double lT, double ye);
double EOS_SC_entropy(double lrho, double lT, double ye);
double EOS_SC_gamma(double lrho, double lT, double ye);
double EOS_SC_temperature(double lT);
double EOS_SC_get_u_of_T(double rho, double T, double ye);
double EOS_SC_u_press(double press, double rho, double ye, double *lTold);
void EOS_SC_mass_fractions(double Xi[NUM_MASS_FRACTIONS], const double* extra);
void EOS_SC_avg_ions(double* Abar, double* Zbar, const double* extra);
void EOS_SC_set_floors(double scale, double rho, double u, double ye,
  double bsqr, double* rhoflr, double* uflr);
double EOS_SC_rho_floor(double scale, double bsq);
double EOS_SC_u_floor(double scale, double bsq, double ye);
double EOS_SC_get_min_lrho();
double EOS_SC_get_min_rho();
double EOS_SC_get_min_lT();
double EOS_SC_get_minu(double rho, double ye);
void EOS_SC_get_polytrope(double lrho, double lT, double ye,
  double* poly_K, double* poly_gamma);
double EOS_SC_hm1_min_adiabat(const struct of_adiabat *a);
int EOS_SC_find_adiabat_1d(double s, double ye,
  double lrho_min, double lrho_max, struct of_adiabat *a);
void EOS_SC_print_adiabat(const struct of_adiabat *a);
void EOS_SC_adiabat_free(struct of_adiabat *a);
void EOS_SC_isoentropy_hm1(double hm1, const struct of_adiabat *a,
  double* lrho_guess, double *rho, double *u);
void EOS_SC_get_bounds(struct of_tablebounds *b);
#endif // EOS_TYPE_TABLE


// eos_gamma.c
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
double EOS_Gamma_pressure_rho0_u(double rho, double u);
double EOS_Gamma_pressure_rho0_w(double rho, double w);
double EOS_Gamma_entropy_rho0_u(double rho, double u);
double EOS_Gamma_enthalpy_rho0_u(double rho, double u);
double EOS_Gamma_adiabatic_constant_rho0_u(double rho, double u);
double EOS_Gamma_sound_speed_rho0_u(double rho, double u);
void EOS_Gamma_set_floors(double scale, double rho, double u, double bsq,
  double* rhoflr, double* uflr);
double EOS_Gamma_rho_floor(double scale, double bsq);
double EOS_Gamma_u_floor(double scale, double bsq);
double EOS_Gamma_u_scale(double rho);
double EOS_Gamma_u_press(double press);
double EOS_Gamma_temp(double rho, double u);
#endif // EOS_TYPE_GAMMA

// eos_stellar_collapse.c
#if EOS == EOS_TYPE_TABLE
void EOS_SC_init(char *name);
void EOS_SC_fill(double* rhoIn, double* uIn, double* yeIn, double* restrict eos);
double EOS_SC_pressure_rho0_u(double lrho, double lT, double ye);
double EOS_SC_pressure_rho0_w(double rho, double w, double ye, double *lTold);
double EOS_SC_specific_enthalpy_rho0_u(double lrho, double lT, double ye);
double EOS_SC_sound_speed(double lrho, double lT, double ye);
double EOS_SC_entropy(double lrho, double lT, double ye);
double EOS_SC_gamma(double lrho, double lT, double ye);
double EOS_SC_temperature(double lT);
double EOS_SC_get_u_of_T(double rho, double T, double ye);
double EOS_SC_u_press(double press, double rho, double ye, double *lTold);
void EOS_SC_mass_fractions(double Xi[NUM_MASS_FRACTIONS], const double* extra);
void EOS_SC_avg_ions(double* Abar, double* Zbar, const double* extra);
void EOS_SC_set_floors(double scale, double rho, double u, double ye,
  double bsqr, double* rhoflr, double* uflr);
double EOS_SC_rho_floor(double scale, double bsq);
double EOS_SC_u_floor(double scale, double bsq, double ye);
double EOS_SC_get_min_lrho();
double EOS_SC_get_min_rho();
double EOS_SC_get_min_lT();
double EOS_SC_get_minu(double rho, double ye);
void EOS_SC_get_polytrope(double lrho, double lT, double ye,
  double* poly_K, double* poly_gamma);
double EOS_SC_hm1_min_adiabat(const struct of_adiabat *a);
int EOS_SC_find_adiabat_1d(double s, double ye,
  double lrho_min, double lrho_max, struct of_adiabat *a);
void EOS_SC_print_adiabat(const struct of_adiabat *a);
void EOS_SC_adiabat_free(struct of_adiabat *a);
void EOS_SC_isoentropy_hm1(double hm1, const struct of_adiabat *a,
  double* lrho_guess, double *rho, double *u);
void EOS_SC_get_bounds(struct of_tablebounds *b);
#endif // EOS_TYPE_TABLE
#endif

//Utilities
void *safe_malloc(int size)
{
  // malloc(0) may or may not return NULL, depending on compiler.
  if (size == 0) return NULL;

  void *A = malloc(size);
  if (A == NULL) {
    fprintf(stderr, "Failed to malloc\n");
    exit(-1);
  }
  return A;
}

// Error-handling wrappers for standard C functions
void safe_system(const char *command)
{
  int systemReturn = system(command);
  if (systemReturn == -1) {
    fprintf(stderr, "system() call %s failed! Exiting!\n", command);
    exit(-1);
  }
}

void safe_fscanf(FILE *stream, const char *format, ...)
{
  va_list args;
  va_start(args, format);
  int vfscanfReturn = vfscanf(stream, format, args);
  va_end(args);
  if (vfscanfReturn == -1) {
    fprintf(stderr, "fscanf() call failed! Exiting!\n");
    exit(-1);
  }
}

double find_min(const double* array, int size) {
  double min = INFINITY;
  for (int i = 0; i < size; i++) {
    if (array[i] < min) min = array[i];
  }
  return min;
}

double find_max(const double* array, int size) {
  double max = -INFINITY;
  for (int i = 0; i < size; i++) {
    if (array[i] > max) max = array[i];
  }
  return max;
}

int find_index(double value, const double* array, int size)
{
  for (int i = 0; i < size; i++) {
    if (array[i] >= value) return i;
  }
  return -1;
}

// A very general 1D linear interpolator
double interp_1d(double x,
                 const double xmin, const double xmax,
                 const int imin, const int imax,
                 const double* restrict tab_x,
                 const double* restrict tab_y)
{
  if (x < xmin) x = xmin;
  if (x > xmax) x = xmax-SMALL;
  const int nx = imax - imin;
  const double dx = (xmax - xmin)/(nx-1);
  const int ix = imin + (x - xmin)/dx;
  const double delx = (x - tab_x[ix])/dx;
  const double out = (1-delx)*tab_y[ix] + delx*tab_y[ix+1];
  /*
  // DEBUGGING
  if (isnan(out) || x > xmax || x < xmin) {
    fprintf(stderr,"[interp_1d]: out is NaN!\n");
    fprintf(stderr,"\t\tx           = %g\n",x);
    fprintf(stderr,"\t\timin        = %d\n",imin);
    fprintf(stderr,"\t\timax        = %d\n",imax);
    fprintf(stderr,"\t\txmin        = %e\n",xmin);
    fprintf(stderr,"\t\txmax        = %e\n",xmax);
    fprintf(stderr,"\t\tnx          = %d\n",nx);
    fprintf(stderr,"\t\tdx          = %e\n",dx);
    fprintf(stderr,"\t\tix          = %d\n",ix);
    fprintf(stderr,"\t\tdelx        = %e\n",delx);
    fprintf(stderr,"\t\ttab_x[ix]   = %e\n",tab_x[ix]);
    fprintf(stderr,"\t\ttab_x[ix+1] = %e\n",tab_x[ix+1]);
    fprintf(stderr,"\t\ttab_y[ix]   = %e\n",tab_y[ix]);
    fprintf(stderr,"\t\ttab_y[ix+1] = %e\n",tab_y[ix+1]);
    fprintf(stderr,"\n");
    exit(1);
  }
  */
  return out;
}
#if NEED_UNITS
void set_units()
{
  #if METRIC == MKS
  L_unit = GNEWT*Mbh/(CL*CL);
  #endif
  T_unit = L_unit/CL;
  RHO_unit = M_unit*pow(L_unit,-3.);
  U_unit = RHO_unit*CL*CL;
  B_unit = CL*sqrt(4.*M_PI*RHO_unit);
  #if EOS == EOS_TYPE_TABLE
  TEMP_unit = MEV;
  #endif

  #if RADIATION
  Ne_unit = RHO_unit/(MP + ME);
  kphys_to_num = ME/M_unit;
  //(RADIATION == RADTYPE_LIGHT) ? ME/M_unit : MP/M_unit;
  // kphys_to_num = MBary/M_unit;
  #if ELECTRONS
  Thetae_unit = MP/ME;
  #else
  Thetae_unit = EOS_Theta_unit();
  #endif // ELECTRONS
  #endif // RADIATION

  #if EOS == EOS_TYPE_TABLE && POLYTROPE_FALLBACK
  rho_poly_thresh = EOS_SC_get_min_rho();
  #endif
}
#endif // NEED UNITS


// ----------------------------------------------------------------------
