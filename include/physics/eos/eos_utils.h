/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * EOS_UTILS.H                                                                     *
 *                                                                            *
 * GLOBAL MACROS, UTILITIES, INCLUDES, AND DECLARATIONS            *
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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include "params.h"

// Fundamental constants in CGS
constexpr double EE =          4.80320680e-10;   // Electron charge
constexpr double CL =          2.99792458e10;    // Speed of light
constexpr double ME =          9.1093826e-28;    // Electron mass
constexpr double MP =          1.67262171e-24;   // Proton mass
constexpr double MN =          1.67492728e-24;   // Neutron mass
constexpr double HPL =         6.6260693e-27;    // Planck constant
constexpr double HBAR =        HPL/(2.*M_PI);    // Reduced Planck constant
constexpr double KBOL =        1.3806505e-16;    // Boltzmann constant
constexpr double GNEWT =       6.6742e-8;        // Gravitational constant
constexpr double SIG =         5.670400e-5;      // Stefan-Boltzmann constant
constexpr double AR =          4*SIG/CL;         // Radiation constant
constexpr double THOMSON =     0.665245873e-24;  // Thomson cross section
constexpr double COULOMB_LOG = 20.;              // Coulomb logarithm
constexpr double ALPHAFS =     0.007299270073;   // Fine structure constant ~ 1./137.
constexpr double GFERM =       1.435850814e-49;  // Fermi constant
constexpr double GA =          -1.272323;        // Axial-vector coupling
constexpr double GA2 =         GA*GA;
constexpr double S2THW =       0.222321;         // sin^2(Theta_W), Theta_W = Weinberg angle
constexpr double S4THW =       S2THW*S2THW;
constexpr double NUSIGMA0 =    1.7611737037e-44; // Fundamental neutrino cross section

// Unit Conversion factors
constexpr double EV =   1.60217653e-12;   // Electron-volt
constexpr double MEV =  1.0e6*EV;         // Mega-Electron-Volt
constexpr double GEV =  1.0e9*EV;         // Giga-Electron-Volt
constexpr double JY =   1.e-23;           // Jansky
constexpr double PC =   3.085678e18;      // Parsec
constexpr double AU =   1.49597870691e13; // Astronomical unit
constexpr double YEAR = 31536000.;
constexpr double DAY =  86400.;
constexpr double HOUR = 3600.;
constexpr double MSUN = 1.989e33;         // Solar mass


// Macros
// ----------------------------------------------------------------------

// Primitive and conserved variables
constexpr int RHO = 0;
constexpr int UU =  1;
constexpr int U1 =  2;
constexpr int U2 =  3;
constexpr int U3 =  4;
constexpr int B1 =  5;
constexpr int B2 =  6;
constexpr int B3 =  7;
constexpr int NVAR_BASE = B3 + 1;

// Passive variables (if present)
#define PASSIVE_START (NVAR_BASE)
#define PASSIVE_STOP (NVAR_BASE + NVAR_PASSIVE)
#define PASSTYPE_INTRINSIC (0)
#define PASSTYPE_NUMBER    (1)
#define YE (PASSIVE_START)

// EOS
#define EOS_TYPE_GAMMA     (0)
#define EOS_TYPE_POLYTROPE (1)
#define EOS_TYPE_TABLE     (2)
#define EOS_NUM_EXTRA (0)
#define EOS_LRHO (0)
#define EOS_LT   (1)
#define EOS_YE   (2)
// mass fractions
#define NUM_MASS_FRACTIONS (4)
#define MF_XA    (0)
#define MF_XH    (1)
#define MF_XN    (2)
#define MF_XP    (3)

// Fixup parameters
// may only apply for EOS GAMMA
constexpr double RHOMINLIMIT = 1.e-17;
constexpr double UUMINLIMIT =  1.e-20;
constexpr double RHOMIN =      1.e-5;
constexpr double UUMIN =       1.e-8;
constexpr double BSQORHOMAX =  50.;
constexpr double BSQOUMAX =    2500.;
constexpr double RHOEPS =      2.0;
constexpr double UORHOMAX =    50.;

// Root finding
constexpr bool ROOT_SUCCESS = true;
constexpr bool ROOT_FAIL = false;
#define FCOUNT_NBINS (6)
#define FCOUNT_MORE  (FCOUNT_NBINS-1)

// Numerical convenience to represent a small (<< 1) non-zero quantity
constexpr double SMALL = 1.e-20;

// Loop over primitive variables
#define PLOOP for(int ip = 0; ip < NVAR; ip++)
#define BASELOOP for (int ip = 0; ip < NVAR_BASE; ip++)

// Loop over extra variables
// TODO: Figure out how to make this conditionally defined.
#define EOS_ELOOP for (int e = 0; e < EOS_NUM_EXTRA; e++)


// ----------------------------------------------------------------------
// Function defs
// TODO : Make it correct order. We define above here because of ordering
double EOS_Poly_pressure_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_pressure_rho0_w(double rho, double w, double K, double Gam);
double EOS_Poly_enthalpy_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_entropy_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_sound_speed_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_rho_floor(double scale, double bsq);
double EOS_Poly_u_floor(double scale, double bsq);
void EOS_Poly_set_floors(double scale, double rho, double u, double bsq,
                         double* rhoflr, double* uflr);
double EOS_Poly_adiabatic_constant(double rho, double u, double K, double Gam);
double EOS_Poly_temperature(double rho, double u, double K, double Gam);
double EOS_Poly_u_press(double press, double rho, double K, double Gam);
double EOS_Poly_Theta_unit();
// HL : this is mostly wrapper with C functions
// TODO : Need to corporate with C++
void *safe_malloc(int size);
void safe_system(const char *command);
void safe_fscanf(FILE *stream, const char *format, ...);
double find_min(const double* array, int size);
double find_max(const double* array, int size);
int find_index(double value, const double* array, int size);
double interp_1d(double x,
                 const double xmin, const double xmax,
                 const int imin, const int imax,
                 const double* tab_x,
                 const double* tab_y);
void set_units();
static int root_secant(double (*f)(const double, const void*),
                       const void* params,
                       const double ytarget, const double xguess,
                       const double xmin,    const double xmax,
                       const double xtol,    const double ytol,
                       double* xroot);
static int root_bisect(double (*f)(const double, const void*),
                       const void* params,
                       const double ytarget, const double xguess,
                       const double xmin,    const double xmax,
                       const double xtol,    const double ytol,
                       double* xroot);
void initialize_root_fcounts();
void print_root_fcounts();
int find_root(double (*f)(const double, const void*),
              const void* params,
              const double ytarget, double xguess,
              const double xmin,    const double xmax,
              const double xtol,    const double ytol,
              double* xroot);






// Structs
// ----------------------------------------------------------------------
struct of_adiabat {
  double s, ye;
  double lrho_min, lrho_max;
  int imin, imax;
  double hm1_min, hm1_max;
  double* lT;
};

struct of_tablebounds {
  int Nrho, NT, NYe;
  double lrho_min, lrho_max, dlrho;
  double lT_min, lT_max, dlT;
  double Ye_min, Ye_max, dYe;
};

// ----------------------------------------------------------------------


// Important global variables
// ----------------------------------------------------------------------
struct GV{
  static double rho_poly_thresh;
  static double poly_K, poly_gam;
  static double Reh;
  static double Risco;
  static double mbh, Mbh, L_unit, T_unit, M_unit, RHO_unit, U_unit, B_unit;
  static double TEMP_unit;
};
// Default values for Global Variables
double GV::B_unit = 1;
double GV::L_unit = 1;
double GV::M_unit = 1;
double GV::T_unit = 1;
double GV::U_unit = 1;
double GV::RHO_unit = 1;
double GV::TEMP_unit = 1;

// Public function APIs
// ----------------------------------------------------------------------

//Adding polytrope eos for fallbakcing
//HL : Since gamma EOS in eos.h is parametrize, I add here for
//     additional polytrope eos but this might be redundant.
//     Possibly, we may want to use gamma law EOS rather than
//     polytrope EOS.
double EOS_Poly_pressure_rho0_u(double rho, double u, double K, double Gam)
{
  rho = fabs(rho+SMALL);
  return K*pow(rho,Gam);
}

double EOS_Poly_pressure_rho0_w(double rho, double w, double K, double Gam)
{
  rho = fabs(rho+SMALL);
  return K*pow(rho,Gam);
}
double EOS_Poly_enthalpy_rho0_u(double rho, double u, double K, double Gam)
{
  double P = EOS_Poly_pressure_rho0_u(rho,u,K,Gam);
  return rho + u + P;
}

double EOS_Poly_entropy_rho0_u(double rho, double u, double K, double Gam)
{
  // TODO: units?
  return K;
}

double EOS_Poly_sound_speed_rho0_u(double rho, double u, double K, double Gam)
{
  rho = std::max(0.0,rho);
  u = std::max(0.0,u);
  double rhogam = K*pow(rho,Gam-1);
  double cs2_num = Gam*rhogam;
  double cs2_den = 1 + u + rhogam;
  double cs2  = cs2_num / cs2_den;
  return sqrt(cs2);
}

double EOS_Poly_rho_floor(double scale, double bsq)
{
  double rhoflr = RHOMIN*scale;
  rhoflr = std::max(rhoflr, RHOMINLIMIT);
  rhoflr = std::max(rhoflr, bsq/BSQORHOMAX);

  return rhoflr;
}

double EOS_Poly_u_floor(double scale, double bsq)
{
  double uflr = UUMIN*scale;
  uflr = std::max(uflr,UUMINLIMIT);
  uflr = std::max(uflr, bsq/BSQOUMAX);

  return uflr;
}

void EOS_Poly_set_floors(double scale, double rho, double u, double bsq,
                         double* rhoflr, double* uflr)
{
  *rhoflr = EOS_Poly_rho_floor(scale,bsq);
  *uflr = EOS_Poly_u_floor(scale,bsq);
}

double EOS_Poly_adiabatic_constant(double rho, double u, double K, double Gam)
{
  return K/(Gam-1.);
}

double EOS_Poly_temperature(double rho, double u, double K, double Gam)
{
  // TODO: is this right?
  return 0.0; // Polytrope is at zero internal energy/temperature
}

double EOS_Poly_u_press(double press, double rho, double K, double Gam)
{
  return 0.0; // pressure and internal energy independent
}

double EOS_Poly_Theta_unit()
{
  return MP/ME;
}

//Utilities
template<typename T>
T *safe_malloc(int size)
{
  size *= sizeof(T);
  // malloc(0) may or may not return NULL, depending on compiler.
  if (size == 0) return NULL;

  T *A = static_cast<T*>(malloc(size));
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
                 const double* tab_x,
                 const double* tab_y)
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

void set_units()
{
  GV::T_unit = GV::L_unit/CL;
  GV::RHO_unit = GV::M_unit*pow(GV::L_unit,-3.);
  GV::U_unit = GV::RHO_unit*CL*CL;
  GV::B_unit = CL*sqrt(4.*M_PI*GV::RHO_unit);
  GV::TEMP_unit = MEV;
}

// Root-finder based on gsl root finder API
//
#define ROOT_DEBUG       (0)
#define ROOT_VERBOSE     (0)
#define ROOT_NAN_OK      (0)
#define SECANT_NITER_MAX (10)

// Define secant and bisection methods in case of simpler
// one fails

static int root_secant(double (*f)(const double, const void*),
                       const void* params,
                       const double ytarget, const double xguess,
                       const double xmin,    const double xmax,
                       const double xtol,    const double ytol,
                       double* xroot)
{
  double dx;
  double x_last, y, yp, ym, dyNum, dyDen, dy;

  double root_fcount[FCOUNT_NBINS]; // TODO : Check this

  double x = xguess;
  unsigned int iter = 0;
  do {
    x_last = x;
    dx = fabs(1.e-7*x) + xtol;
    y = (*f)(x,params) - ytarget;
    yp = (*f)(x+dx,params);
    ym = (*f)(x-dx,params);
    dyNum = yp - ym;
    dyDen = (2.*dx);
    dy = dyNum/dyDen;
    x -= y/dy;
    iter++;
    if ( isnan(x) || isinf(x) ) {
      // can't recover from this
      #if ROOT_DEBUG
      fprintf(stderr,
              "\n\n[root_secant]: NAN or out-of-bounds detected!\n"
              "\txguess  = %.10e\n"
              "\tytarget = %.10e\n"
              "\tx       = %.10e\n"
              "\tx_last  = %.10e\n"
              "\tx_min   = %.10e\n"
              "\tx_max   = %.10e\n"
              "\ty       = %.10e\n"
              "\tdx      = %.10e\n"
              "\typ      = %.10e\n"
              "\tym      = %.10e\n"
              "\tdyNum   = %.10e\n"
              "\tdyDen   = %.10e\n"
              "\tdy      = %.10e\n"
              "\titer    = %d\n"
              "\tsign x  = %d\n",
              xguess,ytarget,
              x,x_last,xmin,xmax,
              y,dx,yp,ym,dyNum,dyDen,dy,
              iter, (int)MY_SIGN(x));
      #endif
      #if ROOT_NAN_OK
      if (isinf(x)) {
        if (x < xmin) x = xmin;
        if (x > xmax) x = xmax;
      } else {
        root_fcount[FCOUNT_MORE]++;
        return ROOT_FAIL;
      }
      #else
      root_fcount[FCOUNT_MORE]++;
      return ROOT_FAIL;
      #endif
      if (x < xmin) x = xmin;
      if (x > xmax) x = xmax;
    }
  } while (iter < SECANT_NITER_MAX && fabs(x-x_last)/(fabs(x)+xtol) > xtol);

  if (iter < FCOUNT_NBINS) root_fcount[iter]++;
  else root_fcount[FCOUNT_MORE]++;

  *xroot = x;

  y = (*f)(x,params);
  const double frac_error = fabs(y-ytarget)/(fabs(y)+ytol);
  #if ROOT_DEBUG
  if ( frac_error > ytol ) {
    fprintf(stderr,
            "\n\n[root_secant]: Failed via too large yerror.\n"
            "\tfractional error = %.10e\n"
            "\tx                = %.10e\n"
            "\ty                = %.10e\n"
            "\tytarget          = %.10e\n"
            "\typ               = %.10e\n"
            "\tym               = %.10e\n"
            "\tdy               = %.10e\n"
            "\tdx               = %.10e\n"
            "\titer             = %d\n",
            frac_error,
            x,y,ytarget,
            yp,ym,dy,dx,
            iter);
  }
  if (fabs(x - x_last) > xtol) {
    fprintf(stderr,
            "\n\n[root_secant]: failed vial dx too big.\n"
            "\tfractional error = %.10e\n"
            "\tx                = %.10e\n"
            "\tx_last           = %.10e\n"
            "\tdx               = %.10e\n"
            "\ty                = %.10e\n"
            "\tytarget          = %.10e\n"
            "\typ               = %.10e\n"
            "\tym               = %.10e\n"
            "\tdy               = %.10e\n"
            "\tdx               = %.10e\n"
            "\titer             = %d\n",
            frac_error,
            x,x_last,
            fabs(x-x_last),
            y,ytarget,
            yp,ym,dy,dx,
            iter);
  }
  #endif

  const int secant_failed = (
                             (fabs(x-x_last) > xtol
                              && fabs(frac_error) > ytol)
                             || isnan(x)
                             || isinf(x)
                             );
  return secant_failed ? ROOT_FAIL : ROOT_SUCCESS;
}

static int root_bisect(double (*f)(const double, const void*),
                       const void* params,
                       const double ytarget, const double xguess,
                       const double xmin,    const double xmax,
                       const double xtol,    const double ytol,
                       double* xroot)
{
  double xl, xr, fl, fr, dx;

  double grow = 0.01;
  double x = xguess;
  if (fabs(x) < xtol) x += xtol;
  do { // Try to find reasonable region for bisection
    dx = fabs(grow*x);
    xl = x - dx;
    xr = x + dx;
    fl = (*f)(xl,params) - ytarget;
    fr = (*f)(xr,params) - ytarget;
    grow *= 1.1;
  } while (fl*fr > 0 && xl >= xmin && xr <= xmax);

  // force back onto the bisection region
  if (xr > xmax) {
    xr = xmax;
    fr = (*f)(xr,params) - ytarget;
  }
  if (xl < xmin) {
    xl = xmin;
    fl = (*f)(xl,params) - ytarget;
  }

  // if they have the same sign, change that.
  // if we can't fix it, fail.
  if (fl*fr > 0) {
    xl = xmin;
    fl = (*f)(xl,params) - ytarget;
    if (fl*fr > 0) {
      xr = xmax;
      fr = (*f)(xr,params) - ytarget;
      if (fl*fr > 0) {
        #if ROOT_DEBUG
        double il = (*f)(xl,params);
        double ir = (*f)(xr,params);
        fprintf(stderr,
                "\n\n[root_bisect]: fl*fr > 0!\n"
                "\txguess  = %.10e\n"
                "\tytarget = %.10e\n"
                "\txl      = %.10e\n"
                "\txr      = %.10e\n"
                "\tfl      = %.10e\n"
                "\tfr      = %.10e\n"
                "\til      = %.10e\n"
                "\tir      = %.10e\n",
                xguess,ytarget,
                xl,xr,
                fl,fr,
                il,ir);
        int nx = 300;
        double dx = (xmax - xmin)/(nx-1);
        fprintf(stderr,"Area map:\nx\ty\n");
        for (int i = 0; i < nx; i++) {
          fprintf(stderr,"%.4f\t%.4e\n",
                  x+i*dx,(*f)(x+i*dx,params));
        }
        #endif
        return ROOT_FAIL;
      }
    }
  }

  do { // bisection algorithm
    double xm = 0.5*(xl + xr);
    double fm = (*f)(xm,params) - ytarget;
    if (fl*fm <= 0) {
      xr = xm;
      fr = fm;
    } else {
      xl = xm;
      fl = fm;
    }
  } while (xr - xl > xtol);

  *xroot = 0.5*(xl+xr);

  #if ROOT_DEBUG
  if (isnan(*xroot)) {
    double il = (*f)(xl,params);
    double ir = (*f)(xr,params);
    fprintf(stderr,
            "\n\n[root_bisect]: NAN DETECTED!\n"
            "\txguess  = %.10e\n"
            "\tytarget = %.10e\n"
            "\txl      = %.10e\n"
            "\txr      = %.10e\n"
            "\tdx      = %.10e\n"
            "\tgrow    = %.10e\n"
            "\txtol    = %.10e\n"
            "\tfl      = %.10e\n"
            "\tfr      = %.10e\n"
            "\til      = %.10e\n"
            "\tir      = %.10e\n"
            "\txmin    = %.10e\n"
            "\txmax    = %.10e\n",
            xguess,ytarget,
            xl,xr,dx,
            grow,xtol,
            fl,fr,
            il,ir,
            xmin,xmax);
  }
  #endif

  return ROOT_SUCCESS;
}

void initialize_root_fcounts()
{
  double root_fcount[FCOUNT_NBINS]; //TODO :Check this
  for (int i = 0; i < FCOUNT_NBINS; i++) root_fcount[i] = 0.0;
}

void print_root_fcounts()
{
  double root_fcount[FCOUNT_NBINS]; //TODO : Check this explicit declaration

  double fcount_tot = 0.0;
  double global_fcount[FCOUNT_NBINS];
  double fcount_percs[FCOUNT_NBINS];

  for (int i = 0; i < FCOUNT_NBINS; i++) {
    global_fcount[i] = root_fcount[i];
    fcount_tot += global_fcount[i];
  }
  for (int i = 0; i < FCOUNT_NBINS; i++) {
    fcount_percs[i] = (100.*global_fcount[i])/fcount_tot;
  }
  #if 0
  if (mpi_io_proc()) {
    fprintf(stdout, "\n********** ROOT FINDING *********\n");
    fprintf(stdout, "   ITERATIONS          PERCENTAGE\n"  );
    for (int i = 0; i < FCOUNT_NBINS; i++) {
      if (i == FCOUNT_NBINS - 1) {
        fprintf(stdout, "         more          %.2e %%\n",
                fcount_percs[i]);
      } else {
        fprintf(stdout, "          %3d          %.2e %%\n",
                i,fcount_percs[i]);
      }
    }
    fprintf(stdout, "*********************************\n\n");
  }
  #endif
  // reset counts
  initialize_root_fcounts();
}

//Root finder routine stars

int find_root(double (*f)(const double, const void*),
              const void* params,
              const double ytarget, double xguess,
              const double xmin,    const double xmax,
              const double xtol,    const double ytol,
              double* xroot)
{
  int status;

  // HL : Where this is originally called? Define here but need to check
  double root_fcount[FCOUNT_NBINS];

  // first check if we're at the max or min values
  const double fmax = (*f)(xmax,params);
  const double errmax = fabs(fmax-ytarget)/(fabs(fmax)+ytol);
  if ( errmax < ytol ) {
    *xroot = xmax;
    root_fcount[0]++;
    return ROOT_SUCCESS;
  }
  const double fmin = (*f)(xmin,params);
  const double errmin = fabs(fmin-ytarget)/(fabs(fmin)+ytol);
  if ( errmin < ytol ) {
    *xroot = xmin;
    root_fcount[0]++;
    return ROOT_SUCCESS;
  }

  if (xguess >= xmax) xguess = xmax-xtol;
  if (xguess < xmin) xguess = xmin;

  // Next try Secant
  status = root_secant(f,params,
                       ytarget,xguess,
                       xmin,xmax,
                       xtol,ytol,
                       xroot);
  if ( status == ROOT_SUCCESS ) return ROOT_SUCCESS;

  #if ROOT_DEBUG
  if (isnan(*xroot)) {
    fprintf(stderr,"xroot is nan after secant\n");
  }
  #endif

  // Secant failed. Try bisection.
  #if ROOT_VERBOSE
  fprintf(stderr,
          "\n\nRoot finding. Secant failed. Trying bisection.\n"
          "\txguess  = %.10g\n"
          "\tytarget = %.10g\n"
          "\txmin    = %.10g\n"
          "\txmax    = %.10g\n",
          xguess,ytarget,xmin,xmax);
  #endif
  status = root_bisect(f, params, ytarget,
                       xguess,
                       xmin, xmax,
                       xtol, ytol,
                       xroot);
  // Check for something horrible happening
  #if ROOT_DEBUG
  if ( isnan(*xroot) || isinf(*xroot) ) {
    fprintf(stderr,"xroot is nan after bisection\n");
    return ROOT_FAIL;
  }
  #endif

  return status;
}
