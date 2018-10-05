/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * DECS.H                                                                     *
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

#include "constants.h"
//#include "params.h" //From Original

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
// TODO: Figure out how to make this conditionally defined. ~JMM
#define EOS_ELOOP for (int e = 0; e < EOS_NUM_EXTRA; e++)
// ----------------------------------------------------------------------


// Structs
// ----------------------------------------------------------------------
// Set global variables that indicate current local metric, etc.
// HL : May not need our code
#if 0
struct of_geom {
  double gcon[NDIM][NDIM];
  double gcov[NDIM][NDIM];
  double g;
  double alpha;
};
#endif

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
#if EOS == EOS_TYPE_POLYTROPE
extern double poly_K,poly_gam;
#endif
#if POLYTROPE_FALLBACK
extern double rho_poly_thresh;
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
// ----------------------------------------------------------------------


// Public function APIs
// ----------------------------------------------------------------------

// eos.c
void init_EOS();
double EOS_bad_eos_error();
double EOS_get_gamma(const double* extra); 
double EOS_pressure_rho0_u(double rho, double u, const double* extra);
//HL : Disable this
#if 0
double EOS_pressure_rho0_w(double rho, double w, double gamma,
			   const struct of_geom *geom,
			   double* extra);
#endif
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
//May not need
#if 0
void do_ye_fixup(int i, int j, int k,
  double pv[NVAR], double pv_prefloor[NVAR]);
#endif
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
#if RADIATION
double EOS_Gamma_Theta_unit();
#endif // RADIATION
#endif // EOS_TYPE_GAMMA

// eos_poly.c
#if EOS == EOS_TYPE_POLYTROPE || POLYTROPE_FALLBACK
double EOS_Poly_pressure_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_pressure_rho0_w(double rho, double w, double K, double Gam);
double EOS_Poly_enthalpy_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_entropy_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_sound_speed_rho0_u(double rho, double u, double K, double Gam);
void EOS_Poly_set_floors(double scale, double rho, double u, double bsq,
  double* rhoflr, double* uflr);
double EOS_Poly_rho_floor(double scale, double bsq);
double EOS_Poly_u_floor(double scale, double bsq);
double EOS_Poly_adiabatic_constant(double rho, double u, double K, double Gam);
#endif // EOS_TYPE_POLY

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
#if 0
void do_ye_fixup(int i, int j, int k,
  double pv[NVAR], double pv_prefloor[NVAR]);
#endif
void EOS_SC_get_bounds(struct of_tablebounds *b);
#endif // EOS_TYPE_TABLE

// util.c
double find_min(const double* array, int size);
double find_max(const double* array, int size);
double interp_1d(double x,
  const double xmin, const double xmax,
  const int imin, const int imax,
  const double* restrict tab_x,
  const double* restrict tab_y);
int find_index(double value, const double* array, int size);
void *safe_malloc(int size);
void safe_system(const char *command);
void safe_fscanf(FILE *stream, const char *format, ...);
#if NEED_UNITS
void set_units();
#endif


// ----------------------------------------------------------------------
