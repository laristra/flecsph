/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * EOS_STELLAR_COLLAPSE.h                                                     *
 *                                                                            *
 * PROTOTYPES FOR READING EOS TABLES PROVIDED ON STELLARCOLLAPSE.ORG          *
 *                                                                            *
 ******************************************************************************/

#pragma once

// Save all utilities
#include "eos_utils.h"

#define TABLE_TOL        (1.e-10)
#define TABLE_FTOL       (1.e-10)
#define SC_DEBUG         (0)
#define SC_MONOTONE_SAFE (1)
#define SC_THROTTLE_CS   (0)

#define EOS_ELEM(irho,iT,iY) (Nrho*((iY)*NT + (iT)) + (irho))
#define MMA_ELEM(irho,iY)    (Nrho*iY + irho)

static int Nrho,NT,NYe;
static double *tab_lrho;
static double *tab_lT;
static double *tab_Ye;
static double *tab_lP;
static double *tab_ent;
static double *tab_dpderho;
static double *tab_dpdrhoe;
static double *tab_cs2;
static double *tab_le;
static double *tab_Xa;
static double *tab_Xh;
static double *tab_Xn;
static double *tab_Xp;
static double *tab_Abar;
static double *tab_Zbar;
static double *tab_lwmrho; // log enthalpy - rho, by volume
static double *tab_hm1;    // enthalpy - 1, by mass
static double* tab_poly_gamma; // Polytrope gamma
static double* tab_poly_K; // polytrope K

// min and max of wmrho given fixed ilrho and iY
static double* tab_le_min_2d;
static double* tab_le_max_2d;
static double* tab_lP_min_2d;
static double* tab_lP_max_2d;
static double* tab_lwmrho_min_2d;
static double* tab_lwmrho_max_2d;
static double* tab_hm1_min_1d;

static double tab_lrho_min,tab_lrho_max;
static double tab_lT_min,  tab_lT_max;
static double tab_Ye_min,  tab_Ye_max;
static double tab_dlrho,   tab_dlT,tab_dYe;

static double tab_lP_min,      tab_lP_max;
static double tab_ent_min,     tab_ent_max;
static double tab_cs2_min,     tab_cs2_max;
static double tab_le_min,      tab_le_max;
static double tab_Xa_min,      tab_Xa_max;
static double tab_Xh_min,      tab_Xh_max;
static double tab_Xn_min,      tab_Xn_max;
static double tab_Xp_min,      tab_Xp_max;
static double tab_Abar_min,    tab_Abar_max;
static double tab_Zbar_min,    tab_Zbar_max;
static double tab_dpderho_min, tab_dpderho_max;
static double tab_dpdrhoe_min, tab_dpdrhoe_max;
static double tab_lwmrho_min,  tab_lwmrho_max;

static double tab_rho_min,   tab_rho_max;
static double tab_T_min,     tab_T_max;
static double tab_e_min,     tab_e_max;
static double tab_P_min,     tab_P_max;
static double tab_wmrho_min, tab_wmrho_max;
static double tab_hm1_min,   tab_hm1_max;

static double energy_shift;
static double enthalpy_shift;
// from eos.c
//
double EOS_bad_eos_error()
{
  fprintf(stderr, "ERROR! UNKNOWN EOS TYPE! TYPE = %i\n",EOS);
  exit(1);
}


/*******************************************************************************
      Wrappers
*******************************************************************************/
void init_EOS()
{
  #if EOS == EOS_TYPE_GAMMA || EOS == EOS_TYPE_POLYTROPE
  return;
  // HL : disalbe
  #if 0
  #elif EOS == EOS_TYPE_TABLE
  EOS_SC_init(eospath);
  #else
  #endif
  EOS_bad_eos_error();
  #endif
}

double EOS_pressure_rho0_u(double rho, double u,
                           const double* extra)
{
  double press;
  #if EOS == EOS_TYPE_POLYTROPE
  press = EOS_Poly_pressure_rho0_u(rho,u,poly_K,poly_gam);
  #elif EOS == EOS_TYPE_TABLE
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double ye   = extra[EOS_YE];
  #if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    double K,Gam;
    lrho = EOS_SC_get_min_lrho();
    lT = EOS_SC_get_min_lT();
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    press = EOS_Poly_pressure_rho0_u(rho,u,K,Gam);
  } else {
    press = EOS_SC_pressure_rho0_u(lrho,lT,ye);
  }
  #else
  press = EOS_SC_pressure_rho0_u(lrho,lT,ye);
  #endif // POLYTROPE_FALLBACK
  #else
  EOS_bad_eos_error();
  #endif
  return press;
}

double EOS_enthalpy_rho0_u(double rho, double u, const double* extra)
{
  double enth;
  #if EOS == EOS_TYPE_POLYTROPE
  enth = EOS_Poly_enthalpy_rho0_u(rho,u,poly_K,poly_gam);
  #elif EOS == EOS_TYPE_TABLE
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double ye   = extra[EOS_YE];
  #if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    double K,Gam;
    lrho = EOS_SC_get_min_lrho();
    lT = EOS_SC_get_min_lT();
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    enth = EOS_Poly_enthalpy_rho0_u(rho,MY_MAX(u,0.0),K,Gam);
  } else {
    double h = EOS_SC_specific_enthalpy_rho0_u(lrho,lT,ye);
    enth = h*rho;
  }
  #else
  double h = EOS_SC_specific_enthalpy_rho0_u(lrho,lT,ye);
  enth = h*rho;
  #endif // POLYTROPE_FALLBACK
  #else
  EOS_bad_eos_error();
  #endif
  return enth;
}

double EOS_entropy_rho0_u(double rho, double u, const double* extra)
{
  double ent;
  #if EOS == EOS_TYPE_POLYTROPE
  ent = EOS_Poly_entropy_rho0_u(rho,u,poly_K,poly_gam);
  #elif EOS == EOS_TYPE_TABLE
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double ye   = extra[EOS_YE];
  #if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    double K,Gam;
    lrho = EOS_SC_get_min_lrho();
    lT = EOS_SC_get_min_lT();
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    ent = EOS_Poly_entropy_rho0_u(rho,u,K,Gam);
  } else {
    ent = EOS_SC_entropy(lrho,lT,ye);
  }
  #else
  ent = EOS_SC_entropy(lrho,lT,ye);
  #endif // POLYTROPE_FALLBACK
  #else
  EOS_bad_eos_error();
  #endif
  return ent;
}

double EOS_sound_speed_rho0_u(double rho, double u, const double* extra)
{
  double cs;
  #if EOS == EOS_TYPE_POLYTROPE
  cs = EOS_Poly_sound_speed_rho0_u(rho,u,poly_K,poly_gam);
  #elif EOS == EOS_TYPE_TABLE
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double ye   = extra[EOS_YE];
  #if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    double K,Gam;
    lrho = EOS_SC_get_min_lrho();
    lT = EOS_SC_get_min_lT();
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    cs = EOS_Poly_sound_speed_rho0_u(rho,u,K,Gam);
  } else {
    cs = EOS_SC_sound_speed(lrho,lT,ye);
  }
  #else
  cs = EOS_SC_sound_speed(lrho,lT,ye);
  #endif // POLYTROPE_FALLBACK
  #else
  EOS_bad_eos_error();
  #endif
  return cs;
}

void EOS_set_floors(double scale, double rho, double u, double bsq,
  double* rhoflr, double* uflr, const double* extra)
{
  #if EOS == EOS_TYPE_POLYTROPE
  EOS_Poly_set_floors(scale, rho, u, bsq, rhoflr, uflr);
  #elif EOS == EOS_TYPE_TABLE
  double ye = extra[EOS_YE];
  EOS_SC_set_floors(scale, rho, u, ye, bsq, rhoflr, uflr);
  #else
  EOS_bad_eos_error();
  #endif
}

double EOS_adiabatic_constant(double rho, double u, const double* extra)
{
  double cad;
  #if EOS == EOS_TYPE_POLYTROPE
  EOS_Poly_adiabatic_constant(rho,u,poly_K,poly_gam);
  #elif EOS == EOS_TYPE_TABLE
  double gam = EOS_get_gamma(extra);
  cad = u*pow(rho,-gam);
  #else
  EOS_bad_eos_error();
  #endif
  return cad;
}

double EOS_get_gamma(const double* extra)
{
  #if EOS == EOS_TYPE_POLYTROPE
  return poly_gam;
  #elif EOS == EOS_TYPE_TABLE
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double ye   = extra[EOS_YE];
  double gam  = EOS_SC_gamma(lrho,lT,ye);
  return gam;
  #else
  EOS_bad_eos_error();
  #endif
}

double EOS_temperature(double rho, double u, const double* extra)
{
  #if EOS == EOS_TYPE_POLYTROPE
  return EOS_Poly_temperature(rho,u,poly_K,poly_gam);
  #elif EOS == EOS_TYPE_TABLE
  double lT = extra[EOS_LT];
  return EOS_SC_temperature(lT);
  #else
  EOS_bad_eos_error();
  #endif
}

double EOS_u_press(double press, double rho, double* extra)
{
  double u;
  #if EOS == EOS_TYPE_PLYTROPE
  u = EOS_Poly_u_press(press,rho,poly_K,poly_Gam);
  #elif EOS == EOS_TYPE_TABLE
  double ye = extra[EOS_YE];
  //double yedens = extra[EOS_YE];
  //double ye     = fabs(yedens) / (fabs(rho) + SMALL);
  double lTold  = extra[EOS_LT];
  #if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    u = EOS_SC_get_minu(rho,ye);
    lTold = EOS_SC_get_min_lT();
  } else {
    u = EOS_SC_u_press(press, rho, ye, &lTold);
  }
  #else
  u = EOS_SC_u_press(press, rho, ye, &lTold);
  #endif // POLYTROPE_FALLBACK
  extra[EOS_LT] = lTold;
  #else
  EOS_bad_eos_error();
  #endif // EOS
  return u;
}

#if 0
//HL : func prototype is not needed because we change everything as
//     header but I keep this for reference anyways. This will be 
//     cleaned later
//
// core
static double EOS_SC_interp(const double lrho, const double lT, const double Ye,
			    double* restrict tab);

// max, min, etc.
static void fill_min_1d(double* tab_min_1d, double* tab);
static void fill_max_min_2d(double* tab_max_2d,
			    double* tab_min_2d, double* tab);
static double interp_2d(const double lrho, const double Ye,
			const double* tab_2d);
static void temp_map(double lrho, double Ye, const double* tab);
static double catch_var_2d(const double lrho, const double Ye,
			   const double var,
			   const double* tab_min_2d,
			   const double* tab_max_2d);
static double catch_var_2d_monotone(const double lrho, const double Ye,
				    const double var,
				    double* restrict tab);

// Root finding
static int find_lT(const double lrho, double lTguess, const double ye,
		   double* restrict tab, const double val,
		   double* lT);
static int find_adiabat_0d(double lrho, double lTguess, double ye,
			   double s, double* lT);
static double lT_f(const double lT, const void* params);
struct of_lT_params {
  double lrho, ye;
  double* restrict tab;
};
static double lT_f_adiabat(double lT, const void* params);
struct of_lT_adiabat_params {
  double* restrict tab;
  const struct of_adiabat *a;
};

// utilities
static double le2e(const double le);
static double e2le(const double e);
//static double lw2w(const double lw);
static double w2lw(const double w);
static double catch_rho(const double rho);
static double catch_lrho(const double lrho);
static double catch_e(const double e);
static double catch_press(const double press);
static double catch_w(const double w);
//static double catch_h(const double h);
static double catch_ye(const double ye);
//static double catch_temp(const double temp);
static double catch_lT(const double lT);
static double catch_s(const double s);
static double catch_hm1(const double hm1);
#endif
