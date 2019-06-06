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
// HL : Func defs from eos.c
// TODO : Remove this with correct order
void init_EOS();
double EOS_pressure_rho0_u(body& b);
double EOS_enthalpy_rho0_u(double rho, double u, const double* extra);
double EOS_entropy_rho0_u(double rho, double u, const double* extra);
double EOS_sound_speed_rho0_u(body& b);
void EOS_set_floors(double scale, double rho, double u, double bsq,
  double* rhoflr, double* uflr, const double* extra);
double EOS_adiabatic_constant(double rho, double u, const double* extra);
double EOS_get_gamma(const double* extra);
double EOS_temperature(double rho, double u, const double* extra);
double EOS_u_press(double press, double rho, double* extra);

//HL : func defs for eos_stellar_collapse.h
//TODO : Re-ordering funcs to remove this
// core
static double EOS_SC_interp(const double lrho, const double lT, const double Ye,
                            double* tab);

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
                                    double* tab);

// Root finding
static int find_lT(const double lrho, double lTguess, const double ye,
                   double* tab, const double val,
                   double* lT);
static int find_adiabat_0d(double lrho, double lTguess, double ye,
                           double s, double* lT);
static double lT_f(const double lT, const void* params);
struct of_lT_params {
  double lrho, ye;
  double* tab;
};
static double lT_f_adiabat(double lT, const void* params);
struct of_lT_adiabat_params {
  double* tab;
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

// SC related function
void EOS_SC_init(char *name);
void EOS_SC_fill(body& b, double* eos);
//void EOS_SC_fill(double* rhoIn, double* uIn, double* yeIn, double* eos);
double EOS_SC_pressure_rho0_u(double lrho, double lT, double ye);
double EOS_SC_specific_enthalpy_rho0_u(double lrho, double lT, double ye);
double EOS_SC_sound_speed(double lrho, double lT, double ye);
double EOS_SC_entropy(double lrho, double lT, double ye);
double EOS_SC_gamma(double lrho, double lT, double ye);
double EOS_SC_temperature(double lT);
double EOS_SC_get_u_of_T(double rho, double T, double ye);
double EOS_SC_pressure_rho0_w(double rho, double w, double ye, double *lTold);
double EOS_SC_u_press(double press, double rho, double ye, double *lTold);
void EOS_SC_mass_fractions(double Xi[NUM_MASS_FRACTIONS], const double* extra);
void EOS_SC_avg_ions(double* Abar, double* Zbar, const double* extra);
void EOS_SC_set_floors(double scale, double rho, double u, double ye,
                       double bsq, double* rhoflr, double* uflr);
double EOS_SC_rho_floor(double scale, double bsq);
double EOS_SC_get_min_lrho();
double EOS_SC_get_min_rho();
double EOS_SC_get_min_lT();
double EOS_SC_get_minu(double rho, double ye);
double EOS_SC_u_floor(double scale, double bsq, double ye);
void EOS_SC_get_polytrope(double lrho, double lT, double ye,
                          double* poly_K, double* poly_gamma);


