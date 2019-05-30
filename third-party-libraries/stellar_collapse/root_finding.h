/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/
/******************************************************************************
 *                                                                            *
 * ROOT_FINDING.h                                                             *
 *                                                                            *
 ******************************************************************************/

/* Implementation based on gsl root finder API */

#include "decs.h"

#define ROOT_DEBUG       (0)
#define ROOT_VERBOSE     (0)
#define ROOT_NAN_OK      (0)
#define SECANT_NITER_MAX (10)

// PROTOTYPES FOR SHARE
// ----------------------------------------------------------------------
int find_root(double (*f)(const double, const void*),
              const void* params,
              const double ytarget, double xguess,
              const double xmin,    const double xmax,
              const double xtol,    const double ytol,
              double* xroot);
// ----------------------------------------------------------------------
