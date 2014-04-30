/* interpolation/hermite.c
 * 
 * Copyright (C) 2014 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * This module contains routines for piecewise cubic Hermite
 * interpolation, following the reference:
 *
 * [1] Burden and Faires, Numerical Analysis, 9th ed
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interpd.h>
#include <gsl/gsl_poly.h>

typedef struct
{
  double *q;         /* Hermite polynomial coefficients */
  double *z;         /* Hermite polynomial z values */
  double *coeff;     /* Taylor coefficients */
  double *work;      /* additional workspace */
  size_t degree;     /* degree of polynomial */
  size_t ncoeff;     /* number of polynomial coefficients */
  size_t npts;       /* number of points needed for piecewise interpolation */
  size_t max_degree; /* maximum possible polynomial degree */
} hermite_state_t;

static void *hermite_alloc(const size_t size, size_t degree);
static int hermite_calc_coeffs(hermite_state_t *state, const double xa[],
                               const double ya[], const double dya[],
                               const size_t size, const double x,
                               gsl_interp_accel *acc);

static void *
hermite_alloc(const size_t size, size_t degree)
{
  hermite_state_t *state;
  const size_t max_degree = 2 * size - 1;

  state = calloc(1, sizeof(hermite_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for state", GSL_ENOMEM);
    }

  if (degree == 0)
    degree = max_degree; /* use maximum possible degree */

  state->degree = degree;
  state->ncoeff = degree + 1;
  state->max_degree = max_degree;

  /* piecewise interpolation window */
  state->npts = (degree - 1) / 2 + 1;

  state->z = malloc(state->ncoeff * sizeof(double));
  if (state->z == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for z", GSL_ENOMEM);
    }

  state->q = malloc(state->ncoeff * sizeof(double));
  if (state->q == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for q", GSL_ENOMEM);
    }

  state->coeff = malloc(state->ncoeff * sizeof(double));
  if (state->coeff == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for coeff", GSL_ENOMEM);
    }

  state->work = malloc(state->ncoeff * sizeof(double));
  if (state->work == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for work", GSL_ENOMEM);
    }

  return state;
}

/* piecewise cubic Hermite polynomials */
static void *
hermite_alloc_cubic(const size_t size)
{
  return hermite_alloc(size, 3);
}

/* maximum degree Hermite polynomial */
static void *
hermite_alloc_all(const size_t size)
{
  return hermite_alloc(size, 0);
}

static void
hermite_free(void * vstate)
{
  hermite_state_t *state = (hermite_state_t *) vstate;

  if (state->z)
    free(state->z);

  if (state->q)
    free(state->q);

  if (state->coeff)
    free(state->coeff);

  if (state->work)
    free(state->work);

  free(state);
} /* hermite_free() */

/*
hermite_init()
  Initialize Hermite polynomial interpolator

If we are using the maximum degree polynomial (degree 2*n - 1),
then the coefficients only need to be computed once for each
dataset and this is done here. For piecewise cubic interpolation
the coefficients need to be computed for each interpolated
point x, and that is done in the evaluation functions.
*/

static int
hermite_init(void * vstate, const double xa[], const double ya[],
             const double dya[], const size_t size)
{
  hermite_state_t *state = (hermite_state_t *) vstate;
  int status = GSL_SUCCESS;

  /*
   * If we're fitting the maximum degree polynomial, we
   * only need to calculate the divided difference representation
   * once. For piecewise cubic Hermite polynomials we need to
   * compute it for each interpolated value in the _eval functions
   */
  if (state->degree == state->max_degree)
    {
      status = gsl_poly_dd_hermite_init(state->q, state->z,
                                        xa, ya, dya, size);
    }

  return status;
} /* hermite_init() */

static int
hermite_eval(const void * vstate, const double xa[], const double ya[], 
             const double dya[], const size_t size, const double x,
             gsl_interp_accel *acc, double *y)
{
  hermite_state_t *state = (hermite_state_t *) vstate;

  if (state->degree != state->max_degree)
    {
      /* calculate cubic Hermite coefficients for this x */
      int status = hermite_calc_coeffs(state, xa, ya, dya, size, x, acc);
      if (status)
        return status;
    }

  /* now evaluate polynomial at the point x */
  *y = gsl_poly_dd_eval(state->q, state->z, state->ncoeff, x);

  return GSL_SUCCESS;
} /* hermite_eval() */

static int
hermite_eval_deriv(const void * vstate, const double xa[], const double ya[], 
                   const double dya[], const size_t size, const double x,
                   gsl_interp_accel *acc, double *dydx)
{
  hermite_state_t *state = (hermite_state_t *) vstate;

  if (state->degree != state->max_degree)
    {
      /* calculate cubic Hermite coefficients for this x */
      int status = hermite_calc_coeffs(state, xa, ya, dya, size, x, acc);
      if (status)
        return status;
    }

  /* now evaluate polynomial derivative at the point x */
  gsl_poly_dd_taylor(state->coeff, x, state->q, state->z,
                     state->ncoeff, state->work);
  *dydx = state->coeff[1];

  return GSL_SUCCESS;
} /* hermite_eval_deriv() */

static int
hermite_eval_deriv2(const void * vstate, const double xa[], const double ya[], 
                    const double dya[], const size_t size, const double x,
                    gsl_interp_accel *acc, double * y_pp)
{
  hermite_state_t *state = (hermite_state_t *) vstate;

  if (state->degree != state->max_degree)
    {
      /* calculate cubic Hermite coefficients for this x */
      int status = hermite_calc_coeffs(state, xa, ya, dya, size, x, acc);
      if (status)
        return status;
    }

  /* now evaluate polynomial derivative at the point x */
  gsl_poly_dd_taylor(state->coeff, x, state->q, state->z,
                     state->ncoeff, state->work);
  *y_pp = 2.0 * state->coeff[2];

  return GSL_SUCCESS;
} /* hermite_eval_deriv2() */

/*
hermite_calc_coeffs()
  This function is called for the case of piecewise cubic
Hermite interpolation. First locate the neighboring points
to the desired point x, then construct coefficients of
a cubic Hermite polynomial matching function values and
derivatives at grid points

Inputs: state - workspace
        xa    - x array data
        ya    - y array data
        dya   - dy/dx array data
        size  - size of arrays
        x     - desired interpolation point
        acc   - accelerator

Notes:

1) On output, the array w->q contains the divided difference
polynomial coefficients, and w->z contains the z array values
(see [1])

2) only cubic interpolation is currently supported
(using the two neighboring grid points surrounding the point x).
Higher order polynomials could easily be implemented by including
more grid points around the point x, but this is not currently
done to keep the code cleaner.
*/

static int
hermite_calc_coeffs(hermite_state_t *state, const double xa[],
                    const double ya[], const double dya[],
                    const size_t size, const double x,
                    gsl_interp_accel *acc)
{
  int s;
  size_t idx;

  /* find idx so that xa[idx] <= x < xa[idx+1] */
  idx = gsl_interp_accel_find(acc, xa, size, x);
  assert((xa[idx] <= x && x < xa[idx + 1]) || idx == size - 2);

  /* compute cubic Hermite polynomial coefficients */
  s = gsl_poly_dd_hermite_init(state->q, state->z,
                               &xa[idx], &ya[idx], &dya[idx],
                               state->npts);

  return s;
} /* hermite_calc_coeffs() */

static const gsl_interpd_type chermite_type =
{
  "chermite",
  3,
  &hermite_alloc_cubic,
  &hermite_init,
  &hermite_eval,
  &hermite_eval_deriv,
  &hermite_eval_deriv2,
  NULL,
  &hermite_free
};

const gsl_interpd_type * gsl_interpd_chermite = &chermite_type;

static const gsl_interpd_type hermite_type =
{
  "hermite",
  2,
  &hermite_alloc_all,
  &hermite_init,
  &hermite_eval,
  &hermite_eval_deriv,
  &hermite_eval_deriv2,
  NULL,
  &hermite_free
};

const gsl_interpd_type * gsl_interpd_hermite = &hermite_type;
