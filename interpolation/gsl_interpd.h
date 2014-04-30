/* interpolation/gsl_interpd.h
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

#ifndef __GSL_INTERPD_H__
#define __GSL_INTERPD_H__

#include <stdlib.h>
#include <gsl/gsl_inline.h>
#include <gsl/gsl_types.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* interpolation object type */
typedef struct {
  const char * name;
  unsigned int min_size;
  void *  (*alloc) (const size_t size);
  int     (*init)  (void *, const double xa[], const double ya[], const double dya[], const size_t size);
  int     (*eval)  (const void *, const double xa[], const double ya[], const double dya[], const size_t size, const double x, gsl_interp_accel *, double * y);
  int     (*eval_deriv)  (const void *, const double xa[], const double ya[], const double dya[], const size_t size, const double x, gsl_interp_accel *, double * y_p);
  int     (*eval_deriv2) (const void *, const double xa[], const double ya[], const double dya[], const size_t size, const double x, gsl_interp_accel *, double * y_pp);
  int     (*eval_integ)  (const void *, const double xa[], const double ya[], const double dya[], const size_t size, gsl_interp_accel *, const double a, const double b, double * result);
  void    (*free)         (void *);

} gsl_interpd_type;

/* general interpolation object */
typedef struct {
  const gsl_interpd_type * type;
  double  xmin;
  double  xmax;
  size_t  size;
  void * state;
} gsl_interpd;

/* available types */
GSL_VAR const gsl_interpd_type * gsl_interpd_hermite;

gsl_interpd *gsl_interpd_alloc(const gsl_interpd_type * T, const size_t n);
     
int gsl_interpd_init(gsl_interpd * obj, const double xa[],
                     const double ya[], const double dya[],
                     const size_t size);

const char * gsl_interpd_name(const gsl_interpd * interp);
unsigned int gsl_interpd_min_size(const gsl_interpd * interp);
unsigned int gsl_interpd_type_min_size(const gsl_interpd_type * T);

int gsl_interpd_eval_e(const gsl_interpd * obj,
                       const double xa[], const double ya[],
                       const double dya[], const double x,
                       gsl_interp_accel * a, double * y);

double gsl_interpd_eval(const gsl_interpd * obj,
                        const double xa[], const double ya[],
                        const double dya[], const double x,
                        gsl_interp_accel * a);

int gsl_interpd_eval_deriv_e(const gsl_interpd * obj,
                             const double xa[], const double ya[],
                             const double dya[], const double x,
                             gsl_interp_accel * a,
                             double * d);

double gsl_interpd_eval_deriv(const gsl_interpd * obj,
                              const double xa[], const double ya[],
                              const double dya[], const double x,
                              gsl_interp_accel * a);

int gsl_interpd_eval_deriv2_e(const gsl_interpd * obj,
                              const double xa[], const double ya[],
                              const double dya[], const double x,
                              gsl_interp_accel * a,
                              double * d2);

double gsl_interpd_eval_deriv2(const gsl_interpd * obj,
                               const double xa[], const double ya[],
                               const double dya[], const double x,
                               gsl_interp_accel * a);

int gsl_interpd_eval_integ_e(const gsl_interpd * obj,
                             const double xa[], const double ya[],
                             const double dya[],
                             const double a, const double b,
                             gsl_interp_accel * acc,
                             double * result);

double gsl_interpd_eval_integ(const gsl_interpd * obj,
                              const double xa[], const double ya[],
                              const double dya[],
                              const double a, const double b,
                              gsl_interp_accel * acc);

void gsl_interpd_free(gsl_interpd * interp);

__END_DECLS

#endif /* __GSL_INTERPD_H__ */
