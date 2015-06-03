#ifndef TRACE_H_INCLUDED
#define TRACE_H_INCLUDED

#include "utils.h"

void trace ( double *RESTRICT q, double *RESTRICT dq,
             double *RESTRICT c, double *RESTRICT qxm,
             double *RESTRICT qxp,
             const double dtdx, const long n, const long Hscheme,
             const long Hnvar, const long Hnxyt );

#endif // TRACE_H_INCLUDED
