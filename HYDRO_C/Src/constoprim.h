#ifndef CONSTOPRIM_H_INCLUDED
#define CONSTOPRIM_H_INCLUDED

#include "utils.h"

void constoprim(double *RESTRICT u, double *RESTRICT q, double *RESTRICT e,
                const long n, const long Hnxyt, const long Hnvar, const double Hsmallr);

#endif // CONSTOPRIM_H_INCLUDED
