#ifndef CMPFLX_H_INCLUDED
#define CMPFLX_H_INCLUDED

#include "utils.h"
void cmpflx(double *RESTRICT qgdnv, double *RESTRICT flux, const long narray,
            const long Hnxyt, const long Hnvar, const double Hgamma);

#endif // CMPFLX_H_INCLUDED
