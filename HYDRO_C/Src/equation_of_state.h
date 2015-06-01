#ifndef EQUATION_OF_STATE_H_INCLUDED
#define EQUATION_OF_STATE_H_INCLUDED

#include "utils.h"
#include "parametres.h"

void equation_of_state ( double *RESTRICT rho, double *RESTRICT eint,
                         double *RESTRICT p, double *RESTRICT c, long imin,
                         long imax, const double Hsmallc, const double Hgamma );

#endif // EQUATION_OF_STATE_H_INCLUDED
