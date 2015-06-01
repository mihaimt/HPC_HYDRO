/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

// #include <stdlib.h>
// #include <unistd.h>
#include <math.h>
#include <stdio.h>

#include "equation_of_state.h"
// #include "parametres.h"
#include "utils.h"
void
equation_of_state ( double *RESTRICT rho, double *RESTRICT eint,
                    double *RESTRICT p, double *RESTRICT c, long imin,
                    long imax, const double Hsmallc, const double Hgamma ) {
    long k;
    double smallp;
    WHERE ( "equation_of_state" );
    smallp = Square ( Hsmallc ) / Hgamma;
    MFLOPS ( 0, 1, 0, 0 );

    for ( k = imin; k < imax; k++ ) {
        p[k] = ( Hgamma - one ) * rho[k] * eint[k];
        p[k] = MAX ( p[k], ( double ) ( rho[k] * smallp ) );
        c[k] = sqrt ( Hgamma * p[k] / rho[k] );
        MFLOPS ( 5, 2, 1, 0 );
    }
}                               // equation_of_state


// EOF
