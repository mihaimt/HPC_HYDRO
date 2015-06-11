/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>

#include "parametres.h"
#include "constoprim.h"
#include "utils.h"




void constoprim ( double *RESTRICT u,
                  double *RESTRICT q,
                  double *RESTRICT e,
                  const long n,
                  const long Hnxyt,
                  const long Hnvar,
                  const double Hsmallr ) {

    long ijmin, ijmax, IN, i;
    double eken;
    const long nxyt = Hnxyt;

    WHERE ( "constoprim" );

    ijmin = 0;
    ijmax = n;

#define IHVW(i,v) ((i) + (v) * nxyt)
    // //pragma omp for
    for ( i = ijmin; i < ijmax; i++ ) {
        q[IHVW ( i, ID )] = MAX ( u[IHVW ( i, ID )], Hsmallr );
        q[IHVW ( i, IU )] = u[IHVW ( i, IU )] / q[IHVW ( i, ID )];
        q[IHVW ( i, IV )] = u[IHVW ( i, IV )] / q[IHVW ( i, ID )];
        eken = half * ( Square ( q[IHVW ( i, IU )] ) + Square ( q[IHVW ( i, IV )] ) );
        q[IHVW ( i, IP )] = u[IHVW ( i, IP )] / q[IHVW ( i, ID )] - eken;
    }
    if ( Hnvar > IP+1 ) {
        // //pragma omp for private(i)
        for ( IN = IP + 1; IN < Hnvar; IN++ ) {
            for ( i = ijmin; i < ijmax; i++ ) {
                q[IHVW ( i, IN )] = u[IHVW ( i, IN )] / q[IHVW ( i, IN )];
            }
        }
    }
    // //pragma omp for
    for ( i = ijmin; i < ijmax; i++ ) {
        e[i] = q[IHVW ( i, IP )];
    }
}                               // constoprim


#undef IHVW
//EOF
