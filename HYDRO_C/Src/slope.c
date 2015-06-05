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
#include "utils.h"
#include "slope.h"

//#ifdef HMPP
// #define DABS(x) (double) ((x)>0?(x):-(x))
// #else
#define DABS(x) (double) fabs((x))
// #endif



void slope ( double *RESTRICT q,
             double *RESTRICT dq,
             const long narray,
             const long Hnvar,
             const long Hnxyt,
             const double slope_type ) {

    long n, i, ijmin, ijmax;

#define IHVW(i, v) ((i) + (v) * Hnxyt)

    WHERE ( "slope" );

    ijmin = 0;
    ijmax = narray;

    #pragma omp for private(i)
    for ( n = 0; n < Hnvar; n++ ) {
        for ( i = ijmin + 1; i < ijmax - 1; i++ ) {
            double dlft, drgt, dcen, dsgn, slop, dlim;
            long ihvwin, ihvwimn, ihvwipn;
            ihvwin = IHVW ( i, n );
            ihvwimn = IHVW ( i - 1, n );
            ihvwipn = IHVW ( i + 1, n );
            dlft = slope_type * ( q[ihvwin] - q[ihvwimn] );
            drgt = slope_type * ( q[ihvwipn] - q[ihvwin] );
            dcen = half * ( dlft + drgt ) / slope_type;
            dsgn = ( dcen > 0 ) ? ( double ) 1.0 : ( double ) -1.0; // sign(one, dcen);
            slop = DABS ( drgt );
            slop = ( double ) MIN ( DABS ( dlft ), slop );
            dlim = slop;
            if ( ( dlft * drgt ) <= zero ) {
                dlim = zero;
            }
            dq[ihvwin] = dsgn * ( double ) MIN ( dlim, DABS ( dcen ) );
        }
    }
} // slope

#undef IHVW











