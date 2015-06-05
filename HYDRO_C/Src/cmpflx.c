/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

#include <math.h>
#include <malloc.h>
// #include <unistd.h>
// #include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "parametres.h"
#include "utils.h"
#include "cmpflx.h"



void cmpflx ( double *RESTRICT qgdnv,
              double *RESTRICT flux,
              const long narray,
              const long Hnxyt,
              const long Hnvar,
              const double Hgamma ) {

#define IHVW(i,v) ((i) + (v) * nxyt)

    WHERE ( "cmpflx" );

    long nface, i, IN;
    double entho, ekin, etot;
    const long nxyt = Hnxyt;

    nface = narray;
    entho = one / ( Hgamma - one );

    // Compute fluxes
    #pragma omp for
    for ( i = 0; i < nface; i++ ) {
        // Mass density
        flux[IHVW ( i, ID )] = qgdnv[IHVW ( i, ID )] * qgdnv[IHVW ( i, IU )];
        // Normal momentum
        flux[IHVW ( i, IU )] = flux[IHVW ( i, ID )] * qgdnv[IHVW ( i, IU )] + qgdnv[IHVW ( i, IP )];
        // Transverse momentum 1
        flux[IHVW ( i, IV )] = flux[IHVW ( i, ID )] * qgdnv[IHVW ( i, IV )];
        // Total energy
        ekin =
            half * qgdnv[IHVW ( i, ID )] * ( Square ( qgdnv[IHVW ( i, IU )] ) +
                                             Square ( qgdnv[IHVW ( i, IV )] ) );
        etot = qgdnv[IHVW ( i, IP )] * entho + ekin;
        flux[IHVW ( i, IP )] = qgdnv[IHVW ( i, IU )] * ( etot + qgdnv[IHVW ( i, IP )] );
    }

    // Other advected quantities
    if ( Hnvar > IP+1 ) {
        #pragma omp for private(i)
        for ( IN = IP + 1; IN < Hnvar; IN++ ) {
            for ( i = 0; i < nface; i++ ) {
                flux[IHVW ( i, IN )] = flux[IHVW ( i, IN )] * qgdnv[IHVW ( i, IN )];
            }
        }
    }
}                               // cmpflx


#undef IHVW

//EOF
