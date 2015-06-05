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
#include "conservar.h"


/**
 * @brief gatherConservativeVars
 * 
 * @param idim IN
 * @param rowcol IN
 * @param uold IN
 * @param u OUT
 * @param Himin IN
 * @param Himax IN
 * @param Hjmin IN
 * @param Hjmax IN
 * @param Hnvar IN
 * @param Hnxt IN
 * @param Hnyt IN
 * @param Hnxyt IN
 * @return void
 */
void gatherConservativeVars ( const long idim,
                              const long rowcol,
                              double *RESTRICT uold,
                              double *RESTRICT u,
                              const long Himin,
                              const long Himax,
                              const long Hjmin,
                              const long Hjmax,
                              const long Hnvar,
                              const long Hnxt,
                              const long Hnyt,
                              const long Hnxyt ) {

#define IHU(i, j, v)  ((i) + Hnxt * ((j) + Hnyt * (v)))
#define IHVW(i, v) ((i) + (v) * Hnxyt)

    WHERE ( "gatherConservativeVars" );

    long i, j, ivar;

    if ( idim == 1 ) {
        // Gather conservative variables
        //pragma omp for
        for ( i = Himin; i < Himax; i++ ) {
            u[IHVW ( i, ID )] = uold[IHU ( i, rowcol, ID )];
            u[IHVW ( i, IU )] = uold[IHU ( i, rowcol, IU )];
            u[IHVW ( i, IV )] = uold[IHU ( i, rowcol, IV )];
            u[IHVW ( i, IP )] = uold[IHU ( i, rowcol, IP )];
        }

        if ( Hnvar > IP + 1 ) {
            //pragma omp for private(i)
            for ( ivar = IP + 1; ivar < Hnvar; ivar++ ) {
                for ( i = Himin; i < Himax; i++ ) {
                    u[IHVW ( i, ivar )] = uold[IHU ( i, rowcol, ivar )];
                }
            }
        }
    }
    else {
        // Gather conservative variables
        //pragma omp for
        for ( j = Hjmin; j < Hjmax; j++ ) {
            u[IHVW ( j, ID )] = uold[IHU ( rowcol, j, ID )];
            u[IHVW ( j, IU )] = uold[IHU ( rowcol, j, IV )];
            u[IHVW ( j, IV )] = uold[IHU ( rowcol, j, IU )];
            u[IHVW ( j, IP )] = uold[IHU ( rowcol, j, IP )];
        }
        if ( Hnvar > IP + 1 ) {
            //pragma omp for private(j)
            for ( ivar = IP + 1; ivar < Hnvar; ivar++ ) {
                for ( j = Hjmin; j < Hjmax; j++ ) {
                    u[IHVW ( j, ivar )] = uold[IHU ( rowcol, j, ivar )];
                }
            }
        }
    }
#undef IHVW
#undef IHU
}


void updateConservativeVars ( const long idim,
                              const long rowcol,
                              const double dtdx,
                              double *RESTRICT uold,
                              double *RESTRICT u,
                              double *RESTRICT flux,
                              const long Himin,
                              const long Himax,
                              const long Hjmin,
                              const long Hjmax,
                              const long Hnvar,
                              const long Hnxt,
                              const long Hnyt,
                              const long Hnxyt ) {

#define IHU(i, j, v)  ((i) + Hnxt * ((j) + Hnyt * (v)))
#define IHVW(i, v) ((i) + (v) * Hnxyt)

    WHERE ( "updateConservativeVars" );

    long i, j, ivar;

    if ( idim == 1 ) {
        // Update conservative variables
        //pragma omp for
        for ( i = Himin + ExtraLayer; i < Himax - ExtraLayer; i++ ) {
            uold[IHU ( i, rowcol, ID )] =
                u[IHVW ( i, ID )] + ( flux[IHVW ( i - 2, ID )] - flux[IHVW ( i - 1, ID )] ) * dtdx;
            uold[IHU ( i, rowcol, IU )] =
                u[IHVW ( i, IU )] + ( flux[IHVW ( i - 2, IU )] - flux[IHVW ( i - 1, IU )] ) * dtdx;
            uold[IHU ( i, rowcol, IV )] =
                u[IHVW ( i, IV )] + ( flux[IHVW ( i - 2, IV )] - flux[IHVW ( i - 1, IV )] ) * dtdx;
            uold[IHU ( i, rowcol, IP )] =
                u[IHVW ( i, IP )] + ( flux[IHVW ( i - 2, IP )] - flux[IHVW ( i - 1, IP )] ) * dtdx;
            MFLOPS ( 12, 0, 0, 0 );
        }

        if ( Hnvar > IP + 1 ) {
            //pragma omp for private(i)
            for ( ivar = IP + 1; ivar < Hnvar; ivar++ ) {
                for ( i = Himin + ExtraLayer; i < Himax - ExtraLayer; i++ ) {
                    uold[IHU ( i, rowcol, ivar )] =
                        u[IHVW ( i, ivar )] + ( flux[IHVW ( i - 2, ivar )] -
                                                flux[IHVW ( i - 1, ivar )] ) * dtdx;
                    MFLOPS ( 3, 0, 0, 0 );
                }
            }
        }
    } else {
        // Update conservative variables
        //pragma omp for
        for ( j = Hjmin + ExtraLayer; j < Hjmax - ExtraLayer; j++ ) {
            uold[IHU ( rowcol, j, ID )] =
                u[IHVW ( j, ID )] + ( flux[IHVW ( j - 2, ID )] - flux[IHVW ( j - 1, ID )] ) * dtdx;
            uold[IHU ( rowcol, j, IP )] =
                u[IHVW ( j, IP )] + ( flux[IHVW ( j - 2, IP )] - flux[IHVW ( j - 1, IP )] ) * dtdx;
            uold[IHU ( rowcol, j, IV )] =
                u[IHVW ( j, IU )] + ( flux[IHVW ( j - 2, IU )] - flux[IHVW ( j - 1, IU )] ) * dtdx;
            uold[IHU ( rowcol, j, IU )] =
                u[IHVW ( j, IV )] + ( flux[IHVW ( j - 2, IV )] - flux[IHVW ( j - 1, IV )] ) * dtdx;
            MFLOPS ( 12, 0, 0, 0 );
        }
        if ( Hnvar > IP + 1 ) {
            //pragma omp for private(j)
            for ( ivar = IP + 1; ivar < Hnvar; ivar++ ) {
                for ( j = Hjmin + ExtraLayer; j < Hjmax - ExtraLayer; j++ ) {
                    uold[IHU ( rowcol, j, ivar )] =
                        u[IHVW ( j, ivar )] + ( flux[IHVW ( j - 2, ivar )] -
                                                flux[IHVW ( j - 1, ivar )] ) * dtdx;
                    MFLOPS ( 3, 0, 0, 0 );
                }
            }
        }
    }
#undef IHVW
#undef IHU
}















