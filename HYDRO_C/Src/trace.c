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
#include "trace.h"


void trace ( double *RESTRICT q,
             double *RESTRICT dq,
             double *RESTRICT c,
             double *RESTRICT qxm,
             double *RESTRICT qxp,
             const double dtdx,
             const long n,
             const long Hscheme,
             const long Hnvar,
             const long Hnxyt ) {

#define IHVW(i, v) ((i) + (v) * Hnxyt)

    WHERE ( "trace" );

    long ijmin, ijmax;
    long i, IN;
    double zerol = 0.0, zeror = 0.0, project = 0.;

    ijmin = 0;
    ijmax = n;

    // if (strcmp(Hscheme, "muscl") == 0) {       // MUSCL-Hancock method
    if ( Hscheme == HSCHEME_MUSCL ) {   // MUSCL-Hancock method
        zerol = -hundred / dtdx;
        zeror = hundred / dtdx;
        project = one;
        MFLOPS ( 0, 2, 0, 0 );
    }
    // if (strcmp(Hscheme, "plmde") == 0) {       // standard PLMDE
    if ( Hscheme == HSCHEME_PLMDE ) {   // standard PLMDE
        zerol = zero;
        zeror = zero;
        project = one;
    }
    // if (strcmp(Hscheme, "collela") == 0) {     // Collela's method
    if ( Hscheme == HSCHEME_COLLELA ) { // Collela's method
        zerol = zero;
        zeror = zero;
        project = zero;
    }

    // //pragma omp for
    for ( i = ijmin + 1; i < ijmax - 1; i++ ) {
        double cc, csq, r, u, v, p;
        double dr, du, dv, dp;
        double alpham, alphap, alpha0r, alpha0v;
        double spminus, spplus, spzero;
        double apright, amright, azrright, azv1right;
        double apleft, amleft, azrleft, azv1left;
        cc = c[i];
        csq = Square ( cc );
        r = q[IHVW ( i, ID )];
        u = q[IHVW ( i, IU )];
        v = q[IHVW ( i, IV )];
        p = q[IHVW ( i, IP )];
        dr = dq[IHVW ( i, ID )];
        du = dq[IHVW ( i, IU )];
        dv = dq[IHVW ( i, IV )];
        dp = dq[IHVW ( i, IP )];
        alpham = half * ( dp / ( r * cc ) - du ) * r / cc;
        alphap = half * ( dp / ( r * cc ) + du ) * r / cc;
        alpha0r = dr - dp / csq;
        alpha0v = dv;

        // Right state
        spminus = ( u - cc ) * dtdx + one;
        spplus = ( u + cc ) * dtdx + one;
        spzero = u * dtdx + one;
        if ( ( u - cc ) >= zeror ) {
            spminus = project;
        }
        if ( ( u + cc ) >= zeror ) {
            spplus = project;
        }
        if ( u >= zeror ) {
            spzero = project;
        }
        apright = -half * spplus * alphap;
        amright = -half * spminus * alpham;
        azrright = -half * spzero * alpha0r;
        azv1right = -half * spzero * alpha0v;
        qxp[IHVW ( i, ID )] = r + ( apright + amright + azrright );
        qxp[IHVW ( i, IU )] = u + ( apright - amright ) * cc / r;
        qxp[IHVW ( i, IV )] = v + ( azv1right );
        qxp[IHVW ( i, IP )] = p + ( apright + amright ) * csq;

        // Left state
        spminus = ( u - cc ) * dtdx - one;
        spplus = ( u + cc ) * dtdx - one;
        spzero = u * dtdx - one;
        if ( ( u - cc ) <= zerol ) {
            spminus = -project;
        }
        if ( ( u + cc ) <= zerol ) {
            spplus = -project;
        }
        if ( u <= zerol ) {
            spzero = -project;
        }
        apleft = -half * spplus * alphap;
        amleft = -half * spminus * alpham;
        azrleft = -half * spzero * alpha0r;
        azv1left = -half * spzero * alpha0v;
        qxm[IHVW ( i, ID )] = r + ( apleft + amleft + azrleft );
        qxm[IHVW ( i, IU )] = u + ( apleft - amleft ) * cc / r;
        qxm[IHVW ( i, IV )] = v + ( azv1left );
        qxm[IHVW ( i, IP )] = p + ( apleft + amleft ) * csq;
    }
    
    
    if ( Hnvar > IP+1 ) {
        
        // //pragma omp for private(i)
        for ( IN = IP + 1; IN < Hnvar; IN++ ) {
            for ( i = ijmin + 1; i < ijmax - 1; i++ ) {
                double spzero;
                double acmpleft, u, a , da;
                double acmpright;
                u = q[IHVW ( i, IU )];
                a = q[IHVW ( i, IN )];
                da = dq[IHVW ( i, IN )];

                // Right state
                spzero = u * dtdx + one;
                if ( u >= zeror ) {
                    spzero = project;
                }
                acmpright = -half * spzero * da;
                qxp[IHVW ( i, IN )] = a + acmpright;

                // Left state
                spzero = u * dtdx - one;
                if ( u <= zerol ) {
                    spzero = -project;
                }
                acmpleft = -half * spzero * da;
                qxm[IHVW ( i, IN )] = a + acmpleft;
                MFLOPS ( 10, 0, 0, 0 );
            }
        }
    }
}                               // trace

#undef IHVW

//EOF
