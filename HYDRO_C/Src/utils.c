/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>

#include "utils.h"
// #include "parametres.h"
double **
allocate ( long imin, long imax, long nvar ) {
    long i;

#ifdef FAST
    double **r = ( double ** ) malloc ( nvar * sizeof ( double * ) );

#else /*  */
    double **r = ( double ** ) calloc ( nvar, sizeof ( double * ) );

#endif /*  */
    assert ( r != NULL );
    for ( i = 0; i < nvar; i++ ) {
        r[i] = DMalloc ( imax - imin + 1 + MallocGuard );
    }
    return r;
}

double *
DMalloc ( long n ) {

#ifdef FAST
    double *r = ( double * ) malloc ( ( n + MallocGuard ) * sizeof ( double ) );

#else /*  */
    double *r = ( double * ) calloc ( ( n + MallocGuard ), sizeof ( double ) );

#endif /*  */
    assert ( r != NULL );
    return r;
}

long *
IMalloc ( long n ) {

#ifdef FAST
    long *r = ( long * ) malloc ( ( n + MallocGuard ) * sizeof ( long ) );

#else /*  */
    long *r = ( long * ) calloc ( ( n + MallocGuard ), sizeof ( long ) );

#endif /*  */
    assert ( r != NULL );
    return r;
}


#include "parametres.h"
#define VALPERLINE 11
void
printuold ( const hydroparam_t H, hydrovar_t * Hv ) {
    long i, j, nvar;
    for ( nvar = 0; nvar < H.nvar; nvar++ ) {
        fprintf ( stdout, "=uold %ld >\n", nvar );
        for ( j = 0; j < H.nyt; j++ ) {
            long nbr = 1;
            for ( i = 0; i < H.nxt; i++ ) {
                // fprintf(stdout, "%13.6e ", Hv->uold[IHv(i, j, nvar)]);
                fprintf ( stdout, "%10.3e ", Hv->uold[IHv ( i, j, nvar )] );
                nbr++;
                if ( nbr == VALPERLINE ) {
                    fprintf ( stdout, "\n" );
                    nbr = 1;
                }
            }
            if ( nbr != 1 )
                fprintf ( stdout, "\n" );
            fprintf ( stdout, "%%\n" );
        }
    }
}
void
printarray ( double *a, long n, const char *nom ) {
    long i, nbr = 1;
    fprintf ( stdout, "=%s >\n", nom );
    for ( i = 0; i < n; i++ ) {
        fprintf ( stdout, "%13.6e ", a[i] );
        nbr++;
        if ( nbr == VALPERLINE ) {
            fprintf ( stdout, "\n" );
            nbr = 1;
        }
    }
    if ( nbr != 1 )
        fprintf ( stdout, "\n" );
}

void
printarrayi ( long *a, long n, const char *nom ) {
    long i, nbr = 1;
    fprintf ( stdout, "=%s >\n", nom );
    for ( i = 0; i < n; i++ ) {
        fprintf ( stdout, "%4ld ", a[i] );
        nbr++;
        if ( nbr == VALPERLINE ) {
            fprintf ( stdout, "\n" );
            nbr = 1;
        }
    }
    if ( nbr != 1 )
        fprintf ( stdout, "\n" );
}

void
printarrayv ( double *a, long n, const char *nom, const hydroparam_t H ) {
    long i, nbr = 1;
    long nvar;
    fprintf ( stdout, "=%s >\n", nom );
    for ( nvar = 0; nvar < H.nvar; nvar++ ) {
        nbr = 1;
        for ( i = 0; i < n; i++ ) {
            fprintf ( stdout, "%13.6e ", a[IHvw ( i, nvar )] );
            nbr++;
            if ( nbr == VALPERLINE ) {
                fprintf ( stdout, "\n" );
                nbr = 1;
            }
        }
        if ( nbr != 1 )
            fprintf ( stdout, "\n" );
        fprintf ( stdout, "---\n" );
    }
}
void
timeToString ( char *buf, const double timeInS ) {
    char ctenth[10];
    long hour = timeInS / 3600;
    long minute = ( timeInS - hour * 3600 ) / 60;
    long second = timeInS - hour * 3600 - minute * 60;
    double tenth = timeInS - hour * 3600 - minute * 60 - second;
    sprintf ( ctenth, "%.3lf", tenth );
    sprintf ( buf, "%02ld:%02ld:%02ld%s", hour, minute, second, &ctenth[1] );
}
double
cclock ( void ) {
    const double micro = 1.0e-06;       /* Conversion constant */
    static long start = 0L, startu;
    struct timeval tp;          /* Structure used by gettimeofday */
    double wall_time;           /* To hold the result */
    if ( gettimeofday ( &tp, NULL ) == -1 )
        wall_time = -1.0e0;

    else if ( !start ) {
        start = tp.tv_sec;
        startu = tp.tv_usec;
        wall_time = 0.0e0;
    } else
        wall_time = ( double ) ( tp.tv_sec - start ) + micro * ( tp.tv_usec - startu );
    return wall_time;
}


//EOF
