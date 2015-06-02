/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "parametres.h"


static void usage ( void ) {

    fprintf ( stderr, "options possibles du programme hydro" );
    fprintf ( stderr, "--help" );
    fprintf ( stderr, "-i input" );
    fprintf ( stderr, "-v :: pour avoir les impressions internes" );
    fprintf ( stderr, "------------------------------------" );
    exit ( 1 );

}


static void default_values ( hydroparam_t * H ) {

    fprintf ( stderr,"Setting default values ...\n" );

    // Default values should be given
    H->prt = 0;                 // no printing of internal arrays

    // Global
    H->nxdomain = 2;
    H->nydomain = 2;
    
    // local (is set in MPI_domain_decomp)
    H->nx = 0;
    H->ny = 0;
    H->nvar = IP + 1;
    H->dx = 1.0;
    H->t = 0.0;
    H->nstep = 0;
    H->tend = 0.0;
    H->gamma = 1.4;
    H->courant_factor = one / two;
    H->smallc = 1e-10;
    H->smallr = 1e-10;
    H->niter_riemann = 10;
    H->iorder = 2;
    H->slope_type = 1.;

    // strcpy(H->scheme, "muscl");
    H->scheme = HSCHEME_MUSCL;
    H->boundary_right = 1;
    H->boundary_left = 1;
    H->boundary_up = 1;
    H->boundary_down = 1;
    H->noutput = 1000000;
    H->nstepmax = 1000000;
    H->dtoutput = 0.0;
}


static void keyval ( char *buffer, char **pkey, char **pval ) {

    char *ptr;
    *pkey = buffer;
    *pval = buffer;

    // kill the newline
    *pval = strchr ( buffer, '\n' );
    if ( *pval ) {
        **pval = 0;
    }

    // suppress lead whites or tabs
    while ( ( **pkey == ' ' ) || ( **pkey == '\t' ) ) {
        ( *pkey ) ++;
    }
    *pval = strchr ( buffer, '=' );
    if ( *pval ) {
        **pval = 0;
        ( *pval ) ++;
    }

    // strip key from white or tab
    while ( ( ptr = strchr ( *pkey, ' ' ) ) != NULL ) {
        *ptr = 0;
    }
    while ( ( ptr = strchr ( *pkey, '\t' ) ) != NULL ) {
        *ptr = 0;
    }
}


static void process_input ( char *datafile, hydroparam_t * H ) {

    FILE *fd = NULL;
    char buffer[1024];
    char *pval, *pkey;
    fd = fopen ( datafile, "r" );
    
    if ( fd == NULL ) {
        fprintf ( stderr, "Provided input file (-i) is not readable\n" );
        exit ( 1 );
    }

    while ( fgets ( buffer, 1024, fd ) == buffer ) {
        keyval ( buffer, &pkey, &pval );

        // int parameters
        if ( strcmp ( pkey, "nstepmax" ) == 0 ) {
            sscanf ( pval, "%ld", &H->nstepmax );
            continue;
        }
        if ( strcmp ( pkey, "nx" ) == 0 ) {
            sscanf ( pval, "%ld", &H->nxdomain );
            continue;
        }
        if ( strcmp ( pkey, "ny" ) == 0 ) {
            sscanf ( pval, "%ld", &H->nydomain );
            continue;
        }
        if ( strcmp ( pkey, "boundary_left" ) == 0 ) {
            sscanf ( pval, "%ld", &H->boundary_left );
            continue;
        }
        if ( strcmp ( pkey, "boundary_right" ) == 0 ) {
            sscanf ( pval, "%ld", &H->boundary_right );
            continue;
        }
        if ( strcmp ( pkey, "boundary_up" ) == 0 ) {
            sscanf ( pval, "%ld", &H->boundary_up );
            continue;
        }
        if ( strcmp ( pkey, "boundary_down" ) == 0 ) {
            sscanf ( pval, "%ld", &H->boundary_down );
            continue;
        }
        if ( strcmp ( pkey, "niter_riemann" ) == 0 ) {
            sscanf ( pval, "%ld", &H->niter_riemann );
            continue;
        }
        if ( strcmp ( pkey, "noutput" ) == 0 ) {
            sscanf ( pval, "%ld", &H->noutput );
            continue;
        }
        if ( strcmp ( pkey, "iorder" ) == 0 ) {
            sscanf ( pval, "%ld", &H->iorder );
            continue;
        }
        // float parameters
        if ( strcmp ( pkey, "slope_type" ) == 0 ) {
            sscanf ( pval, "%lf", &H->slope_type );
            continue;
        }
        if ( strcmp ( pkey, "tend" ) == 0 ) {
            sscanf ( pval, "%lf", &H->tend );
            continue;
        }
        if ( strcmp ( pkey, "dx" ) == 0 ) {
            sscanf ( pval, "%lf", &H->dx );
            continue;
        }
        if ( strcmp ( pkey, "courant_factor" ) == 0 ) {
            sscanf ( pval, "%lf", &H->courant_factor );
            continue;
        }
        if ( strcmp ( pkey, "smallr" ) == 0 ) {
            sscanf ( pval, "%lf", &H->smallr );
            continue;
        }
        if ( strcmp ( pkey, "smallc" ) == 0 ) {
            sscanf ( pval, "%lf", &H->smallc );
            continue;
        }
        if ( strcmp ( pkey, "dtoutput" ) == 0 ) {
            sscanf ( pval, "%lf", &H->dtoutput );
            continue;
        }
        // string parameter
        if ( strcmp ( pkey, "scheme" ) == 0 ) {
            if ( strcmp ( pval, "muscl" ) == 0 ) {
                H->scheme = HSCHEME_MUSCL;
            } else if ( strcmp ( pval, "plmde" ) == 0 ) {
                H->scheme = HSCHEME_PLMDE;
            } else if ( strcmp ( pval, "collela" ) == 0 ) {
                H->scheme = HSCHEME_COLLELA;
            } else {
                fprintf ( stderr,
                          "Nom de schema <%s> inconnu, devrait etre l'un de [muscl,plmde,collela]\n",
                          pval );
                exit ( 1 );
            }
            continue;
        }
    }
    fclose ( fd );
    // petit resume de la situation

    printf ( "+-------------------+\n" );
    printf ( "|nx=%-7ld         |\n", H->nx );
    printf ( "|ny=%-7ld         |\n", H->ny );
    printf ( "|tend=%-10.3f    |\n", H->tend );
    printf ( "|nstepmax=%-7ld   |\n", H->nstepmax );
    printf ( "|noutput=%-7ld    |\n", H->noutput );
    printf ( "|dtoutput=%-10.3f|\n", H->dtoutput );
    printf ( "+-------------------+\n" );

    //exit(0);
}


/**
 * @brief processes the commandline arguments and the input file
 * 
 * @param argc ...
 * @param argv ...
 * @param H ...
 * @return void
 */
void process_args ( long argc, char **argv, hydroparam_t * H ) {

    long n = 1;
    char donnees[512];
    default_values ( H );
    while ( n < argc ) {
        if ( strcmp ( argv[n], "--help" ) == 0 ) {
            usage();
            n++;
            continue;
        }
        if ( strcmp ( argv[n], "-v" ) == 0 ) {
            n++;
            H->prt++;
            continue;
        }
        if ( strcmp ( argv[n], "-i" ) == 0 ) {
            n++;
            strncpy ( donnees, argv[n], 512 );
            donnees[511] = 0;   // security
            n++;
            continue;
        }
        fprintf ( stderr, "Argument %s unknown\n", argv[n] );
        n++;
    }
    if ( donnees != NULL ) {
        process_input ( donnees, H );
    } else {
        fprintf ( stderr, "Option -i is missing (it's mandatory!)\n" );
        exit ( 1 );
    }
}








