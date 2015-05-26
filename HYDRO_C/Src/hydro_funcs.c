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
#include <mpi.h>
#include <assert.h>

#include "utils.h"
#include "hydro_funcs.h"


void mpi_init ( hydroparam_t* H, int* argc, char*** argv ) {


    int status;

    H->mpi_is_init = 0;

    status = MPI_Init ( argc, argv );
    MPI_Comm_size ( MPI_COMM_WORLD, &H->n_procs );
    MPI_Comm_rank ( MPI_COMM_WORLD, &H->rank );

    if ( status != 0 ) {
        printf ( "MPI_Init: Error %i\n", status );
        exit ( 1 );
    } else {
        dbg_print ( "MPI_Init: successful\n" );
        dbg_print ( "          found %i n_procs\n", H->n_procs );
        dbg_sprint ( "process %03i reporting: ok\n", H->rank );
    }

    H->mpi_is_init = 1;


    // debug: test if mpi is really working:

    if ( DEBUG && H->n_procs > 1) {
        dbg_print ( "testing mpi: simple rend recv\n" );

        int num = 5434; //random int to test sending
        int recv = 0;
        MPI_Status stat;

        if ( H->rank == 0 ) {
            // rank0 sends sone random int

            MPI_Send ( &num, 1, MPI_INTEGER, 1, 0, MPI_COMM_WORLD );
            printf ( "sent: %i\n", num );

        } else if ( H->rank == 1 ) {

            MPI_Recv ( &recv, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &stat );

            printf ( "recv: %i\n", recv );
            assert ( num == recv );
        }
        MPI_Barrier ( MPI_COMM_WORLD );
        dbg_print ( "test complete\n" );
    }

}


void mpi_finish ( hydroparam_t* H ) {

    MPI_Finalize ( );

}


/*
 * Does the domain decompostition
 */
void domain_decomp ( hydroparam_t* H ) {

    dbg_print ( "Domain decomp\n" );

    float frac;
    int lo, hi;

    H->ny = H->ny_domain;
    frac = 1.0 * H->nx_domain / H->n_procs;
    lo = ( int ) ( frac * H->rank );
    hi = ( int ) ( frac * ( H->rank+1 ) );
    H->nx = hi-lo;

    // adjust the domain border properties:
    // remember:    0: cyclic border
    //              1: solid border (mirror)
    //              2: soft border (sink / damped)

    // default for all domains in between
    H->boundary_left = 2;
    H->boundary_right = 2;

    if ( H->rank == 0 ) { // left most domain
        H->boundary_left = 1;
    }
    if ( H->rank+1 == H->n_procs ) { // right most domain (dont use else if in case we have only 1 proc)
        H->boundary_right = 1;
    }

    dbg_print ( "whole domain: %i x %i\n", H->nx_domain, H->ny_domain );
    dbg_sprint ( "rank %03i: nx:%06i ny:%06i hi:%06i lo:%06i frac:%6.4f\n", H->rank, H->nx, H->ny, hi, lo, frac );
}



void hydro_init ( hydroparam_t * H, hydrovar_t * Hv ) {

    long i, j;
    long x, y;

    // *WARNING* : we will use 0 based arrays everywhere since it is C code!
    H->imin = H->jmin = 0;

    // We add two extra layers left/right/top/bottom
    H->imax = H->nx + ExtraLayerTot;
    H->jmax = H->ny + ExtraLayerTot;
    H->nxt = H->imax - H->imin; // column size in the array
    H->nyt = H->jmax - H->jmin; // row size in the array
    // maximum direction size
    H->nxyt = ( H->nxt > H->nyt ) ? H->nxt : H->nyt;

    H->arSz = ( H->nxyt + 2 );
    H->arVarSz = ( H->nxyt + 2 ) * H->nvar;

    // allocate uold for each conservative variable
    Hv->uold = ( double * ) calloc ( H->nvar * H->nxt * H->nyt, sizeof ( double ) );

    // wind tunnel with point explosion
    for ( j = H->jmin + ExtraLayer; j < H->jmax - ExtraLayer; j++ ) {
        for ( i = H->imin + ExtraLayer; i < H->imax - ExtraLayer; i++ ) {
            Hv->uold[IHvP ( i, j, ID )] = one;
            Hv->uold[IHvP ( i, j, IU )] = zero;
            Hv->uold[IHvP ( i, j, IV )] = zero;
            Hv->uold[IHvP ( i, j, IP )] = 1e-5;
        }
    }
    // point explosion at middle of the domian
    /*    x = (H->imax - H->imin) / 2 + ExtraLayer * 0;
    y = (H->jmax - H->jmin) / 2 + ExtraLayer * 0;

     printf("PFL %d %d\n", x, y);
     Hv->uold[IHvP(x, y, IP)] = one / H->dx / H->dx;*/
    // point explosion at corner (top,left)
    // Hv->uold[IHvP ( H->imin+ExtraLayer, H->jmin+ExtraLayer, IP )] = one / H->dx / H->dx;

    // point explosion at the most left domain, top/left corner
    if ( H->rank == 0 ) {
        Hv->uold[IHvP ( H->imin+ExtraLayer, H->jmin+ExtraLayer, IP )] = one / H->dx / H->dx;
    }
}                               // hydro_init



void hydro_finish ( const hydroparam_t H, hydrovar_t * Hv ) {

    Free ( Hv->uold );

}                               // hydro_finish



void allocate_work_space ( const hydroparam_t H, hydrowork_t * Hw, hydrovarwork_t * Hvw ) {

    WHERE ( "allocate_work_space" );

    Hvw->u = DMalloc ( H.arVarSz );
    Hvw->q = DMalloc ( H.arVarSz );
    Hvw->dq = DMalloc ( H.arVarSz );
    Hvw->qxm = DMalloc ( H.arVarSz );
    Hvw->qxp = DMalloc ( H.arVarSz );
    Hvw->qleft = DMalloc ( H.arVarSz );
    Hvw->qright = DMalloc ( H.arVarSz );
    Hvw->qgdnv = DMalloc ( H.arVarSz );
    Hvw->flux = DMalloc ( H.arVarSz );
    Hw->e = DMalloc ( H.arSz );
    Hw->c = DMalloc ( H.arSz );
    Hw->rl = DMalloc ( H.arSz );
    Hw->ul = DMalloc ( H.arSz );
    Hw->pl = DMalloc ( H.arSz );
    Hw->cl = DMalloc ( H.arSz );
    Hw->rr = DMalloc ( H.arSz );
    Hw->ur = DMalloc ( H.arSz );
    Hw->pr = DMalloc ( H.arSz );
    Hw->cr = DMalloc ( H.arSz );
    Hw->ro = DMalloc ( H.arSz );
    Hw->uo = DMalloc ( H.arSz );
    Hw->po = DMalloc ( H.arSz );
    Hw->co = DMalloc ( H.arSz );
    Hw->rstar = DMalloc ( H.arSz );
    Hw->ustar = DMalloc ( H.arSz );
    Hw->pstar = DMalloc ( H.arSz );
    Hw->cstar = DMalloc ( H.arSz );
    Hw->wl = DMalloc ( H.arSz );
    Hw->wr = DMalloc ( H.arSz );
    Hw->wo = DMalloc ( ( H.arSz ) );
    Hw->sgnm = IMalloc ( H.arSz );
    Hw->spin = DMalloc ( H.arSz );
    Hw->spout = DMalloc ( H.arSz );
    Hw->ushock = DMalloc ( H.arSz );
    Hw->frac = DMalloc ( H.arSz );
    Hw->scr = DMalloc ( H.arSz );
    Hw->delp = DMalloc ( H.arSz );
    Hw->pold = DMalloc ( H.arSz );
    Hw->ind = IMalloc ( H.arSz );
    Hw->ind2 = IMalloc ( H.arSz );
} // allocate_work_space


/*
static void
VFree(double **v, const hydroparam_t H)
{
    long i;
    for (i = 0; i < H.nvar; i++) {
        Free(v[i]);
    }
    Free(v);
} // VFree
*/


void deallocate_work_space ( const hydroparam_t H, hydrowork_t * Hw, hydrovarwork_t * Hvw ) {

    WHERE ( "deallocate_work_space" );

    //
    Free ( Hw->e );
    //
    Free ( Hvw->u );
    Free ( Hvw->q );
    Free ( Hvw->dq );
    Free ( Hvw->qxm );
    Free ( Hvw->qxp );
    Free ( Hvw->qleft );
    Free ( Hvw->qright );
    Free ( Hvw->qgdnv );
    Free ( Hvw->flux );

    //
    Free ( Hw->c );
    Free ( Hw->rl );
    Free ( Hw->ul );
    Free ( Hw->pl );
    Free ( Hw->cl );
    Free ( Hw->rr );
    Free ( Hw->ur );
    Free ( Hw->pr );
    Free ( Hw->cr );
    Free ( Hw->ro );
    Free ( Hw->uo );
    Free ( Hw->po );
    Free ( Hw->co );
    Free ( Hw->rstar );
    Free ( Hw->ustar );
    Free ( Hw->pstar );
    Free ( Hw->cstar );
    Free ( Hw->wl );
    Free ( Hw->wr );
    Free ( Hw->wo );
    Free ( Hw->sgnm );
    Free ( Hw->spin );
    Free ( Hw->spout );
    Free ( Hw->ushock );
    Free ( Hw->frac );
    Free ( Hw->scr );
    Free ( Hw->delp );
    Free ( Hw->pold );
    Free ( Hw->ind );
    Free ( Hw->ind2 );

} // deallocate_work_space



