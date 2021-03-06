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
#include <omp.h>
#include <assert.h>

#include "debug.h"

#include "utils.h"
#include "hydro_funcs.h"


/**
 * @brief init mpi and according variables
 * 
 * MPI_init() initializes the MPI library and allocates memory for the
 * MPI variables. We need to do this before anything else so that MPI
 * is initialized and ready to use.
 * 
 * Don't do much if no mpi is used, but init variables never the less
 * and make sure to set n_procs and rank correctly!
 * 
 * Expects: -none-
 * Sets:    mpi_is_init, n_procs, rank
 * 
 * @param H ...
 * @param argc ...
 * @param argv ...
 * @return void
 */
void MPI_init ( hydroparam_t * H, int * argc, char *** argv ) {

    LOC ( H->rank );

    H->mpi_is_init = 0;

/*
    // Allocate Status (use one for all, overwrite old ones)
    H->mpi_status = malloc ( sizeof ( MPI_Status ) );

    // Allocate Requests (need 4 because of the 4 parallel MPI calls:
    // 1 doublerow on each side x 2 sides (left & right) x 2 (send & recv)
    H->mpi_req = malloc ( 4 * sizeof ( MPI_Request ) );
*/

    if ( USE_MPI ) {

        // Initialize MPI library
        H->mpi_error = MPI_Init ( argc, argv );

        // Get the props of the MPI world
        MPI_Comm_size ( MPI_COMM_WORLD, &H->n_procs );
        MPI_Comm_rank ( MPI_COMM_WORLD, &H->rank );

        if ( H->mpi_error != 0 ) {
            ERR ( "MPI_Init: Error %i\n", H->mpi_error );
            exit ( 1 );
        }
        H->mpi_is_init = 1;

        TRC ( H->rank, "Yep, I'm alive" );
        INF_if ( H->rank==0, "mpi init successful, using %i procs\n", H->n_procs );
    }
    else {
        // the rest of the code relies on those, so make sure they are initialized
        // to sensible values!

        H->n_procs = 1;
        H->rank = 0;
        H->mpi_is_init = 1; // yes, that doesn't make much sense.. but anyways..

        INF_if ( H->rank==0, "mpi not running, using %i procs\n", H->n_procs );
    }


}



/**
 * @brief Shutdown MPI
 * 
 * @param H ...
 * @return void
 */
void MPI_finish ( hydroparam_t *H ) {
    //TODO (CR) Dont we need a hydroparam_t *H rather than a const
    // hydroparam_t H here??

    LOC ( H->rank );

/*
    Free ( H->mpi_status );
    Free ( H->mpi_req );
*/

    if ( USE_MPI ) {

        // Free MPI data type
        MPI_Type_free ( &H->mpi_hydro_vector_type );

        H->mpi_error = MPI_Finalize();

        if ( H->mpi_error != 0 ) {
            ERR ( "MPI_Finalize: Error %i\n",H->mpi_error );
            exit ( 1 );
        }
    }

    H->mpi_is_init = 0;

}



/**
 * @brief Initializes OPENMP
 * 
 * @param H ...
 * @return void
 */
void OPENMP_init ( hydroparam_t* H ) {

    if ( USE_OPENMP ) {
        omp_set_dynamic ( 0 ); // dont adjust number of threads, always use the same
        omp_set_num_threads ( N_OMP_THREADS );

        // test openmp, abort of for some reason can't spawn enough threads
        int tid, n_threads;
// //pragma omp parallel private(tid, n_threads)
        {
            tid = omp_get_thread_num();
            H->n_threads = omp_get_num_threads();
            TRC ( H->rank, "rank %i thread %i of %i: Alive, but only for short time", H->rank, tid+1, H->n_threads );
            assert ( H->n_threads == N_OMP_THREADS );
        }

        INF_if ( H->rank==0, "OpenMP init successful, using %i threads\n", H->n_threads );

    }
    else {
        H->n_threads = 1;
        INF_if ( H->rank==0, "OpenMP not running\n");
    }

}


/**
 * @brief Shutdown OPENMP
 * 
 * This is probably empty
 * 
 * @param H ...
 * @return void
 */
void OPENMP_finish ( hydroparam_t* H ) {

}






/**
 * @brief Simple domain decomposition
 * 
 * We use simple vertical slides as sub-domains, assuming a wind tunnel
 * domain. Each process gets one slide of the whole domain.
 * Use some slight integer magic here.
 * (remember: casting to int rounds DOWN)
 * 
 * Expects: n_proc, rank, nxdomain, nydomain
 * Sets   : nx, ny
 * 
 * @param H ...
 * @return void
 */
void MPI_domain_decomp ( hydroparam_t *H ) {

    LOC ( H->rank );

    float frac; // size of one sub domain
    int lo, up; // lower and upper boundary of this proc in the whole domain

    // X DIRECTION: Slicing
    frac = 1.0 * H->nxdomain / H->n_procs;
    lo = ( int ) ( frac * H->rank );
    up = ( int ) ( frac * ( H->rank+1 ) );

    H->nx = up - lo;

    // Y DIRECTION: Dont change ny (only split along the x direction)
    H->ny = H->nydomain;

    // (CR,RK) debug test if domain decomposition makes sense
    if ( DEBUG && USE_MPI ) {
        int nxtotal = 0.0;
        MPI_Allreduce ( &H->nx, &nxtotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
        TRC ( H->rank, "nx=%i ny=%i nxdomain=%i nydomain=%i nxtotal=%i", H->nx, H->ny, H->nxdomain, H->nydomain, nxtotal);
        assert ( nxtotal == H->nxdomain );
    }

} // MPI_domain_decomp




//TODO delete this func
void
hydro_init ( hydroparam_t * H, hydrovar_t * Hv ) {
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
    Hv->uold[IHvP ( H->imin+ExtraLayer, H->jmin+ExtraLayer, IP )] = one / H->dx / H->dx;


}                               // hydro_init






/**
 * @brief Init the simulation
 * 
 * MPI_hydro_init() is basically the same function as hydro_init() but
 * it only allocates the local computational domain and sets the initial
 * conditions depending on the rank of the current MPI process.
 * 
 * Note: the rank is also set (==0) if MPI is disabled!
 * 
 * Expects: mpi_is_init, nx, ny, nvar
 * Sets   : imin, imax, nxt, nyt, nxyt, arSz, arVarSz, 
 *          mpi_hydro_vector_type
 *          uold
 * 
 * @param H ...
 * @param Hv ...
 * @return void
 */
void MPI_hydro_init ( hydroparam_t * H, hydrovar_t * Hv ) {

    LOC ( H->rank );

    long i, j;
    long x, y;

    // Make sure that MPI is initialized before we use it.
    assert ( H->mpi_is_init );

    // Make sure that we did the domain decomposition.
    assert ( H->nx > 0 && H->ny > 0 );

    // WARNING: we will use 0 based arrays everywhere since it is C code!
    H->imin = H->jmin = 0;

    // We add two extra layers left/right/top/bottom
    H->imax = H->nx + ExtraLayerTot;
    H->jmax = H->ny + ExtraLayerTot;
    H->nxt = H->imax - H->imin; // column size in the array
    H->nyt = H->jmax - H->jmin; // row size in the array
    H->nxyt = ( H->nxt > H->nyt ) ? H->nxt : H->nyt; // maximum direction size

    H->arSz = ( H->nxyt + 2 );
    H->arVarSz = ( H->nxyt + 2 ) * H->nvar;


    // Define a new MPI data type vector for transmitting ghost cells
    if ( USE_MPI ) {
        // Just define one column at the moment
        MPI_Type_vector ( H->nvar*H->nyt, ExtraLayer, H->nxt, MPI_DOUBLE,
                          &H->mpi_hydro_vector_type );
        MPI_Type_commit ( &H->mpi_hydro_vector_type );
    }


    // allocate uold for each conservative variable
    Hv->uold = ( double* ) calloc ( H->nvar * H->nxt * H->nyt, sizeof ( double ) );

    // wind tunnel initial condition: silence
    for ( j = H->jmin + ExtraLayer; j < H->jmax - ExtraLayer; j++ ) {
        for ( i = H->imin + ExtraLayer; i < H->imax - ExtraLayer; i++ ) {
            Hv->uold[IHvP ( i, j, ID )] = one;
            Hv->uold[IHvP ( i, j, IU )] = zero;
            Hv->uold[IHvP ( i, j, IV )] = zero;
            Hv->uold[IHvP ( i, j, IP )] = 1e-5;
        }
    }

    // point explosion at corner (bottom, left)
    if ( H->rank == 0 ) {
        Hv->uold[IHvP ( H->imin+ExtraLayer, H->jmin+ExtraLayer, IP )] = one / H->dx / H->dx;
    }

} // MPI_hydro_init





//TODO delete this func
void hydro_finish ( const hydroparam_t H, hydrovar_t * Hv ) {

    Free ( Hv->uold );

} // hydro_finish




/**
 * @brief Cleans up the hydro vars
 * 
 * @param H ...
 * @param Hv ...
 * @return void
 */
void MPI_hydro_finish ( hydroparam_t *H, hydrovar_t * Hv ) {

    LOC ( H->rank );

    Free ( Hv->uold );

} // MPI_hydro_finish







/**
 * @brief Allocates all needed vars for the actual work
 * 
 * @param H ...
 * @param Hw ...
 * @param Hvw ...
 * @return void
 */
void allocate_work_space ( const hydroparam_t H, hydrowork_t * Hw, hydrovarwork_t * Hvw ) {

    LOC ( H.rank );

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




/**
 * @brief Free work variables
 * 
 * @param H ...
 * @param Hw ...
 * @param Hvw ...
 * @return void
 */
void deallocate_work_space ( const hydroparam_t H, hydrowork_t * Hw, hydrovarwork_t * Hvw ) {

    LOC ( H.rank );
    
    //TRC ( H.rank, "thread %i", omp_get_thread_num() );

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














