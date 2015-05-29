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

void
MPI_init(hydroparam_t * H, int * argc, char *** argv)
{
	/*
	** MPI_init() initializes the MPI library and allocates memory for the
	** MPI variables. We need to do this before anything else so that MPI
	** is initialized and ready to use.
	*/
	H->bInit = 0; 

   	/* We need to allocate memory for this variable! */
	H->MPIStatus = malloc(sizeof(MPI_Status));

	// Allocate H->MPI_req !!!
	H->MPI_req = malloc(8*sizeof(MPI_Request));	

	/* Initialize MPI library */
	H->iMPIError = MPI_Init(argc,argv);
	
	MPI_Comm_size(MPI_COMM_WORLD,&H->iNProc);
	MPI_Comm_rank(MPI_COMM_WORLD,&H->iProc);
	
	if (H->iMPIError != 0)
	{
		printf("MPI_Init: Error %i\n",H->iMPIError);
		exit(1);
	}
	H->bInit = 1;
}                               // MPI_init

void
MPI_finish(hydroparam_t *H)
{
	/* (CR) Dont we need a hydroparam_t *H rather than a const hydroparam_t H here?? */
	/* Finalize MPI library */
	H->iMPIError = MPI_Finalize();

	Free(H->MPIStatus);
	// Free MPI_req
	Free(H->MPI_req);

	if (H->iMPIError != 0)
	{
		printf("MPI_Finalize: Error %i\n",H->iMPIError);
		exit(1);
	}

	H->bInit = 0;
}                               // MPI_finish


/*
** Do a domain decomposition. Very simple in our case.
*/
void
MPI_domain_decomp(hydroparam_t *H)
{
	float frac;
	int lo, up;
//	int i,j;

	// Dont change ny (only split along the x direction)
	H->ny = H->nydomain;

	frac = 1.0 * H->nxdomain / H->iNProc;
	lo = (int)(frac * H->iProc);
	up = (int)(frac * (H->iProc+1));

	H->nx = up-lo;
	// (CR) Debug
	fprintf(stderr,"Rank %i: nx=%i ny=%i nxdomain=%i nydomain=%i\n", H->iProc, H->nx, H->ny, H->nxdomain, H->nydomain);
}                               // MPI_domain_decomp			

void
hydro_init(hydroparam_t * H, hydrovar_t * Hv)
{
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
    H->nxyt = (H->nxt > H->nyt) ? H->nxt : H->nyt;

    H->arSz = (H->nxyt + 2);
    H->arVarSz = (H->nxyt + 2) * H->nvar;

    // allocate uold for each conservative variable
    Hv->uold = (double *) calloc(H->nvar * H->nxt * H->nyt, sizeof(double));

    // wind tunnel with point explosion
    for (j = H->jmin + ExtraLayer; j < H->jmax - ExtraLayer; j++) {
        for (i = H->imin + ExtraLayer; i < H->imax - ExtraLayer; i++) {
            Hv->uold[IHvP(i, j, ID)] = one;
            Hv->uold[IHvP(i, j, IU)] = zero;
            Hv->uold[IHvP(i, j, IV)] = zero;
            Hv->uold[IHvP(i, j, IP)] = 1e-5;
        }
    }
    // point explosion at middle of the domian
    /*    x = (H->imax - H->imin) / 2 + ExtraLayer * 0;
    y = (H->jmax - H->jmin) / 2 + ExtraLayer * 0;

     printf("PFL %d %d\n", x, y);
     Hv->uold[IHvP(x, y, IP)] = one / H->dx / H->dx;*/
    // point explosion at corner (top,left)
    Hv->uold[IHvP(H->imin+ExtraLayer, H->jmin+ExtraLayer, IP)] = one / H->dx / H->dx;
}                               // hydro_init

/*
** MPI_hydro_init() is basically the same function as hydro_init() but
** it only allocates the local computational domain and sets the initial
** conditions depending on the rank of the current MPI process.
*/
void
MPI_hydro_init(hydroparam_t * H, hydrovar_t * Hv)
{
    long i, j;
    long x, y;

	/* Make sure that MPI is initialized before we use it. */
	assert(H->bInit);
	
	/* Make sure that we did the domain decomposition. */
	assert( H->nx > 0 && H->ny > 0);

    // *WARNING* : we will use 0 based arrays everywhere since it is C code!
    H->imin = H->jmin = 0;

    // We add two extra layers left/right/top/bottom
    H->imax = H->nx + ExtraLayerTot;
    H->jmax = H->ny + ExtraLayerTot;
    H->nxt = H->imax - H->imin; // column size in the array
    H->nyt = H->jmax - H->jmin; // row size in the array
    // maximum direction size
    H->nxyt = (H->nxt > H->nyt) ? H->nxt : H->nyt;

    H->arSz = (H->nxyt + 2);
    H->arVarSz = (H->nxyt + 2) * H->nvar;

	// Define a new MPI data type
//	MPI_Type_vector( H->nvar*H->nyt*ExtraLayer, ExtraLayer, H->nxt, MPI_DOUBLE, &H->MPI_Hydro_vars );
//	Just define one column at the moment
	MPI_Type_vector( H->nvar*H->nyt, 1, H->nxt, MPI_DOUBLE, &H->MPI_Hydro_vars );
	MPI_Type_commit( & H->MPI_Hydro_vars );
		
    // allocate uold for each conservative variable
    Hv->uold = (double *) calloc(H->nvar * H->nxt * H->nyt, sizeof(double));

    // wind tunnel with point explosion
	for (j = H->jmin + ExtraLayer; j < H->jmax - ExtraLayer; j++) {
		for (i = H->imin + ExtraLayer; i < H->imax - ExtraLayer; i++) {
			Hv->uold[IHvP(i, j, ID)] = one;
            Hv->uold[IHvP(i, j, IU)] = zero;
			Hv->uold[IHvP(i, j, IV)] = zero;
			Hv->uold[IHvP(i, j, IP)] = 1e-5;
        }
    }
    // point explosion at middle of the domian
    /*    x = (H->imax - H->imin) / 2 + ExtraLayer * 0;
    y = (H->jmax - H->jmin) / 2 + ExtraLayer * 0;

     printf("PFL %d %d\n", x, y);
     Hv->uold[IHvP(x, y, IP)] = one / H->dx / H->dx;*/
    // point explosion at corner (top,left)
	// (CR) Bottom left??
	if (H->iProc == 0)
	{
		/* Only in the domain of process zero we have the initial explosion. */
	    Hv->uold[IHvP(H->imin+ExtraLayer, H->jmin+ExtraLayer, IP)] = one / H->dx / H->dx;
	}
}                               // MPI_hydro_init

void
hydro_finish(const hydroparam_t H, hydrovar_t * Hv)
{
    Free(Hv->uold);
}                               // hydro_finish

void
MPI_hydro_finish(hydroparam_t *H, hydrovar_t * Hv)
{
	/* (CR) Dont we need a hydroparam_t *H rather than a const hydroparam_t H here?? */

	// Free MPI data type
	MPI_Type_free( &H->MPI_Hydro_vars );

	/* Finalize MPI library */
	if (H->MPIStatus != NULL) MPI_finish(H);

    Free(Hv->uold);
}
                               // hydro_finish
void
allocate_work_space(const hydroparam_t H, hydrowork_t * Hw, hydrovarwork_t * Hvw)
{
    WHERE("allocate_work_space");
    Hvw->u = DMalloc(H.arVarSz);
    Hvw->q = DMalloc(H.arVarSz);
    Hvw->dq = DMalloc(H.arVarSz);
    Hvw->qxm = DMalloc(H.arVarSz);
    Hvw->qxp = DMalloc(H.arVarSz);
    Hvw->qleft = DMalloc(H.arVarSz);
    Hvw->qright = DMalloc(H.arVarSz);
    Hvw->qgdnv = DMalloc(H.arVarSz);
    Hvw->flux = DMalloc(H.arVarSz);
    Hw->e = DMalloc(H.arSz);
    Hw->c = DMalloc(H.arSz);
    Hw->rl = DMalloc(H.arSz);
    Hw->ul = DMalloc(H.arSz);
    Hw->pl = DMalloc(H.arSz);
    Hw->cl = DMalloc(H.arSz);
    Hw->rr = DMalloc(H.arSz);
    Hw->ur = DMalloc(H.arSz);
    Hw->pr = DMalloc(H.arSz);
    Hw->cr = DMalloc(H.arSz);
    Hw->ro = DMalloc(H.arSz);
    Hw->uo = DMalloc(H.arSz);
    Hw->po = DMalloc(H.arSz);
    Hw->co = DMalloc(H.arSz);
    Hw->rstar = DMalloc(H.arSz);
    Hw->ustar = DMalloc(H.arSz);
    Hw->pstar = DMalloc(H.arSz);
    Hw->cstar = DMalloc(H.arSz);
    Hw->wl = DMalloc(H.arSz);
    Hw->wr = DMalloc(H.arSz);
    Hw->wo = DMalloc((H.arSz));
    Hw->sgnm = IMalloc(H.arSz);
    Hw->spin = DMalloc(H.arSz);
    Hw->spout = DMalloc(H.arSz);
    Hw->ushock = DMalloc(H.arSz);
    Hw->frac = DMalloc(H.arSz);
    Hw->scr = DMalloc(H.arSz);
    Hw->delp = DMalloc(H.arSz);
    Hw->pold = DMalloc(H.arSz);
    Hw->ind = IMalloc(H.arSz);
    Hw->ind2 = IMalloc(H.arSz);
}                               // allocate_work_space


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
void
deallocate_work_space(const hydroparam_t H, hydrowork_t * Hw, hydrovarwork_t * Hvw)
{
    WHERE("deallocate_work_space");

    //
    Free(Hw->e);
    //
    Free(Hvw->u);
    Free(Hvw->q);
    Free(Hvw->dq);
    Free(Hvw->qxm);
    Free(Hvw->qxp);
    Free(Hvw->qleft);
    Free(Hvw->qright);
    Free(Hvw->qgdnv);
    Free(Hvw->flux);

    //
    Free(Hw->c);
    Free(Hw->rl);
    Free(Hw->ul);
    Free(Hw->pl);
    Free(Hw->cl);
    Free(Hw->rr);
    Free(Hw->ur);
    Free(Hw->pr);
    Free(Hw->cr);
    Free(Hw->ro);
    Free(Hw->uo);
    Free(Hw->po);
    Free(Hw->co);
    Free(Hw->rstar);
    Free(Hw->ustar);
    Free(Hw->pstar);
    Free(Hw->cstar);
    Free(Hw->wl);
    Free(Hw->wr);
    Free(Hw->wo);
    Free(Hw->sgnm);
    Free(Hw->spin);
    Free(Hw->spout);
    Free(Hw->ushock);
    Free(Hw->frac);
    Free(Hw->scr);
    Free(Hw->delp);
    Free(Hw->pold);
    Free(Hw->ind);
    Free(Hw->ind2);

}                               // deallocate_work_space


// EOF
