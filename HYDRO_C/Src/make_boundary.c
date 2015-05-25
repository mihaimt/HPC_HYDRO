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

#include "parametres.h"
#include "make_boundary.h"
#include "utils.h"
void
make_boundary(long idim, const hydroparam_t H, hydrovar_t * Hv)
{

    // - - - - - - - - - - - - - - - - - - -
    // Cette portion de code est à vérifier
    // détail. J'ai des doutes sur la conversion
    // des index depuis fortran.
    // - - - - - - - - - - - - - - - - - - -
    long i, ivar, i0, j, j0;
    double sign;
    WHERE("make_boundary");
    
    //dbg_Print("make boundary\n");

 
    if (idim == 1) {

        // Left boundary
        //dbg_sPrint("rk %03i: left boundary\n", H.rank);

        if ( H.rank == 0 ) { // this is the left most domain
            printf("rk %03i: most left\n", H.rank);
            for (ivar = 0; ivar < H.nvar; ivar++) {
                for (i = 0; i < ExtraLayer; i++) {
                    sign = 1.0;
                    if (H.boundary_left == 1) {
                        i0 = ExtraLayerTot - i - 1;
                        if (ivar == IU) {
                            sign = -1.0;
                        }
                    } else if (H.boundary_left == 2) {
                        i0 = 2;
                    } else {
                        i0 = H.nx + i;
                    }
                    for (j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++) {
                        Hv->uold[IHv(i, j, ivar)] = Hv->uold[IHv(i0, j, ivar)] * sign;
                        MFLOPS(1, 0, 0, 0);
                    }
                }
            }
        }
        else {
            printf("rk %03i: not most left\n", H.rank);
            // get the values from left neighboring domain
            // send values to left
            MPI_Status stat;
            MPI_Request req;
            double buf = 600.0;
            
            for (ivar = 0; ivar < H.nvar; ivar++) {
                for (i = 0; i < ExtraLayer; i++) {
                    for (j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++) {
                        
                        
                        MPI_Recv( &buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
                        Hv->uold[IHv(i, j, ivar)] = buf;
                        
                        //MPI_Send( &Hv->uold[IHv(i+ExtraLayer, j, ivar)], 1, MPI_DOUBLE, H.rank-1, 0, MPI_COMM_WORLD );

                        //Hv->uold[IHv(i, j, ivar)] = Hv->uold[IHv(i0, j, ivar)] * sign;
                    }
                }
            }
        }

	/* fprintf(stderr,"PFL H.nvar %d H.nx %d\n",H.nvar,H.nx);
	fprintf(stderr,"PFL ExtraLayer %d ExtraLayerTot %d\n",ExtraLayer,ExtraLayerTot);
	fprintf(stderr,"PFL H.jmin %d H.jmax %d\n",H.jmin,H.jmax); */

        // Right boundary
        //dbg_sPrint("rk %03i: right boundary\n", H.rank);
        
        if ( H.rank == H.n_procs-1 ) { // this is the right most domain
            printf("rk %03i: right most\n", H.rank);

            for (ivar = 0; ivar < H.nvar; ivar++) {
                for (i = H.nx + ExtraLayer; i < H.nx + ExtraLayerTot; i++) {
                    sign = 1.0;
                    if (H.boundary_right == 1) {
                        i0 = 2 * H.nx + ExtraLayerTot - i - 1;
                        if (ivar == IU) {
                            sign = -1.0;
                        }
                    } else if (H.boundary_right == 2) {
                        i0 = H.nx + ExtraLayer;
                    } else {
                        i0 = i - H.nx;
                    }
                    for (j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++) {
                        /* fprintf(stderr,"PFL %d %d\n",i,j); */ 
                        Hv->uold[IHv(i, j, ivar)] = Hv->uold[IHv(i0, j, ivar)] * sign;
                        /* fprintf(stderr,"PFL \n"); */

                        MFLOPS(1, 0, 0, 0);
                    }
                }
            }
        }
        else { // get the values from right neigbouring domain
            printf("rk %03i: not right most\n", H.rank);

            MPI_Status stat;
            MPI_Request req;
            double buf = 0.0;

            for (ivar = 0; ivar < H.nvar; ivar++) {
                for (i = H.nx + ExtraLayer; i < H.nx + ExtraLayerTot; i++) {
                    for (j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++) {

                        MPI_Send( &Hv->uold[IHv(i-ExtraLayer, j, ivar)], 1, MPI_DOUBLE, H.rank+1, 0, MPI_COMM_WORLD);

                        //MPI_Recv( &buf, 1, MPI_DOUBLE, H.rank + 1, 0, MPI_COMM_WORLD, &stat);
                        Hv->uold[IHv(i, j, ivar)] = buf;
                        
                        
                    }
                }
            }
        }
            

        
        
        
        
    } else {

        // Lower boundary
        j0 = 0;
        for (ivar = 0; ivar < H.nvar; ivar++) {
            for (j = 0; j < ExtraLayer; j++) {
                sign = 1.0;
                if (H.boundary_down == 1) {
                    j0 = ExtraLayerTot - j - 1;
                    if (ivar == IV) {
                        sign = -1.0;
                    }
                } else if (H.boundary_down == 2) {
                    j0 = ExtraLayerTot;
                } else {
                    j0 = H.ny + j;
                }
                for (i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++) {
                    Hv->uold[IHv(i, j, ivar)] = Hv->uold[IHv(i, j0, ivar)] * sign;
                    MFLOPS(1, 0, 0, 0);
                }
            }
        }

        // Upper boundary
        for (ivar = 0; ivar < H.nvar; ivar++) {
            for (j = H.ny + ExtraLayer; j < H.ny + ExtraLayerTot; j++) {
                sign = 1.0;
                if (H.boundary_up == 1) {
                    j0 = 2 * H.ny + ExtraLayerTot - j - 1;
                    if (ivar == IV) {
                        sign = -1.0;
                    }
                } else if (H.boundary_up == 2) {
                    j0 = H.ny + 1;
                } else {
                    j0 = j - H.ny;
                }
                for (i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++) {
                    Hv->uold[IHv(i, j, ivar)] = Hv->uold[IHv(i, j0, ivar)] * sign;
                    MFLOPS(1, 0, 0, 0);
                }
            }
        }
    }
}                               // make_boundary


//EOF
