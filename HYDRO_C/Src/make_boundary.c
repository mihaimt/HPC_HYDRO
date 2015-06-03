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
#include <assert.h>

#include "debug.h"

#include "parametres.h"
#include "make_boundary.h"
#include "utils.h"







void
make_boundary ( long idim, const hydroparam_t H, hydrovar_t * Hv ) {

    long i, ivar, i0, j, j0;
    double sign;
    WHERE ( "make_boundary" );


    if ( idim == 1 ) {

        // Left boundary
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( i = 0; i < ExtraLayer; i++ ) {
                sign = 1.0;
                if ( H.boundary_left == 1 ) {
                    i0 = ExtraLayerTot - i - 1;
                    if ( ivar == IU ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_left == 2 ) {
                    i0 = 2;
                } else {
                    i0 = H.nx + i;
                }
                for ( j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++ ) {
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i0, j, ivar )] * sign;
                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }

        /* fprintf(stderr,"PFL H.nvar %d H.nx %d\n",H.nvar,H.nx);
        fprintf(stderr,"PFL ExtraLayer %d ExtraLayerTot %d\n",ExtraLayer,ExtraLayerTot);
        fprintf(stderr,"PFL H.jmin %d H.jmax %d\n",H.jmin,H.jmax); */

        // Right boundary
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( i = H.nx + ExtraLayer; i < H.nx + ExtraLayerTot; i++ ) {
                sign = 1.0;
                if ( H.boundary_right == 1 ) {
                    i0 = 2 * H.nx + ExtraLayerTot - i - 1;
                    if ( ivar == IU ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_right == 2 ) {
                    i0 = H.nx + ExtraLayer;
                } else {
                    i0 = i - H.nx;
                }
                for ( j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++ ) {
                    /* fprintf(stderr,"PFL %d %d\n",i,j); */
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i0, j, ivar )] * sign;
                    /*		  fprintf(stderr,"PFL \n"); */

                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }
    } else {

        // Lower boundary
        j0 = 0;
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( j = 0; j < ExtraLayer; j++ ) {
                sign = 1.0;
                if ( H.boundary_down == 1 ) {
                    j0 = ExtraLayerTot - j - 1;
                    if ( ivar == IV ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_down == 2 ) {
                    j0 = ExtraLayerTot;
                } else {
                    j0 = H.ny + j;
                }
                for ( i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++ ) {
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i, j0, ivar )] * sign;
                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }

        // Upper boundary
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( j = H.ny + ExtraLayer; j < H.ny + ExtraLayerTot; j++ ) {
                sign = 1.0;
                if ( H.boundary_up == 1 ) {
                    j0 = 2 * H.ny + ExtraLayerTot - j - 1;
                    if ( ivar == IV ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_up == 2 ) {
                    j0 = H.ny + 1;
                } else {
                    j0 = j - H.ny;
                }
                for ( i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++ ) {
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i, j0, ivar )] * sign;
                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }
    }
}                               // make_boundary









/**
 * @brief Defines the boundary cells
 * 
 * Exchange the boundary conditions with neighboring domains that are on
 * different processes. We use MPI_IRecv() and MPI_ISend() to exchange the
 * data. This means we have to make sure that all data was transvered before
 * we can do the computation for the most outer layer of the grid using
 * MPI_get_boundary_end().
 * 
 * Note: NO special handling for non MPI mode is required throuout the file!
 * The MPI calls to send/recv will never happen anyways if rank=0 and n_procs=1
 * 
 * @param idim ...
 * @param H ...
 * @param Hv ...
 * @param MPI_req ...
 * @return void
 */
void MPI_get_boundary_start ( long idim, const hydroparam_t H, hydrovar_t* Hv,
                              MPI_Request* MPI_req ) {

    LOC ( H.rank );

    long i, ivar, i0, j, j0, k;
    double sign;

    // In which direction are we scanning?
    if ( idim == 1 ) {
        // scanning horizontally (left / right)


        // Make sure MPI_req is allocated
        assert ( MPI_req != NULL );

        // Set physical boundary conditions for the left ghost layer of this domain

        if ( H.rank == 0 ) {
            // this is the left most domain, so enforce physical boundary conditions here

            for ( ivar = 0; ivar < H.nvar; ivar++ ) {
                for ( i = 0; i < ExtraLayer; i++ ) {
                    sign = 1.0;
                    if ( H.boundary_left == 1 ) {   // 1 is fixed / reflecting border
                        i0 = ExtraLayerTot - i - 1;
                        if ( ivar == IU ) {
                            sign = -1.0;
                        }
                    }
                    else if ( H.boundary_left == 2 ) { // 2 is absorbing border
                        i0 = 2;
                    }
                    else {   // everything else is circular connected domain
                        // TODO This is not implemented! additional MPI calls between first and last domain
                        // would be required..
                        // SO we break if anyone wants to do this and let him know...
                        i0 = H.nx + i;
                        ERR ( "circular boundary conditions are setup\nThis is NOT supported at the moment, sorry!" );
                        exit ( 1 );
                    }
                    // Copy the data into the ghost cells
                    for ( j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++ ) {
                        Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i0, j, ivar )] * sign;
                        //TODO remove this //MFLOPS ( 1, 0, 0, 0 );
                    }
                }
            }
        }
        else {
            // this is NOT the left most domain, but somewhere in the middle or right

#if COM_METHOD == _CM_VECTOR
            // Exchange MPI boundary conditions with proc / domain to the left (send two columns at once)
            // get 2 rows from the actual domain of the left neighbor and use those as ghost cells in the current cell.
            // send 2 rows of physical data to the left, such that it can use it as ghost cells
            MPI_Irecv ( &Hv->uold[IHv ( 0, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 2, MPI_COMM_WORLD, MPI_req );
            MPI_Isend ( &Hv->uold[IHv ( 2, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 0, MPI_COMM_WORLD, MPI_req+1 );
#elif COM_METHOD == _CM_SINGLE
            MPI_Status stat;
            MPI_Request req;
            double buf = 0.0;
            for (ivar = 0; ivar < H.nvar; ivar++) {
                for (i = 0; i < ExtraLayer; i++) {
                    for (j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++) {
                        MPI_Recv( &buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
                        Hv->uold[IHv(i, j, ivar)] = buf;
                        
                        MPI_Send( &Hv->uold[IHv(i+ExtraLayer, j, ivar)], 1, MPI_DOUBLE, H.rank-1, 0, MPI_COMM_WORLD );
                    }
                }
            }
#endif
        }

        // Set physical boundary conditions for the right ghost layer of this domain

        if ( H.rank == H.n_procs-1 ) {
            // this is the right most domain, so enforce physical boundary conditions here

            for ( ivar = 0; ivar < H.nvar; ivar++ ) {  // explanation see above (left hand part)
                for ( i = H.nx + ExtraLayer; i < H.nx + ExtraLayerTot; i++ ) {
                    sign = 1.0;
                    if ( H.boundary_right == 1 ) {
                        i0 = 2 * H.nx + ExtraLayerTot - i - 1;
                        if ( ivar == IU ) {
                            sign = -1.0;
                        }
                    } else if ( H.boundary_right == 2 ) {
                        i0 = H.nx + ExtraLayer;
                    } else {
                        i0 = i - H.nx;
                    }
                    for ( j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++ ) {
                        Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i0, j, ivar )] * sign;
                        //TODO delete this //MFLOPS ( 1, 0, 0, 0 );
                    }
                }
            }
        }
        else {
            // this is NOT the right most domain, but somewhere in the middle or to the left

#if COM_METHOD == _CM_VECTOR
            // Exchange MPI boundary conditions with proc / domain to the right (send two columns at once)
            // get 2 rows from the actual domain of the right neighbor and use those as ghost cells in the current cell.
            // send 2 rows of physical data to the right, such that it can use it as ghost cells
            MPI_Irecv ( &Hv->uold[IHv ( H.nx+ExtraLayer, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 0, MPI_COMM_WORLD, MPI_req+2 );
            MPI_Isend ( &Hv->uold[IHv ( H.nx+ExtraLayer-2, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 2, MPI_COMM_WORLD, MPI_req+3 );
#elif COM_METHOD == _CM_SINGLE
            MPI_Status stat;
            MPI_Request req;
            double buf = 0.0;
            for (ivar = 0; ivar < H.nvar; ivar++) {
                for (i = 0; i < ExtraLayer; i++) {
                    for (j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++) {
                        MPI_Send( &Hv->uold[IHv(i+ExtraLayer, j, ivar)], 1, MPI_DOUBLE, H.rank+1, 0, MPI_COMM_WORLD );

                        MPI_Recv( &buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
                        Hv->uold[IHv(i, j, ivar)] = buf;
                    }
                }
            }
#endif

        }

    }
    else {
        // scanning vertically (top / down)
        // this part is still like the original code

        // Lower boundary
        j0 = 0;
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( j = 0; j < ExtraLayer; j++ ) {
                sign = 1.0;
                if ( H.boundary_down == 1 ) {
                    j0 = ExtraLayerTot - j - 1;
                    if ( ivar == IV ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_down == 2 ) {
                    j0 = ExtraLayerTot;
                } else {
                    j0 = H.ny + j;
                }
                for ( i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++ ) {
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i, j0, ivar )] * sign;
                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }

        // Upper boundary
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( j = H.ny + ExtraLayer; j < H.ny + ExtraLayerTot; j++ ) {
                sign = 1.0;
                if ( H.boundary_up == 1 ) {
                    j0 = 2 * H.ny + ExtraLayerTot - j - 1;
                    if ( ivar == IV ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_up == 2 ) {
                    j0 = H.ny + 1;
                } else {
                    j0 = j - H.ny;
                }
                for ( i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++ ) {
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i, j0, ivar )] * sign;
                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }
    }
} // MPI_get_boundary_start










/**
 * @brief Finish up domain transfers
 * 
 * Make sure that all boundary cells have been successfully exchanged between
 * the processes before we continue.
 * 
 * @param idim ...
 * @param H ...
 * @param Hv ...
 * @param MPI_req ...
 * @return void
 */
void MPI_get_boundary_end ( long idim, const hydroparam_t H, hydrovar_t * Hv, MPI_Request *MPI_req ) {

    LOC ( H.rank );

    int count, offset;
    int error;
    MPI_Status *status;

    status = malloc( 4*sizeof( MPI_Status ));

    offset = 0;

    if ( idim == 1) {
        // Don't do this if we sweep the grid in ny direction.

#if COM_METHOD == _CM_VEKTOR
        if ( USE_MPI && H.n_procs>1 ) {
            if (H.rank == 0 || H.rank == H.n_procs-1) {
                // For the most outer domains we have less b.c. to exchange.
                count = 1;
                if ( H.rank == 0 ) {
                    offset = 2;
                }
            } else {
                count = 2;
            }
            error = MPI_Waitall( 2*count, MPI_req + offset, status);
            assert( error == 0 );
        }
#elif COM_METHOD == _CM_SINGLE
#endif
    }

    
/* 
    // (CR) I have no idea why MPI_Waitall() doesnt work
    // (RK) I did somehow fix it.. now idea how..

    // Here we send two columns at once
    if ( idim == 1 && H.n_procs > 1 ) {
        // Dont do this if we sweep the grid in ny direction.
        if ( H.rank == 0 ) {
            // Only exchanged cells with the right layer
            // (CR) Debug
            fprintf ( stderr,"Rank %i: MPI_Wait 2\n", H.rank );
            MPI_Wait ( MPI_req+2, &status );
            fprintf ( stderr,"Rank %i: MPI_Wait 3\n", H.rank );
            MPI_Wait ( MPI_req+3, &status );
        } else if ( H.rank == H.n_procs-1 ) {
            // Only exchanged cells with the left layer
            // (CR) Debug
            fprintf ( stderr,"Rank %i: MPI_Wait 0\n", H.rank );
            MPI_Wait ( MPI_req, &status );
            fprintf ( stderr,"Rank %i: MPI_Wait 1\n", H.rank );
            MPI_Wait ( MPI_req+1, &status );
        } else {
            MPI_Wait ( MPI_req, &status );
            MPI_Wait ( MPI_req+1, &status );
            MPI_Wait ( MPI_req+2, &status );
            MPI_Wait ( MPI_req+3, &status );
        }
    }
*/

} // MPI_get_boundary_end
















/*
** Here we use MPI_sendrecv rather than more complicated stuff.
*/
void
MPI_get_boundary_sendrecv ( long idim, const hydroparam_t H, hydrovar_t * Hv ) {
    long i, ivar, i0, j, j0, k;
    double sign;
    MPI_Status status;
    WHERE ( "MPI_get_boundary" );

    // Use:
    // MPI_Sendrecv( sendbuff, sendcount, sendtype, dest, sendtag,
    //				 recvbuff, recvcount, recvtype, source, recvtag, comm, status);
    // MPI_Send( buf, count, datatype, dest, tag, comm);
    // MPI_Recv( buf, count, datatype, source, tag, comm, status);
    if ( idim == 1 ) {
        // (CR) Debug
        fprintf ( stderr,"Rank %i: Send copy layer to the right and receive ghost layer from the left\n", H.rank );

        //
        // (CR) Print debug information.
        //

        /*		for (ivar = 0; ivar < H.nvar; ivar++) {
        			for (j = 0; j < H.nyt; j++) {
        				for (i = 0; i < H.nxt; i++) {
        					fprintf(stdout, "%10.3e ", Hv->uold[IHv(i, j, nvar)]);
                       	}
                    }
                }*/

        hydrovar_t *Hold;
        Hold->uold = ( double * ) calloc ( H.nvar * H.nxt * H.nyt, sizeof ( double ) );

        if ( H.rank == 1 ) {
            i = H.nx+ExtraLayer;
            fprintf ( stdout, "Rank %i: i = %i\n", H.rank, i );
            for ( j = 0; j < H.nyt; j++ ) {
                for ( i = 0; i < H.nxt; i++ ) {
//					fprintf(stdout, "%1.1e ", Hv->uold[IHv(i, j, ID)]);
                    fprintf ( stdout, "%i ", ( int ) Hv->uold[IHv ( i, j, ID )] );
                }
                fprintf ( stdout, "\n" );
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////
        // Send right copy layer to the right domain and receive left ghost layer from the left domain.
        ///////////////////////////////////////////////////////////////////////////////////////////////
        if ( H.rank == 0 ) {
            // Rank 0 (most left domain)
            // Set physical boundary conditions for the left ghost cells
            for ( ivar = 0; ivar < H.nvar; ivar++ ) {
                for ( i = 0; i < ExtraLayer; i++ ) {
                    sign = 1.0;
                    if ( H.boundary_left == 1 ) {
                        i0 = ExtraLayerTot - i - 1;
                        if ( ivar == IU ) {
                            sign = -1.0;
                        }
                    } else if ( H.boundary_left == 2 ) {
                        i0 = 2;
                    } else {
                        i0 = H.nx + i;
                    }
                    for ( j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++ ) {
                        Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i0, j, ivar )] * sign;
                        MFLOPS ( 1, 0, 0, 0 );
                    }
                }
            }
            // Send right copy layer to the right domain (but dont receive the left ghost layer from the left).
            // (CR) Debug
            fprintf ( stderr, "iProc %i: Sending right copy layer to iProc %i\n", H.rank, H.rank+1 );

            MPI_Send ( &Hv->uold[IHv ( H.nx+ExtraLayer-2, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 0, MPI_COMM_WORLD );
            MPI_Send ( &Hv->uold[IHv ( H.nx+ExtraLayer-1, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 1, MPI_COMM_WORLD );

            fprintf ( stderr, "iProc %i: Sending right copy layer to iProc %i. Done\n", H.rank, H.rank+1 );
        } else if ( H.rank == H.n_procs-1 ) {
            // Rank iNProc-1 (most right domain)
            // Receive left ghost layer from the left domain (but dont send the right copy layer to the right).
            // (CR) Debug
            fprintf ( stderr, "iProc %i: Receiving left ghost layer from iProc %i\n", H.rank, H.rank-1 );

            MPI_Recv ( &Hv->uold[IHv ( 0, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 0, MPI_COMM_WORLD, &status );
            MPI_Recv ( &Hv->uold[IHv ( 1, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 1, MPI_COMM_WORLD, &status );

            fprintf ( stderr, "iProc %i: Receiving left ghost layer from iProc %i. Done\n", H.rank, H.rank-1 );
        } else {
            // Exchange MPI boundary conditions
            // (CR) Debug
            fprintf ( stderr, "iProc %i: Sendrecv() 1\n", H.rank );
            // Send the first right copy layer (nx+ExtraLayer-2) to the right domain (iProc+1) and
            // receive the first left ghost layer (0) from the left domain (iProc-1)
            MPI_Sendrecv ( &Hv->uold[IHv ( H.nx+ExtraLayer-2, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 4,
                           &Hv->uold[IHv ( 0, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 5, MPI_COMM_WORLD, &status );

            fprintf ( stderr, "iProc %i: Sendrecv() 2\n", H.rank );

            // Send the second right copy layer (nx+ExtraLayer-2) to the right domain (iProc+1) and
            // receive the second left ghost layer (0) from the left domain (iProc-1)
            MPI_Sendrecv ( &Hv->uold[IHv ( H.nx+ExtraLayer-1, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 6,
                           &Hv->uold[IHv ( 1, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 7, MPI_COMM_WORLD, &status );
        }

        // (CR) Debug
        fprintf ( stderr,"Rank %i: Synchronizing\n", H.rank );
        MPI_Barrier ( MPI_COMM_WORLD );
        fprintf ( stderr,"Rank %i: Synchronizing. Done\n", H.rank );

        fprintf ( stderr,"Rank %i: Send copy layer to the left and receive ghost layer from the right\n", H.rank );

        ///////////////////////////////////////////////////////////////////////////////////////////////
        // Send left copy layer to the left domain and receive right ghost layer from the right domain.
        ///////////////////////////////////////////////////////////////////////////////////////////////
        if ( H.rank == 0 ) {
            // Rank 0 (most left domain)
            // Receive right ghost cells from the right domain (but dont send the left copy layer to the left).
//			fprintf(stderr, "Sending right copy layer for iProc %i\n", H.iProc);
            // MPI_Recv( buf, count, datatype, source, tag, comm, status)
            MPI_Recv ( &Hv->uold[IHv ( H.nx+ExtraLayer+1, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 0, MPI_COMM_WORLD, &status );
            MPI_Recv ( &Hv->uold[IHv ( H.nx+ExtraLayer, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 1, MPI_COMM_WORLD, &status );
        } else if ( H.rank == H.n_procs - 1 ) {
            // Rank: iNProc-1 (most right domain)
            // Set physical boundary conditions for the right ghost cells
            for ( ivar = 0; ivar < H.nvar; ivar++ ) {
                for ( i = H.nx + ExtraLayer; i < H.nx + ExtraLayerTot; i++ ) {
                    sign = 1.0;
                    if ( H.boundary_right == 1 ) {
                        i0 = 2 * H.nx + ExtraLayerTot - i - 1;
                        if ( ivar == IU ) {
                            sign = -1.0;
                        }
                    } else if ( H.boundary_right == 2 ) {
                        i0 = H.nx + ExtraLayer;
                    } else {
                        i0 = i - H.nx;
                    }
                    for ( j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++ ) {
                        /* fprintf(stderr,"PFL %d %d\n",i,j); */
                        Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i0, j, ivar )] * sign;
                        /* fprintf(stderr,"PFL \n"); */
                        MFLOPS ( 1, 0, 0, 0 );
                    }
                }
            }
            // Send left copy layer to the left domain (but dont receive the right ghost layer from the right).
            MPI_Send ( &Hv->uold[IHv ( 2, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 0, MPI_COMM_WORLD );
            MPI_Send ( &Hv->uold[IHv ( 3, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 1, MPI_COMM_WORLD );
        } else {
            // Exchange MPI boundary conditions

            // Send the first left copy layer (2) and receive the first right ghost layer (nx+ExtraLayer)
            MPI_Sendrecv ( &Hv->uold[IHv ( 2, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 0,
                           &Hv->uold[IHv ( H.nx+ExtraLayer, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 1, MPI_COMM_WORLD, &status );

            // Send the second left copy layer (3) and receive the second right ghost layer (nx+ExtraLayer+1)
            MPI_Sendrecv ( &Hv->uold[IHv ( 3, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 0,
                           &Hv->uold[IHv ( H.nx+ExtraLayer+1, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 1, MPI_COMM_WORLD, &status );
        }
        Free ( Hold );

    } else {
        // Lower boundary
        j0 = 0;
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( j = 0; j < ExtraLayer; j++ ) {
                sign = 1.0;
                if ( H.boundary_down == 1 ) {
                    j0 = ExtraLayerTot - j - 1;
                    if ( ivar == IV ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_down == 2 ) {
                    j0 = ExtraLayerTot;
                } else {
                    j0 = H.ny + j;
                }
                for ( i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++ ) {
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i, j0, ivar )] * sign;
                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }

        // Upper boundary
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( j = H.ny + ExtraLayer; j < H.ny + ExtraLayerTot; j++ ) {
                sign = 1.0;
                if ( H.boundary_up == 1 ) {
                    j0 = 2 * H.ny + ExtraLayerTot - j - 1;
                    if ( ivar == IV ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_up == 2 ) {
                    j0 = H.ny + 1;
                } else {
                    j0 = j - H.ny;
                }
                for ( i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++ ) {
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i, j0, ivar )] * sign;
                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }
    }
}                               // MPI_get_boundary











/*
 * The most simple version possible. Only two processes and I implemented every send and receive
 * manually.
 */
void
MPI_get_boundary_simple ( long idim, const hydroparam_t H, hydrovar_t * Hv ) {
    long i, ivar, i0, j, j0, k;
    double sign;
    MPI_Status status;
    int source, dest;
    WHERE ( "MPI_get_boundary_simple" );

    // (CR) Debug
    assert ( H.n_procs == 2 );

    fprintf ( stdout,"Rank %i: MPI_get_boundary_simple()\n", H.rank );
    // Use:
    // MPI_Sendrecv( sendbuff, sendcount, sendtype, dest, sendtag,
    //				 recvbuff, recvcount, recvtype, source, recvtag, comm, status);
    // MPI_Send( buf, count, datatype, dest, tag, comm);
    // MPI_Recv( buf, count, datatype, source, tag, comm, status);
    if ( idim == 1 ) {

        //
        // (CR) Print debug information.
        //

        /*		for (ivar = 0; ivar < H.nvar; ivar++) {
        			for (j = 0; j < H.nyt; j++) {
        				for (i = 0; i < H.nxt; i++) {
        					fprintf(stdout, "%10.3e ", Hv->uold[IHv(i, j, nvar)]);
                       	}
                    }
                }*/

//		hydrovar_t *Hold;
//		Hold->uold = (double *) calloc(H.nvar * H.nxt * H.nyt, sizeof(double));

        /*		if ( H.iProc == 1 )
        		{
        			i = H.nx+ExtraLayer;
        			fprintf(stdout, "Rank %i: i = %i\n", H.iProc, i);
        			for (j = 0; j < H.nyt; j++) {
        				for (i = 0; i < H.nxt; i++) {
        //					fprintf(stdout, "%1.1e ", Hv->uold[IHv(i, j, ID)]);
        					fprintf(stdout, "%i ", (int) Hv->uold[IHv(i, j, ID)]);
        				}
        				fprintf(stdout, "\n");
        			}
        		}*/

        ///////////////////////////////////////////////////////////////////
        // Make physical b.c. for the left ghost layer of the first domain.
        ///////////////////////////////////////////////////////////////////
        if ( H.rank == 0 ) {
            ///////////////////////////////////////////////////////////////////
            // Make physical b.c. for the left ghost layer of the first domain.
            ///////////////////////////////////////////////////////////////////
            for ( ivar = 0; ivar < H.nvar; ivar++ ) {
                for ( i = 0; i < ExtraLayer; i++ ) {
                    sign = 1.0;
                    if ( H.boundary_left == 1 ) {
                        i0 = ExtraLayerTot - i - 1;
                        if ( ivar == IU ) {
                            sign = -1.0;
                        }
                    } else if ( H.boundary_left == 2 ) {
                        i0 = 2;
                    } else {
                        i0 = H.nx + i;
                    }
                    for ( j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++ ) {
                        Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i0, j, ivar )] * sign;
                        MFLOPS ( 1, 0, 0, 0 );
                    }
                }
            }
            //////////////////////////////////////////////
            // Send right copy layer to the second domain.
            //////////////////////////////////////////////
            fprintf ( stderr, "iProc %i: Sending right copy layer to iProc %i\n", H.rank, H.rank+1 );

            // Send two columns at once
            MPI_Send ( &Hv->uold[IHv ( H.nx+ExtraLayer-2, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 0, MPI_COMM_WORLD );
//			MPI_Send( &Hv->uold[IHv(H.nx+ExtraLayer-2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 0, MPI_COMM_WORLD);
//			MPI_Send( &Hv->uold[IHv(H.nx+ExtraLayer-1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 1, MPI_COMM_WORLD);

            fprintf ( stderr, "iProc %i: Sending right copy layer to iProc %i. Done\n", H.rank, H.rank+1 );
            ////////////////////////////////////////////////////
            // Receive right ghost layer from the second domain.
            ////////////////////////////////////////////////////
            fprintf ( stderr, "iProc %i: Receiving right ghost layer from iProc %i\n", H.rank, H.rank+1 );

            // Send two columns at once

            MPI_Recv ( &Hv->uold[IHv ( H.nx+ExtraLayer, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank+1, 2, MPI_COMM_WORLD, &status );
//			MPI_Recv( &Hv->uold[IHv(H.nx+ExtraLayer, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 2, MPI_COMM_WORLD, &status);
//			MPI_Recv( &Hv->uold[IHv(H.nx+ExtraLayer+1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 3, MPI_COMM_WORLD, &status);

            fprintf ( stderr, "iProc %i: Receiving right ghost layer from iProc %i. Done\n", H.rank, H.rank+1 );
        } else  {
            /////////////////////////////////////////////////////////////////////
            // Make physical b.c. for the right ghost layer of the second domain.
            /////////////////////////////////////////////////////////////////////
            for ( ivar = 0; ivar < H.nvar; ivar++ ) {
                for ( i = H.nx + ExtraLayer; i < H.nx + ExtraLayerTot; i++ ) {
                    sign = 1.0;
                    if ( H.boundary_right == 1 ) {
                        i0 = 2 * H.nx + ExtraLayerTot - i - 1;
                        if ( ivar == IU ) {
                            sign = -1.0;
                        }
                    } else if ( H.boundary_right == 2 ) {
                        i0 = H.nx + ExtraLayer;
                    } else {
                        i0 = i - H.nx;
                    }
                    for ( j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++ ) {
                        /* fprintf(stderr,"PFL %d %d\n",i,j); */
                        Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i0, j, ivar )] * sign;
                        /* fprintf(stderr,"PFL \n"); */
                        MFLOPS ( 1, 0, 0, 0 );
                    }
                }
            }
            //////////////////////////////////////////////////
            // Receive left ghost layer from the first domain.
            //////////////////////////////////////////////////
            fprintf ( stderr, "iProc %i: Receiving left ghost layer from iProc %i\n", H.rank, H.rank-1 );

            // Send two columns at once
            MPI_Recv ( &Hv->uold[IHv ( 0, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 0, MPI_COMM_WORLD, &status );
//			MPI_Recv( &Hv->uold[IHv(0, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 0, MPI_COMM_WORLD, &status);
//			MPI_Recv( &Hv->uold[IHv(1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 1, MPI_COMM_WORLD, &status);

            fprintf ( stderr, "iProc %i: Receiving left ghost layer from iProc %i. Done\n", H.rank, H.rank-1 );

            //////////////////////////////////////////////
            // Send right copy layer to the second domain.
            //////////////////////////////////////////////
            fprintf ( stderr, "iProc %i: Sending left copy layer to iProc %i\n", H.rank, H.rank-1 );

            // Send two columns at once

            MPI_Send ( &Hv->uold[IHv ( 2, 0, ID )], 1, H.mpi_hydro_vector_type, H.rank-1, 2, MPI_COMM_WORLD );
//			MPI_Send( &Hv->uold[IHv(2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 2, MPI_COMM_WORLD);
//			MPI_Send( &Hv->uold[IHv(3, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 3, MPI_COMM_WORLD);

            fprintf ( stderr, "iProc %i: Sending left copy layer to iProc %i. Done\n", H.rank, H.rank-1 );
        }

        fprintf ( stdout,"All layers exchanged\n" );
        //	Free( Hold );

    } else {
        // Lower boundary
        j0 = 0;
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( j = 0; j < ExtraLayer; j++ ) {
                sign = 1.0;
                if ( H.boundary_down == 1 ) {
                    j0 = ExtraLayerTot - j - 1;
                    if ( ivar == IV ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_down == 2 ) {
                    j0 = ExtraLayerTot;
                } else {
                    j0 = H.ny + j;
                }
                for ( i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++ ) {
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i, j0, ivar )] * sign;
                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }

        // Upper boundary
        for ( ivar = 0; ivar < H.nvar; ivar++ ) {
            for ( j = H.ny + ExtraLayer; j < H.ny + ExtraLayerTot; j++ ) {
                sign = 1.0;
                if ( H.boundary_up == 1 ) {
                    j0 = 2 * H.ny + ExtraLayerTot - j - 1;
                    if ( ivar == IV ) {
                        sign = -1.0;
                    }
                } else if ( H.boundary_up == 2 ) {
                    j0 = H.ny + 1;
                } else {
                    j0 = j - H.ny;
                }
                for ( i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++ ) {
                    Hv->uold[IHv ( i, j, ivar )] = Hv->uold[IHv ( i, j0, ivar )] * sign;
                    MFLOPS ( 1, 0, 0, 0 );
                }
            }
        }
    }
}                               // MPI_get_boundary_simple









/**
 * @brief Exchange the boundary conditions
 * 
 * Exchange the boundary conditions with neighboring domains that are on
 * different processes.
 * 
 * @param idim ...
 * @param H ...
 * @param Hv ...
 * @return void
 */
void MPI_make_boundary ( long idim, const hydroparam_t H, hydrovar_t * Hv ) {

    LOC ( H.rank );
    
    // Allocate MPI_req !!!
    MPI_Request *MPI_req;
    MPI_req = malloc ( 4 * sizeof ( MPI_Request ) );

    // Initiate send and receive requests
    MPI_get_boundary_start ( idim, H, Hv, MPI_req );

    // Make sure the data was successfully exchanged before we continue.
    MPI_get_boundary_end ( idim, H, Hv, MPI_req );

    // (CR) Debug
    // MPI_Barrier( MPI_COMM_WORLD );
    Free ( MPI_req );

} // MPI_make_boundary











