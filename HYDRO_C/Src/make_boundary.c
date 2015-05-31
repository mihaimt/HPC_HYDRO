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

 
    if (idim == 1) {

        // Left boundary
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

	/* fprintf(stderr,"PFL H.nvar %d H.nx %d\n",H.nvar,H.nx);
	fprintf(stderr,"PFL ExtraLayer %d ExtraLayerTot %d\n",ExtraLayer,ExtraLayerTot);
	fprintf(stderr,"PFL H.jmin %d H.jmax %d\n",H.jmin,H.jmax); */

        // Right boundary
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
		    /*		  fprintf(stderr,"PFL \n"); */

                    MFLOPS(1, 0, 0, 0);
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

/*
** Exchange the boundary conditions with neighboring domains that are on
** different processes. We use MPI_IRecv() and MPI_ISend() to exchange the
** data. This means we have to make sure that all data was transvered before
** we can do the computation for the most outer layer of the grid using
** MPI_get_boundary_end().
*/
void
MPI_get_boundary_start(long idim, const hydroparam_t H, hydrovar_t * Hv, MPI_Request *MPI_req)
{
    long i, ivar, i0, j, j0, k;
    double sign;
    WHERE("MPI_get_boundary_start");
		
//	MPI_Irecv(values, count, type, source, tag, comm, req)
//	H.iMPIError = MPI_Irecv(Hv, 1, columntype, H.iProc-1, tag, MPI_COMM_WORLD,req);
//	MPI_Isend(values, count, datatype, dest, tag, comm, req)
//	H.iMPIError = MPI_Isend(Hv, 1, columntype, H.iProc-1, tag, MPI_COMM_WORLD,MPI_Request *request);
 
	if (idim == 1) {
		// Make sure MPI_req is allocated
		assert(MPI_req != NULL);

		// Get values from the left domain.
		if (H.iProc == 0)
		{
			// Set physical boundary conditions for the left ghost layer
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
		} else 	{

			// Exchange MPI boundary conditions to the left (send two columns at once)
			MPI_Irecv( &Hv->uold[IHv(0, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 2, MPI_COMM_WORLD, MPI_req );
			MPI_Isend( &Hv->uold[IHv(2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 0, MPI_COMM_WORLD, MPI_req+1 );

			// Exchange MPI boundary conditions to the left (one column per call)
//			MPI_Irecv( &Hv->uold[IHv(0, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 2, MPI_COMM_WORLD, MPI_req );
//			MPI_Irecv( &Hv->uold[IHv(1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 3, MPI_COMM_WORLD, MPI_req+1 );

//			MPI_Isend( &Hv->uold[IHv(2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 0, MPI_COMM_WORLD, MPI_req+2 );
//			MPI_Isend( &Hv->uold[IHv(3, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 1, MPI_COMM_WORLD, MPI_req+3 );
		}

		// Get values from the right domain.
		if (H.iProc == H.iNProc-1)
		{
			// Set physical boundary conditions for the right ghost layer 
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
		} else {
			// Exchange MPI boundary conditions to the right (send two columns at once)
			MPI_Irecv( &Hv->uold[IHv(H.nx+ExtraLayer, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 0, MPI_COMM_WORLD, MPI_req+2 );
			MPI_Isend( &Hv->uold[IHv(H.nx+ExtraLayer-2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 2, MPI_COMM_WORLD, MPI_req+3 );

			// Exchange MPI boundary conditions to the left (one column per call)
//			MPI_Irecv( &Hv->uold[IHv(H.nx+ExtraLayer, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 0, MPI_COMM_WORLD, MPI_req+4 );
//			MPI_Irecv( &Hv->uold[IHv(H.nx+ExtraLayer+1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 1, MPI_COMM_WORLD, MPI_req+5 );

//			MPI_Isend( &Hv->uold[IHv(H.nx+ExtraLayer-2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 2, MPI_COMM_WORLD, MPI_req+6 );
//			MPI_Isend( &Hv->uold[IHv(H.nx+ExtraLayer-1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 3, MPI_COMM_WORLD, MPI_req+7 );
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
}                               // MPI_get_boundary_start

/*
** Make sure that all boundary cells have been successfully exchanged between
** the processes before we continue.
*/
void
MPI_get_boundary_end(long idim, const hydroparam_t H, hydrovar_t * Hv, MPI_Request *MPI_req)
{
	int count, offset, MPIError;
	//MPI_Status *status;

	//status = malloc(8*sizeof(status));

	MPI_Status status;
	offset = 0;
/*
	if ( idim == 1) {
		// Dont do this if we sweep the grid in ny direction.
		if (H.iProc == 0 || H.iProc == H.iNProc-1)
		{
			// For the most outer domains we have less b.c. to exchange.
			count = 1;
			if ( H.iProc == 0 )
			{
				offset = 4;
			}
		} else {
			count = 2;
		}
		MPIError = MPI_Waitall(4*count, MPI_req+offset, status);
		assert( MPIError == 0 );
	}
*/

	// I have no idea why MPI_Waitall() doesnt work
	// Here we send two columns at once
	if ( idim == 1 && H.iNProc > 1) {
		// Dont do this if we sweep the grid in ny direction.
		if (H.iProc == 0 )
		{
			// Only exchanged cells with the right layer
			// (CR) Debug
			fprintf(stderr,"Rank %i: MPI_Wait 2\n", H.iProc);
			MPI_Wait( MPI_req+2, &status);
			fprintf(stderr,"Rank %i: MPI_Wait 3\n", H.iProc);
			MPI_Wait( MPI_req+3, &status);
		} else if ( H.iProc == H.iNProc-1 ) {
			// Only exchanged cells with the left layer
			// (CR) Debug
			fprintf(stderr,"Rank %i: MPI_Wait 0\n", H.iProc);
			MPI_Wait( MPI_req, &status);
			fprintf(stderr,"Rank %i: MPI_Wait 1\n", H.iProc);
			MPI_Wait( MPI_req+1, &status);
		} else {
			MPI_Wait( MPI_req, &status);
			MPI_Wait( MPI_req+1, &status);
			MPI_Wait( MPI_req+2, &status);
			MPI_Wait( MPI_req+3, &status);
		}
	}
/*
	if ( idim == 1) {
		// Dont do this if we sweep the grid in ny direction.
		if (H.iProc == 0 )
		{
			// Only exchanged cells with the right layer
			MPI_Wait( MPI_req+4, status);
			MPI_Wait( MPI_req+5, status);
			MPI_Wait( MPI_req+6, status);
			MPI_Wait( MPI_req+7, status);
		} else if ( H.iProc == H.iNProc-1 ) {
			// Only exchanged cells with the left layer
			MPI_Wait( MPI_req, status);
			MPI_Wait( MPI_req+1, status);
			MPI_Wait( MPI_req+2, status);
			MPI_Wait( MPI_req+3, status);
		} else {
			MPI_Wait( MPI_req, status);
			MPI_Wait( MPI_req+1, status);
			MPI_Wait( MPI_req+2, status);
			MPI_Wait( MPI_req+3, status);
			MPI_Wait( MPI_req+4, status);
			MPI_Wait( MPI_req+5, status);
			MPI_Wait( MPI_req+6, status);
			MPI_Wait( MPI_req+7, status);
		}
	}
*/
//	Free( status );
}                               // MPI_get_boundary_end

/*
** Here we use MPI_sendrecv rather than more complicated stuff.
*/
void
MPI_get_boundary(long idim, const hydroparam_t H, hydrovar_t * Hv)
{
    long i, ivar, i0, j, j0, k;
    double sign;
	MPI_Status status;
    WHERE("MPI_get_boundary");

	// Use:		
	// MPI_Sendrecv( sendbuff, sendcount, sendtype, dest, sendtag,
	//				 recvbuff, recvcount, recvtype, source, recvtag, comm, status);
	// MPI_Send( buf, count, datatype, dest, tag, comm);
	// MPI_Recv( buf, count, datatype, source, tag, comm, status);
	if (idim == 1) {
		// (CR) Debug
		fprintf(stderr,"Rank %i: Send copy layer to the right and receive ghost layer from the left\n", H.iProc);
		
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
		Hold->uold = (double *) calloc(H.nvar * H.nxt * H.nyt, sizeof(double));

		if ( H.iProc == 1 )
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
		}
        
		///////////////////////////////////////////////////////////////////////////////////////////////
		// Send right copy layer to the right domain and receive left ghost layer from the left domain.
		///////////////////////////////////////////////////////////////////////////////////////////////
		if (H.iProc == 0)
		{
			// Rank 0 (most left domain)
			// Set physical boundary conditions for the left ghost cells
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
			// Send right copy layer to the right domain (but dont receive the left ghost layer from the left).
			// (CR) Debug
			fprintf(stderr, "iProc %i: Sending right copy layer to iProc %i\n", H.iProc, H.iProc+1);

			MPI_Send( &Hv->uold[IHv(H.nx+ExtraLayer-2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 0, MPI_COMM_WORLD);
			MPI_Send( &Hv->uold[IHv(H.nx+ExtraLayer-1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 1, MPI_COMM_WORLD);

			fprintf(stderr, "iProc %i: Sending right copy layer to iProc %i. Done\n", H.iProc, H.iProc+1);
		} else if ( H.iProc == H.iNProc-1) {
			// Rank iNProc-1 (most right domain)
			// Receive left ghost layer from the left domain (but dont send the right copy layer to the right).
			// (CR) Debug
			fprintf(stderr, "iProc %i: Receiving left ghost layer from iProc %i\n", H.iProc, H.iProc-1);

			MPI_Recv( &Hv->uold[IHv(0, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 0, MPI_COMM_WORLD, &status);
			MPI_Recv( &Hv->uold[IHv(1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 1, MPI_COMM_WORLD, &status);

			fprintf(stderr, "iProc %i: Receiving left ghost layer from iProc %i. Done\n", H.iProc, H.iProc-1);
		} else {
			// Exchange MPI boundary conditions
			// (CR) Debug
			fprintf(stderr, "iProc %i: Sendrecv() 1\n", H.iProc);
			// Send the first right copy layer (nx+ExtraLayer-2) to the right domain (iProc+1) and
			// receive the first left ghost layer (0) from the left domain (iProc-1)
			MPI_Sendrecv( &Hv->uold[IHv(H.nx+ExtraLayer-2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 4,
					      &Hv->uold[IHv(0, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 5, MPI_COMM_WORLD, &status);

			fprintf(stderr, "iProc %i: Sendrecv() 2\n", H.iProc);

			// Send the second right copy layer (nx+ExtraLayer-2) to the right domain (iProc+1) and
			// receive the second left ghost layer (0) from the left domain (iProc-1)
			MPI_Sendrecv( &Hv->uold[IHv(H.nx+ExtraLayer-1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 6,
					      &Hv->uold[IHv(1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 7, MPI_COMM_WORLD, &status);
		}
		
		// (CR) Debug
		fprintf(stderr,"Rank %i: Synchronizing\n", H.iProc);
		MPI_Barrier( MPI_COMM_WORLD );
		fprintf(stderr,"Rank %i: Synchronizing. Done\n", H.iProc);

		fprintf(stderr,"Rank %i: Send copy layer to the left and receive ghost layer from the right\n", H.iProc);

		///////////////////////////////////////////////////////////////////////////////////////////////
		// Send left copy layer to the left domain and receive right ghost layer from the right domain.
		///////////////////////////////////////////////////////////////////////////////////////////////
		if (H.iProc == 0)
		{
			// Rank 0 (most left domain)
			// Receive right ghost cells from the right domain (but dont send the left copy layer to the left).
//			fprintf(stderr, "Sending right copy layer for iProc %i\n", H.iProc);
			// MPI_Recv( buf, count, datatype, source, tag, comm, status)
			MPI_Recv( &Hv->uold[IHv(H.nx+ExtraLayer+1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 0, MPI_COMM_WORLD, &status);
			MPI_Recv( &Hv->uold[IHv(H.nx+ExtraLayer, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 1, MPI_COMM_WORLD, &status);
		} else if (H.iProc == H.iNProc - 1) {
			// Rank: iNProc-1 (most right domain)
			// Set physical boundary conditions for the right ghost cells 
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
			// Send left copy layer to the left domain (but dont receive the right ghost layer from the right).
			MPI_Send( &Hv->uold[IHv(2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 0, MPI_COMM_WORLD);
			MPI_Send( &Hv->uold[IHv(3, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 1, MPI_COMM_WORLD);
		} else {
			// Exchange MPI boundary conditions

			// Send the first left copy layer (2) and receive the first right ghost layer (nx+ExtraLayer)
			MPI_Sendrecv( &Hv->uold[IHv(2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 0,
					      &Hv->uold[IHv(H.nx+ExtraLayer, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 1, MPI_COMM_WORLD, &status);

			// Send the second left copy layer (3) and receive the second right ghost layer (nx+ExtraLayer+1)
			MPI_Sendrecv( &Hv->uold[IHv(3, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 0,
					      &Hv->uold[IHv(H.nx+ExtraLayer+1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 1, MPI_COMM_WORLD, &status);
		}
		Free( Hold );

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
}                               // MPI_get_boundary

/*
 * The most simple version possible. Only two processes and I implemented every send and receive
 * manually.
 */
void
MPI_get_boundary_simple(long idim, const hydroparam_t H, hydrovar_t * Hv)
{
    long i, ivar, i0, j, j0, k;
    double sign;
	MPI_Status status;
	int source, dest;
    WHERE("MPI_get_boundary_simple");

	// (CR) Debug
	assert( H.iNProc == 2 );

	fprintf(stdout,"Rank %i: MPI_get_boundary_simple()\n", H.iProc);
	// Use:		
	// MPI_Sendrecv( sendbuff, sendcount, sendtype, dest, sendtag,
	//				 recvbuff, recvcount, recvtype, source, recvtag, comm, status);
	// MPI_Send( buf, count, datatype, dest, tag, comm);
	// MPI_Recv( buf, count, datatype, source, tag, comm, status);
	if (idim == 1) {
		
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
		if (H.iProc == 0)
		{
			///////////////////////////////////////////////////////////////////
			// Make physical b.c. for the left ghost layer of the first domain.
			///////////////////////////////////////////////////////////////////
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
			//////////////////////////////////////////////
			// Send right copy layer to the second domain.
			//////////////////////////////////////////////
			fprintf(stderr, "iProc %i: Sending right copy layer to iProc %i\n", H.iProc, H.iProc+1);

			// Send two columns at once
			MPI_Send( &Hv->uold[IHv(H.nx+ExtraLayer-2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 0, MPI_COMM_WORLD);
//			MPI_Send( &Hv->uold[IHv(H.nx+ExtraLayer-2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 0, MPI_COMM_WORLD);
//			MPI_Send( &Hv->uold[IHv(H.nx+ExtraLayer-1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 1, MPI_COMM_WORLD);

			fprintf(stderr, "iProc %i: Sending right copy layer to iProc %i. Done\n", H.iProc, H.iProc+1);
			////////////////////////////////////////////////////
			// Receive right ghost layer from the second domain.
			////////////////////////////////////////////////////
			fprintf(stderr, "iProc %i: Receiving right ghost layer from iProc %i\n", H.iProc, H.iProc+1);

			// Send two columns at once

			MPI_Recv( &Hv->uold[IHv(H.nx+ExtraLayer, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 2, MPI_COMM_WORLD, &status);
//			MPI_Recv( &Hv->uold[IHv(H.nx+ExtraLayer, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 2, MPI_COMM_WORLD, &status);
//			MPI_Recv( &Hv->uold[IHv(H.nx+ExtraLayer+1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc+1, 3, MPI_COMM_WORLD, &status);

			fprintf(stderr, "iProc %i: Receiving right ghost layer from iProc %i. Done\n", H.iProc, H.iProc+1);
		} else  {
			/////////////////////////////////////////////////////////////////////
			// Make physical b.c. for the right ghost layer of the second domain.
			/////////////////////////////////////////////////////////////////////
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
			//////////////////////////////////////////////////
			// Receive left ghost layer from the first domain.
			//////////////////////////////////////////////////
			fprintf(stderr, "iProc %i: Receiving left ghost layer from iProc %i\n", H.iProc, H.iProc-1);

			// Send two columns at once
			MPI_Recv( &Hv->uold[IHv(0, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 0, MPI_COMM_WORLD, &status);
//			MPI_Recv( &Hv->uold[IHv(0, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 0, MPI_COMM_WORLD, &status);
//			MPI_Recv( &Hv->uold[IHv(1, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 1, MPI_COMM_WORLD, &status);

			fprintf(stderr, "iProc %i: Receiving left ghost layer from iProc %i. Done\n", H.iProc, H.iProc-1);

			//////////////////////////////////////////////
			// Send right copy layer to the second domain.
			//////////////////////////////////////////////
			fprintf(stderr, "iProc %i: Sending left copy layer to iProc %i\n", H.iProc, H.iProc-1);

			// Send two columns at once

			MPI_Send( &Hv->uold[IHv(2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 2, MPI_COMM_WORLD);
//			MPI_Send( &Hv->uold[IHv(2, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 2, MPI_COMM_WORLD);
//			MPI_Send( &Hv->uold[IHv(3, 0, ID)], 1, H.MPI_Hydro_vars, H.iProc-1, 3, MPI_COMM_WORLD);

			fprintf(stderr, "iProc %i: Sending left copy layer to iProc %i. Done\n", H.iProc, H.iProc-1);
		}

		fprintf(stdout,"All layers exchanged\n");
	//	Free( Hold );

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
}                               // MPI_get_boundary_simple
/*
** Exchange the boundary conditions with neighboring domains that are on
** different processes. 
*/
void
MPI_make_boundary(long idim, const hydroparam_t H, hydrovar_t * Hv)
{
	// Allocate MPI_req !!!
	MPI_Request *MPI_req;
	MPI_req = malloc(8*sizeof(MPI_Request));

	// Initiate send and receive requests
	MPI_get_boundary_start(idim, H, Hv, MPI_req);
	
	// Make sure the data was successfully exchanged before we continue.
	MPI_get_boundary_end(idim, H, Hv, MPI_req);


	// (CR) Debug
//	MPI_Barrier( MPI_COMM_WORLD );
	// Free MPI_req
	Free(MPI_req);
}                               // MPI_make_boundary

//EOF
