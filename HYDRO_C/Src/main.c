/*
   A simple 2D hydro code
   (C) Romain Teyssier : CEA/IRFU           -- original F90 code
   (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
   (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
 */

#include <stdio.h>
#include <time.h>

#include "parametres.h"
#include "hydro_funcs.h"
#include "vtkfile.h"
#include "compute_deltat.h"
#include "hydro_godunov.h"
#include "utils.h"

hydroparam_t   H;
hydrovar_t     Hv;
hydrovarwork_t Hvw;
hydrowork_t    Hw;
unsigned long flops = 0;



int main ( int argc, char **argv ) {

    // DECLARE VARS
    //------------------------------------------------------------------------
    
    int nb_th = 1;
    double dt = 0;
    long nvtk = 0;
    char outnum[80];
    long time_output = 0;

    // double output_time = 0.0;
    double next_output_time = 0;
    double start_time = 0, end_time = 0;
    double start_iter = 0, end_iter = 0;
    double elaps = 0;

    // The smallest possible time step for the entire computational domain.
    double dtmin = 0.0;

    // The maximal time used to run
    double tmax = 0.0;

    
    // INIT EVERYTHING
    //------------------------------------------------------------------------

    start_time = cclock();

    // Initialize MPI library (and allocate memory for the MPI variables).
    MPI_init ( &H, &argc, &argv );

    // Parse command line variables and input file(s)
    process_args ( argc, argv, &H );

    // Domain decomposition
    MPI_domain_decomp ( &H );

    int nxtotal = 0.0;
    // (CR) Debug
    MPI_Allreduce ( &H.nx,&nxtotal,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD );
    fprintf ( stderr,"Process %i: nx = %i nxtotal = %i\n",H.rank,H.nx,nxtotal );

    // Initialize the hydro variables and set initial conditions.
    MPI_hydro_init ( &H, &Hv );
    PRINTUOLD ( H, &Hv );

    fprintf ( stderr,"Process %i: nxt=%i nyt=%i\n",H.rank,H.nxt,H.nyt );
    if ( H.rank == 0 ) {
        printf ( "Hydro starts - MPI version \n" );
        printf ( "Running on %i processes\n", H.n_proc );
#ifdef MPI_NO_OUTPUT
        printf ( "NOT WRITING ANY OUTPUT\n" );
#endif
    }

    // vtkfile(nvtk, H, &Hv);
    if ( H.dtoutput > 0 ) {
        // outputs are in physical time not in time steps
        time_output = 1;
        next_output_time = next_output_time + H.dtoutput;
    }

    /*
     ** The main loop.
     */
    while ( ( H.t < H.tend ) && ( H.nstep < H.nstepmax ) ) {
//		fprintf(stderr,"Main loop: nstep = %i \n",H.nstep);
        start_iter = cclock();
        outnum[0] = 0;
        flops = 0;
        if ( ( H.nstep % 2 ) == 0 ) {
            // We calculate the new time step for every even step.
            compute_deltat ( &dt, H, &Hw, &Hv, &Hvw );
            if ( H.nstep == 0 ) {
                dt = dt / 2.0;
            }

            // Get the smallest possible time step for all processes.
            H.mpi_error = MPI_Allreduce ( &dt,&dtmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD );
//			fprintf(stderr,"Process %i: dt = %g dtmin = %g\n",H.iProc,dt,dtmin);
            dt = dtmin;
        }

        // This is the acutal calculation
        if ( ( H.nstep % 2 ) == 0 ) {
            MPI_hydro_godunov ( 1, dt, H, &Hv, &Hw, &Hvw );
            MPI_hydro_godunov ( 2, dt, H, &Hv, &Hw, &Hvw );
        } else {
            MPI_hydro_godunov ( 2, dt, H, &Hv, &Hw, &Hvw );
            MPI_hydro_godunov ( 1, dt, H, &Hv, &Hw, &Hvw );
        }

        end_iter = cclock();
        H.nstep++;
        H.t += dt;

        if ( flops > 0 ) {
            double iter_time = ( double ) ( end_iter - start_iter );
            if ( iter_time > 1.e-9 ) {
                double mflops = ( double ) flops / ( double ) 1.e+6 / iter_time;
                sprintf ( outnum, "%s {%.3f Mflops} (%.3fs)", outnum, mflops, iter_time );
            }
        } else {
            double iter_time = ( double ) ( end_iter - start_iter );
            sprintf ( outnum, "%s (%.3fs)", outnum, iter_time );
        }

        if ( time_output == 0 ) {
            if ( ( H.nstep % H.noutput ) == 0 ) {
                vtkfile ( ++nvtk, H, &Hv );
                sprintf ( outnum, "%s [%04ld]", outnum, nvtk );
            }
        } else {
            if ( H.t >= next_output_time ) {
                vtkfile ( ++nvtk, H, &Hv );
                next_output_time = next_output_time + H.dtoutput;
                sprintf ( outnum, "%s [%04ld]", outnum, nvtk );
            }
        }

        // Synchronize all processes
        MPI_Barrier ( MPI_COMM_WORLD );

        if ( H.rank == 0 ) {
            fprintf ( stdout, "--> step=%-4ld %12.5e, %10.5e %s\n", H.nstep, H.t, dt, outnum );
        }
    }   // end while loop

    end_time = cclock();
    elaps = ( double ) ( end_time - start_time );



    // Get the largest time for all processes.
    H.mpi_error = MPI_Allreduce ( &elaps,&tmax,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD );

    if ( H.rank == 0 ) {
        if ( elaps < tmax ) {
            elaps = tmax;
        }
        timeToString ( outnum, elaps );
        fprintf ( stdout, "Hydro ends in %ss (%.3lf).\n", outnum, elaps );
    }

    // Finalize MPI and free memory.
    MPI_hydro_finish ( &H, &Hv );

    return 0;
}

