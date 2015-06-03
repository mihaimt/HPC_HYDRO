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

#include "debug.h"

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
    long nvtk = 0; // counter for steps / states written to file
    char outnum[160];
    long time_output = 0;

    // double output_time = 0.0;
    double next_output_time = 0;
    double start_time = 0, end_time = 0;
    double start_iter = 0, end_iter = 0;
    double iter_time = 0; // time for one iteration
    double it_min = 0, it_max = 0; // used for global stats about iter times
    double elaps = 0;

    // The smallest possible time step for the entire computational domain.
    double dtmin = 0.0;

    // The maximal time used to run
    double tmax = 0.0;



    // INIT EVERYTHING
    //------------------------------------------------------------------------

    // global program timer
    start_time = cclock();

    // Initialize MPI library (and allocate memory for the MPI variables).
    MPI_init ( &H, &argc, &argv );

    // Initialise OpenMP
    OPENMP_init ( &H );

    // Parse command line variables and input file(s)
    process_args ( argc, argv, &H );

    // Domain decomposition
    MPI_domain_decomp ( &H );

    // Initialize the hydro variables and set initial conditions.
    MPI_hydro_init ( &H, &Hv );
    
    // initialize the output file to write the iteration timings into
    if ( WRITE_TIMING ) {
        timingfile_init ( &H );
    }


    TRC ( H.rank, "nxt=%i nyt=%i", H.nxt, H.nyt );
    INF_if ( H.rank==0, "Hydro setup successful\n" );
    INF_if ( H.rank==0, "   use MPI:      %s\n", __str(USE_MPI) );
    INF_if ( H.rank==0, "   use OPENMP:   %s\n", __str(USE_OPENMP) );
    INF_if ( H.rank==0, "   use DEBUG:    %s\n", __str(DEBUG) );
    INF_if ( H.rank==0, "   use ASSSERTS: %s\n", __str(DO_ASSERTS) );
    INF_if ( H.rank==0, "   use COLOR:    %s\n\n", __str(USE_COLOR) );

    // give every proc time to start up and read from filesys
    if ( USE_MPI ) { MPI_Barrier( MPI_COMM_WORLD ); }

    // write the initial state to file
    if ( WRITE_INIT_STATE ) {
        vtkfile ( nvtk, H, &Hv );
    }

    // check in which mode to write the state to the vtk files
    // * if {dtoutput} is defined in the input file, then write state after this
    //   amount of physical time has passed
    // * else (dtoutput==0) simply write output after {noutput} steps happened
    //   (default value is VERY big! -> no output )
    if ( H.dtoutput > 0 ) {
        // outputs are in physical time not in time steps
        time_output = 1;
        next_output_time = next_output_time + H.dtoutput;
    }

    INF_if ( H.rank==0, "starting mainloop...\n" );

    //-------------------------------------------------------------------------
    // The main loop
    //-------------------------------------------------------------------------

    while ( ( H.t < H.tend ) && ( H.nstep < H.nstepmax ) ) {

        DBG_if ( H.rank==0, "main loop | start: nstep = %i \n", H.nstep);

        start_iter = cclock();
        outnum[0] = 0; // delete string by setting first char to 0 byte
        flops = 0;

        // Calculate new time step for every even step.
        if ( ( H.nstep % 2 ) == 0 ) {

            compute_deltat ( &dt, H, &Hw, &Hv, &Hvw );
            if ( H.nstep == 0 ) {
                dt = dt / 2.0;
            }

            // Get the smallest possible time step for all processes.
            if ( USE_MPI ) {
                H.mpi_error = MPI_Allreduce ( &dt, &dtmin, 1, MPI_DOUBLE,
                                              MPI_MIN, MPI_COMM_WORLD );
                TRC ( H.rank, "time sync: dt=%.2e dtmin=%.2e", dt, dtmin);
                dt = dtmin;
            }

        }


        // This is the actual calculation
        if ( ( H.nstep % 2 ) == 0 ) {
            MPI_hydro_godunov ( 1, dt, H, &Hv, &Hw, &Hvw );
            MPI_hydro_godunov ( 2, dt, H, &Hv, &Hw, &Hvw );
        }
        else {
            MPI_hydro_godunov ( 2, dt, H, &Hv, &Hw, &Hvw );
            MPI_hydro_godunov ( 1, dt, H, &Hv, &Hw, &Hvw );
        }


        end_iter = cclock();
        H.nstep++;
        H.t += dt;


        // some outdated code that calculates the flops with a rather ad hoc method..
        // I (RK) don't see the point in this...
        // The else part is below thou
/*
        if ( flops > 0 ) {
            double iter_time = ( double ) ( end_iter - start_iter );
            if ( iter_time > 1.e-9 ) {
                double mflops = ( double ) flops / ( double ) 1.e+6 / iter_time;
                sprintf ( outnum, "%.3f Mflops t_iter=%.3fs) ", mflops, iter_time );
            }
        }
        else {
            double iter_time = ( double ) ( end_iter - start_iter );
            sprintf ( outnum, "t_iter=%.3es ", iter_time );
        }
*/


        // get time elapsed in one iteration
        if ( GET_LOCAL_ITER_TIME || GET_GLOBAL_ITER_TIME ) {
            iter_time = end_iter - start_iter;
            sprintf ( outnum, "t_iter=%.3es ", iter_time );
        }

        // get iter timings from all procs to rank 0
        if ( USE_MPI && GET_GLOBAL_ITER_TIME ) {
            it_min = 0.0;
            it_max = 0.0;
            H.mpi_error = MPI_Reduce ( &iter_time, &it_min, 1, MPI_DOUBLE,
                                       MPI_MIN, 0, MPI_COMM_WORLD );
            H.mpi_error = MPI_Reduce ( &iter_time, &it_max, 1, MPI_DOUBLE,
                                       MPI_MAX, 0, MPI_COMM_WORLD );
            if ( H.rank == 0 ) {
                sprintf ( outnum, "%s (global %.3e .. %.3e)", outnum, it_min, it_max );
            }
        }

        // write the time for this step to the output file
        if ( GET_LOCAL_ITER_TIME && WRITE_TIMING ) {
            timingfile_write ( H.nstep, iter_time, H );
        }

        // write the current state to vtk file
        if ( WRITE_INTER_STATE ) {

            // write each nth step
            if ( time_output == 0 ) {
                if ( ( H.nstep % H.noutput ) == 0 ) {
                    vtkfile ( ++nvtk, H, &Hv );
                    sprintf ( outnum, "%s [filenr=%04ld]", outnum, nvtk );
                }
            }
            // or write after dtoutput physical time has passed
            else {
                if ( H.t >= next_output_time ) {
                    vtkfile ( ++nvtk, H, &Hv );
                    next_output_time = next_output_time + H.dtoutput;
                    sprintf ( outnum, "%s [filenr=%04ld]", outnum, nvtk );
                }
            }
        }



        // Synchronize all processes if debugging for nicer output
        if ( DEBUG && USE_MPI ) { MPI_Barrier ( MPI_COMM_WORLD ); }

        // and some infos about the steps on rank 0 (step is already increased!)
        DBG_if ( H.rank == 0, "main loop | end: step=%04li t=%.4e dt=%.4e %s\n", H.nstep-1, H.t, dt, outnum );

    }   // end main loop
    
    INF_if ( H.rank==0, "main loop finished, cleaning up...\n" );


    // write the final state to file
    // be aware that this last step could have another time spacing as the
    // intermediate steps above!
    if ( WRITE_FINAL_STATE ) {
        vtkfile ( ++nvtk, H, &Hv );
    }

    end_time = cclock();
    elaps = ( double ) ( end_time - start_time );


    // Get the largest time for all processes.
    if ( USE_MPI ) {
        H.mpi_error = MPI_Allreduce ( &elaps, &tmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    }


    if ( H.rank == 0 ) {
        if ( elaps < tmax ) {
            elaps = tmax;
        }
        timeToString ( outnum, elaps );
        INF ( "Hydro ends in %ss (%.3lf).\n", outnum, elaps );

        // write a status report
        write_stat ( elaps, H.nstep, nvtk, H );
    }


    // Finalize MPI and free memory.
    MPI_finish ( &H );
    MPI_hydro_finish ( &H, &Hv );
    if ( WRITE_TIMING ) { timingfile_finish( &H ); }


    // test message printing
    if (DEBUG && OFF) {
        TRC ( H.rank, "this is a trace" );
        TRC ( H.rank, "%s", "this is a trace" );
        DBG ( "a debug\n" );
        DBG ( "%s", "a debug\n" );
        ERR ( "an error\n" );
        ERR ( "%s", "an error\n" );
        ERR_at ( -1, "%s", "an error\n" );
        WRN ( "warn\n" );
        WRN ( "%s", "warn\n" );
        INF ( "info\n" );
        INF ( "%s", "info\n" );
        INF_if ( 1==1, "%s", "cond. info, true\n" );
        INF_if ( 1==2, "%s", "cond. info, false\n" );
    }


    return 0;
}

