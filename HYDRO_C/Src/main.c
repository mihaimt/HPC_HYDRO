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

hydroparam_t H;
hydrovar_t Hv;                  // nvar
hydrovarwork_t Hvw;             // nvar
hydrowork_t Hw;
unsigned long flops = 0;

int
main(int argc, char **argv)
{
	int nb_th=1;
	double dt = 0;
	long nvtk = 0;
	char outnum[80];
	long time_output = 0;

	// double output_time = 0.0;
	double next_output_time = 0;
	double start_time = 0, end_time = 0;
	double start_iter = 0, end_iter = 0;
	double elaps = 0;

	/*
	** The smallest possible time step for the entire computational
	** domain.
	*/
	double dtmin = 0.0;

	start_time = cclock();
	process_args(argc, argv, &H);
	
	/* Initialize MPI library and allocate memory for the variables. */
	MPI_hydro_init(&H, &Hv,&argc,&argv);
	PRINTUOLD(H, &Hv);

	printf("Hydro starts - MPI version \n");

	// vtkfile(nvtk, H, &Hv);
	if (H.dtoutput > 0) 
	{	
		// outputs are in physical time not in time steps
		time_output = 1;
		next_output_time = next_output_time + H.dtoutput;
	}

	/*
	 ** The main loop.
	 */
	while ((H.t < H.tend) && (H.nstep < H.nstepmax)) 
	{	
		start_iter = cclock();
		outnum[0] = 0;
		flops = 0;
		if ((H.nstep % 2) == 0)
		{
			/* We calculate the new time step for every even step. */
			compute_deltat(&dt, H, &Hw, &Hv, &Hvw);
			if (H.nstep == 0) {
				dt = dt / 2.0;
			}
				
			/* Get the smallest possible time step for all processes. */
			H.iMPIError = MPI_Allreduce(&dt,&dtmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
			dt = dtmin;
		}
		
		/* This is the acutal calculation */
		if ((H.nstep % 2) == 0) {
			MPI_hydro_godunov(1, dt, H, &Hv, &Hw, &Hvw);
			MPI_hydro_godunov(2, dt, H, &Hv, &Hw, &Hvw); 
		} else {
			MPI_hydro_godunov(2, dt, H, &Hv, &Hw, &Hvw);
			MPI_hydro_godunov(1, dt, H, &Hv, &Hw, &Hvw);
		}

		end_iter = cclock();
		H.nstep++;
		H.t += dt;

		if (flops > 0) {
			double iter_time = (double) (end_iter - start_iter);
			if (iter_time > 1.e-9) {
				double mflops = (double) flops / (double) 1.e+6 / iter_time;
				sprintf(outnum, "%s {%.3f Mflops} (%.3fs)", outnum, mflops, iter_time);
			}
		} else {
			double iter_time = (double) (end_iter - start_iter);
			sprintf(outnum, "%s (%.3fs)", outnum, iter_time);
		}
			
		if (time_output == 0) {
			if ((H.nstep % H.noutput) == 0) {
				vtkfile(++nvtk, H, &Hv);
				sprintf(outnum, "%s [%04ld]", outnum, nvtk);
			}
		} else {
			if (H.t >= next_output_time) {
				vtkfile(++nvtk, H, &Hv);
				next_output_time = next_output_time + H.dtoutput;
				sprintf(outnum, "%s [%04ld]", outnum, nvtk);
			}
		}

		fprintf(stdout, "--> step=%-4ld %12.5e, %10.5e %s\n", H.nstep, H.t, dt, outnum);
	}   // end while loop
		
	/* Finalize MPI and free memory. */
	MPI_hydro_finish(&H, &Hv);
	end_time = cclock();
	elaps = (double) (end_time - start_time);
	timeToString(outnum, elaps); 

	fprintf(stdout, "Hydro ends in %ss (%.3lf).\n", outnum, elaps);

	return 0;
}

