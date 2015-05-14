#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
	int iProc;
	int iNumProc;
        int i;

	/* Initialize MPI library */
	i = MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD,&iNumProc);
	MPI_Comm_rank(MPI_COMM_WORLD,&iProc);

	printf("MPI_Init: %i process %i of %i\n",i,iProc,iNumProc-1);

	i = MPI_Finalize();

	printf("MPI_Finalize: %i process %i of %i\n",i,iProc,iNumProc-1);

	return 0;
}

