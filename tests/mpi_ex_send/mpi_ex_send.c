/*
** This program is a very simple example of how to use MPI_Send() and MPI_Recv()
** to do point-to-point communication between two nodes.
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
	int iProc = 0;
	int iNumProc = 0;
        int i = 0;
	int iSource = 0;
	int iDest = 0;
	int N = 0;
	int tag = 100;
	MPI_Status *status;

	/* We need to allocate memory for this variable! */
	status = malloc(sizeof(MPI_Status));
	
	/* Initialize MPI library */
	i = MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD,&iNumProc);
	MPI_Comm_rank(MPI_COMM_WORLD,&iProc);
	
	if (i != 0)
	{
		printf("MPI_Init: Error %i\n",i);
		exit(1);
	}
	
	/* Make sure we have more than one process */
	if (iNumProc == 1)
	{
		printf("Only %i processes initialized.\n",iNumProc);
		exit(1);
	}

	iDest = iNumProc-1;

	if (iProc == iSource)
	{
		/* We only want one initialization message */	
		printf("MPI: %i processes initialized.\n",iNumProc);

		N = 1000;
		printf("Process %i: Sending data to process %i (N = %i)\n",iProc,iDest,N);
		i = MPI_Send(&N,1,MPI_INT,iDest,tag,MPI_COMM_WORLD);
	}

	if (iProc == iDest)
	{
		i = MPI_Recv(&N,1,MPI_INT,iSource,tag,MPI_COMM_WORLD,status);
		printf("Process %i: Receiving data from process %i (N = %i)\n",iProc,iDest,N);
	}

	/* Finalize MPI library */
	i = MPI_Finalize();

	return 0;
}

