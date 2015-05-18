/*
** This program is a very simple example of how to use MPI_Isend() and MPI_Irecv()
** to do point-to-point communication between two nodes.
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define NX 6
#define NY 6

#define LEFT_DOMAIN 0
#define MIDDLE_DOMAIN 1
#define RIGHT_DOMAIN 2

int main(int argc, char *argv[]) {
	// MPI variables
	int iProc = 0;
	int iNumProc = 0;
	int tag = 0;
	MPI_Status status;
	MPI_Datatype MPI_COLUMN_TYPE;
	MPI_Request MPI_Send_req_left, MPI_Send_req_right;
	MPI_Request MPI_Recv_req_left, MPI_Recv_req_right;
	// delay
	//lont wait = 1000; // 1sec
	struct timespec t1, t2;
	int iError = 0;
	int N = 0;
	// Array that we want to exchange
	double a[NY][NX] = {{0}};	

	/* Initialize MPI library */
	iError = MPI_Init(&argc,&argv);
	if (iError != 0)
	{
		printf("MPI_Init: Error %i\n",iError);
		exit(1);
	}

	MPI_Comm_size(MPI_COMM_WORLD,&iNumProc);
	MPI_Comm_rank(MPI_COMM_WORLD,&iProc);

	printf("Process %i: MPI successfully initialized.\n",iProc);
	
	/* Make sure we have at least three processes */
	if (iNumProc < 3)
	{
		printf("Only %i processes initialized.\n",iNumProc);
		exit(1);
	}
	
	// Define a new MPI data type
	MPI_Type_vector(NY, 1, NX, MPI_DOUBLE, &MPI_COLUMN_TYPE);
	MPI_Type_commit(&MPI_COLUMN_TYPE);

	// 1 sec delay
	t1.tv_sec = 1;
	t1.tv_nsec = 0;
		
	// random delay to avoid problems with output
	t2.tv_sec = 2*iProc+1;
	
	// Fill the array
	for (int i = 0; i < NY; i++)
	{
		for (int j = 0; j < NX; j++)
		{
			a[i][j] = iProc;
		}
	}
	
	if (iProc == MIDDLE_DOMAIN)
	{
		// delay
//		nanosleep(&t2, NULL);
		printf("iProc: %i \n      ",iProc);
		for (int i = 0; i < NY; i++)
		{
			for (int j = 0; j < NX; j++)
			{
				printf("%4.0g ",a[i][j]);
			}
			printf("\n      ");
		}
		printf("\n");
	}
	// synchronize
	MPI_Barrier(MPI_COMM_WORLD);

	if (iProc == 0)
	{
		printf("Start exchanging MPI data.\n");
	}


	if (iProc == MIDDLE_DOMAIN)
	{
		// Receive column 0 and 1 from process 0
//		MPI_Irecv(&a[0][0], 2, MPI_COLUMN_TYPE, LEFT_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Recv_req_left);
		// Receive column 4 and 5 from process 2
//		MPI_Irecv(&a[0][4], 2, MPI_COLUMN_TYPE, RIGHT_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Recv_req_right);
		// Just one column to debug
		MPI_Irecv(&a[0][0], 1, MPI_COLUMN_TYPE, LEFT_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Recv_req_left);
		MPI_Irecv(&a[0][5], 1, MPI_COLUMN_TYPE, RIGHT_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Recv_req_right);


		// Send column 2 and 3 to the other processes
//		MPI_Isend(&a[0][2], 2, MPI_COLUMN_TYPE, LEFT_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Send_req_left);
//		MPI_Isend(&a[0][2], 2, MPI_COLUMN_TYPE, RIGHT_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Send_req_right);
	}

	if (iProc == LEFT_DOMAIN)
	{
		// Receive column 2 and 3 from process 1
//		MPI_Irecv(&a[0][2], 2, MPI_COLUMN_TYPE, 1, tag, MPI_COMM_WORLD, &MPI_Recv_req_right);

		// Send column 4 and 5 to the process 1
//		MPI_Isend(&a[0][4], 2, MPI_COLUMN_TYPE, MIDDLE_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Send_req_right);
		MPI_Isend(&a[0][4], 1, MPI_COLUMN_TYPE, MIDDLE_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Send_req_right);

	}

	if (iProc == RIGHT_DOMAIN)
	{
		// Receive column 2 and 3 from process 1
//		MPI_Irecv(&a[0][2], 2, MPI_COLUMN_TYPE, 1, tag, MPI_COMM_WORLD, &MPI_Recv_req_left);

		// Send column 0 and 1 to the process 1
//		MPI_Isend(&a[0][4], 2, MPI_COLUMN_TYPE, MIDDLE_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Send_req_left);
		// Send column 0 and 1 to the process 1
		MPI_Isend(&a[0][4], 1, MPI_COLUMN_TYPE, MIDDLE_DOMAIN, tag, MPI_COMM_WORLD, &MPI_Send_req_left);
	}
/*
	// synchronize
	MPI_Barrier(MPI_COMM_WORLD);

	if (iProc == 0)
	{
		printf("Requests sent.\n");
	}
*/
	nanosleep(&t1, NULL);

	// wait until all communication is done
	if (iProc == LEFT_DOMAIN)
	{
//		MPI_Wait(&MPI_Send_req_left, &status);
		MPI_Wait(&MPI_Send_req_right, &status);
//		MPI_Wait(&MPI_Recv_req_left, &status);
//		MPI_Wait(&MPI_Recv_req_right, &status);
		printf("Process %i: Done.\n",iProc);
	} else if (iProc == MIDDLE_DOMAIN) {
//		MPI_Wait(&MPI_Send_req_left, &status);
//		MPI_Wait(&MPI_Send_req_right, &status);
		MPI_Wait(&MPI_Recv_req_left, &status);
		MPI_Wait(&MPI_Recv_req_right, &status);
		printf("Process %i: Done.\n",iProc);
	} else if (iProc == RIGHT_DOMAIN) {
		MPI_Wait(&MPI_Send_req_left, &status);
//		MPI_Wait(&MPI_Send_req_right, &status);
//		MPI_Wait(&MPI_Recv_req_left, &status);
//		MPI_Wait(&MPI_Recv_req_right, &status);
		printf("Process %i: Done.\n",iProc);
	}

	if (iProc == 0)
	{
		printf("Done.\n");
	}

	// synchronize
	MPI_Barrier(MPI_COMM_WORLD);

	if (iProc == MIDDLE_DOMAIN)
	{
		// delay
//		nanosleep(&t2, NULL);

		printf("iProc: %i \n      ",iProc);
		for (int i = 0; i < NY; i++)
		{
			for (int j = 0; j < NX; j++)
			{
				printf("%4.0g ",a[i][j]);
			}
			printf("\n      ");
		}
		printf("\n");
	}
	printf("\n");

	/* Finalize MPI library */
	iError = MPI_Finalize();

	return 0;
}

