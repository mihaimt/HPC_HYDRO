#include "mpi.h" 
#include <stdio.h> 
#define SIZE 4 
int main(int argc; char **argv)
{ 

int numtasks, rank, source=0, dest, tag=1, i; 
float a[SIZE][SIZE] =  
  {1.0, 2.0, 3.0, 4.0,   
   5.0, 6.0, 7.0, 8.0,  
   9.0, 10.0, 11.0, 12.0, 
  13.0, 14.0, 15.0, 16.0}; 

float b[SIZE][SIZE] =  
  {0.0, 0.0, 0.0, 0.0,   
   0.0, 0.0, 0.0, 0.0,  
   0.0, 0.0, 0.0, 0.0, 
   0.0, 0.0, 0.0, 0.0}; 
  
MPI_Status stat; 
MPI_Datatype columntype; 

MPI_Init(NULL,NULL); 
MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
MPI_Comm_size(MPI_COMM_WORLD, &numtasks); 
    
MPI_Type_vector(SIZE, 1, SIZE, MPI_FLOAT, &columntype); 
MPI_Type_commit(&columntype); 
if (numtasks == SIZE) { 
  if (rank == 0) { 
    for (i=0; i < numtasks; i++)  
     	MPI_Send(&a[0][i], 1, columntype, i, tag, MPI_COMM_WORLD); 
	}
	
	MPI_Recv(&b[0][i], SIZE, columntype, source, tag, MPI_COMM_WORLD, &stat); 

	for (i=0;i< numtasks; i++)
		printf("rank= %d  b= %3.1f %3.1f %3.1f %3.1f\n", 
        rank,b[0][i],b[1][i],b[2][i],b[3][i]); 
} 
else 
  printf("Must specify %d processors. Terminating.\n",SIZE); 
 MPI_Finalize(); 
} 
