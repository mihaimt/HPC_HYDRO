#include "mpi.h" 
#include <stdio.h>
#include <assert.h>
#define SIZE 4

int main(int argc, char *argv[])
{ 
    int numtasks, rank; 
    int i=0;

    float a[SIZE][SIZE] =  
    {
        1.0, 2.0, 3.0, 4.0,   
        5.0, 6.0, 7.0, 8.0,  
        9.0, 10.0, 11.0, 12.0, 
        13.0, 14.0, 15.0, 16.0
    }; 

    float b[SIZE][SIZE] =  
    {
        0.0, 0.0, 0.0, 0.0,   
        0.0, 0.0, 0.0, 0.0,  
        0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0
    }; 

    MPI_Status stat; 
    MPI_Datatype columntype; 

    MPI_Init(NULL,NULL); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks); 

    MPI_Type_vector(SIZE, 1, SIZE, MPI_FLOAT, &columntype); 
    MPI_Type_commit(&columntype); 

    assert ( numtasks == 2 );
    
    if (rank == 0)
    {        
        MPI_Send(&a[0][2], 1, columntype, 1, 0, MPI_COMM_WORLD);

        printf("rank= %d; a=\n", rank);
        
        for ( i=0; i<4; i++ )
        {
            printf("   %4.1f %4.1f %4.1f %4.1f\n", a[i][0],a[i][1],a[i][2],a[i][3]); 
        }
    }
    else
    {
        MPI_Recv(&b[0][3], SIZE, columntype, 0, 0, MPI_COMM_WORLD, &stat);

        printf("rank= %d; b=\n", rank);
        
        for ( i=0; i<4; i++ )
        {
            printf("   %4.1f %4.1f %4.1f %4.1f\n", b[i][0],b[i][1],b[i][2],b[i][3]); 
        }
    }

    MPI_Finalize(); 
} 
