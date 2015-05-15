/*
 * demonstrates the inplace swap of two coloumns using Isend / MPI_Irecv
 * 
 * compile: mpicc -o dtype_inplace dtype_inplace.c
 * run:     mpirun -n 3 ./dtype_inplace
 * 
 *  reminder 2d array refinition and x/y convention
 *   -----> x
 *  a[0][0]  a[0][1] ...
 *  a[1][0]  a[1][1] ...
 *  ...
 * 
 *  flat: a[0][0] a[0][1] a[0][2] ... a[1][0]...
 * 
 * 
 *  that means:
 *  a[y][x] !!!
 * 
 *  and loops like:
 *  for (y in SIZEy) {
 *      for (x in SIZEx) {
 *          bla
 *      }
 *  }
 */

#include "mpi.h" 
#include <stdio.h>
#include <assert.h>
#include <time.h>

// nr of rows
#define SIZEx 6 
// nr of cols
#define SIZEy 4


void printArr (int rank, float a[][SIZEx]){
    int i, j;
    
    printf("r:%d;    ", rank);        
    for ( i=0; i<SIZEy; ++i )
    {
        for (j=0; j<SIZEx; ++j) {
            printf("%4.0f ", a[i][j]);
        }
        printf("\n        ");
    }
    printf("\n");
}


int main(int argc, char *argv[])
{ 
    int numtasks, rank;
    
    assert ( SIZEx>=4 ); // for this to work the array need at least 2 cols
    
    MPI_Status stat; 
    MPI_Datatype columntype;

    MPI_Request send_req_left, send_req_right,
                recv_req_left, recv_req_right;

    MPI_Init(NULL,NULL); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks); 

    MPI_Type_vector(SIZEy, 1, SIZEx, MPI_FLOAT, &columntype); 
    MPI_Type_commit(&columntype); 

    assert ( numtasks >= 2 ); // only works with >=2 procs


    long wait = 100; // in ms
    struct timespec tim, tim2;
    // some artificial delay
    tim.tv_sec = 0;
    tim.tv_nsec = wait * 1000000; //sleep in ns
    
    // small rnd wait bevor printing the final results (should prevent interleaving of printed lines)
    srand(time(NULL)+rank);
    tim2.tv_sec = 0;
    tim2.tv_nsec = (rand()%10) * 1000000; //sleep in ns



    
    // define and fill the array see comment below how it looks like
    float a[SIZEy][SIZEx] =  {{0}};
    int i, j;
    for ( i=0; i<SIZEy; ++i) {
        for ( j=0; j<SIZEx; ++j ) {
            //a[i][j] = (1+j+SIZEy*i) * (-1*(rank*2-1)); //last part assigns the sign
            a[i][j] = rank;
        }
    }
    /*
    if ( rank == 1 ) {
        float a[SIZEy][SIZEx] =  {
            1.0,  2.0,  3.0,  4.0,   
            5.0,  6.0,  7.0,  8.0,  
            9.0, 10.0, 11.0, 12.0, 
            13.0, 14.0, 15.0, 16.0
        };
        target = 1;
    } else {
        float a[SIZEy][SIZEx] =  {
            - 1.0, - 2.0, - 3.0, - 4.0,   
            - 5.0, - 6.0, - 7.0, - 8.0,  
            - 9.0, -10.0, -11.0, -12.0, 
            -13.0, -14.0, -15.0, -16.0
        };
        target = 0;
    }
    */

    // tartet rank of comm.. 
    //int target = 1-rank; // it rank 2: for 0 its 1 and for 1 its 0..
    // this wraps around for demonstration..
    const int LEFT_NODE  = (rank-1<0) ? numtasks-1 : rank-1;
    const int RIGHT_NODE = (rank+1>=numtasks) ? 0 : rank+1;
    
    printf("\n/----------------- at start (%d)\n", rank);
    printArr(rank, a);

    // NOTE: we can use tags!
    // each cell is sending 2 vectors and expecting 2 in an undefined order!

    // define tags for the code to be better readable
    // in the end, better make those defines.. 
    const int TO_LEFT=1,  FROM_RIGHT=1,
              TO_RIGHT=2, FROM_LEFT=2;
    
    printf("r:%d; before the action\n", rank);
    printf("r:%d; ln:%d rn:%d\n", rank, LEFT_NODE, RIGHT_NODE);

/* USING TAGS TO MATCH SEND / RECV    
    MPI_Isend ( &a[0][1],       1, columntype, LEFT_NODE,  TO_LEFT,  MPI_COMM_WORLD, &send_req_left);
    MPI_Isend ( &a[0][SIZEx-2], 1, columntype, RIGHT_NODE, TO_RIGHT, MPI_COMM_WORLD, &send_req_right);
    printf("r:%d; started sending\n", rank);

    MPI_Irecv ( &a[0][0],       SIZEy, columntype, MPI_ANY_SOURCE, FROM_LEFT,  MPI_COMM_WORLD, &recv_req_left);
    MPI_Irecv ( &a[0][SIZEx-1], SIZEy, columntype, MPI_ANY_SOURCE, FROM_RIGHT, MPI_COMM_WORLD, &recv_req_right);
    printf("r:%d; started recv\n", rank);
*/

    // reconsidering: instead of using tags, we can just set the targets properly
    // sorry was already late
    MPI_Isend ( &a[0][1],       1, columntype, LEFT_NODE,  0, MPI_COMM_WORLD, &send_req_left);
    MPI_Isend ( &a[0][SIZEx-2], 1, columntype, RIGHT_NODE, 0, MPI_COMM_WORLD, &send_req_right);
    MPI_Irecv ( &a[0][0],       SIZEy, columntype, LEFT_NODE, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_req_left);
    MPI_Irecv ( &a[0][SIZEx-1], SIZEy, columntype, RIGHT_NODE,  MPI_ANY_TAG,  MPI_COMM_WORLD, &recv_req_right);


    // NOTE: instead of using an MPI_Irecv, followed by a wait
    // one uses BETTER simply a MPI_Recv! (as its the same)


    
/*
    printf("\nr:%d; interm state\n", rank);
    printArr(rank, a);
    
    printf("\nr:%d; going to sleep for %d msec\n", rank, wait);
    
    nanosleep(&tim, NULL); 
    
    printf("\nr:%d; after sleep\n", rank);
    printArr(rank, a);
*/    

    // probably one can only wait for the recv requests !?!?
    // or just use blocking receives?
    MPI_Wait(&send_req_left, &stat);
    MPI_Wait(&send_req_right, &stat);
    MPI_Wait(&recv_req_left, &stat);
    MPI_Wait(&recv_req_right, &stat);

    nanosleep(&tim2, NULL); 
    printf("\nr:%d; final state -------------------\n", rank);
    printArr(rank, a);
    
    
    MPI_Finalize(); 
} 
