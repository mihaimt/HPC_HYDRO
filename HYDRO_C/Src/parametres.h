#ifndef PARAMETRES_H_INCLUDED
#define PARAMETRES_H_INCLUDED
#include <mpi.h>


extern unsigned long flops;


typedef struct _hydroparam {
    long prt;

    // time control
    double t, tend;
    long nstep, nstepmax;
    long noutput;
    double dtoutput;

    // dimensions (local)
    long imin, imax, jmin, jmax, nx, ny;
    long nxt, nyt, nxyt;

    // Global domain
    long nxdomain, nydomain;

    long arSz, arVarSz;          // taille des buffers alloues

    /*
       nx, ny: real useful size of the domain
       i/j min, max: the total domain index range
       nxt, nyt: the total domain size (includes 2 extra layers around the domain)
       nxyt: maximum of nxt and nyt to minimize allocations
     */
    // physics
    long nvar;
    double dx;
    double gamma;
    double courant_factor;
    double smallc, smallr;

    // numerical scheme
    long niter_riemann;
    long iorder;
    double slope_type;

    // char scheme[20];
    long scheme;
    long boundary_right, boundary_left, boundary_down, boundary_up;

    /*
     * Variables needed for the MPI.
     */
    int rank;           // rank of the current process
    int n_procs;        // total number of processes in MPI_COMM_WORLD
    int mpi_error;      // error returned from last MPI call (reusable)
    int mpi_is_init;    // (bool) is the MPI library initialized?
    MPI_Datatype mpi_hydro_vector_type; // data vector for transport of ghost cells
/*
    MPI_Request* mpi_req;   // pointer to array of currently hanging MPI requests
    MPI_Status* mpi_status; // pointer to array with statuses of requests
*/

    // variables for openmp
    int n_threads;      /// how many threads we are using (basically the macro N_OMP_THREADS)

} hydroparam_t;


#define HSCHEME_MUSCL 1
#define HSCHEME_PLMDE 2
#define HSCHEME_COLLELA 3


#ifndef IDX3D
#define IDX3D(x, y, z, nx, ny) ( (x) + (nx) * ( (y) + (ny) * (z) ) )
#define IDX2D(x, y, nx)        ( (x) + (nx) * ( (y) ) )
#endif


// Hydrovar holds the whole 2D problem for all variables
typedef struct _hydrovar {
    double* uold;               // nxt, nyt, nvar allocated as (nxt * nyt), nvar
} hydrovar_t;                   // 1:nvar

#ifndef IHv
// #define IHv(i,j,v) ((i) + (j) * H.nxt + (H.nxt * H.nyt) * (v))
#define IHv(i,j,v) ((i) + (H.nxt * (H.nyt * (v)+ (j))))
#define IHvP(i,j,v) ((i) + (j) * H->nxt + (H->nxt * H->nyt) * (v))
#endif /*  */


// work arrays along one direction for all variables
typedef struct _hydrovarwork {
    double *u, *q, *qxm, *qxp, *dq;     // (nxt or nyt), nvar
    double *qleft, *qright, *qgdnv, *flux;      // (nx+1 or ny+1), nvar
} hydrovarwork_t;               // 1:nvar

#ifndef IHvw
#define IHvw(i,v) ((i) + (v) * H.nxyt)
#define IHvwP(i,v) ((i) + (v) * H->nxyt)
#endif /*  */


// works arrays along one direction
typedef struct _hydrowork {
    double *c;                  // nxt or nyt
    double *e;                  // nxt or nyt
    // all others nx+1 or ny+1
    double *rl, *ul, *pl, *cl, *wl;
    double *rr, *ur, *pr, *cr, *wr;
    double *ro, *uo, *po, *co, *wo;
    double *rstar, *ustar, *pstar, *cstar;
    long *sgnm;
    double *spin, *spout, *ushock;
    double *frac, *scr, *delp, *pold;
    long *ind, *ind2;
} hydrowork_t;


// All variables are grouped in structs for clarity sake.
/*
	Warning : no global variables are declared.
	They are passed as arguments.
*/


// useful constants to force double promotion
#ifdef pourfairepropre
static const double zero = ( double ) 0.0;
static const double one = ( double ) 1.0;
static const double two = ( double ) 2.0;
static const double three = ( double ) 3.0;
static const double hundred = ( double ) 100.0;
static const double two3rd = ( double ) 2.0 / ( double ) 3.0;
static const double half = ( double ) 1.0 / ( double ) 2.0;
static const double third = ( double ) 1.0 / ( double ) 3.0;
static const double forth = ( double ) 1.0 / ( double ) 4.0;
static const double sixth = ( double ) 1.0 / ( double ) 6.0;

// conservative variables with C indexing
static const long ID = 1 - 1;
static const long IU = 2 - 1;
static const long IV = 3 - 1;
static const long IP = 4 - 1;

// The current scheme has two extra layers around the domain.
static const long ExtraLayer = 2;
static const long ExtraLayerTot = 2 * 2;

#else /*  */

#define zero   ((double) 0.0)
#define one    ((double) 1.0)
#define two    ((double) 2.0)
#define three  ((double) 3.0)
#define hundred  ((double) 100.0)
#define two3rd ((double) 2.0 / (double) 3.0)
#define half   ((double) 1.0 / (double) 2.0)
#define third  ((double) 1.0 / (double) 3.0)
#define forth  ((double) 1.0 / (double) 4.0)
#define sixth  ((double) 1.0 / (double) 6.0)
#define ID     (0)
#define IU     (1)
#define IV     (2)
#define IP     (3)
#define ExtraLayer    (2)
#define ExtraLayerTot (2 * 2)
#endif /*  */


void process_args ( long argc, char **argv, hydroparam_t * H );


#ifndef MFLOPS
#if defined(FLOPS) && !defined(HMPP)
/*
	Pour le calcul des flops on utilise le bareme CRAY-Los Alamos
	1 flop pour +-*
	1 flop pour ABS MAX MIN SIGN (quelque soit la forme de l'expression)
	4 flops pour / et SQRT
	8 flops pour SIN, COS, EXP, ...
*/
#define MFLOPS(simple, reciproque, maxmin, transcendant) do { flops += ((simple) + (maxmin) + 4 * (reciproque) + 8 * (transcendant)); } while (0)
#else
#define MFLOPS(simple, reciproque, maxmin, transcendant)
// do { } while (0)
#endif
#endif

#endif // PARAMETRES_H_INCLUDED






