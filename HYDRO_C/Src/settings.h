/**
 * settings.h
 * 
 * adjust your hydro experience here
 * 
 * @author: Rafael Kueng <rafi.kueng@gmx.ch>
 */




// SETTINGS  (remember 1-on; 0-off; or the defs 4 lines above ;) )
//-----------------------------------------------------------------------------

// compile with mpi
#define USE_MPI YES

// compile with OPENMP
#define USE_OPENMP NO

// compile tests and other debug code
#define DEBUG ON

// run the asserts
#define DO_ASSERTS ON

// measure run time of one iteration locally / globally (min/max)
// GLOBAL results in overhead of two additional MPI_Reduce per step!
#define GET_LOCAL_ITER_TIME  ON
#define GET_GLOBAL_ITER_TIME ON

// write initital | intermediate | final states to vtk file
// attention: the last / final step write will probably write a state
//     with a different time-spacing to the previous one.
//     So dt (or d_steps) between two frames is always the same, expect for
//     the last write! (default: OFF)
#define WRITE_INIT_STATE  ON
#define WRITE_INTER_STATE ON
#define WRITE_FINAL_STATE OFF

// use color output
#define USE_COLOR YES

// --- [ COM METHOD ] ---------------------------------------------------------
// which method should be used to progagate the domain boundaries?
// this is only for testing purposes..
// for real application use _CM_VEKTOR

// definition of possibilities
#define _CM_VEKTOR 0 // use the vector method (one send / recv per domain border)
#define _CM_SINGLE 1 // use many single, blocking send / recv (only works for n_procs ==2 so far)

// select a defined method here
#define COM_METHOD _CM_VEKTOR

// --- [ DEBUG PRINT OUTPUT ] -------------------------------------------------------

// print locations in code. at start of each function, do a print.
// used to be the WHERE macro in the original code
// I (RK) think they are pretty useless and very verbose, turn OFF
#define LOCATION_PRINT OFF

// do traces output (lowest level, print where in the code it is as well,
// usually OFF)
#define TRACE_PRINT ON

// do debug output (usually OFF)
#define DEBUG_PRINT ON

// do run time error output (stuff for the user to see if something is
// massively screwed up and the program will shut down, usually ON)
#define ERROR_PRINT ON

// do run time warn output (stuff for the user to see if something is screwed
// up, usually ON)
#define WARN_PRINT ON

// do run time info output (stuff for the user to see in normal operation,
// usually ON)
#define INFO_PRINT ON

