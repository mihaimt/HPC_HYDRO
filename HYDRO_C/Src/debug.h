/**
 * debug.h
 * 
 * one place to set up everything for debug.
 * This might only work for C99 using GNU C compiler!
 * Make sure to include stdio in the file where you use this!
 * 
 * For any prints to the screen, use the functions below instead prints
 * - TRC (trace)
 * - DBG (debug)
 * - WRN (warn)
 * - INF (info)
 * 
 * - CHECK(test, text if false)
 * 
 * And adjust the settings in the section before compiling to decive, how
 * verbose the prog will be.
 * 
 * based on this SO question:
 * http://stackoverflow.com/questions/1644868/c-define-macro-for-debug-printing
 * 
 * @author: Rafael Kueng <rafi.kueng@gmx.ch>
 */

#ifndef __DEBUG_H__
#define __DEBUG_H__

// just some handy preliminary defs for the confused ones
#define ON 1
#define OFF 0
#define YES 1
#define NO 0
#define TRUE 1
#define FALSE 0




// SETTINGS  (remember 1-on; 0-off; or the defs 4 lines above ;) )
//-----------------------------------------------------------------------------

// compile with mpi
#define USE_MPI YES

// compile with OPENMP
#define USE_OPENMP NO

// compile tests and other debug code
#define DEBUG ON

// do traces output (lowest level, print where in the code it is as well, usually off)
#define TRACE_PRINT ON

// do debug output (usually off)
#define DEBUG_PRINT ON

// do runtime warn output (stuff for the user to see if something is screwed up, usually on)
#define WARN_PRINT ON

// do runtime info output (stuff for the user to see in normal operation, usually on)
#define INFO_PRINT ON







// ACTUAL DEFINITIONS (no changes here plz)
//-----------------------------------------------------------------------------

// Note to myself: there's a good reason this uses NO ifdef / else clauses! read the SO question!!
// (--> compiler can check the print code)
// the optimiser should remove while{if{0}} never the less

#define TRC(FORMAT, ...) \
    do { if (TRACE_PRINT) fprintf(stderr, "! %s() in %s, line %i: " FORMAT "\n", __func__, __FILE__, __LINE__, __VA_ARGS__); } while (0)

#define DBG(...) \
    do { if (DEBUG_PRINT) fprintf(stderr, __VA_ARGS__); } while (0)

#define WRN(...) \
    do { if (WARN_PRINT) fprintf(stderr, __VA_ARGS__); } while (0)

#define INF(...) \
    do { if (INFO_PRINT) fprintf(stdout, __VA_ARGS__); } while (0)





#endif // DEBUG_H









