/**
 * debug.h
 * 
 * one place to set up everything for debug.
 * This might only work for C99 using GNU C compiler!
 * Make sure to include stdio in the file where you use this!
 * 
 * For any prints to the screen, use the functions below instead prints
 * - TRC (trace)  TRC(rank, format, args)
 *      prints a nice trace, with leading file:function:line rank and ending newline
 *      for debugging
 * - DBG (debug)  
 *      same purpose as TRC, but without nice output, just basic print.
 * - ERR (error)  
 *      output if an critical error occured and the program will exit
 * - WRN (warn)
 *      warings about errors, but the program tries to continue
 * - INF (info)
 *      just some nice information..
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

// use color output
#define USE_COLOR YES

// do traces output (lowest level, print where in the code it is as well, usually off)
#define TRACE_PRINT ON

// do debug output (usually off)
#define DEBUG_PRINT ON

// do runtime error output (stuff for the user to see if something is massivly screwed up and the programm will shut down, usually on)
#define ERROR_PRINT ON

// do runtime warn output (stuff for the user to see if something is screwed up, usually on)
#define WARN_PRINT ON

// do runtime info output (stuff for the user to see in normal operation, usually on)
#define INFO_PRINT ON







// ACTUAL DEFINITIONS (no changes here plz)
//-----------------------------------------------------------------------------

// color definition according to:
// http://stackoverflow.com/questions/3585846/color-text-in-terminal-aplications-in-unix


#if USE_COLOR==TRUE
    #define KNRM  "\x1B[0m"
    #define KRED  "\x1B[31m"
    #define KGRN  "\x1B[32m"
    #define KYEL  "\x1B[33m"
    #define KBLU  "\x1B[34m"
    #define KMAG  "\x1B[35m"
    #define KCYN  "\x1B[36m"
    #define KWHT  "\x1B[37m"
    #define RESET "\033[0m"
#else
    #define KNRM  ""
    #define KRED  ""
    #define KGRN  ""
    #define KYEL  ""
    #define KBLU  ""
    #define KMAG  ""
    #define KCYN  ""
    #define KWHT  ""
    #define RESET ""
#endif

// Note to myself: there's a good reason this uses NO ifdef / else clauses! read the SO question!!
// (--> compiler can check the print code)
// the optimizer should remove while{if{0}} never the less

#define TRC(_RANK, _FORMAT, ...) \
    do { if (TRACE_PRINT) fprintf(stderr, KMAG "rank%04i " KMAG "%s:%s():%i > " KCYN _FORMAT RESET "\n", _RANK, __FILE__, __func__, __LINE__, ##__VA_ARGS__); } while (0)

#define DBG(_FORMAT, ...) \
    do { if (DEBUG_PRINT) fprintf(stderr, KCYN _FORMAT RESET, ##__VA_ARGS__); } while (0)

#define ERR(_FORMAT, ...) \
    do { if (ERROR_PRINT) fprintf(stderr, KRED _FORMAT RESET, ##__VA_ARGS__); } while (0)

#define WRN(_FORMAT, ...) \
    do { if (WARN_PRINT) fprintf(stderr, KYEL _FORMAT RESET, ##__VA_ARGS__); } while (0)

#define INF(_FORMAT, ...) \
    do { if (INFO_PRINT) fprintf(stdout, _FORMAT, ##__VA_ARGS__); } while (0)





#endif // DEBUG_H









