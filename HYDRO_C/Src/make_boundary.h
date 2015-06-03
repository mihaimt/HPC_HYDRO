#ifndef MAKE_BOUNDARY_H_INCLUDED
#define MAKE_BOUNDARY_H_INCLUDED

// Define values for MPI_get_boundary_start
#define LEFT_GHOST_CELLS 0
#define RIGHT_GHOST_CELLS 1
void make_boundary ( long idim, const hydroparam_t H, hydrovar_t * Hv );
void MPI_get_boundary_start ( long idim, const hydroparam_t H, hydrovar_t * Hv, MPI_Request *MPI_req );
void MPI_get_boundary_end ( long idim, const hydroparam_t H, hydrovar_t * Hv, MPI_Request *MPI_req );

//void MPI_get_boundary_sendrecv( long idim, const hydroparam_t H, hydrovar_t * Hv );
//void MPI_get_boundary_simple ( long idim, const hydroparam_t H, hydrovar_t * Hv );

void MPI_make_boundary ( long idim, const hydroparam_t H, hydrovar_t * Hv );

#endif // MAKE_BOUNDARY_H_INCLUDED
