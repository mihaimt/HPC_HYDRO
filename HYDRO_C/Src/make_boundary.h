#ifndef MAKE_BOUNDARY_H_INCLUDED
#define MAKE_BOUNDARY_H_INCLUDED

void make_boundary ( long idim, const hydroparam_t H, hydrovar_t * Hv );
void MPI_get_boundary_start ( long idim, const hydroparam_t H, hydrovar_t * Hv );
void MPI_get_boundary_end ( long idim, const hydroparam_t H, hydrovar_t * Hv );
void MPI_get_boundary ( long idim, const hydroparam_t H, hydrovar_t * Hv );

#endif // MAKE_BOUNDARY_H_INCLUDED
