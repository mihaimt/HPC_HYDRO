#ifndef HYDRO_FUNCS_H_INCLUDED
#define HYDRO_FUNCS_H_INCLUDED

#include "parametres.h"
void MPI_init(hydroparam_t * H, int * argc, char *** argv);
void MPI_finish(hydroparam_t * H);
void MPI_domain_decomp(hydroparam_t * H);
void hydro_init(hydroparam_t * H, hydrovar_t * Hv);
void MPI_hydro_init(hydroparam_t * H, hydrovar_t * Hv);
void hydro_finish(const hydroparam_t H, hydrovar_t * Hv);
void MPI_hydro_finish(hydroparam_t *H, hydrovar_t * Hv);

void allocate_work_space(const hydroparam_t H, hydrowork_t * Hw, hydrovarwork_t * Hvw);
void deallocate_work_space(const hydroparam_t H, hydrowork_t * Hw, hydrovarwork_t * Hvw);

#endif // HYDRO_FUNCS_H_INCLUDED
