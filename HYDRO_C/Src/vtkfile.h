#ifndef VTKFILE_H_INCLUDED
#define VTKFILE_H_INCLUDED


void vtkfile ( long step, const hydroparam_t H, hydrovar_t * Hv );

void timingfile_init ( hydroparam_t* H );
void timingfile_finish ( hydroparam_t* H );
inline void timingfile_write ( const hydroparam_t H, const TIMINGS T );

void write_stat  ( const hydroparam_t H, const TIMINGS T, int nstates, double tmax, double tmin );

#endif // VTKFILE_H_INCLUDED
