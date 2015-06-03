#ifndef VTKFILE_H_INCLUDED
#define VTKFILE_H_INCLUDED
void vtkfile ( long step, const hydroparam_t H, hydrovar_t * Hv );

void timingfile_init ( hydroparam_t* H );
void timingfile_finish ( hydroparam_t* H );
inline void timingfile_write ( long step, double time, const hydroparam_t H);

void write_stat  ( double elapsed, long nsteps, const hydroparam_t H );

#endif // VTKFILE_H_INCLUDED
