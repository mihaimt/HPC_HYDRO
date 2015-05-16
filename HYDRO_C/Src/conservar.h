#ifndef CONSERVAR_H_INCLUDED
#define CONSERVAR_H_INCLUDED


void gatherConservativeVars(const long idim, const long rowcol,
                            double *RESTRICT uold,
                            double *RESTRICT u,
                            const long Himin,
                            const long Himax,
                            const long Hjmin,
                            const long Hjmax,
                            const long Hnvar,
                            const long Hnxt, const long Hnyt, const long Hnxyt);

void updateConservativeVars(const long idim, const long rowcol,
                            const double dtdx,
                            double *RESTRICT uold,
                            double *RESTRICT u,
                            double *RESTRICT flux,
                            const long Himin,
                            const long Himax,
                            const long Hjmin,
                            const long Hjmax,
                            const long Hnvar,
                            const long Hnxt, const long Hnyt, const long Hnxyt);

#endif // CONSERVAR_H_INCLUDED
