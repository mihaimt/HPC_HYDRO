/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>

#include "parametres.h"
#include "utils.h"
#include "qleftright.h"
void
// qleftright(const long idim, const hydroparam_t H, hydrovarwork_t * Hvw)
qleftright(const long idim, const long Hnx, const long Hny, const long Hnxyt,
           const long Hnvar,
           double *RESTRICT qxm, double *RESTRICT qxp,
           double *RESTRICT qleft, double *RESTRICT qright)
{
#define IHVW(i,v) ((i) + (v) * Hnxyt)
    long nvar, i;
    long bmax;
    WHERE("qleftright");
    if (idim == 1) {
        bmax = Hnx + 1;
    } else {
        bmax = Hny + 1;
    }
    for (nvar = 0; nvar < Hnvar; nvar++) {
        for (i = 0; i < bmax; i++) {
            qleft[IHVW(i, nvar)] = qxm[IHVW(i + 1, nvar)];
            qright[IHVW(i, nvar)] = qxp[IHVW(i + 2, nvar)];
        }
    }
}

#undef IHVW

// EOF
