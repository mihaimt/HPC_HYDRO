/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

#include <stdio.h>
// #include <stdlib.h>
#include <malloc.h>
// #include <unistd.h>
#include <math.h>

#include "parametres.h"
#include "compute_deltat.h"
#include "utils.h"
#include "equation_of_state.h"

#define DABS(x) (double) fabs((x))

static void
ComputeQEforRow(const long j, double *uold, double *q, double *e,
                const double Hsmallr, const long Hnx, const long Hnxt,
                const long Hnyt, const long nxyt)
{
    long i;
    double eken;

#define IHV(i, j, v)  ((i) + Hnxt * ((j) + Hnyt * (v)))
#define IHVW(i, v) ((i) + (v) * nxyt)

    for (i = 0; i < Hnx; i++) {
        long idxuID = IHV(i + ExtraLayer, j, ID);
        long idxuIU = IHV(i + ExtraLayer, j, IU);
        long idxuIV = IHV(i + ExtraLayer, j, IV);
        long idxuIP = IHV(i + ExtraLayer, j, IP);
        q[IHVW(i, ID)] = MAX(uold[idxuID], Hsmallr);
        q[IHVW(i, IU)] = uold[idxuIU] / q[IHVW(i, ID)];
        q[IHVW(i, IV)] = uold[idxuIV] / q[IHVW(i, ID)];
        eken = half * (Square(q[IHVW(i, IU)]) + Square(q[IHVW(i, IV)]));
        q[IHVW(i, IP)] = uold[idxuIP] / q[IHVW(i, ID)] - eken;
        e[i] = q[IHVW(i, IP)];
    }
#undef IHV
#undef IHVW
} 

static void
courantOnXY(double *cournox_loc, double *cournoy_loc, const long Hnx, const long nxyt,
            double *c, double *q)
{
    long i;
    double maxValC = zero;
    double tmp1, tmp2;

#define IHVW(i,v) ((i) + (v) * nxyt)
    for (i = 0; i < Hnx; i++) {
      tmp1 = c[i] + DABS(q[IHVW(i, IU)]);
      tmp2 = c[i] + DABS(q[IHVW(i, IV)]);
      *cournox_loc = MAX(*cournox_loc, tmp1);
      *cournoy_loc = MAX(*cournoy_loc, tmp2);
    }


#undef IHVW
}

void
compute_deltat(double *dt, const hydroparam_t H, hydrowork_t * Hw,
               hydrovar_t * Hv, hydrovarwork_t * Hvw)
{
  double cournox, cournoy;

  long j;
  WHERE("compute_deltat");

#define IHVW(i,v) ((i) + (v) * nxyt)

    // compute time step on grid interior
    cournox = zero;
    cournoy = zero;

    Hvw->q = (double *) calloc(H.nvar * H.nxyt, sizeof(double));
    Hw->e = (double *) malloc(H.nx * sizeof(double));
    Hw->c = (double *) malloc(H.nx * sizeof(double));

    for (j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++) {
      ComputeQEforRow(j, Hv->uold, Hvw->q, Hw->e, H.smallr, H.nx, H.nxt, H.nyt, H.nxyt);
      equation_of_state(&Hvw->q[IHvw(0, ID)], Hw->e,
			&Hvw->q[IHvw(0, IP)], Hw->c, 0, H.nx, H.smallc, H.gamma);
      courantOnXY(&cournox, &cournoy, H.nx, H.nxyt, Hw->c, Hvw->q);
    }

    Free(Hvw->q);
    Free(Hw->e);
    Free(Hw->c);

    *dt = H.courant_factor * H.dx / MAX(cournox, MAX(cournoy, H.smallc));

#undef IHVW
}                               // compute_deltat


//EOF
