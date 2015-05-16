#ifndef QLEFTRIGHT_H_INCLUDED
#define QLEFTRIGHT_H_INCLUDED

void
  qleftright(const long idim, const long Hnx, const long Hny, const long Hnxyt,
             const long Hnvar,
             double *RESTRICT qxm, double *RESTRICT qxp,
             double *RESTRICT qleft, double *RESTRICT qright);

#endif // QLEFTRIGHT_H_INCLUDED
