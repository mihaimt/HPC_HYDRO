#ifndef RIEMANN_H_INCLUDED
#define RIEMANN_H_INCLUDED


void riemann(double *RESTRICT qleft, double *RESTRICT qright,
             double *RESTRICT qgdnv, double *RESTRICT rl,
             double *RESTRICT ul, double *RESTRICT pl,
             double *RESTRICT cl, double *RESTRICT wl,
             double *RESTRICT rr, double *RESTRICT ur,
             double *RESTRICT pr, double *RESTRICT cr,
             double *RESTRICT wr, double *RESTRICT ro,
             double *RESTRICT uo, double *RESTRICT po,
             double *RESTRICT co, double *RESTRICT wo,
             double *RESTRICT rstar, double *RESTRICT ustar,
             double *RESTRICT pstar, double *RESTRICT cstar,
             long *RESTRICT sgnm, double *RESTRICT spin,
             double *RESTRICT spout, double *RESTRICT ushock,
             double *RESTRICT frac, double *RESTRICT scr,
             double *RESTRICT delp, double *RESTRICT pold,
             long *RESTRICT ind, long *RESTRICT ind2,
             const long narray,
             const double Hsmallr,
             const double Hsmallc,
             const double Hgamma, const long Hniter_riemann, const long Hnvar,
             const long Hnxyt);

#endif // RIEMANN_H_INCLUDED
