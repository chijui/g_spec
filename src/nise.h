#ifndef NISE_H
#define NISE_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

// Propagator
int gen_UMat(GNCOMP **UMat, GNREAL *Ham, GNREAL *Evals, int N, int xfr, double expfac, double wo); 
int gen_pert_prop( GNCOMP *U1Q, GNCOMP *U2Q, int nosc, int n2Q, double expfac, double delta);

// FTIR subroutines
int initialize_psi(GNCOMP ***psi, GNREAL **Dip1QMat[3], int tinit, int xfr, int nosc, int nelem);
int calc_CorrFunc(GNCOMP *CorrFunc, GNCOMP **U1QMat, GNCOMP **psi_1d[3], GNREAL ***Dip1QMat, GNCOMP **cDip1Q, int nosc, int window, int fr, int nbuffer, double expfac, double wo, BLAS_gemv_opt *gemv_opt, int nthreads);

// 2D IR subroutines
int initialize_psi2Q(GNCOMP ***psi, GNREAL ***Dip2QMat, GNCOMP **psi_2Q, int nosc, int n2Q, int tau1, int tau2, int nelem, int nthreads);

int propagate_psi1Q(GNCOMP **psi[3], GNREAL **Dip1QMat[3], GNCOMP **U1QMat, int nosc, int fr, int tau0, int xfr0, int nbuffer, int winsize, BLAS_gemv_opt *gemv_opt, int nthreads);

int propagate_psi_2Q(GNCOMP **U2QMat, GNCOMP **psi_ca, GNCOMP **psi_cb, GNCOMP *psi2Q, int n2Q, int xfr3, BLAS_gemv_opt *gemv_opt, int nthreads);
int propagate_psi_2Q_pert(GNCOMP **U1QMat, GNCOMP *U2Qs, GNCOMP **psi_ca, GNCOMP **psi_cb, GNCOMP *psi2Q, int nosc, int n2Q, real delta, double expfac, int xfr3, int nthreads);

int calc_2dir_nise(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, int *POL, int nosc, int n2Q, int npol, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, int nthreads);

int pathway_single_side(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val);
int pathway_double_sides(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val);

int pathway_REPH_2Q(GNCOMP *psi2Qcb, GNCOMP *psi2Qa_aa, GNCOMP *psi2Qa_ba, int n2Q, int p, GNREAL popfac2Q, GNCOMP *REPH_private);
int pathway_NREPH_2Q(GNCOMP *psi2Qca, GNCOMP *psi2Qb1_aa, GNCOMP *psi2Qb1_ba, GNCOMP *psi2Qb1_ab, GNCOMP *psi2Qb1_bb, int n2Q, int p, GNREAL popfac2Q, GNCOMP *NREPH_private);

#endif
