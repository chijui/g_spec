#ifndef SPEC_MOD_H
#define SPEC_MOD_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

#define nosc_switch 0 // Switch between home-built gemv and optimized BLAS gemv routines.

int gen_ExDip2Q(GNREAL ***ExDip2QAr, GNREAL ***Dip2QAr, GNREAL **Ham1QAr, GNREAL **Ham2QAr, int nosc, int n2Q, int tid);

int gen_ham_2Q(GNREAL *Ham1Q, int nosc, GNREAL *Ham2Q, int n2Q, real delta);
int gen_dip_2Q(GNREAL *Dip1Q, GNREAL *Dip2Q, int nosc, int n2Q);

int gen_avg_ham(GNREAL **ham, GNREAL *avg_ham, GNREAL *hann, int nosc, int fr, int whann, int window, int nbuffer);
int gen_avg_dip(GNREAL ***dip1Q, GNREAL ***dip2Q, GNREAL ***avg_dip1Q, GNREAL ***avg_dip2Q, int nosc, int n2Q, int fr, int nAvgd, int do2d, int pert, int pertvec, int whann, int window, int nbuffer);

GNREAL orient(GNREAL **Dip1[3], GNREAL **Dip2[3], int tid, int ndxA, int ndxB, int ndxC, int ndxD, int pol);
// int gen_pert_prop( GNCOMP *U1Q, GNCOMP *U2Q, int nosc, int n2Q, double expfac, double delta);
int gen_perturb_2Q_energies(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evals2Q, int nosc, real delta );
int gen_perturb_2Q_vectors(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evecs2Q, GNREAL *Evals2Q, int nosc, int n2Q, real delta );
int gen_perturb_2Q_matrix(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evecs2Q, GNREAL *Evals2Q, GNREAL *Temp2Q, int nosc, int n2Q, real delta );
// cjfeng 04/24/2017
int gen_popdecay(GNREAL *popdecay, GNREAL tstep, GNREAL TauP, int window, int win2d);
// cjfeng 07/20/2017
int gen_BLAS_gemv_opt(BLAS_gemv_opt *gemv_opt);

// int calc_CorrFunc(GNCOMP *CorrFunc, GNCOMP **U1QMat, GNCOMP **psi_1d[3], GNREAL ***Dip1QMat, GNCOMP **cDip1Q, int nosc, int window, int fr, int nbuffer, double expfac, double wo, BLAS_gemv_opt *gemv_opt, int nthreads);

// int gen_UMat(GNCOMP **UMat, GNREAL *Ham, GNREAL *Evals, int N, int xfr, double expfac, double wo); 

int calc_2dir(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], GNREAL **ExDip2QAr[3], int tid, int nosc, int n2Q, int npts, GNREAL res, GNREAL start, GNREAL stop, GNCOMP *REPH[4], GNCOMP *NREPH[4], int POL[4], int npol, int reph, int nreph);
int calc_2dir_pert(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], int tid, int nosc, int n2Q, int npts, GNREAL res, GNREAL start, GNREAL stop, GNCOMP *REPH[4], GNCOMP *NREPH[4], int POL[4], int npol, int reph, int nreph);

// int initialize_psi(GNCOMP ***psi, GNREAL **Dip1QMat[3], int tinit, int xfr, int nosc, int nelem);
// int initialize_psi2Q(GNCOMP ***psi, GNREAL ***Dip2QMat, GNCOMP **psi_2Q, int nosc, int n2Q, int tau1, int tau2, int nelem, int nthreads);

// int propagate_psi1Q(GNCOMP **psi[3], GNREAL **Dip1QMat[3], GNCOMP **U1QMat, int nosc, int fr, int tau0, int xfr0, int nbuffer, int winsize, BLAS_gemv_opt *gemv_opt, int nthreads);

// int propagate_psi_2Q(GNCOMP **U2QMat, GNCOMP **psi_ca, GNCOMP **psi_cb, GNCOMP *psi2Q, int n2Q, int xfr3, BLAS_gemv_opt *gemv_opt, int nthreads);
// int propagate_psi_2Q_pert(GNCOMP **U1QMat, GNCOMP *U2Qs, GNCOMP **psi_ca, GNCOMP **psi_cb, GNCOMP *psi2Q, int nosc, int n2Q, real delta, double expfac, int xfr3, int nthreads);

// int calc_2dir_nise(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, int *POL, int nosc, int n2Q, int npol, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, int nthreads);

// int pathway_single_side(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val);
// int pathway_double_sides(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val);
// 
// int pathway_REPH_2Q(GNCOMP *psi2Qcb, GNCOMP *psi2Qa_aa, GNCOMP *psi2Qa_ba, int n2Q, int p, GNREAL popfac2Q, GNCOMP *REPH_private);
// int pathway_NREPH_2Q(GNCOMP *psi2Qca, GNCOMP *psi2Qb1_aa, GNCOMP *psi2Qb1_ba, GNCOMP *psi2Qb1_ab, GNCOMP *psi2Qb1_bb, int n2Q, int p, GNREAL popfac2Q, GNCOMP *NREPH_private);

#endif
