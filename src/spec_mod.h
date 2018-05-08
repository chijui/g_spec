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
int gen_perturb_2Q_energies(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evals2Q, int nosc, real delta );
int gen_perturb_2Q_vectors(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evecs2Q, GNREAL *Evals2Q, int nosc, int n2Q, real delta );
int gen_perturb_2Q_matrix(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evecs2Q, GNREAL *Evals2Q, GNREAL *Temp2Q, int nosc, int n2Q, real delta );
// cjfeng 04/24/2017
int gen_popdecay(GNREAL *popdecay, GNREAL tstep, GNREAL TauP, int window, int win2d);
// cjfeng 07/20/2017
int gen_BLAS_gemv_opt(BLAS_gemv_opt *gemv_opt);

int calc_2dir(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], GNREAL **ExDip2QAr[3], int tid, int nosc, int n2Q, int npts, GNREAL res, GNREAL start, GNREAL stop, GNCOMP *REPH[4], GNCOMP *NREPH[4], int POL[4], int npol, int reph, int nreph);
int calc_2dir_pert(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], int tid, int nosc, int n2Q, int npts, GNREAL res, GNREAL start, GNREAL stop, GNCOMP *REPH[4], GNCOMP *NREPH[4], int POL[4], int npol, int reph, int nreph);

#endif
