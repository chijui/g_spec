#ifndef SPEC_MOD_H
#define SPEC_MOD_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

#define nosc_switch 0 // Switch between home-built gemv and optimized BLAS gemv routines.

int gen_ExDip2Q(GNREAL ***ExDip2QAr, GNREAL ***Dip2QAr, GNREAL **Ham1QAr, GNREAL **Ham2QAr, traj_param *traj_param, int tid);

int gen_ham_2Q(GNREAL *Ham1Q, GNREAL *Ham2Q, traj_param *traj_param, real delta);
int gen_dip_2Q(GNREAL **Dip1Q, GNREAL **Dip2Q, traj_param *traj_param);
int gen_dip_Sp2Q_ind(int *Dip2Qrow, int *Dip2Qcol, traj_param *traj_param);
int gen_dip_Sp2Q(GNREAL **Dip1Q, GNREAL **Dip2Q, traj_param *traj_param);

int gen_avg_ham(GNREAL **ham, GNREAL *avg_ham, GNREAL *hann, int fr, w_param *w_param, traj_param *traj_param);
int gen_avg_dip(GNREAL ***dip1Q, GNREAL ***dip2Q, GNREAL ***avg_dip1Q, GNREAL ***avg_dip2Q, int fr, int nAvgd, traj_param *traj_param, spec_param *spec_param, w_param *w_param);

GNREAL orient(GNREAL **Dip1[3], GNREAL **Dip2[3], int tid, int ndxA, int ndxB, int ndxC, int ndxD, int pol, POL_info *pol_info);
int gen_perturb_2Q_energies(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evals2Q, int nosc, real delta );
int gen_popdecay(GNREAL *popdecay, spec_param *spec_param, w_param *w_param, int if2Q);
int gen_hann(GNREAL *hann, int window, phys_const *Const);
int apply_shift(GNREAL *Ham1Q, shift_info *shift_info, int nosc);

int gen_U1Q(GNREAL **Ham1QMat, GNREAL **Evals1QAr, GNCOMP **U1QMat, spec_param *spec_param, w_param *w_param, traj_param *traj_param, int justread, int readframe, double expfac);
int gen_U2Q(GNREAL **Ham1QMat, GNREAL **Ham2QMat, GNREAL **Evals2QAr, GNCOMP **U2QMat, spec_param *spec_param, w_param *w_param, traj_param *traj_param, int justread, int readframe, real delta, double expfac);
int calc_ftir(GNREAL *ftir, GNREAL **Evals1QAr, GNREAL **ExDip1QAr[3], w_param *w_param, int nosc, int tid);
void eig_avg_Ham2Q(GNREAL **Ham1QAr, GNREAL **Ham2QAr, GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip2QAr[3], int *isuppz2Q, spec_param *spec_param, w_param *w_param, traj_param *traj_param, int tid, real delta);
int calc_2dir(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], GNREAL **ExDip2QAr[3], int tid, traj_param *traj_param, w_param *w_param, GNCOMP *REPH[4], GNCOMP *NREPH[4], POL_info *pol_info, spec_param *spec_param);

#endif
