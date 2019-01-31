#ifndef NISE_H
#define NISE_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

#ifndef TROTTER_H
#include "trotter.h"
#endif

// Propagator
int gen_UMat(GNCOMP **UMat, GNREAL *Ham, GNREAL *Evals, int N, int xfr, double expfac); 
int gen_pert_prop( GNCOMP *U1Q, GNCOMP *U2Q, traj_param *traj_param, double expfac, double delta);

// FTIR subroutines
int initialize_psi1Q(GNCOMP **psi, GNREAL **Dip1QMat, int nosc);
int calc_CorrFunc(GNCOMP *CorrFunc, GNCOMP **U1QMat, GNCOMP **psi_1d[3], GNREAL ***Dip1QMat, traj_param *traj_param, int window, int fr, double expfac, BLAS_gemv_opt *gemv_opt, int nthreads);

// 2D IR subroutines
int initialize_psi2Q_Sp(GNCOMP **psi, GNREAL **Dip2QSpMat, int *Dip2Qrow, int *Dip2Qcol, GNCOMP **psi_2Q, traj_param *traj_param, int nthreads);
int initialize_psi2Q(GNCOMP **psi, GNREAL **Dip2QMat, GNCOMP **psi_2Q, traj_param *traj_param, int nthreads);

int propagate_psi1Q(GNCOMP **psi[3], GNREAL **Dip1QMat[3], GNCOMP **U1QMat, traj_param *traj_param, int fr, int tau0, int xfr0, int winsize, BLAS_gemv_opt *gemv_opt, int nthreads);
int propagate_psi2Q(GNCOMP **U2QMat, GNCOMP **psi, GNCOMP *psi2Q, int n2Q, int xfr3, BLAS_gemv_opt *gemv_opt, int nthreads);

int calc_2dir_nise_Sp_old(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param, int nthreads);
int calc_2dir_nise_Sp(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param, int nthreads);
int calc_2dir_nise(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param, int nthreads);

int pathway_single_side(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val);
int pathway_double_sides(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val);

int pathway_REPH_2Q(GNCOMP *psi2Qcb, GNCOMP *psi2Qa_aa, GNCOMP *psi2Qa_ba, int n2Q, int p, GNREAL popfac2Q, GNCOMP *REPH_private);
int pathway_NREPH_2Q(GNCOMP *psi2Qca, GNCOMP *psi2Qb1_aa, GNCOMP *psi2Qb1_ba, GNCOMP *psi2Qb1_ab, GNCOMP *psi2Qb1_bb, int n2Q, int p, GNREAL popfac2Q, GNCOMP *NREPH_private);
int pathway_2Q(GNCOMP *psi2Qc, GNCOMP *psi2Q, int n2Q, GNREAL popfac2Q, GNCOMP *val);

#endif
