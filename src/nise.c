#include "nise.h"

// Propagator
int gen_UMat(GNCOMP **UMat, GNREAL *Ham, GNREAL *Evals, int N, int xfr, double expfac) {
  int i, j;
  for(i=0; i<N; i++) {
    int ni = i*N;
    GNREAL cre, cim;
    for(j=0; j<N; j++) {
      int k, nj = j*N;
      cre = 0.0;
      cim = 0.0;
      for(k=0; k<N; k++) {
        GNREAL val = expfac * Evals[k];
        GNREAL str = Ham[ni+k] * Ham[nj+k];
        cre += str * cos( val );
        cim += str * sin( val );
      }
      UMat[xfr][ni+j] = cre + I*cim;
    }
  }
  return 0;
}

// cjfeng 01/04/2019
// Added traj_param
int gen_pert_prop( GNCOMP *U1Q, GNCOMP *U2Q, traj_param *traj_param, double expfac, double delta) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int i,j,ip,jp;
  GNCOMP fac;
  double sqrt2inv = 0.707106781186547;
  // negative sign on sin() term since delta is provided as a positive argument.
  GNREAL val = expfac*delta/2;
  GNCOMP facnew = sqrt2inv * ( cos(val) - I*sin(val) );
  for(i=0; i<nosc; i++) {
    int ni = i*nosc;
    for(j=0; j<=i; j++) {
      int nj = j*nosc;
      int NDX2Qij = NDX2Q(i,j)*n2Q;
      for(ip=0; ip<nosc; ip++) {
        GNCOMP U1Qniip = U1Q[ni+ip];
        GNCOMP U1Qnjip = U1Q[nj+ip];
        for(jp=0; jp<=ip; jp++) {
          fac = 1.0;
          if(i==j) fac *= facnew;
          if(ip==jp) fac *= facnew;
          U2Q[NDX2Qij+NDX2Q(ip,jp)] = fac*( U1Qniip*U1Q[nj+jp] + U1Qnjip*U1Q[ni+jp] );
        }
      }
    }
  }
  return 1;
}

// FTIR subroutines
int initialize_psi1Q(GNCOMP **psi, GNREAL **Dip1QMat, int nosc) {
  int i, j;
  for(i=0; i<3; i++) for(j=0; j<nosc; j++) psi[i][j] = Dip1QMat[i][j];
  return 0;
}

// cjfeng 01/03/2019
// Using cval2 instead of cDip1Q
// cjfeng 01/04/2019
// Added traj_param
int calc_CorrFunc(GNCOMP *CorrFunc, GNCOMP **U1QMat, GNCOMP ***psi_1d, GNREAL ***Dip1QMat, traj_param *traj_param, int window, int fr, double expfac, BLAS_gemv_opt *gemv_opt, int nthreads) {
  int nosc = traj_param->nosc;
  int nbuffer = traj_param->nbuffer;
  int noscsq = nosc * nosc;
  int i, j, n;
  int xfr;
  GNCOMP cval;
  initialize_psi1Q(psi_1d[0], Dip1QMat[fr], nosc);
  
  for(n=0; n<window; n++) {
    xfr = (fr+n) % nbuffer;
    cval = 0.0;
    #if OMP_PARALLEL
    #pragma omp parallel if(nosc >= NBCOMP && nthreads>1) shared(nosc, Dip1QMat, xfr, fr) private(i,j)
    #pragma omp for schedule(guided) collapse(2) reduction(+:cval)
    #endif
    for(i=0; i<3; i++) {
      for(j=0; j<nosc; j++) {
        GNCOMP cval2 = psi_1d[0][i][j];
        cval += Dip1QMat[xfr][i][j] * cval2;
      }
    }
    CorrFunc[n] += cval;
  
    // Propagation  
    #if BLAS_subroutines
    if(nosc > nosc_switch) for(i=0; i<3; i++) cgemv_(&(gemv_opt->ta), &nosc, &nosc, &(gemv_opt->alpha), U1QMat[xfr], &nosc, psi_1d[0][i], &(gemv_opt->inc), &(gemv_opt->beta), psi_1d[1][i], &(gemv_opt->inc));  
    else for(i=0; i<3; i++) mvmult_comp(U1QMat[xfr], psi_1d[0][i], psi_1d[1][i], nosc, nthreads);
    #else
    for(i=0; i<3; i++) mvmult_comp(U1QMat[xfr], psi_1d[0][i], psi_1d[1][i], nosc, nthreads);
    #endif
    for(i=0; i<3; i++) copy_gncomp(psi_1d[1][i], psi_1d[0][i], nosc, 1);
  }
  return 0;
}

// 2D IR subroutines
// cjfeng 12/19/2018
// Using sparse matrix
// cjfeng 01/04/2019
// Added traj_param
int initialize_psi2Q_Sp(GNCOMP **psi, GNREAL **Dip2QSpMat, int *Dip2Qrow, int *Dip2Qcol, GNCOMP **psi_2Q, traj_param *traj_param, int nthreads) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int i, j;
  for(j=0; j<3; j++) for(i=0; i<3; i++) mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2QSpMat[j], psi[i], psi_2Q[i*3+j], nosc, n2Q, nthreads);
  return 0;
}

int initialize_psi2Q(GNCOMP **psi, GNREAL **Dip2QMat, GNCOMP **psi_2Q, traj_param *traj_param, int nthreads) {
  int i, j;
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  for(j=0; j<3; j++) for(i=0; i<3; i++) mvmult_comp_trans_x(Dip2QMat[j], psi[i], psi_2Q[i*3+j], nosc, n2Q, nthreads);
  return 0;
}

// cjfeng 01/04/2019
// Added traj_param
int propagate_psi1Q(GNCOMP ***psi, GNREAL ***Dip1QMat, GNCOMP **U1QMat, traj_param *traj_param, int fr, int tau0, int xfr0, int winsize, BLAS_gemv_opt *gemv_opt, int nthreads) {
  int i, tau;
  int nosc = traj_param->nosc;
  int nbuffer = traj_param->nbuffer;
  // Initial state
  initialize_psi1Q(psi[tau0], Dip1QMat[xfr0], nosc);
  #if BLAS_subroutines
  if(nosc > nosc_switch) {
    for(tau=tau0; tau<tau0+winsize-1; tau++) for(i=0; i<3; i++) cgemv_(&(gemv_opt->ta), &nosc, &nosc, &(gemv_opt->alpha), U1QMat[(fr+tau)%nbuffer], &nosc, psi[tau][i], &(gemv_opt->inc), &(gemv_opt->beta), psi[tau+1][i], &(gemv_opt->inc));
  }
  else {
    for(tau=tau0; tau<tau0+winsize-1; tau++) for(i=0; i<3; i++) mvmult_comp(U1QMat[(fr+tau)%nbuffer], psi[tau][i], psi[tau+1][i], nosc, nthreads );
  }
  #else
  for(tau=tau0; tau<tau0+winsize-1; tau++) for(i=0; i<3; i++) mvmult_comp(U1QMat[(fr+tau)%nbuffer], psi[tau][i], psi[tau+1][i], nosc, nthreads );
  #endif
  return 0;
}

int propagate_psi2Q(GNCOMP **U2QMat, GNCOMP **psi, GNCOMP *psi2Q, int n2Q, int xfr3, BLAS_gemv_opt *gemv_opt, int nthreads) {
  int i, j;
  for(i=0; i<9; i++) {
    // First initialize psi and propagate it.
    for(j=0; j<n2Q; j++) psi2Q[j] = psi[i][j];
    #if BLAS_subroutines
    if (n2Q > nosc_switch) cgemv_(&(gemv_opt->ta), &n2Q, &n2Q, &(gemv_opt->alpha), U2QMat[xfr3], &n2Q, psi2Q, &(gemv_opt->inc), &(gemv_opt->beta), psi[i], &(gemv_opt->inc));
    else mvmult_comp(U2QMat[xfr3], psi2Q, psi[i], n2Q, nthreads);
    #else
    mvmult_comp(U2QMat[xfr3], psi2Q, psi[i], n2Q, nthreads);
    #endif
  }
  return 0;
}

// cjfeng 03/23/2017
// This routine aims for computing the nonlinear response function in nise algorithm.
// cjfeng 12/07/2018
// Using pol_info, and spec_param. Added reph and nreph flags
// cjfeng 12/19/2018
// Using sparse matrix
// cjfeng 01/04/2019
// Added traj_param
int calc_2dir_nise_Sp_old(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param, int nthreads) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int tau1window = tau1 * window;
  int P, p, i, j, k, l, n, a, b;
  int npol = pol_info->npol;
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;

  // cjfeng 07/02/2017
  GNCOMP *psi2Qa_aa, *psi2Qa_ba, *psi2Qb1_aa, *psi2Qb1_ab, *psi2Qb1_ba, *psi2Qb1_bb; 
  set_gncomp_1d_array(&psi2Qa_aa, n2Q, 1);
  set_gncomp_1d_array(&psi2Qa_ba, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_aa, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_ab, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_ba, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_bb, n2Q, 1);
  // cjfeng 07/03/2017
  // Change the order of loops
  // cjfeng 07/02/2017
  // Further optimzations to reduce number of 2Q computations.
  // Recognizing that i is always equal to a, and 2Q transitions are only related to i, and j,
  // we can pull out the i terms first and reduce operations.
  // p sums over the forms iiii, iijj, ijji, and ijij.
  for(a=0; a<3; a++) {
    GNREAL *Dip2Qa = Dip2QSpMat[xfr3][a];
    mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qa, psi_a[tau1+tau2+tau3][a], psi2Qa_aa, nosc, n2Q, nthreads);
    mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qa, psi_b1[tau1+tau2+tau3][a], psi2Qb1_aa, nosc, n2Q, nthreads);
    for(b=0; b<3; b++) {
      // cjfeng 07/02/2017
      GNREAL *Dip2Qb = Dip2QSpMat[xfr3][b];
      if ( b != a ) {
        mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qb, psi_a[tau1+tau2+tau3][a], psi2Qa_ba, nosc, n2Q, nthreads);
        mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qa, psi_b1[tau1+tau2+tau3][b], psi2Qb1_ab, nosc, n2Q, nthreads);
        mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qb, psi_b1[tau1+tau2+tau3][a], psi2Qb1_ba, nosc, n2Q, nthreads);
        mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qb, psi_b1[tau1+tau2+tau3][b], psi2Qb1_bb, nosc, n2Q, nthreads);
      }
      else {
        for(n=0; n<n2Q; n++) psi2Qa_ba[n] = psi2Qa_aa[n];
        for(n=0; n<n2Q; n++) psi2Qb1_ab[n] = psi2Qb1_aa[n];
        for(n=0; n<n2Q; n++) psi2Qb1_ba[n] = psi2Qb1_aa[n];
        for(n=0; n<n2Q; n++) psi2Qb1_bb[n] = psi2Qb1_aa[n];
      }
      for(p=0; p<4; p++) {
        if(p==0) { i=a; j=a; k=a; l=a; }
        else if(p==1) { i=a; j=a; k=b; l=b; }
        else if(p==2) { i=a; j=b; k=b; l=a; }
        else { i=a; j=b; k=a; l=b; } // p==3
        // For p==0, we only add the signal when a==b.
        if( (p!=0) ) {
          if ( a != b) {
            GNCOMP REPH_private = 0.0;
            GNCOMP NREPH_private = 0.0;
            GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
            Dipj = Dip1QMat[xfr1][j];
            Dipk = Dip1QMat[xfr2][k];
            Dipl = Dip1QMat[xfr3][l];
            Dip2Q = Dip2QSpMat[xfr3][l];
            int ik = i*3 + k;
            int jk = j*3 + k;
            // cjfeng 07/03/2017
            // Replaced the original code by sub-routines
            if (reph) {
              pathway_double_sides(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &REPH_private);
              pathway_double_sides(Dipl, Dipk, psi_b1[tau1+tau2+tau3][j], psi_a[tau1+tau2][i], nosc, popfac1Q, &REPH_private);
              pathway_REPH_2Q(psi_cb[jk], psi2Qa_aa, psi2Qa_ba, n2Q, p, popfac2Q, &REPH_private);
            }
            if (nreph) {
              pathway_single_side(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &NREPH_private);
              pathway_double_sides(Dipl, Dipk, psi_a[tau1+tau2+tau3][i], psi_b1[tau1+tau2][j], nosc, popfac1Q, &NREPH_private);
              pathway_NREPH_2Q(psi_ca[ik], psi2Qb1_aa, psi2Qb1_ba, psi2Qb1_ab, psi2Qb1_bb, n2Q, p, popfac2Q, &NREPH_private);
            }
            // cjfeng 07/03/2017
            // Computing different polarized spectrum
            for(P=0; P<npol; P++) {
              int nP = pol_info->POL[P];
              GNREAL orient_fac = pol_info->M_ijkl_IJKL[P][p];
              if (reph) REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
              if (nreph) NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
            }
          }
        }
        else if( (a==b) ) {
          GNCOMP REPH_private = 0.0;
          GNCOMP NREPH_private = 0.0;
          GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
          Dipj = Dip1QMat[xfr1][j];
          Dipk = Dip1QMat[xfr2][k];
          Dipl = Dip1QMat[xfr3][l];
          Dip2Q = Dip2QSpMat[xfr3][l];
          int ik = i*3 + k;
          int jk = j*3 + k;
          // cjfeng 07/03/2017
          // Replaced the original code by sub-routines
          if (reph) {
            pathway_double_sides(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &REPH_private);
            pathway_double_sides(Dipl, Dipk, psi_b1[tau1+tau2+tau3][j], psi_a[tau1+tau2][i], nosc, popfac1Q, &REPH_private);
            pathway_REPH_2Q(psi_cb[jk], psi2Qa_aa, psi2Qa_ba, n2Q, p, popfac2Q, &REPH_private);
          }
          if (nreph) {
            pathway_single_side(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &NREPH_private);
            pathway_double_sides(Dipl, Dipk, psi_a[tau1+tau2+tau3][i], psi_b1[tau1+tau2][j], nosc, popfac1Q, &NREPH_private);
            pathway_NREPH_2Q(psi_ca[ik], psi2Qb1_aa, psi2Qb1_ba, psi2Qb1_ab, psi2Qb1_bb, n2Q, p, popfac2Q, &NREPH_private);
          }
          // cjfeng 07/03/2017
          // Computing different polarized spectrum
          for(P=0; P<npol; P++) {
            int nP = pol_info->POL[P];
            GNREAL orient_fac = pol_info->M_ijkl_IJKL[P][p];
            if (reph) REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
            if (nreph) NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
          }
        }
      }
    }
  }
  unset_gncomp_1d_array(psi2Qa_aa);
  unset_gncomp_1d_array(psi2Qa_ba);
  unset_gncomp_1d_array(psi2Qb1_aa);
  unset_gncomp_1d_array(psi2Qb1_ab);
  unset_gncomp_1d_array(psi2Qb1_ba);
  unset_gncomp_1d_array(psi2Qb1_bb);
  return 0;
}

void calc_R3_iiii_Sp(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param, int nthreads) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int tau1window = tau1 * window;
  int P, i, j, k, l, n, a, b;
  int npol = pol_info->npol;
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;

  GNCOMP *psi2Qa_aa, *psi2Qb1_aa; 
  set_gncomp_1d_array(&psi2Qa_aa, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_aa, n2Q, 1);
  GNCOMP REPH_private = 0.0;
  GNCOMP NREPH_private = 0.0;
  for(a=0; a<3; a++) {
    GNREAL *Dip2Qa = Dip2QSpMat[xfr3][a];
    mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qa, psi_a[tau1+tau2+tau3][a], psi2Qa_aa, nosc, n2Q, nthreads);
    mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qa, psi_b1[tau1+tau2+tau3][a], psi2Qb1_aa, nosc, n2Q, nthreads);
    i=a; j=a; k=a; l=a;
    GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
    Dipj = Dip1QMat[xfr1][j];
    Dipk = Dip1QMat[xfr2][k];
    Dipl = Dip1QMat[xfr3][l];
    Dip2Q = Dip2QSpMat[xfr3][l];
    int ik = i*3 + k;
    int jk = j*3 + k;
    if (reph) {
      pathway_double_sides(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &REPH_private);
      pathway_double_sides(Dipl, Dipk, psi_b1[tau1+tau2+tau3][j], psi_a[tau1+tau2][i], nosc, popfac1Q, &REPH_private);
      pathway_2Q(psi_cb[jk], psi2Qa_aa, n2Q, popfac2Q, &REPH_private);
    }
    if (nreph) {
      pathway_single_side(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &NREPH_private);
      pathway_double_sides(Dipl, Dipk, psi_a[tau1+tau2+tau3][i], psi_b1[tau1+tau2][j], nosc, popfac1Q, &NREPH_private);
      pathway_2Q(psi_ca[ik], psi2Qb1_aa, n2Q, popfac2Q, &NREPH_private);
    }
  }
  // cjfeng 07/03/2017
  // Computing different polarized spectrum
  for(P=0; P<npol; P++) {
    int nP = pol_info->POL[P];
    GNREAL orient_fac = pol_info->M_ijkl_IJKL[P][0];
    if (reph) REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
    if (nreph) NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
  }
  unset_gncomp_1d_array(psi2Qa_aa);
  unset_gncomp_1d_array(psi2Qb1_aa);
}

void calc_R3_iijj_Sp(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param, int nthreads) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int tau1window = tau1 * window;
  int P, i, j, k, l, n, a, b;
  int npol = pol_info->npol;
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;

  GNCOMP *psi2Qa_ba, *psi2Qb1_ba; 
  set_gncomp_1d_array(&psi2Qa_ba, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_ba, n2Q, 1);
  GNCOMP REPH_private = 0.0;
  GNCOMP NREPH_private = 0.0;
  for(a=0; a<3; a++) {
    for(b=0; b<3; b++) {
      if ( b != a ) {
        GNREAL *Dip2Qb = Dip2QSpMat[xfr3][b];
        mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qb, psi_a[tau1+tau2+tau3][a], psi2Qa_ba, nosc, n2Q, nthreads);
        mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qb, psi_b1[tau1+tau2+tau3][a], psi2Qb1_ba, nosc, n2Q, nthreads);
        i=a; j=a; k=b; l=b;
        GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
        Dipj = Dip1QMat[xfr1][j];
        Dipk = Dip1QMat[xfr2][k];
        Dipl = Dip1QMat[xfr3][l];
        Dip2Q = Dip2QSpMat[xfr3][l];
        int ik = i*3 + k;
        int jk = j*3 + k;
        if (reph) {
          pathway_double_sides(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &REPH_private);
          pathway_double_sides(Dipl, Dipk, psi_b1[tau1+tau2+tau3][j], psi_a[tau1+tau2][i], nosc, popfac1Q, &REPH_private);
          pathway_2Q(psi_cb[jk], psi2Qa_ba, n2Q, popfac2Q, &REPH_private);
        }
        if (nreph) {
          pathway_single_side(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &NREPH_private);
          pathway_double_sides(Dipl, Dipk, psi_a[tau1+tau2+tau3][i], psi_b1[tau1+tau2][j], nosc, popfac1Q, &NREPH_private);
          pathway_2Q(psi_ca[ik], psi2Qb1_ba, n2Q, popfac2Q, &NREPH_private);
        }
      }
    }
  }
  // cjfeng 07/03/2017
  // Computing different polarized spectrum
  for(P=0; P<npol; P++) {
    int nP = pol_info->POL[P];
    GNREAL orient_fac = pol_info->M_ijkl_IJKL[P][1];
    if (reph) REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
    if (nreph) NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
  }
  unset_gncomp_1d_array(psi2Qa_ba);
  unset_gncomp_1d_array(psi2Qb1_ba);
}

void calc_R3_ijji_Sp(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param, int nthreads) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int tau1window = tau1 * window;
  int P, i, j, k, l, n, a, b;
  int npol = pol_info->npol;
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;

  // cjfeng 07/02/2017
  GNCOMP *psi2Qa_aa, *psi2Qb1_ab; 
  set_gncomp_1d_array(&psi2Qa_aa, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_ab, n2Q, 1);
  GNCOMP REPH_private = 0.0;
  GNCOMP NREPH_private = 0.0;
  for(a=0; a<3; a++) {
    GNREAL *Dip2Qa = Dip2QSpMat[xfr3][a];
    mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qa, psi_a[tau1+tau2+tau3][a], psi2Qa_aa, nosc, n2Q, nthreads);
    for(b=0; b<3; b++) {
      if ( b != a ) {
        GNREAL *Dip2Qb = Dip2QSpMat[xfr3][b];
        mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qa, psi_b1[tau1+tau2+tau3][b], psi2Qb1_ab, nosc, n2Q, nthreads);
        i=a; j=b; k=b; l=a;
        GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
        Dipj = Dip1QMat[xfr1][j];
        Dipk = Dip1QMat[xfr2][k];
        Dipl = Dip1QMat[xfr3][l];
        Dip2Q = Dip2QSpMat[xfr3][l];
        int ik = i*3 + k;
        int jk = j*3 + k;
        if (reph) {
          pathway_double_sides(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &REPH_private);
          pathway_double_sides(Dipl, Dipk, psi_b1[tau1+tau2+tau3][j], psi_a[tau1+tau2][i], nosc, popfac1Q, &REPH_private);
          pathway_2Q(psi_cb[jk], psi2Qa_aa, n2Q, popfac2Q, &REPH_private);
        }
        if (nreph) {
          pathway_single_side(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &NREPH_private);
          pathway_double_sides(Dipl, Dipk, psi_a[tau1+tau2+tau3][i], psi_b1[tau1+tau2][j], nosc, popfac1Q, &NREPH_private);
          pathway_2Q(psi_ca[ik], psi2Qb1_ab, n2Q, popfac2Q, &NREPH_private);
        }
      }
    }
  }
  // cjfeng 07/03/2017
  // Computing different polarized spectrum
  for(P=0; P<npol; P++) {
    int nP = pol_info->POL[P];
    GNREAL orient_fac = pol_info->M_ijkl_IJKL[P][2];
    if (reph) REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
    if (nreph) NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
  }
  unset_gncomp_1d_array(psi2Qa_aa);
  unset_gncomp_1d_array(psi2Qb1_ab);
}

void calc_R3_ijij_Sp(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param, int nthreads) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int tau1window = tau1 * window;
  int P, i, j, k, l, n, a, b;
  int npol = pol_info->npol;
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;

  GNCOMP *psi2Qa_ba, *psi2Qb1_bb; 
  set_gncomp_1d_array(&psi2Qa_ba, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_bb, n2Q, 1);
  GNCOMP REPH_private = 0.0;
  GNCOMP NREPH_private = 0.0;
  for(a=0; a<3; a++) {
    for(b=0; b<3; b++) {
      if ( b != a ) {
        GNREAL *Dip2Qb = Dip2QSpMat[xfr3][b];
        mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qb, psi_a[tau1+tau2+tau3][a], psi2Qa_ba, nosc, n2Q, nthreads);
        mvmult_comp_trans_Sp(Dip2Qrow, Dip2Qcol, Dip2Qb, psi_b1[tau1+tau2+tau3][b], psi2Qb1_bb, nosc, n2Q, nthreads);
        i=a; j=b; k=a; l=b; 
        GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
        Dipj = Dip1QMat[xfr1][j];
        Dipk = Dip1QMat[xfr2][k];
        Dipl = Dip1QMat[xfr3][l];
        Dip2Q = Dip2QSpMat[xfr3][l];
        int ik = i*3 + k;
        int jk = j*3 + k;
        if (reph) {
          pathway_double_sides(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &REPH_private);
          pathway_double_sides(Dipl, Dipk, psi_b1[tau1+tau2+tau3][j], psi_a[tau1+tau2][i], nosc, popfac1Q, &REPH_private);
          pathway_2Q(psi_cb[jk], psi2Qa_ba, n2Q, popfac2Q, &REPH_private);
        }
        if (nreph) {
          pathway_single_side(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &NREPH_private);
          pathway_double_sides(Dipl, Dipk, psi_a[tau1+tau2+tau3][i], psi_b1[tau1+tau2][j], nosc, popfac1Q, &NREPH_private);
          pathway_2Q(psi_ca[ik], psi2Qb1_bb, n2Q, popfac2Q, &NREPH_private);
        }
      }
    }
  }
  // cjfeng 07/03/2017
  // Computing different polarized spectrum
  for(P=0; P<npol; P++) {
    int nP = pol_info->POL[P];
    GNREAL orient_fac = pol_info->M_ijkl_IJKL[P][3];
    if (reph) REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
    if (nreph) NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
  }
  unset_gncomp_1d_array(psi2Qa_ba);
  unset_gncomp_1d_array(psi2Qb1_bb);
}

// cjfeng 01/10/2019
// New calculation algorithm by explicitly calculating iiii, iijj, ijji, and ijij
int calc_2dir_nise_Sp(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param, int nthreads) {
  calc_R3_iiii_Sp(psi_a, psi_b1, psi_b12, pol_info, traj_param, popfac1Q, popfac2Q, tau1, tau2, tau3, xfr1, xfr2, xfr3, window, REPH, NREPH, spec_param, nthreads);
  calc_R3_iijj_Sp(psi_a, psi_b1, psi_b12, pol_info, traj_param, popfac1Q, popfac2Q, tau1, tau2, tau3, xfr1, xfr2, xfr3, window, REPH, NREPH, spec_param, nthreads);
  calc_R3_ijji_Sp(psi_a, psi_b1, psi_b12, pol_info, traj_param, popfac1Q, popfac2Q, tau1, tau2, tau3, xfr1, xfr2, xfr3, window, REPH, NREPH, spec_param, nthreads);
  calc_R3_ijij_Sp(psi_a, psi_b1, psi_b12, pol_info, traj_param, popfac1Q, popfac2Q, tau1, tau2, tau3, xfr1, xfr2, xfr3, window, REPH, NREPH, spec_param, nthreads);
  return 0;
}

// cjfeng 03/23/2017
// This routine aims for computing the nonlinear response function in nise algorithm.
// cjfeng 12/07/2018
// Using pol_info, and spec_param. Added reph and nreph flags
// cjfeng 01/04/2019
// Added traj_param
// cjfeng 01/10/2019
// Swapped xyz and window
int calc_2dir_nise(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, POL_info *pol_info, traj_param *traj_param, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, spec_param *spec_param,int nthreads) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int tau1window = tau1 * window;
  int P, p, i, j, k, l, n, a, b;
  int npol = pol_info->npol;
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;

  // cjfeng 07/02/2017
  GNCOMP *psi2Qa_aa, *psi2Qa_ba, *psi2Qb1_aa, *psi2Qb1_ab, *psi2Qb1_ba, *psi2Qb1_bb; 
  set_gncomp_1d_array(&psi2Qa_aa, n2Q, 1);
  set_gncomp_1d_array(&psi2Qa_ba, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_aa, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_ab, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_ba, n2Q, 1);
  set_gncomp_1d_array(&psi2Qb1_bb, n2Q, 1);
  // cjfeng 07/03/2017
  // Change the order of loops
  // cjfeng 07/02/2017
  // Further optimzations to reduce number of 2Q computations.
  // Recognizing that i is always equal to a, and 2Q transitions are only related to i, and j,
  // we can pull out the i terms first and reduce operations.
  // p sums over the forms iiii, iijj, ijji, and ijij.
  for(a=0; a<3; a++) {
    // cjfeng 07/02/2017
    GNREAL *Dip2Qa = Dip2QMat[xfr3][a];
    mvmult_comp_trans_x(Dip2Qa, psi_a[tau1+tau2+tau3][a], psi2Qa_aa, nosc, n2Q, nthreads);
    mvmult_comp_trans_x(Dip2Qa, psi_b1[tau1+tau2+tau3][a], psi2Qb1_aa, nosc, n2Q, nthreads);
    for(b=0; b<3; b++) {
      // cjfeng 07/02/2017
      GNREAL *Dip2Qb = Dip2QMat[xfr3][b];
      if ( b != a ) {
        mvmult_comp_trans_x(Dip2Qb, psi_a[tau1+tau2+tau3][a], psi2Qa_ba, nosc, n2Q, nthreads);
        mvmult_comp_trans_x(Dip2Qa, psi_b1[tau1+tau2+tau3][b], psi2Qb1_ab, nosc, n2Q, nthreads);
        mvmult_comp_trans_x(Dip2Qb, psi_b1[tau1+tau2+tau3][a], psi2Qb1_ba, nosc, n2Q, nthreads);
        mvmult_comp_trans_x(Dip2Qb, psi_b1[tau1+tau2+tau3][b], psi2Qb1_bb, nosc, n2Q, nthreads);
      }
      else {
        for(n=0; n<n2Q; n++) psi2Qa_ba[n] = psi2Qa_aa[n];
        for(n=0; n<n2Q; n++) psi2Qb1_ab[n] = psi2Qb1_aa[n];
        for(n=0; n<n2Q; n++) psi2Qb1_ba[n] = psi2Qb1_aa[n];
        for(n=0; n<n2Q; n++) psi2Qb1_bb[n] = psi2Qb1_aa[n];
      }
      for(p=0; p<4; p++) {
        if(p==0) { i=a; j=a; k=a; l=a; }
        else if(p==1) { i=a; j=a; k=b; l=b; }
        else if(p==2) { i=a; j=b; k=b; l=a; }
        else { i=a; j=b; k=a; l=b; } // p==3
        // For p==0, we only add the signal when a==b.
        if( (p!=0) ) {
          if ( (a!=b) ) {
            GNCOMP REPH_private = 0.0;
            GNCOMP NREPH_private = 0.0;
            GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
            Dipj = Dip1QMat[xfr1][j];
            Dipk = Dip1QMat[xfr2][k];
            Dipl = Dip1QMat[xfr3][l];
            Dip2Q = Dip2QMat[xfr3][l];
            int ik = i*3 + k;
            int jk = j*3 + k;
            // cjfeng 07/03/2017
            // Replaced the original code by sub-routines
            if (reph) {
              pathway_double_sides(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &REPH_private);
              pathway_double_sides(Dipl, Dipk, psi_b1[tau1+tau2+tau3][j], psi_a[tau1+tau2][i], nosc, popfac1Q, &REPH_private);
              pathway_REPH_2Q(psi_cb[jk], psi2Qa_aa, psi2Qa_ba, n2Q, p, popfac2Q, &REPH_private);
            }
            if (nreph) {
              pathway_single_side(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &NREPH_private);
              // cjfeng 12/18/2018
              // psi_a[i] instead of psi_a[j], and psi_b1[j] instead of psi_b1[i]
              pathway_double_sides(Dipl, Dipk, psi_a[tau1+tau2+tau3][i], psi_b1[tau1+tau2][j], nosc, popfac1Q, &NREPH_private);
              pathway_NREPH_2Q(psi_ca[ik], psi2Qb1_aa, psi2Qb1_ba, psi2Qb1_ab, psi2Qb1_bb, n2Q, p, popfac2Q, &NREPH_private);
            }
            // cjfeng 07/03/2017
            // Computing different polarized spectrum
            for(P=0; P<npol; P++) {
              int nP = pol_info->POL[P];
              GNREAL orient_fac = pol_info->M_ijkl_IJKL[P][p];
              if (reph) REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
              if (nreph) NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
            }
          }
        }
        else if( (a==b) ) {
          GNCOMP REPH_private = 0.0;
          GNCOMP NREPH_private = 0.0;
          GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
          Dipj = Dip1QMat[xfr1][j];
          Dipk = Dip1QMat[xfr2][k];
          Dipl = Dip1QMat[xfr3][l];
          Dip2Q = Dip2QMat[xfr3][l];
          int ik = i*3 + k;
          int jk = j*3 + k;
          // cjfeng 07/03/2017
          // Replaced the original code by sub-routines
          if (reph) {
            pathway_double_sides(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &REPH_private);
            pathway_double_sides(Dipl, Dipk, psi_b1[tau1+tau2+tau3][j], psi_a[tau1+tau2][i], nosc, popfac1Q, &REPH_private);
            pathway_REPH_2Q(psi_cb[jk], psi2Qa_aa, psi2Qa_ba, n2Q, p, popfac2Q, &REPH_private);
          }
          if (nreph) {
            pathway_single_side(Dipl, Dipj, psi_b12[tau1+tau2+tau3][k], psi_a[tau1][i], nosc, popfac1Q, &NREPH_private);
            // cjfeng 12/18/2018
            // psi_a[i] instead of psi_a[j], and psi_b1[j] instead of psi_b1[i]
            pathway_double_sides(Dipl, Dipk, psi_a[tau1+tau2+tau3][i], psi_b1[tau1+tau2][j], nosc, popfac1Q, &NREPH_private);
            pathway_NREPH_2Q(psi_ca[ik], psi2Qb1_aa, psi2Qb1_ba, psi2Qb1_ab, psi2Qb1_bb, n2Q, p, popfac2Q, &NREPH_private);
          }
          // cjfeng 07/03/2017
          // Computing different polarized spectrum
          for(P=0; P<npol; P++) {
            int nP = pol_info->POL[P];
            GNREAL orient_fac = pol_info->M_ijkl_IJKL[P][p];
            if (reph) REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
            if (nreph) NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
          }
        }
      }
    }
  }
  unset_gncomp_1d_array(psi2Qa_aa);
  unset_gncomp_1d_array(psi2Qa_ba);
  unset_gncomp_1d_array(psi2Qb1_aa);
  unset_gncomp_1d_array(psi2Qb1_ab);
  unset_gncomp_1d_array(psi2Qb1_ba);
  unset_gncomp_1d_array(psi2Qb1_bb);
  return 0;
}

// cjfeng 07/05/2017
// Grouping individual pathways into non-redundant functions
int pathway_single_side(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val) {
  // Pathway involving only the ket side: 
  //
  // Non-rephasing pathway:
  //  | b 0 |       -- tau1+tau2+tau3
  //  | 0 0 |       -- tau1+tau2
  //  | a 0 |       -- tau1
  //  | 0 0 |       -- 0
  //  Dipl * psi_b1
  //  Dipj * psi_a    
  GNCOMP cval1, cval2;
  int n;
  cval1 = 0.0; cval2 = 0.0;
  for(n=0; n<nosc; n++) cval1 += Dip1[n] * psi1[n]; // Dipl and psi_a
  for(n=0; n<nosc; n++) cval2 += Dip2[n] * psi2[n]; // Dipk and psi_b1
  *val += cval1 * cval2 * popfac;
  return 0;
}

int pathway_double_sides(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val) {
  // Pathways involving both sides:
  //
  // Rephasing pathways:
  //  | b 0 |         | b 0 |       -- tau1+tau2+tau3
  //  | 0 0 |         | b a |       -- tau1+tau2
  //  | 0 a |         | 0 a |       -- tau1
  //  | 0 0 |         | 0 0 |       -- 0
  //  Dipl * psi_b12  Dipl * psi_b1
  //  Dipj * psi_a    Dipk * psi_a 
  //
  // Non-rephasing pathway:
  //  | a 0 |       -- tau1+tau2+tau3
  //  | a b |       -- tau1+tau2
  //  | a 0 |       -- tau1
  //  | 0 0 |       --
  //  Dipl * psi_a
  //  Dipk * psi_b1
  // Note that the complex conjugate is required for the R.H.S.
  GNCOMP cval1, cval2;
  int n;
  cval1 = 0.0; cval2 = 0.0;
  for(n=0; n<nosc; n++) cval1 += Dip1[n] * psi1[n];
  for(n=0; n<nosc; n++) cval2 += Dip2[n] * psi2[n];
  *val += cval1 * conj(cval2) * popfac;

  return 0;
}

int pathway_REPH_2Q(GNCOMP *psi2Qcb, GNCOMP *psi2Qa_aa, GNCOMP *psi2Qa_ba, int n2Q, int p, GNREAL popfac2Q, GNCOMP *REPH_private) {
  // Final rephasing pathway:
  //
  //  | c a |
  //  | b a |
  //  | 0 a |
  //  | 0 0 |
  //
  // The contribution to the rephasing spectrum is the orientationally averaged value
  // ( psi_cb[j*3+k] ) * ( Dip2Q[tau1+tau2+tau3] * psi_a[tau1+tau2+tau3] )^\dagger
  GNCOMP cval1 = 0.0;
  int n;
  if ( p == 0 || p == 2 ) for(n=0; n<n2Q; n++) cval1 += psi2Qcb[n] * conj(psi2Qa_aa[n]);
  else for(n=0; n<n2Q; n++) cval1 += psi2Qcb[n] * conj(psi2Qa_ba[n]);
  *REPH_private -= cval1 * popfac2Q;

  return 0;
}

int pathway_NREPH_2Q(GNCOMP *psi2Qca, GNCOMP *psi2Qb1_aa, GNCOMP *psi2Qb1_ba, GNCOMP *psi2Qb1_ab, GNCOMP *psi2Qb1_bb, int n2Q, int p, GNREAL popfac2Q, GNCOMP *NREPH_private) {
  // Final non-rephasing pathway:
  //
  //  | c b |
  //  | a b |
  //  | a 0 |
  //  | 0 0 |
  //
  // The contribution to the non-rephasing spectrum is the orientationally averaged value
  // ( psi_ca[k*3+i] ) * ( Dip2QMat[tau1+tau2+tau3] * psi_b1[tau1+tau2+tau3] )
  GNCOMP cval1 = 0.0;
  int n;
  if ( p == 0 ) for(n=0; n<n2Q; n++) cval1 += psi2Qca[n] * conj(psi2Qb1_aa[n]);
  else if ( p == 1 ) for(n=0; n<n2Q; n++) cval1 += psi2Qca[n] * conj(psi2Qb1_ba[n]);
  else if ( p == 2 ) for(n=0; n<n2Q; n++) cval1 += psi2Qca[n] * conj(psi2Qb1_ab[n]);
  else for(n=0; n<n2Q; n++) cval1 += psi2Qca[n] * conj(psi2Qb1_bb[n]);
  *NREPH_private -= cval1 * popfac2Q;

  return 0;
}

// cjfeng 01/10/2019
// The 2Q pathways can be further reduced to a single subroutine
// Final rephasing pathway:
//
//  | c a |
//  | b a |
//  | 0 a |
//  | 0 0 |
//
// The contribution to the rephasing spectrum is the orientationally averaged value
// ( psi_cb[j*3+k] ) * ( Dip2Q[tau1+tau2+tau3] * psi_a[tau1+tau2+tau3] )^\dagger
// Final non-rephasing pathway:
//
//  | c b |
//  | a b |
//  | a 0 |
//  | 0 0 |
//
// The contribution to the non-rephasing spectrum is the orientationally averaged value
// ( psi_ca[k*3+i] ) * ( Dip2QMat[tau1+tau2+tau3] * psi_b1[tau1+tau2+tau3] )
int pathway_2Q(GNCOMP *psi2Qc, GNCOMP *psi2Q, int n2Q, GNREAL popfac2Q, GNCOMP *val) {
  GNCOMP cval1 = 0.0;
  int n;
  for(n=0; n<n2Q; n++) cval1 += psi2Qc[n] * conj(psi2Q[n]);
  *val -= cval1 * popfac2Q;
  return 0;
}
