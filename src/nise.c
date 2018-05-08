#include "nise.h"

// Propagator
int gen_UMat(GNCOMP **UMat, GNREAL *Ham, GNREAL *Evals, int N, int xfr, double expfac, double wo) {
  int i, j;
  int Nsq = N * N;
  // Initialize UMat
  for(i=0; i<Nsq; i++) UMat[xfr][i] = 0.0;
  // Compute U2Q elements
  for(i=0; i<N; i++) {
    int ni = i*N;
    GNREAL cre, cim;
    for(j=0; j<N; j++) {
      int k, nj = j*N;
      cre = 0.0;
      cim = 0.0;
      for(k=0; k<N; k++) {
        cre += Ham[ni+k] * Ham[nj+k] * cos( expfac * (Evals[k] - wo) );
        cim += Ham[ni+k] * Ham[nj+k] * sin( expfac * (Evals[k] - wo) );
      }
      UMat[xfr][ni+j] = cre + I*cim;
    }
  }
  return 0;
}

int gen_pert_prop( GNCOMP *U1Q, GNCOMP *U2Q, int nosc, int n2Q, double expfac, double delta) {
  int i,j,ip,jp;
  GNCOMP fac;
  double sqrt2inv = 0.707106781186547;
  // negative sign on sin() term since delta is provided as a positive argument.
  GNCOMP facnew = sqrt2inv * ( cos(0.5*expfac*delta) - I*sin(0.5*expfac*delta) );
  // cjfeng 08/01/2016
  // No initialization required
  // for(i=0; i<n2Qsq; i++) U2Q[i] = 0.0;
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
int initialize_psi(GNCOMP ***psi, GNREAL **Dip1QMat[3], int tinit, int xfr, int nosc, int nelem) {
  int i, j;
  for(i=0; i<nelem; i++) for(j=0; j<nosc; j++) psi[i][tinit][j] = Dip1QMat[i][xfr][j];
  return 0;
}

int calc_CorrFunc(GNCOMP *CorrFunc, GNCOMP **U1QMat, GNCOMP **psi_1d[3], GNREAL ***Dip1QMat, GNCOMP *cDip1Q[3], int nosc, int window, int fr, int nbuffer, double expfac, double wo, BLAS_gemv_opt *gemv_opt, int nthreads) {
  int noscsq = nosc * nosc;
  int i, j, n;
  int xfr;
  GNCOMP cval;
  initialize_psi(psi_1d, Dip1QMat, 0, fr, nosc, 3);
  
  for(n=0; n<window; n++) {
    xfr = (fr+n) % nbuffer;
    cval = 0.0;
    #if OMP_PARALLEL
    #pragma omp parallel if(nosc >= NBCOMP && nthreads>1) shared(nosc, cDip1Q, Dip1QMat, xfr, fr) private(i,j)
    #pragma omp for schedule(guided) collapse(2) reduction(+:cval)
    #endif
    for(i=0; i<3; i++) {
      for(j=0; j<nosc; j++) {
        cDip1Q[i][j] = psi_1d[i][0][j];
        cval += Dip1QMat[i][xfr][j] * cDip1Q[i][j];
      }
    }
    CorrFunc[n] += cval;
    
    #if BLAS_subroutines
    if(nosc > nosc_switch) for(i=0; i<3; i++) cgemv_(&(gemv_opt->ta), &nosc, &nosc, &(gemv_opt->alpha), U1QMat[xfr], &nosc, psi_1d[i][0], &(gemv_opt->inc), &(gemv_opt->beta), psi_1d[i][1], &(gemv_opt->inc));  
    else for(i=0; i<3; i++) mvmult_comp(U1QMat[xfr], psi_1d[i][0], psi_1d[i][1], nosc, nthreads);
    #else
    for(i=0; i<3; i++) mvmult_comp(U1QMat[xfr], psi_1d[i][0], psi_1d[i][1], nosc, nthreads);
    #endif
    for(i=0; i<3; i++) copy_gncomp(psi_1d[i][1], psi_1d[i][0], nosc, 1);
  }
  return 0;
}

// 2D IR subroutines
int initialize_psi2Q(GNCOMP ***psi, GNREAL ***Dip2QMat, GNCOMP **psi_2Q, int nosc, int n2Q, int nelem, int tau1, int tau2, int nthreads) {
  int i, j;
  for(i=0; i<nelem; i++) {
    int ii = i*3;
    for(j=0; j<nelem; j++) {
      mvmult_comp_trans_x(Dip2QMat[tau1+tau2][j], psi[i][tau1+tau2], psi_2Q[ii+j], nosc, n2Q, nthreads);
    }
  }
  return 0;
}
int propagate_psi1Q(GNCOMP **psi[3], GNREAL **Dip1QMat[3], GNCOMP **U1QMat, int nosc, int fr, int tau0, int xfr0, int nbuffer, int winsize, BLAS_gemv_opt *gemv_opt, int nthreads) {
  int i, tau;
  // Initial state
  initialize_psi(psi, Dip1QMat, tau0, xfr0, nosc, 3);
  #if BLAS_subroutines
  if(nosc > nosc_switch) {
    for(i=0; i<3; i++) for(tau=tau0; tau<tau0+winsize-1; tau++) cgemv_(&(gemv_opt->ta), &nosc, &nosc, &(gemv_opt->alpha), U1QMat[(fr+tau)%nbuffer], &nosc, psi[i][tau], &(gemv_opt->inc), &(gemv_opt->beta), psi[i][tau+1], &(gemv_opt->inc));
  }
  else {
    for(i=0; i<3; i++) for(tau=tau0; tau<tau0+winsize-1; tau++) mvmult_comp(U1QMat[(fr+tau)%nbuffer], psi[i][tau], psi[i][tau+1], nosc, nthreads );
  }
  #else
  for(i=0; i<3; i++) for(tau=tau0; tau<tau0+winsize-1; tau++) mvmult_comp(U1QMat[(fr+tau)%nbuffer], psi[i][tau], psi[i][tau+1], nosc, nthreads );
  #endif
  return 0;
}

int propagate_psi_2Q(GNCOMP **U2QMat, GNCOMP **psi_ca, GNCOMP **psi_cb, GNCOMP *psi2Q, int n2Q, int xfr3, BLAS_gemv_opt *gemv_opt, int nthreads) {
  int i, j;
  for(i=0; i<9; i++) {
    // First initialize psi_ca and propagate it.
    for(j=0; j<n2Q; j++) psi2Q[j] = psi_ca[i][j];
    #if BLAS_subroutines
    if (n2Q > nosc_switch) cgemv_(&(gemv_opt->ta), &n2Q, &n2Q, &(gemv_opt->alpha), U2QMat[xfr3], &n2Q, psi2Q, &(gemv_opt->inc), &(gemv_opt->beta), psi_ca[i], &(gemv_opt->inc));
    else mvmult_comp(U2QMat[xfr3], psi2Q, psi_ca[i], n2Q, nthreads);
    #else
    mvmult_comp(U2QMat[xfr3], psi2Q, psi_ca[i], n2Q, nthreads);
    #endif
    // Initialize psi_cb and propagate it.
    for(j=0; j<n2Q; j++) psi2Q[j] = psi_cb[i][j];
    #if BLAS_subroutines
    if (n2Q > nosc_switch) cgemv_(&(gemv_opt->ta), &n2Q, &n2Q, &(gemv_opt->alpha), U2QMat[xfr3], &n2Q, psi2Q, &(gemv_opt->inc), &(gemv_opt->beta), psi_cb[i], &(gemv_opt->inc));
    else mvmult_comp(U2QMat[xfr3], psi2Q, psi_cb[i], n2Q, nthreads);
    #else
    mvmult_comp(U2QMat[xfr3], psi2Q, psi_cb[i], n2Q, nthreads);
    #endif
  }
  return 0;
}

int propagate_psi_2Q_pert(GNCOMP **U1QMat, GNCOMP *U2Qs, GNCOMP **psi_ca, GNCOMP **psi_cb, GNCOMP *psi2Q, int nosc, int n2Q, real delta, double expfac, int xfr3, int nthreads) {
  int i, j;
  gen_pert_prop(U1QMat[xfr3], U2Qs, nosc, n2Q, expfac, delta);
  for(i=0; i<9; i++) {
    // First initialize psi_ca and propagate it.
    for(j=0; j<n2Q; j++) psi2Q[j] = psi_ca[i][j];
    mvmult_comp(U2Qs, psi2Q, psi_ca[i], n2Q, nthreads);
    // Initialize psi_cb and propagate it.
    for(j=0; j<n2Q; j++) psi2Q[j] = psi_cb[i][j];
    mvmult_comp(U2Qs, psi2Q, psi_cb[i], n2Q, nthreads);
  }
  return 0;
}

// cjfeng 03/23/2017
// This routine aims for computing the nonlinear response function in nise algorithm.
int calc_2dir_nise(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, int *POL, int nosc, int n2Q, int npol, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, int nthreads) {
  int tau1window = tau1 * window;
  int P, p, i, j, k, l, n, a, b;
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
    mvmult_comp_trans_x(Dip2Qa, psi_a[a][tau1+tau2+tau3], psi2Qa_aa, nosc, n2Q, nthreads);
    mvmult_comp_trans_x(Dip2Qa, psi_b1[a][tau1+tau2+tau3], psi2Qb1_aa, nosc, n2Q, nthreads);
    for(b=0; b<3; b++) {
      // cjfeng 07/02/2017
      GNREAL *Dip2Qb = Dip2QMat[xfr3][b];
      if ( b != a ) {
        mvmult_comp_trans_x(Dip2Qb, psi_a[a][tau1+tau2+tau3], psi2Qa_ba, nosc, n2Q, nthreads);
        mvmult_comp_trans_x(Dip2Qa, psi_b1[b][tau1+tau2+tau3], psi2Qb1_ab, nosc, n2Q, nthreads);
        mvmult_comp_trans_x(Dip2Qb, psi_b1[a][tau1+tau2+tau3], psi2Qb1_ba, nosc, n2Q, nthreads);
        mvmult_comp_trans_x(Dip2Qb, psi_b1[b][tau1+tau2+tau3], psi2Qb1_bb, nosc, n2Q, nthreads);
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
        if( (p!=0) || (a==b) ) {
          GNCOMP REPH_private = 0.0;
          GNCOMP NREPH_private = 0.0;
          GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
          Dipj = Dip1QMat[j][xfr1];
          Dipk = Dip1QMat[k][xfr2];
          Dipl = Dip1QMat[l][xfr3];
          Dip2Q = Dip2QMat[xfr3][l];
          int ik = i*3 + k;
          int jk = j*3 + k;
          // cjfeng 07/03/2017
          // Replaced the original code by sub-routines
          pathway_double_sides(Dipl, Dipj, psi_b12[k][tau1+tau2+tau3], psi_a[i][tau1], nosc, popfac1Q, &REPH_private);
          pathway_double_sides(Dipl, Dipk, psi_b1[j][tau1+tau2+tau3], psi_a[i][tau1+tau2], nosc, popfac1Q, &REPH_private);
          pathway_REPH_2Q(psi_cb[jk], psi2Qa_aa, psi2Qa_ba, n2Q, p, popfac2Q, &REPH_private);
          pathway_single_side(Dipl, Dipj, psi_b12[k][tau1+tau2+tau3], psi_a[i][tau1], nosc, popfac1Q, &NREPH_private);
          pathway_double_sides(Dipl, Dipk, psi_a[j][tau1+tau2+tau3], psi_b1[i][tau1+tau2], nosc, popfac1Q, &NREPH_private);
          pathway_NREPH_2Q(psi_ca[ik], psi2Qb1_aa, psi2Qb1_ba, psi2Qb1_ab, psi2Qb1_bb, n2Q, p, popfac2Q, &NREPH_private);
          // cjfeng 07/03/2017
          // Computing different polarized spectrum
          for(P=0; P<npol; P++) {
            int nP = POL[P];
            GNREAL orient_fac = M_ijkl_IJKL[P][p];
            REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
            NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
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
  // cjfeng 07/02/2017
  if ( p == 0 || p == 2) for(n=0; n<n2Q; n++) cval1 += psi2Qcb[n] * conj(psi2Qa_aa[n]);
  else for(n=0; n<n2Q; n++) cval1 += psi2Qcb[n] * conj(psi2Qa_ba[n]);
  // cjfeng 07/02/2017
  // The original code
  // mvmult_comp_trans_x(Dip2Q, psi_a[i][tau1+tau2+tau3], psi2Q, nosc, n2Q, nthreads);
  // mvmult_comp_trans_x(Dip2QMat[xfr3][l], psi_a[i][tau1+tau2+tau3], psi2Q, nosc, n2Q, nthreads);
  // for(n=0; n<n2Q; n++) cval1 += psi_cb[jk][n] * conj(psi2Q[n]);
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
