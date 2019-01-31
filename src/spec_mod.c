#include "spec_mod.h"

/*********************************
* Spectral parameter generations *
**********************************/

// cjfeng 01/04/2019
// Added traj_param
int gen_ham_2Q(GNREAL *Ham1Q, GNREAL *Ham2Q, traj_param *traj_param, real delta) {
  int m,n,p;
  int ndx;
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int n2Qsq = n2Q * n2Q;
  const GNREAL sqrt2 = sqrt(2.0);
  
  // States are numbered in order 00, 10, 11, 20, 21, 22,...,(nosc-1)(nosc-1)
  // The index of state mn is (m+1)*m/2 + n.
  for(m=0; m<n2Qsq; m++) Ham2Q[m] = 0.0;
  for(m=0; m<nosc; m++) {
    int mm = m*nosc;
    int nn = n*nosc;
    // Double excitation at site m
    Ham2Q[NDX2Q(m,m)*n2Q + NDX2Q(m,m)] = 2*Ham1Q[m*nosc+m]-delta;
    // Single excitations at sites m and n
    for(n=0; n<m; n++) Ham2Q[ NDX2Q(m,n)*n2Q + NDX2Q(m,n) ] = Ham1Q[m*nosc+m]+Ham1Q[n*nosc+n];
    // Coupling between mm and mn states: <mm|H|mn> = sqrt(2)*<m|H|n>
    for(n=0; n<m; n++) Ham2Q[ NDX2Q(m,m)*n2Q + NDX2Q(m,n) ] = sqrt2*Ham1Q[m*nosc+n];
    // Coupling between mm and nm states: <mm|H|nm> = sqrt(2)*<m|H|n>
    for(n=m+1; n<nosc; n++) Ham2Q[ NDX2Q(m,m)*n2Q + NDX2Q(n,m) ] = sqrt2*Ham1Q[m*nosc+n];
    // Coupling between nm and mp states (n<p<m): <mn|H|mp> = <n|H|p>
    for(p=0; p<m; p++) for(n=0; n<p; n++) Ham2Q[ NDX2Q(m,n)*n2Q + NDX2Q(m,p) ] = Ham1Q[n*nosc+p];
    // Coupling between pm and mn states (n<m<p): <mn|H|pm> = <n|H|p>
    for(n=0; n<m; n++) for(p=m+1; p<nosc; p++) Ham2Q[ NDX2Q(m,n)*n2Q + NDX2Q(p,m) ] = Ham1Q[n*nosc+p];
    // Coupling between pm and nm states (m<n<p): <nm|H|pm> = <n|H|p>
    for(n=m+1; n<nosc; n++) for(p=n+1; p<nosc; p++) Ham2Q[ NDX2Q(n,m)*n2Q + NDX2Q(p,m) ] = Ham1Q[n*nosc+p];
  }
  for(m=0; m<n2Q; m++) for(n=0; n<m; n++) Ham2Q[m*n2Q+n] += Ham2Q[n*n2Q+m];
  for(m=0; m<n2Q; m++) for(n=m+1; n<n2Q; n++) Ham2Q[m*n2Q+n] = Ham2Q[n*n2Q+m];
  return 1;
}

// cjfeng 01/04/2019
// Added traj_param
// cjfeng 01/10/2019
// Including xyz loop into this sub-routine
int gen_dip_2Q(GNREAL **Dip1Q, GNREAL **Dip2Q, traj_param *traj_param) {
  int i, m, n, M;
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int n2Qnosc = n2Q * nosc;
  GNREAL fac;
  GNREAL sqrt2 = sqrt(2.0);
  for(i=0; i<3; i++) for(M=0; M<n2Qnosc; M++) Dip2Q[i][M] = 0.0;
 
  for(i=0; i<3; i++) { 
    for(n=0; n<nosc; n++) {
      GNREAL Dipval = Dip1Q[i][n];
      int nn2Q = n*(n+1)/2;
      for(m=0; m<nosc; m++) {
        if(m!=n) fac = 1.0;
        else fac = sqrt2;
        if(m<n) M = (nn2Q+m);
        else M = ((m*(m+1)/2)+n);
        
        Dip2Q[i][M*nosc+m] = fac * Dipval;
      }
    }
  }
  return 1;
}

int gen_dip_Sp2Q_ind(int *Dip2Qrow, int *Dip2Qcol, traj_param *traj_param) {
  int m, n, M;
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int noscsq = nosc * nosc;
  for(n=0; n<nosc; n++) {
    int nn2Q = n*(n+1)/2;
    int nnosc = n*nosc;
    for(m=0; m<nosc; m++) {
      if(m<n) M = (nn2Q+m);
      else M = ((m*(m+1)/2)+n);
      Dip2Qrow[nnosc+m] = M;
      Dip2Qcol[nnosc+m] = m;
    }
  }
  return 1;
}

int gen_dip_Sp2Q(GNREAL **Dip1Q, GNREAL **Dip2Q, traj_param *traj_param) {
  int i, m, n, M;
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int noscsq = nosc * nosc;
  GNREAL fac;
  GNREAL sqrt2 = sqrt(2.0);
  for(i=0; i<3; i++) for(M=0; M<noscsq; M++) Dip2Q[i][M] = 0.0;
  
  for(i=0; i<3; i++) {
    for(n=0; n<nosc; n++) {
      GNREAL Dipval = Dip1Q[i][n];
      int nnosc = n*nosc;
      for(m=0; m<nosc; m++) {
        if(m!=n) fac = 1.0;
        else fac = sqrt2;
        Dip2Q[i][nnosc+m] = fac * Dipval;
      }
    }
  }
  return 1;
}

// cjfeng 01/04/2019
// Added traj_param
int gen_ExDip2Q(GNREAL ***ExDip2QAr, GNREAL ***Dip2QAr, GNREAL **Ham1QAr, GNREAL **Ham2QAr, traj_param *traj_param, int tid) {
  int i, n, N, k;
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  for(i=0; i<3; i++) {
    for(n=0; n<nosc; n++) {
      int nn = n * n2Q;
      for(N=0; N<n2Q; N++) ExDip2QAr[i][tid][nn+N] = 0.0;
    }
    // cjfeng 04/11/2017
    // This part can be potentially faster.
    for(n=0; n<nosc; n++) {
      int nn = n * n2Q;
      for(k=0; k<nosc; k++) {
        int nk = k * nosc;
        int Nk = k * n2Q;
        GNREAL Hamval = Ham1QAr[tid][nk+n];
        for(N=0; N<n2Q; N++) ExDip2QAr[i][tid][nn+N] += Hamval * Dip2QAr[i][tid][Nk+N];  
      }
    }
    for(n=0; n<nosc; n++) {
      int nn = n * n2Q;
      for(N=0; N<n2Q; N++) Dip2QAr[i][tid][nn+N] = ExDip2QAr[i][tid][nn+N];
    }
    for(n=0; n<nosc; n++) {
      int nn = n * n2Q;
      for(N=0; N<n2Q; N++) ExDip2QAr[i][tid][nn+N] = 0.0;
    }
    for(n=0; n<nosc; n++) {
      int nn = n * n2Q;
      for(k=0; k<n2Q; k++) {
        int Nk = k * n2Q;
        for(N=0; N<n2Q; N++) ExDip2QAr[i][tid][nn+N] += Dip2QAr[i][tid][nn+k] * Ham2QAr[tid][Nk+N];
      }
    }
  }
  return 0;
}

// cjfeng 01/04/2019
// Added traj_param
int gen_avg_ham(GNREAL **ham, GNREAL *avg_ham, GNREAL *hann, int fr, w_param *w_param, traj_param *traj_param) {
  int nosc = traj_param->nosc;
  int nbuffer = traj_param->nbuffer;
  int whann = w_param->whann;
  int window = w_param->window;
  int noscsq = nosc * nosc;
  int xfr_half_wind = (fr + ((int) (window/2))) % nbuffer;
  int i, j;
  // Initialize averaged Hamiltonian.
  for(j=0; j<noscsq; j++) avg_ham[j] = 0.0;
  // Compute averaged Hamiltonian
  if(whann) {
    for(i=0; i<window; i++) {
      int xfri = (fr + i) % nbuffer;
      for(j=0; j<noscsq; j++) avg_ham[j] += ham[xfri][j] * hann[i];
    }
  }
  else {
    for(i=0; i<window; i++) {
      int xfri = (fr + i) % nbuffer;
      for(j=0; j<noscsq; j++) avg_ham[j] += ham[xfri][j] / window;
    }
  } 
  return 0;
}

// cjfeng 12/09/2018
// Using spec_param, and w_param
// cjfeng 01/04/2019
// Added traj_param
int gen_avg_dip(GNREAL ***dip1Q, GNREAL ***dip2Q, GNREAL ***avg_dip1Q, GNREAL ***avg_dip2Q, int fr, int nAvgd, traj_param *traj_param, spec_param *spec_param, w_param *w_param) {
  int i, j, k;
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int nbuffer = traj_param->nbuffer;
  int nosc_n2Q = nosc * n2Q;
  int do2d = spec_param->do2d;
  int pert = spec_param->pert;
  int whann = w_param->whann;
  int window = w_param->window;
  int xfr_half_wind = (fr + ((int) (window/2))) % nbuffer; 
  // Initialize averaged dipole moments.
  for(i=0; i<3; i++) for(j=0; j<nosc; j++) avg_dip1Q[i][nAvgd][j] = 0.0;
  if( (do2d) && (!pert) ) for(i=0; i<3; i++) for(j=0; j<nosc_n2Q; j++) avg_dip2Q[i][nAvgd][j] = 0.0;
  if (whann) {
    // one-quantum dipole moment
    for(i=0; i<3; i++) for(j=0; j<nosc; j++) avg_dip1Q[i][nAvgd][j] = dip1Q[xfr_half_wind][i][j];
    // two-quantum dipole moment
    if( (do2d) && (!pert) ) {
      for(i=0; i<3; i++) for(j=0; j<nosc; j++) {
        int jj = j*n2Q;
        for(k=0; k<n2Q; k++) avg_dip2Q[i][nAvgd][jj+k] = dip2Q[xfr_half_wind][i][k*nosc+j];
      }
    }
  }
  else {
    // one-quantum dipole moment
    for(i=0; i<3; i++) for(j=0; j<nosc; j++) avg_dip1Q[i][nAvgd][j] = dip1Q[xfr_half_wind][i][j];
    // two-quantum dipole moment
    if( (do2d) && (!pert) ) {
      for(i=0; i<3; i++) for(j=0; j<nosc; j++) {
        int jj = j * n2Q;
        for(k=0; k<n2Q; k++) avg_dip2Q[i][nAvgd][jj+k] = dip2Q[xfr_half_wind][i][k*nosc+j];
      }
    }
  }  
  return 0;
}

int gen_perturb_2Q_energies(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evals2Q, int nosc, real delta ) {
  int m,n,i;
  GNREAL val, anh;
  for(m=0; m<nosc; m++) {
    for(n=0; n<=m; n++) {
      val = Evals1Q[m] + Evals1Q[n]; // harmonic 
      anh = 0.0;
      for(i=0; i<nosc; i++) {
        int ni = i*nosc;
        anh += Evecs1Q[ni+n]*Evecs1Q[ni+m]*Evecs1Q[ni+n]*Evecs1Q[ni+m];
      }
      if(m!=n) anh = anh * 2.0;
      Evals2Q[NDX2Q(m,n)] = val-delta*anh;
    }
  }
  return 1;
}

// cjfeng 12/08/2018
// Using spec_param and w_param
int gen_popdecay(GNREAL *popdecay, spec_param *spec_param, w_param *w_param, int if2Q) {
  if(popdecay != NULL) {
    int i;
    GNREAL tstep = spec_param->tstep;
    GNREAL TauP;
    if(!if2Q) TauP = spec_param->TauP;
    else TauP = spec_param->TauP2Q;
    int window = w_param->window;
    int win2d = w_param->win2d;

    int npopdec = 2*(win2d-window);
    GNREAL expfac = tstep/(2*TauP);
    for(i=0; i<npopdec; i++) popdecay[i] = exp(-i*expfac);
  }
  return 0;
}

// cjfeng 01/04/2019
// Added phys_const
int gen_hann(GNREAL *hann, int window, phys_const *Const) {
  if(window==1) hann[0] = 1;
  else {
    int i;
    for(i=0; i<window; i++) hann[i] = 0.5*(1+cos((Const->pi*i)/(window-1)));
    // And normalize s.t. sum(hann) = 1.
    GNREAL val = 0.0;
    for(i=0; i<window; i++) val += hann[i];
    for(i=0; i<window; i++) hann[i] /= val;
  }
  return 0;
}

// cjfeng 12/09/2018
// Using subroutine for adding site-specific frequency shift
int apply_shift(GNREAL *Ham1Q, shift_info *shift_info, int nosc) {
  int i;
  int nshift = shift_info->nshift;
  for(i=0; i<nshift; i++) {
    int SHIFTNDX = shift_info->SHIFTNDX[i];
    GNREAL SHIFT = shift_info->SHIFT[i];
    Ham1Q[SHIFTNDX*nosc+SHIFTNDX] += SHIFT;
  }
  return 0;
}

/********************************************************************************
* Spectral calculations.              *
********************************************************************************/

// cjfeng 12/12/2018
// Wrapping parts of generating U1QMat into a single subroutine
// cjfeng 01/04/2019
// Added traj_param
int gen_U1Q(GNREAL **Ham1QMat, GNREAL **Evals1QAr, GNCOMP **U1QMat, spec_param *spec_param, w_param *w_param, traj_param *traj_param, int justread, int readframe, double expfac) {
  int nosc = traj_param->nosc;
  int nbuffer = traj_param->nbuffer;
  // cjfeng 03/31/2016
  // The parallel section won't be executed under single thread to avoid overhead.
  #if OMP_PARALLEL 
  #pragma omp parallel if(spec_param->nthreads>1) shared(justread, Ham1QMat, Evals1QAr, U1QMat) firstprivate(readframe, nosc, expfac)
  #endif
  {
    // cjfeng 07/19/2017
    // Use temporary struct for solving eigenvalue and eigenvectors.
    int fr, i;
    int noscsq = nosc*nosc;
    eig_array eig;
    set_eig_array(&eig, nosc);
    #if OMP_PARALLEL
    #pragma omp for schedule(guided) nowait
    #endif
    for(fr=readframe-justread; fr<readframe; fr++) {
      // cjfeng 06/27/2016
      // Reduce modulus operations.
      int xfr = fr%nbuffer;
      int tid = 0;
      lapack_int info;
      #if OMP_PARALLEL
        tid = omp_get_thread_num();
      #endif
      // Find one-quantum eigenvalues
      // Note that Ham1QMat[xfr] now contains eigenvectors
      // cjfeng 06/27/2016
      // Relatively Robust Representation to compute eigenvalues and eigenvectors.
      for (i=0; i<noscsq; i++) eig.Ham[i] = Ham1QMat[xfr][i];
      // cjfeng 12/08/2018
      // Using w_param
      info = SYEVR (LAPACK_ROW_MAJOR, 'V', 'A', 'U', nosc, eig.Ham, nosc, w_param->wstart, w_param->wstop, 0, nosc-1, 0.00001, &nosc, eig.Evals, eig.Ham, nosc, eig.isuppz);

      // Note that our Hamiltonian is actually H/(h*c)
      // The exponent we calculate is -i*2*pi*tstep*c*Ham1Q
      // cjfeng 04/07/2016
      // Wrapping operations into subroutines
      gen_UMat(U1QMat, eig.Ham, eig.Evals, eig.n, xfr, expfac);
    }
    unset_eig_array(&eig);
  }
  return 0;
}

// cjfeng 01/03/2019
// Putting if2Qfiles into spec_param
// cjfegn 01/04/2019
// Added traj_param
int gen_U2Q(GNREAL **Ham1QMat, GNREAL **Ham2QMat, GNREAL **Evals2QAr, GNCOMP **U2QMat, spec_param *spec_param, w_param *w_param, traj_param *traj_param, int justread, int readframe, real delta, double expfac) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int nbuffer = traj_param->nbuffer;
  #if OMP_PARALLEL 
  #pragma omp parallel if(spec_param->nthreads>1) shared(justread, Ham1QMat, Evals2QAr, U2QMat) firstprivate(readframe, n2Q, spec_param, expfac, delta, nosc)
  #endif
  {
    int fr, i;
    // cjfeng 07/19/2017
    // Use temporary struct.
    int n2Qsq = n2Q*n2Q;
    eig_array eig;
    set_eig_array(&eig, n2Q);
    #if OMP_PARALLEL
    #pragma omp for schedule(guided) nowait
    #endif
    for(fr=readframe-justread; fr<readframe; fr++) {
      // cjfeng 06/27/2016
      // Reduce modulus operations.
      int xfr = fr%nbuffer;
      int tid = 0;
      lapack_int info;
      #if OMP_PARALLEL
        tid = omp_get_thread_num();
      #endif
      // cjfeng 05/07/2018
      // Skip the construction of 2Q Ham if loading from external files
      if ( !(spec_param->if2Qfiles) ) gen_ham_2Q(Ham1QMat[xfr], eig.Ham, traj_param, delta);
      else { // Copy the hamiltonian from the 2Q trajectory
        for (i=0; i<n2Qsq; i++) eig.Ham[i] = Ham2QMat[xfr][i];  
      }
      // printf("Starting eigenvalue calculation for thread %d\n", tid);
      info = SYEVR(LAPACK_ROW_MAJOR, 'V', 'A', 'U', n2Q, eig.Ham, n2Q, w_param->wstart, w_param->wstop, 0, n2Q-1, 0.00001, &n2Q, eig.Evals, eig.Ham, n2Q, eig.isuppz);
      // printf("Finishing eigenvalue calculation for thread %d\n", tid);
      // cjfeng 04/07/2016
      // Wrapping operations into subroutines
      gen_UMat(U2QMat, eig.Ham, eig.Evals, eig.n, xfr, expfac);
    }
    unset_eig_array(&eig);
  }
  return 0;
}

// cjfeng 12/04/2018
// Using pol_info
GNREAL orient(GNREAL **Dip1[3], GNREAL **Dip2[3], int tid, int ndxA, int ndxB, int ndxC, int ndxD, int pol, POL_info *pol_info) {
  GNREAL Val = 0.0;
  int a,b;
  // cjfeng 06/27/2016
  // Improving locality.
  GNREAL M0, M1, M2, M3;
  M0 =  pol_info->M_ijkl_IJKL[pol][0];
  M1 =  pol_info->M_ijkl_IJKL[pol][1];
  M2 =  pol_info->M_ijkl_IJKL[pol][2];
  M3 =  pol_info->M_ijkl_IJKL[pol][3];
  for(a=0; a<3; a++) {
    GNREAL D1aA, D1aB, D1bB, D2aC, D2bC, D2aD, D2bD;
    D1aA = Dip1[a][tid][ndxA];
    D1aB = Dip1[a][tid][ndxB];
    D2aC = Dip2[a][tid][ndxC];
    D2aD = Dip2[a][tid][ndxD];
    Val += M0*D1aA*D1aB*D2aC*D2aD;
    for(b=0; b<3; b++) {
      if(b!=a) {
        D1bB = Dip1[b][tid][ndxB];
        D2bC = Dip2[b][tid][ndxC];
        D2bD = Dip2[b][tid][ndxD];
        Val += M1*D1aA*D1aB*D2bC*D2bD;
        Val += M2*D1aA*D1bB*D2bC*D2aD;
        Val += M3*D1aA*D1bB*D2aC*D2bD;
      }
    }
  }
  return Val; 
}

// cjfeng 01/03/2019
// Using a single subroutine to calculate linear absorption spectrum using sum-over-state approach
int calc_ftir(GNREAL *ftir, GNREAL **Evals1QAr, GNREAL **ExDip1QAr[3], w_param *w_param, int nosc, int tid) {
  int n, d, ndx;
  GNREAL osc;
  for(n=0; n<nosc; n++) {
    osc = 0.0;
    for(d=0; d<3; d++) osc += ExDip1QAr[d][tid][n] * ExDip1QAr[d][tid][n];
    ndx = floor( (Evals1QAr[tid][n] - w_param->wstart ) / w_param->wres + 0.5 );
    if ( ndx < w_param->npts && ndx >= 0 ) ftir[ndx] += osc;
  }
  return 0;
}

// cjfeng 01/04/2019
// Added traj_param
void eig_avg_Ham2Q(GNREAL **Ham1QAr, GNREAL **Ham2QAr, GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip2QAr[3], int *isuppz2Q, spec_param *spec_param, w_param *w_param, traj_param *traj_param, int tid, real delta) {
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q; 
  if( !(spec_param->pert) ) {
    // Two-quantum eigenvalues
    // Note that Ham2QAr[tid] now contains eigenvectors
    SYEVR( LAPACK_ROW_MAJOR, 'V', 'A', 'U', n2Q, Ham2QAr[tid], n2Q, w_param->wstart, w_param->wstop, 0, n2Q-1, 0.00001, &n2Q, Evals2QAr[tid], Ham2QAr[tid], n2Q, isuppz2Q);
    gen_ExDip2Q(ExDip2QAr, Dip2QAr, Ham1QAr, Ham2QAr, traj_param, tid);
  } 
  else {
    // Generate first-order anharmonically perturbed 2Q energies.
    // Remember that Ham1QAr[tid] now contains eigenvectors. 
    gen_perturb_2Q_energies(Ham1QAr[tid], Evals1QAr[tid], Evals2QAr[tid], nosc, delta);
  }
}

// cjfeng 12/10/2018
// Using a single subroutine to combine calc_2dir and calc_2dir_pert
// cjfeng 01/04/2019
// Added traj_param
int calc_2dir(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], GNREAL **ExDip2QAr[3], int tid, traj_param *traj_param, w_param *w_param, GNCOMP *REPH[4], GNCOMP *NREPH[4], POL_info *pol_info, spec_param *spec_param) {
  int i, a, b, c;
  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int ndx, ndx1, ndx3;
  int npol = pol_info->npol;
  int POL[4];
  for(i=0; i<4; i++) POL[i] = pol_info->POL[i];
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;
  int pert = spec_param->pert;
  int npts = w_param->npts;
  GNREAL res = w_param->wres;
  GNREAL start = w_param->wstart;
  GNREAL stop = w_param->wstop;

  GNREAL *Evals1Q = Evals1QAr[tid];
  GNREAL *Evals2Q = Evals2QAr[tid];

  for(a=0; a<nosc; a++) {
    ndx1 = floor( (Evals1Q[a]-start)/res + 0.5 );
    int a2Q = a*n2Q;
    if( (ndx1>0) && (ndx1<npts) ) {
      int ndx1npts = ndx1*npts;
      for(b=0; b<nosc; b++) {
        int b2Q = b*n2Q;
        ndx3 = floor( (Evals1Q[b]-start)/res + 0.5 );
        if( (ndx3>0) && (ndx3<npts) ) {
          if(reph) {
            // Pathway 1 - rephasing
            for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i], pol_info);
            // Pathway 2 - rephasing
            for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, a, b, POL[i], pol_info);
          }
          if(nreph) {
            // Pathway 1 - non-rephasing
            for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i], pol_info);
            // Pathway 2 - non-rephasing
            ndx3 = floor( (Evals1Q[a]-start)/res + 0.5 );
            if( (ndx3>0) && (ndx3<npts) ) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, b, a, POL[i], pol_info);
          }
 
          if(!pert) {
            for(c=0; c<n2Q; c++) {
              if(reph) {
                // Pathway 3 - rephasing
                ndx3 = floor( ( (Evals2Q[c]-Evals1Q[a])-start)/res + 0.5 );
                if( (ndx3>0) && (ndx3<npts) ) for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip2QAr, tid, a, b, b2Q+c, a2Q+c, POL[i], pol_info);
              }
              if(nreph) {
                // Pathway 3 - non-rephasing
                ndx3 = floor( ( (Evals2Q[c]-Evals1Q[b])-start)/res + 0.5 );
                if( (ndx3>0) && (ndx3<npts) ) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip2QAr, tid, a, b, a2Q+c, b2Q+c, POL[i], pol_info);
              }
            }
          }
          else {
            if(a>b) ndx = ((a*(a+1)/2)+b);
            else ndx = ((b*(b+1)/2)+a);
            if (reph) {
              ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[a])-start)/res + 0.5 );
              if( (ndx3>0) && (ndx3<npts) ) {
                // Pathway 3 - rephasing (0,0) -> (0,a) -> (a,a) -> (ab,a)
                for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i], pol_info);
                // Pathway 3 - rephasing (0,0) -> (0,a) -> (b,a) -> (ab,a)
                for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, b, a, b, POL[i], pol_info);
              }
            }
            if (nreph) {
              // Pathway 3 - non-rephasing (0,0) -> (0,a) -> (a,a) -> (a,ab)
              ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[a])-start)/res + 0.5 );
              if( (ndx3>0) && (ndx3<npts) ) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i], pol_info);
              // Pathway 3 - non-rephasing (0,0) -> (0,a) -> (b,a) -> (b,ab)
              ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[b])-start)/res + 0.5 );
              if( (ndx3>0) && (ndx3<npts) ) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, b, b, a, POL[i], pol_info);
            }
          }
        }

      }
    }
  }
  return 1;
}
