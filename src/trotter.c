#include "trotter.h"

int set_trotter_1Qarray(trotter_array *trotter, int nosc, int nbuffer) {
  trotter->nosc = nosc;
  trotter->ncoup = nosc * (nosc-1)/2;
  trotter->nbuffer = nbuffer;
  
  set_gncomp_2d_array(&(trotter->D), trotter->nbuffer, trotter->nosc);
  if (trotter->ncoup != 0) {
    set_gnreal_2d_array(&(trotter->J), trotter->nbuffer, 2*(trotter->ncoup));
    set_int_1d_array(&(trotter->row), trotter->ncoup, 1);
    set_int_1d_array(&(trotter->col), trotter->ncoup, 1);
  }
  return 0;
}

// cjfeng 07/06/2017
// We only take non-zero coupling into account.
int set_trotter_2Qarray(trotter_array *trotter, int nosc, int nbuffer) {
  int n2Q = nosc * (nosc+1) / 2;
  trotter->nosc = n2Q;
  trotter->ncoup = nosc * nosc * (nosc-1) / 2;
  trotter->nbuffer = nbuffer;
  
  set_gncomp_2d_array(&(trotter->D), trotter->nbuffer, trotter->nosc);
  if (trotter->ncoup != 0) {
    set_gnreal_2d_array(&(trotter->J), trotter->nbuffer, 2*(trotter->ncoup));
    set_int_1d_array(&(trotter->row), trotter->ncoup, 1);
    set_int_1d_array(&(trotter->col), trotter->ncoup, 1);
  }
  return 0;
}

int unset_trotter_array(trotter_array *trotter) {
  unset_gncomp_2d_array(trotter->D, trotter->nbuffer);
  if (trotter->ncoup != 0) {
    unset_gnreal_2d_array(trotter->J, trotter->nbuffer);
    unset_int_1d_array(trotter->row);
    unset_int_1d_array(trotter->col);
  }
  return 0;
}

int construct_trotter_mat(trotter_array *trotter, GNREAL *Ham, int xfr, double expfac) {
  construct_trotter_D(trotter, Ham, xfr, expfac);
  if (trotter->ncoup != 0) construct_trotter_J(trotter, Ham, xfr, expfac);
  return 0;
}

int construct_trotter_D(trotter_array *trotter, GNREAL *Ham, int xfr, GNREAL expfac) {
  int n;
  for(n=0; n<(trotter->nosc); n++) {
    int N = trotter->nosc * n;
    GNREAL cre, cim;
    cre = cos(expfac * Ham[N+n]);
    cim = sin(expfac * Ham[N+n]);
    trotter->D[xfr][n] = cre + I * cim;
  }
  return 0;
}

int construct_trotter_index(trotter_array *trotter) {
  int n, m, k = 0;
  for(n=0; n<(trotter->nosc); n++) {
    int N = n * (trotter->nosc);
    for(m=n+1; m<(trotter->nosc); m++) {
      int kk = 2 * k;
        trotter->row[k] = n;
        trotter->col[k] = m;
        k++;
    }
  }
  return 0;
}

int construct_trotter_2Qindex(trotter_array *trotter, int nosc) {
  int i, j, k, l=0;
  // Generating index map for constructing nonzero J array.
  // The index of state mn is (m+1)*m/2 + n.
  // We only count the upper diagonal part of the Hamiltonian.
  for(i=0; i<nosc; i++) {
    int NDXii = NDX2Q(i,i);
    // Coupling between ii and ij states: <ii|H|ij>
    for(j=0; j<i; j++) {
      int NDXij = NDX2Q(i,j);
      trotter->row[l] = NDXii;
      trotter->col[l] = NDXij;
      l++;
    }
    // Coupling between ii and ji states: <ii|H|ji>
    for(j=i+1; j<nosc; j++) {
      int NDXji = NDX2Q(j,i);
      trotter->row[l] = NDXii;
      trotter->col[l] = NDXji;
      l++;
    }
     // Coupling between ij and ik states (j<k<i)
    for(j=0; j<i-1; j++) {
      int NDXij = NDX2Q(i,j);
      for(k=j+1; k<i; k++) {
        int NDXik = NDX2Q(i,k);
        trotter->row[l] = NDXij;
        trotter->col[l] = NDXik;
        l++;
      }
    }
    // Coupling between ki and ij states (j<i<k)
    for(j=0; j<i; j++) {
      int NDXij = NDX2Q(i,j);
      for(k=i+1; k<nosc; k++) {
        int NDXki = NDX2Q(k,i);
        trotter->row[l] = NDXij;
        trotter->col[l] = NDXki;
        l++;
      }
    }
    // Coupling between ki and ji states (i<j<k)
    for(j=i+1; j<nosc; j++) {
      int NDXij = NDX2Q(i,j);
      for(k=j+1; k<nosc; k++) {
        int NDXki = NDX2Q(k,i);
        trotter->row[l] = NDXij;
        trotter->col[l] = NDXki;
        l++;
      }
    }
  }
  // printf("l = %d\t, nosc*nosc*(nosc-1)/2 = %d\n", l, nosc*nosc*(nosc-1)/2);
  return 0;
}

// cjfeng 07/06/2017
// We use row and col vectors to search non-zero coupling values in the Hamiltonian
int construct_trotter_J(trotter_array *trotter, GNREAL *Ham, int xfr, GNREAL expfac) {
  int k;
  // printf("trotter->ncoup = %d\n", trotter->ncoup);
  for(k=0; k<(trotter->ncoup); k++) {
    int kk = 2 * k;
    int ind = (trotter->row[k]) * (trotter->nosc) + (trotter->col[k]);
    GNREAL Hamval = Ham[ind];
    // printf("Hamval = %f\n", Hamval);
    trotter->J[xfr][kk] = cos(expfac * Hamval);
    trotter->J[xfr][kk+1] = sin(expfac * Hamval);
  }
  return 0;  
}

int calc_CorrFunc_trotter(GNCOMP *CorrFunc, trotter_array *trotter, GNCOMP **psi[3], GNREAL ***Dip1QMat, int window, int fr, int nthreads) {
  int i, j, n;
  int k;
  int xfr;
  GNCOMP cval;
  initialize_psi(psi, Dip1QMat, 0, fr, trotter->nosc, 3);
  
  for(n=0; n<window; n++) {
    xfr = (fr+n) % trotter->nbuffer;
    cval = 0.0;
    #if OMP_PARALLEL
    #pragma omp parallel if(trotter->nosc >= NBCOMP && nthreads>1) shared(trotter, Dip1QMat, xfr, fr) private(i,j,k)
    #pragma omp for schedule(guided) collapse(2) reduction(+:cval)
    #endif
    for(i=0; i<3; i++) {
      for(j=0; j<(trotter->nosc); j++) {
        GNCOMP cval2 = psi[i][0][j];
        cval += Dip1QMat[i][xfr][j] * cval2;
      }
    }
    CorrFunc[n] += cval;

    // cjfeng 06/28/2017
    int i1, i2;
    GNCOMP c1, c2;
    // First operation on the diagonal
    for(i=0; i<3; i++) for(j=0; j<(trotter->nosc); j++) psi[i][0][j] = trotter->D[xfr][j] * psi[i][0][j];
    // Off-diagonal terms
    for(k=0; k<(trotter->ncoup); k++) {
      int kk = 2 * k;
      GNREAL cosJ, sinJ;
      cosJ = trotter->J[xfr][kk];
      sinJ = trotter->J[xfr][kk+1];
      int i1, i2;
      i1 = trotter->row[k];
      i2 = trotter->col[k];
      for(i=0; i<3; i++) {
        GNCOMP c1, c2;
        c1 = psi[i][0][i1];
        c2 = psi[i][0][i2];
        psi[i][0][i1] = cosJ * c1 + I * sinJ * c2;
        psi[i][0][i2] = I * sinJ * c1 + cosJ * c2;
      }
    }
  }
  return 0;
}

int propagate_psi1Q_trotter(GNCOMP **psi[3], GNREAL **Dip1QMat[3], trotter_array *trotter, int tau0, int fr, int xfrtau, int winsize, int nthreads) {
  int tau, i, j, k;
  // Initial state
  initialize_psi(psi, Dip1QMat, tau0, xfrtau, trotter->nosc, 3);
  for(tau=tau0; tau<tau0+winsize-1; tau++) {
    int xfr = (fr + tau ) % trotter->nbuffer;
    // First operation on the diagonal
    for(i=0; i<3; i++) for(j=0; j<(trotter->nosc); j++) psi[i][tau+1][j] = trotter->D[xfr][j] * psi[i][tau][j];
    // Off-diagonal terms
    for(k=0; k<(trotter->ncoup); k++) {
      int kk = 2 * k;
      GNREAL cosJ, sinJ;
      cosJ = trotter->J[xfr][kk];
      sinJ = trotter->J[xfr][kk+1];
      int i1, i2;
      i1 = trotter->row[k];
      i2 = trotter->col[k];
      for(i=0; i<3; i++) {
        GNCOMP c1, c2;
        c1 = psi[i][tau+1][i1];
        c2 = psi[i][tau+1][i2];
        psi[i][tau+1][i1] = cosJ * c1 + I * sinJ * c2;
        psi[i][tau+1][i2] = I * sinJ * c1 + cosJ * c2;
      }
    }
  }
  return 0;
}

int propagate_psi2Q_trotter(trotter_array *trotter, GNCOMP **psi_c, GNCOMP *psi2Q, int xfr3, int nthreads) {
  int i, j, k;
  // Initialize psi_ca/psi_cb and propagate it by Trotter propagator
  // Propagation by the diagonal operator
  for(i=0; i<9; i++) for(j=0; j<(trotter->nosc); j++) psi_c[i][j] = trotter->D[xfr3][j] * psi_c[i][j];
  // Propagation by the off-diagonal operator
  for(k=0; k<(trotter->ncoup); k++) {
    int kk = 2 * k;
    GNREAL cosJ, sinJ;
    cosJ = trotter->J[xfr3][kk];
    sinJ = trotter->J[xfr3][kk+1];
    int i1, i2;
    i1 = trotter->row[k];
    i2 = trotter->col[k];
    for(i=0; i<9; i++) {
      GNCOMP c1, c2;
      c1 = psi_c[i][i1];
      c2 = psi_c[i][i2];
      psi_c[i][i1] = cosJ * c1 + I * sinJ * c2;
      psi_c[i][i2] = I * sinJ * c1 + cosJ * c2;
    }
  }
  return 0;
}
