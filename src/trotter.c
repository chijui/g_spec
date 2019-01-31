#include "trotter.h"

void initialize_trotter_array(trotter_array *trotter) {
  trotter->nosc = 0;
  trotter->ncoup = 0;
  trotter->nbuffer = 0;
  trotter->D = NULL;
  trotter->J = NULL;
  trotter->row = NULL;
  trotter->col = NULL;
}

// cjfeng 01/23/2019
// Use a single sub-routine to set trotter arrays
int set_trotter_array(trotter_array *trotter, int nosc, int nbuffer, int if2Q) {
  if(if2Q) {
    int n2Q = nosc * (nosc+1) / 2;
    trotter->nosc = n2Q;
    trotter->ncoup = nosc * nosc * (nosc-1) / 2;
  }
  else {
    trotter->nosc = nosc;
    trotter->ncoup = nosc * (nosc-1)/2;
  }
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
    cre = cos(expfac * Ham[N+n]/2);
    cim = sin(expfac * Ham[N+n]/2);
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
      int NDXji = NDX2Q(j,i);
      for(k=j+1; k<nosc; k++) {
        int NDXki = NDX2Q(k,i);
        trotter->row[l] = NDXji;
        trotter->col[l] = NDXki;
        l++;
      }
    }
  }
  return 0;
}

// cjfeng 07/06/2017
// We use row and col vectors to search non-zero coupling values in the Hamiltonian
int construct_trotter_J(trotter_array *trotter, GNREAL *Ham, int xfr, GNREAL expfac) {
  int k;
  for(k=0; k<(trotter->ncoup); k++) {
    int kk = 2 * k;
    int ind = (trotter->row[k]) * (trotter->nosc) + (trotter->col[k]);
    GNREAL Hamval = Ham[ind];
    trotter->J[xfr][kk] = cos(expfac * Hamval);
    trotter->J[xfr][kk+1] = sin(expfac * Hamval);
  }
  return 0;  
}

// cjfeng 01/10/2019
// Swapped xyz and window/nbuffer
int calc_CorrFunc_trotter(GNCOMP *CorrFunc, trotter_array *trotter, GNCOMP ***psi, GNREAL ***Dip1QMat, int window, int fr, int nthreads) {
  int i, j, n;
  int k;
  int xfr;
  GNCOMP cval;
  initialize_psi1Q(psi[0], Dip1QMat[fr], trotter->nosc);
  
  for(n=0; n<window; n++) {
    xfr = (fr+n) % trotter->nbuffer;
    cval = 0.0;
    for(i=0; i<3; i++) {
      for(j=0; j<(trotter->nosc); j++) {
        GNCOMP cval2 = psi[0][i][j];
        cval += Dip1QMat[xfr][i][j] * cval2;
      }
    }
    CorrFunc[n] += cval;

    propagate_psi_trotter(trotter, psi[0], psi[0], 3, xfr, nthreads);
  }
  return 0;
}

int propagate_psi1Q_trotter(GNCOMP ***psi, GNREAL ***Dip1QMat, trotter_array *trotter, int tau0, int fr, int xfrtau, int winsize, int nthreads) {
  int tau, i, j, k;
  // Initial state
  initialize_psi1Q(psi[tau0], Dip1QMat[xfrtau], trotter->nosc);
  for(tau=tau0; tau<tau0+winsize-1; tau++) {
    int xfr = (fr + tau) % trotter->nbuffer;
    propagate_psi_trotter(trotter, psi[tau], psi[tau+1], 3, xfr, nthreads);
  }
  return 0;
}

int propagate_psi2Q_trotter(trotter_array *trotter, GNCOMP **psi_c, int xfr3, int nthreads) {
  int i, j, k;
  propagate_psi_trotter(trotter, psi_c, psi_c, 9, xfr3, nthreads);
  return 0;
}

// cjfeng 01/23/2019
// Use a single sub-routine for all propagations 
void propagate_psi_trotter(trotter_array *trotter, GNCOMP **psi_in, GNCOMP **psi_out, int nelem, int xfr, int nthreads) {
  int i, j, k;
  // Propagation by the diagonal operator
  for(i=0; i<nelem; i++) for(j=0; j<(trotter->nosc); j++) psi_out[i][j] = trotter->D[xfr][j] * psi_in[i][j];
  // Propagation by the off-diagonal operator
  for(k=0; k<trotter->ncoup; k++) {
    int kk = 2 * k;
    GNREAL cosJ, sinJ;
    cosJ = trotter->J[xfr][kk];
    sinJ = trotter->J[xfr][kk+1];
    int i1, i2;
    i1 = trotter->row[k];
    i2 = trotter->col[k];
    for(i=0; i<nelem; i++) {
      GNCOMP c1, c2;
      c1 = psi_out[i][i1];
      c2 = psi_out[i][i2];
      psi_out[i][i1] = cosJ * c1 + I * sinJ * c2;
      psi_out[i][i2] = I * sinJ * c1 + cosJ * c2;
    }
  }
  // Propagation by the diagonal operator
  for(i=0; i<nelem; i++) for(j=0; j<(trotter->nosc); j++) psi_out[i][j] = trotter->D[xfr][j] * psi_out[i][j];
}
