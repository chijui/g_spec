#include "mat_util.h"

// cjfeng 06/26/2017
// Copy matrix A to B
int copy_gncomp(GNCOMP *A, GNCOMP *B, int nrows, int ncols) {
  if ( (nrows * ncols) ) {
    int i, N = nrows * ncols;
    for(i=0; i<N; i++) B[i] = A[i];
  }
  return 0;
}

int mvmult_real_serial_trans(GNREAL *A, GNREAL *vin, GNREAL *vout, int dim, int nt) {
  int i, j;
  if( dim/nt >= THRESHOLD && nt>1) {
    #if OMP_PARALLEL 
    omp_set_num_threads(nt);
    #pragma omp parallel shared(A, vin, vout) private(i, j) firstprivate(dim)
    #pragma omp for schedule(guided) nowait
    #endif
    for(i=0; i<dim; i++) {
      GNREAL vout_private = 0.0;
      for(j=0; j<dim; j++) vout_private += A[j*dim+i]*vin[j];
      vout[i] = vout_private;
    }
  }
  else {
    for(i=0; i<dim; i++) {
      GNREAL vout_private = 0.0;
      for(j=0; j<dim; j++) vout_private += A[j*dim+i]*vin[j];
      vout[i] = vout_private;
    }
  }
  return 0;
}

int mvmult_comp(GNCOMP *A, GNCOMP *vin, GNCOMP *vout, int dim, int nt) {
  int i, j;
  if(dim >= NBCOMP && nt > 1) {
    #if OMP_PARALLEL 
    omp_set_num_threads(nt);
    #pragma omp parallel shared(A, vin, vout) private(i, j) firstprivate(dim)
    #pragma omp for schedule(guided) nowait
    #endif
    for(i=0; i<dim; i++) {
      int ni = i*dim;
      GNCOMP vout_private= 0.0 + 0.0i;
      for(j=0; j<dim; j++) vout_private+= A[ni+j] * vin[j];
      vout[i] = vout_private;
    }
  }
  else {
    for(i=0; i<dim; i++) {
      int ni = i*dim;
      GNCOMP vout_private= 0.0 + 0.0i;
      for(j=0; j<dim; j++) vout_private += A[ni+j] * vin[j];
      vout[i] = vout_private;
    }
  }
  return 0;
}

// Matrix/vector multiplication for a non-square real matrix and complex vector. 
// Note: Actually the transpose of A. Dip2Q has dimension (nosc)x(n2Q). We want
// to multiply the transpose of that matrix with the input vector (a wavefunction
// of length nosc).

// cjfeng 03/27/2016
// Made change to use the transpose matrix directly, but effectively doing the same 
// multiplication as before while increasing cache memory access count.
int mvmult_comp_trans_x(GNREAL *A_t, GNCOMP *vin, GNCOMP *vout, int dimin, int dimout, int nt) {
  int i, j;
  if( dimout >= NBREAL && nt>1) {
    #if OMP_PARALLEL 
    #pragma omp parallel shared(A_t, vin, vout) private(i, j) firstprivate(dimin, dimout)
    #pragma omp for schedule(guided) nowait
    #endif
    for(i=0; i<dimout; i++) {
      int ni = i*dimin;
      GNCOMP vout_private = 0.0 +0.0i;
      for(j=0; j<dimin; j++) vout_private += A_t[ni+j]*vin[j];
      vout[i] = vout_private;
    }
  }
  else {
    for(i=0; i<dimout; i++) {
      int ni = i*dimin;
      GNCOMP vout_private = 0.0 + 0.0i;
      for(j=0; j<dimin; j++) vout_private += A_t[ni+j]*vin[j];
      vout[i] = vout_private;
    }
  }
  return 0;
}

// cjfeng 12/19/2018
// Using sparse matrix
int mvmult_comp_trans_Sp(int *A_row, int *A_col, GNREAL *A_Sp, GNCOMP *vin, GNCOMP *vout, int dimin, int dimout, int nt) {
  int i, j;
  int diminsq = dimin * dimin;
  for(j=0; j<dimout; j++) vout[j] = 0.0 + 0.0i;
  for(i=0; i<diminsq; i++) vout[A_row[i]] += A_Sp[i] * vin[A_col[i]];
  return 0;
}
