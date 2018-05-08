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

// cjfeng 06/26/2017
// initialization
int zero_gncomp_mat(int nrows, int ncols, GNCOMP *Mat) {
	if( (nrows * ncols) ) {
		int i, N = nrows * ncols;
		for(i=0; i<N; i++) Mat[i] = 0.0;
	}
	return 1;
}

int gncomp_identity_mat(int N, GNCOMP *Mat) {
	if ( N ) {
		int i, Nsq = N * N;
		for(i=0; i<Nsq; i++) Mat[i] = 0.0;
		for(i=0; i<N; i++) Mat[i*N+i] = 1.0;
	}
	return 0;
}

int trans_comp(GNCOMP *A, GNCOMP *A_t, int nrows, int ncols, int nt) {
	int i,j;
	if (nrows >= NBREAL && nt>1) {
		#if OMP_PARALLEL
		#pragma omp parallel shared(A, A_t) private(i,j) firstprivate(nrows, ncols)
		#pragma omp for schedule(guided) nowait
		#endif
		for(i=0; i<nrows; i++) {
			int ni = i*ncols;
			for(j=0; j<ncols; j++) A_t[j*nrows+i] = A[ni+j];
		}
	}
	else {
		for(i=0; i<nrows; i++) {
			int ni = i*ncols;
			for(j=0; j<ncols; j++) A_t[j*nrows+i] = A[ni+j];
		}
	}
	return 0;
}

int trans_real(GNREAL *A, GNREAL *A_t, int nrows, int ncols, int nt) {
	int i,j;
	if (nrows >= NBREAL && nt>1) {
		#if OMP_PARALLEL
		omp_set_num_threads(nt);
		#pragma omp parallel shared(A, A_t) private(i,j) firstprivate(nrows, ncols)
		#pragma omp for schedule(guided) nowait
		#endif
		for(i=0; i<nrows; i++) {
			int ni = i * ncols;
			for(j=0; j<ncols; j++) A_t[j*nrows+i] = A[ni+j];
		}
	}
	else {
		for(i=0; i<nrows; i++) {
			int ni = i * ncols;
			for(j=0; j<ncols; j++) A_t[j*nrows+i] = A[ni+j];
		}
	}
	return 0;
}

// cjfeng 03/23/2016
// This function utilizes the transpose of B to make the original multiplication 
// cache friendly. It is equivalent with C=A*B, but now the matrices are running
// row-major.
int mmult_comp_trans(GNCOMP *A, GNCOMP *B_t, GNCOMP *C, int dim, int nt) {
	int i, j, k;
	if (dim/nt >= THRESHOLD && nt>1) {
		#if OMP_PARALLEL 
		omp_set_num_threads(nt);
		#pragma omp parallel shared(A, B_t, C) private(i, j, k) firstprivate(dim)
		#pragma omp for schedule(guided) nowait //collapse(2) 
		#endif
		for(i=0; i<dim; i++) {
			int ni = i*dim;
			for(j=0; j<dim; j++) {
				int nj = j*dim;
				GNCOMP C_private = 0.0 + 0.0i;
				for(k=0; k<dim; k++) C_private += A[ni+k]*B_t[nj+k];
				C[ni+j] = C_private;
			}
		} 
	}
	else {
		for(i=0; i<dim; i++) {
			int ni = i*dim;
			for(j=0; j<dim; j++) {
				int nj = j*dim;
				GNCOMP C_private = 0.0 + 0.0i;
				for(k=0; k<dim; k++) C_private += A[ni+k]*B_t[nj+k];
				C[ni+j] = C_private;
			}
		}
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
