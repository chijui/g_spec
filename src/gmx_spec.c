#include "fftw3.h"			// fftw for Fourier Transform

#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "futil.h"
#include "tpxio.h"
#include "physics.h"		// Physical constants used in Gromacs
#include "macros.h"
#include "gmx_fatal.h"
#include "index.h"
#include "smalloc.h"
#include "vec.h"
// #include "xvgr.h"
#include "gstat.h"
#include "matio.h"
#include "string2.h"
// #include "pbc.h"

#include "lapacke.h"		// LAPACKE as C interface for solving eigen problems.
#include "complex.h"		// Complex utilities

// cjfeng 07/14/2016
// Tracking execution time and parallel computation
#if OMP_PARALLEL 
	#include "omp.h"
#else
	#include "time.h"
#endif

// cjfeng 06/27/2016
// Tracking memory usage of g_spec
#include <sys/resource.h>

#define THRESHOLD 16
// #define DOUBLE_PRECISION

#if DOUBLE_PRECISION
// #ifdef DOUBLE_PRECISION
	#define GNREAL double
	#define GNCOMP double complex
	#define SYEVR LAPACKE_dsyevr
	#define NBREAL 900
	#define NBCOMP 650
#else
	#define GNREAL float
	#define GNCOMP float complex
	#define SYEVR LAPACKE_ssyevr
	#define NBREAL 1300
	#define NBCOMP 900
#endif

#define NDX2Q(m,n) ( ((m*(m+1)/2)+n) )
#ifndef MIN
	#define MIN(X,Y) ((X) < (Y) ? (X) : (Y)),
#endif

int pertvec = 0;							// first-order perturbative correction on both site energies and vectors

const int maxchar = 1024;

int POL[4];
int* SHIFTNDX = NULL;
real* SHIFT = NULL;

// File pointers
FILE *hfp, *Dfp[3], *sfp, *afp, *ffp, *lfp, *rfp[4], *nrfp[4], *ifp;
// Trajfp is an array of trajectory file pointers: time-stamp, FTIR, rzzzz, rzzyy, rzyyz, rzyzy, nrzzzz, nrzzyy, nrzyyz, nrzyzy
FILE *Trajfp[10];

int *NNums, *CNums;
char **NNames, **CNames;

// cjfeng 07/12/2016
// Time domain responses
FILE *fid_reph[4];
FILE *fid_nreph[4];
char fidbase[8][50] = { "tresponse_reph_zzzz.txt", "tresponse_reph_zzyy.txt", "tresponse_reph_zyyz.txt", "tresponse_reph_zyzy.txt", "treponse_nreph_zzzz.txt", "tresponse_nreph_zzyy.txt", "tresponse_nreph_zyyz.txt", "tresponse_nreph_zyzy.txt" };
char fidname[8][100];

// Now spectra
GNCOMP *CorrFunc;
GNCOMP *NetCorrFunc;
GNREAL *popdecay1Q;
GNREAL *popdecay2Q;
GNREAL *hann;
GNREAL *ftir;
GNREAL *netftir;
fftw_complex *FTin1D, *FTout1D;
fftw_plan FTplan1D;
GNCOMP *REPH[4], *NREPH[4];
GNCOMP *NetREPH[4], *NetNREPH[4];
fftw_complex *FTin2D, *FTout2D;
fftw_plan FTplan2D;

GNCOMP *U1Q;
GNCOMP *U1Qs;
GNCOMP *cDip1Q[3];

// Generic wavefunctions
GNCOMP **psi_a[3];
GNCOMP **psi_b1[3];
GNCOMP **psi_b12[3];
GNCOMP *psi_cb[9];
GNCOMP *psi_ca[9];
GNCOMP *psi2Q;

GNREAL **Ham1QMat;
GNREAL **Dip1QMat[3];
GNCOMP **U1QMat;

GNREAL **Ham1QAr;
GNREAL **Evals1QAr;
GNREAL **Dip1QAr[3];
GNREAL **ExDip1QAr[3];

GNCOMP *U2Qs;

// cjfeng 04/05/2016
// Swap the order of Dip2QMat to be a
// nbuffer*3*n2Q*nosc array.
GNREAL ***Dip2QMat;
GNCOMP **U2QMat;

GNREAL **Ham2QAr;
GNREAL **Evals2QAr;
GNREAL **Dip2QAr[3];
GNREAL **ExDip2QAr[3];

GNREAL *SitesBuffer;

#define f15 ( (1.0)/(5.0) )
#define f115 ( (1.0)/(15.0) )
#define f215 ( (2.0)/(15.0) )
#define f130 (-(1.0)/(30.0) )

// 				zzzz	zzyy	zyyz	zyzy
GNREAL M_ijkl_IJKL[4][4] = { {	f15,	f115,	f115,	f115 }, 	// ZZZZ
			     {	f115,	f215,	f130,	f130 },		// ZZYY
			     {	f115,	f130,	f215,	f130 },		// ZYYZ
			     {	f115,	f130,	f130, f215 } };	// ZYZY


/********************************************************************************
* File IO									*
********************************************************************************/

int count_lines(char* fname) {
	int lines=0,length=0;
	char temp, *c ;
	FILE* fp = fopen(fname, "r");
	if(fp==NULL) {
		printf("Error opening file %s: please check input", fname);
		return -1;
	}
	do {
		temp = fgetc(fp);
		length++;
	} while( temp!='\n' );
	lines++;
	c= (char *) malloc(2*length);		

	while(fgets(c,2*length,fp) != NULL){
		lines++;
	}
	fclose(fp);
	free(c);
	return lines;
}

// MER 04/06/2016
// Modified string tokenizer to include additional tokens beyond '\t'. 
// Also set num to start at zero rather than -1. Previously files whose
// lines did not include a terminating \t character were undercounted. 
int count_entries(char* fname) {
	int num = 0, length=0;
	char *c, temp, *ch,*d="\t, \n";
	FILE* fp = fopen(fname, "r");
	if(fp==NULL) {
		printf("Error openting file %s: please check input", fname);
		return -1;
	}
	do{
		temp = fgetc(fp);
		length++;
	} while( temp!='\n' );
	c = (char *) malloc (2*length);
	fgets(c,2*length,fp);
	ch=strtok(c,d);
	while(ch){
		ch=strtok(NULL,d); 
		num++;
	}
	fclose(fp);
	free(c);
	return num; 
}

int read_line(FILE* fp, int nrows, int ncols, GNREAL *Mat) {
	int i,n=nrows*ncols;
	real val;
	for(i=0; i<n; i++) {
			if(fscanf(fp, "%f%*[ ,]", &val)==0) return 0;
			else Mat[i] = val;
	}
	return 1;
}
/********************\
 * Matrix utilities *
\********************/ 

// cjfeng 04/06/2016
// Added a complex version of print_mat
int print_mat_comp(int nrows, int ncols, GNCOMP* Mat) {
	int i,j;
	for(i=0; i<nrows; i++) {
		for(j=0; j<ncols; j++) {
			printf("\t%6.8f + %6.8fi", creal(Mat[i*ncols+j]), cimag(Mat[i*ncols+j]));
		}
		printf("\n");
	}
	return 1;
}

int print_mat(int nrows, int ncols, GNREAL* Mat) {
	int i,j;
	for(i=0; i<nrows; i++) {
		for(j=0; j<ncols; j++) {
			printf("\t%6.8f", Mat[i*ncols+j]);
		}
		printf("\n");
	}
	return 1;
}

int trans_comp(GNCOMP *A, GNCOMP *A_t, int nrows, int ncols, int nt) {
	int i,j;
	if (nrows >= NBREAL && nt>1) {
		#if OMP_PARALLEL
		omp_set_num_threads(nt);
		#pragma omp parallel shared(A, A_t) private(i,j) firstprivate(nrows, ncols)
		#pragma omp for schedule(guided) nowait
		#endif
		for(i=0; i<nrows; i++) {
			int ni = i*ncols;
			for(j=0; j<ncols; j++){
				A_t[j*nrows+i] = A[ni+j];
			}
		}
	}
	else {
		for(i=0; i<nrows; i++) {
			int ni = i*ncols;
			for(j=0; j<ncols; j++) {
				A_t[j*nrows+i] = A[ni+j];
			}
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
			int ni = i*ncols;
			for(j=0; j<ncols; j++){
				A_t[j*nrows+i] = A[ni+j];
			}
		}
	}
	else {
		for(i=0; i<nrows; i++) {
			int ni = i*ncols;
			for(j=0; j<ncols; j++) {
				A_t[j*nrows+i] = A[ni+j];
			}
		}
	}
	return 0;
}

// cjfeng 03/27/2016
// Loop tiling or blocking:
// This function partitions all of the matrices into smaller blocks to
// make the whole blocks fit into the cache size (accessible in /proc/cpusize).
// Current system has the cache size 20480 kb. Hence for float complex case,
// it can roughly contain 3 * 934-by-934 blocks while it can only accommodate
// 3 * 660-by-660 blocks for double complex case. In addition, the matrix multiplication 
// are treated to be row-major for both A and B (indeed B transpose).
int mmult_comp_block (GNCOMP *A, GNCOMP *B_t, GNCOMP *C, int dim, int nt) {
	int i,j,k,jj,kk;
	if(dim/nt >= THRESHOLD && nt>1) {	// Parallel calculation
		#if OMP_PARALLEL 
		omp_set_num_threads(nt);
		#endif
		for (jj=0; jj<dim; jj+=NBCOMP) {	// NBCOMP: size of a block
			for (kk=0; kk<dim; kk+=NBCOMP) {
				#if OMP_PARALLEL 
				#pragma omp parallel shared(A, B_t, C) private(i,j,k) firstprivate(dim, jj, kk)
				#pragma omp for schedule(guided) nowait
				#endif
				for (i=0; i<dim; i++) {
					int ni = i*dim;
					for (j=jj; j< min(jj+NBCOMP, dim); j++) {
						GNCOMP C_private = 0.0;
						int nj = j*dim;
						for (k=kk; k<min(kk+NBCOMP, dim); k++) {
							C_private += A[ni+k] * B_t[nj+k];
						}
						C[ni+j] = C_private;
					}
				}
			}
		}
	}
	else {	// Serial calculation
		for (jj=0; jj<dim; jj+=NBCOMP) {
			for (kk=0; kk<dim; kk+=NBCOMP) {
				for (i=0; i<dim; i++) {
					int ni = i*dim;
					for (j=jj; j< min(jj+NBCOMP, dim); j++) {
						int nj = j*dim;
						GNCOMP C_private = 0.0;
						for (k=kk; k<min(kk+NBCOMP, dim); k++) {
							C_private += A[ni+k] * B_t[nj+k];
						}
						C[ni+j] = C_private;
					}
				}
			}
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
				GNCOMP C_private = 0.0 + 0.0i;
				int nj = j*dim;
				for(k=0; k<dim; k++) {
					C_private += A[ni+k]*B_t[nj+k];
				}
				C[ni+j] = C_private;
			}
		} 
	}
	else {
		for(i=0; i<dim; i++) {
			int ni = i*dim;
			for(j=0; j<dim; j++) {
				GNCOMP C_private = 0.0 + 0.0i;
				int nj = j*dim;
				for(k=0; k<dim; k++) {
					C_private += A[ni+k]*B_t[nj+k];
				}
				C[ni+j] = C_private;
			}
		}
	}
	return 0;
}

int mmult_comp(GNCOMP *A, GNCOMP *B, GNCOMP *C, int dim, int nt) {
	int i, j, k;
	// cjfeng 07/13/2016
	if (dim/nt >= THRESHOLD && nt > 1) {
		#if OMP_PARALLEL 
		omp_set_num_threads(nt);
		#pragma omp parallel shared(A, B, C) private(i, j, k) firstprivate(dim)
		#pragma omp for schedule(guided) nowait
		#endif
		for(i=0; i<dim; i++) {
			for(j=0; j<dim; j++) {
				GNCOMP C_private = 0.0 + 0.0i;
				for(k=0; k<dim; k++) {
					C_private += A[i*dim+k]*B[k*dim+j];
				}
				C[i*dim+j] = C_private;
			}
		}
	}
	else {
		for(i=0; i<dim; i++) {
			for(j=0; j<dim; j++) {
				GNCOMP C_private = 0.0 + 0.0i;
				for(k=0; k<dim; k++) {
					C_private += A[i*dim+k]*B[k*dim+j];
				}
				C[i*dim+j] = C_private;
			}
		}
	}
	return 0;
}

int mmult_comp_adjointB(GNCOMP *A, GNCOMP *B, GNCOMP *C, int dim, int nt) {
	int i, j, k;
	if( dim/nt >= THRESHOLD && nt>1) {
		#if OMP_PARALLEL 
		omp_set_num_threads(nt);
		#pragma omp parallel shared(A, B, C) private(i, j, k) firstprivate(dim)
		#pragma omp for schedule(guided) nowait
		#endif
		for(i=0; i<dim; i++) {
			int ni = i*dim;
			for(j=0; j<dim; j++) {
				GNCOMP C_private = 0.0 + 0.0i;
				int nj = j*dim;
				for(k=0; k<dim; k++) {
					C_private += A[ni+k]*conj(B[nj+k]);
				}
				C[ni+j] = C_private;
			}
		}
	}
	else {
		for(i=0; i<dim; i++) {
			int ni = i*dim;
			for(j=0; j<dim; j++) {
				GNCOMP C_private = 0.0 + 0.0i;
				int nj = j*dim;
				for(k=0; k<dim; k++) {
					C_private += A[ni+k]*conj(B[nj+k]);
				}
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
			for(j=0; j<dim; j++) {
				vout_private += A[j*dim+i]*vin[j];
			}
			vout[i] = vout_private;
		}
	}
	else {
		for(i=0; i<dim; i++) {
			GNREAL vout_private = 0.0;
			for(j=0; j<dim; j++) {
				vout_private += A[j*dim+i]*vin[j];
			}
			vout[i] = vout_private;
		}
	}
	return 0;
}

// cjfeng 06/27/2016 cjfeng
// Block optimizied version of matrix-vector multiplication
// Putting omp into inner loop to avoid race condition.
int mvmult_comp_block (GNCOMP *A, GNCOMP *vin, GNCOMP *vout, int dimin, int dimout, int nt) {
	int i,j,ii,jj;
	if( dimout >= NBCOMP && nt>1) {	// Parallel calculation
		#if OMP_PARALLEL 
		omp_set_num_threads(nt);
		#endif
		// for (ii=0; ii<dimout; ii+=NBCOMP) {
		// #pragma omp parallel shared(A, vin, vout, ii, jj) private(i, j) firstprivate(dimin, dimout)
		// #pragma omp for schedule(guided) nowait collapse(2)
			for (jj=0; jj<dimin; jj+=NBCOMP) {
				#if OMP_PARALLEL 
				#pragma omp parallel shared(A, vin, vout, ii, jj) private(i, j) firstprivate(dimin, dimout)
				#pragma omp for schedule(guided) nowait 
				#endif
				for (i=0; i< dimout; i++) {
				// for (i=0; i< min(ii+nb,dimout); i++) {
					GNCOMP vout_private = 0.0 + 0.0i;
					int ni = i*dimin;
					for (j=jj; j< min(jj+NBREAL, dimin); j++) {
						vout_private += A[ni+j] * vin[j];
					}
					vout[i] = vout_private;
				}
			}
		// }
	}
	else {		// Serial calculation
		// for (ii=0; ii<dimout; ii+=NBREAL) {
			for (jj=0; jj<dimin; jj+=NBCOMP) {
				for (i=0; i<dimout; i++) {
				// for (i=ii; i< min(ii+NBREAL,dimout); i++) {
					GNCOMP vout_private = 0.0 + 0.0i;
					int ni = i*dimin;
					for (j=jj; j< min(jj+NBCOMP, dimin); j++) {
						vout_private += A[ni+j] * vin[j];
					}
					vout[i] = vout_private;
				}
			}
		// }
	} 
	return 1;
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
			GNCOMP vout_private= 0.0 + 0.0i;
			int ni = i*dim;
			for(j=0; j<dim; j++) {
				vout_private+= A[ni+j]*vin[j];
			}
			vout[i] = vout_private;
		}
	}
	else {
		for(i=0; i<dim; i++) {
			GNCOMP vout_private= 0.0 + 0.0i;
			int ni = i*dim;
			for(j=0; j<dim; j++) {
				vout_private += A[ni+j]*vin[j];
			}
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
// multiplication as before while increasing cache memory access count. In addition,
// the block optimization has been applied here. 
int mvmult_mix_block (GNREAL *A, GNCOMP *vin, GNCOMP *vout, int dimin, int dimout, int nt) {
	int i,j,ii,jj;
	if( dimout >= NBREAL && nt > 1 ) {	// Parallel calculation
		#if OMP_PARALLEL 
		omp_set_num_threads(nt);
		#endif
		// for (ii=0; ii<dimout; ii+=NBREAL) {
		// Recognizing that the only matrices involved here is indeed real.
			for (jj=0; jj<dimin; jj+=NBREAL) {	
				#if OMP_PARALLEL 
				#pragma omp parallel shared(A, vin, vout) private(i, j) firstprivate(dimin, dimout, ii, jj)
				#pragma omp for schedule(guided) nowait
				#endif
				for (i=0; i< dimout; i++) {
				// for (i=0; i<min(ii+NBREAL,dimout); i++) {
					GNCOMP vout_private = 0.0 + 0.0i;
					int ni = i*dimin;
					for (j=jj; j<min(jj+NBREAL, dimin); j++) {
						vout_private += A[ni+j] * vin[j];
					}
				vout[i] = vout_private;
				}
			}
		// }
	}
	else {		// Serial calculation
		// for (ii=0; ii<dimout; ii+=NBREAL) {
			for (jj=0; jj<dimin; jj+=NBREAL) {
				for (i=0; i<dimout; i++) {
				// for (i=ii; i<min(ii+NBREAL,dimout); i++) {
					GNCOMP vout_private = 0.0 + 0.0i;
					int ni = i*dimin;
					for (j=jj; j<min(jj+NBREAL, dimin); j++) {
						vout_private += A[ni+j] * vin[j];
					}
					vout[i] = vout_private;
				}
			}
		// }
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
		omp_set_num_threads(nt);
		#pragma omp parallel shared(A_t, vin, vout) private(i, j) firstprivate(dimin, dimout)
		#pragma omp for schedule(guided) nowait
		#endif
		for(i=0; i<dimout; i++) {
			int ni = i *dimin;
			GNCOMP vout_private = 0.0 +0.0i;
			for(j=0; j<dimin; j++) {
				vout_private += A_t[ni+j]*vin[j];
			}
			vout[i] = vout_private;
		}
	}
	else {
		for(i=0; i<dimout; i++) {
			GNCOMP vout_private = 0.0 + 0.0i;
			int ni = i*dimin;
			for(j=0; j<dimin; j++) {
				vout_private += A_t[ni+j]*vin[j];
			}
			vout[i] = vout_private;
		}
	}
	return 0;
}

int mvmult_comp_x(GNREAL *A, GNCOMP *vin, GNCOMP *vout, int dimin, int dimout, int nt) {
	int i, j;
	if(nt>1) {
		#if OMP_PARALLEL 
		omp_set_num_threads(nt);
		#pragma omp parallel shared(A, vin, vout) private(i, j) firstprivate(dimin, dimout)
		#pragma omp for schedule(guided) nowait
		#endif
		for(i=0; i<dimout; i++) {
			GNCOMP vout_private = 0.0 +0.0i;
			for(j=0; j<dimin; j++) {
				vout_private += A[j*dimout+i]*vin[j];
			}
			vout[i] = vout_private;
		}
	}
	else {
		for(i=0; i<dimout; i++) {
			GNCOMP vout_private = 0.0 + 0.0i;
			for(j=0; j<dimin; j++) {
				vout_private += A[j*dimout+i]*vin[j];
			}
			vout[i] = vout_private;
		}
	}
	return 0;
}

int read_info( FILE *fp, char* Infoname, char*** p_NNames, int** p_NNums, char*** p_CNames, int** p_CNums) {
	
	char line[1024];
	int i,j;
	int nchar = 20;
	char **NNames = NULL;
	char **CNames = NULL; 
	int *NNums = NULL;
	int *CNums = NULL;
	int nbonds = 0;
	int error = 0; 
	fp = fopen(Infoname, "r");
	if(fp==NULL) return 0;
	else {
		while( (!error) && (fgets(line, maxchar, fp)!=NULL)) {
			if(strncmp(line, "BONDS:", 6)==0) {
				if(sscanf(line, "%*s %d", &nbonds)!=1) {
					printf("Error reading number of bonds from info file. \n");
					error = 1;
				} else {
					if(!error) NNames = (char**) malloc(nbonds*sizeof(char*));
					if(NNames==NULL) error = 1; 
					if(!error) CNames = (char**) malloc(nbonds*sizeof(char*));
					if(CNames==NULL) error = 1;
					if(!error) NNums = (int*) malloc(nbonds*sizeof(int));
					if(NNums==NULL) error = 1;
					if(!error) CNums = (int*) malloc(nbonds*sizeof(int));
					for(i=0; i<nbonds; i++) {
						if(!error) {
							NNames[i] = (char*) malloc(nchar*sizeof(char));
							if(NNames[i]==NULL) {
								for(j=0; j<i; j++) free(NNames[j]);
								error = 1;
							}
						} else NNames[i] = NULL;
					}
					for(i=0; i<nbonds; i++) {
						if(!error) {
							CNames[i] = (char*) malloc(nchar*sizeof(char));
							if(CNames[i]==NULL) {
								for(j=0; j<nbonds; j++) free(NNames[j]);
								for(j=0; j<i; j++) free(CNames[j]);
								error = 1;
							} 
						} else CNames[i] = NULL; 
					}
				}
				for(i=0; i<nbonds; i++) {
					if( (fgets(line, maxchar, fp)==NULL) || (sscanf(line, "%s %d %s %d", NNames[i], &NNums[i], CNames[i], &CNums[i])!=4) ) {
						error = 1;
						break;
					}
				}
			}
		}
		if(error) {
			printf("Error reading from info file.\n");
			if(NNames!=NULL) for(i=0; i<nbonds; i++) if(NNames[i]!=NULL) free(NNames[i]);
			if(CNames!=NULL) for(i=0; i<nbonds; i++) if(CNames[i]!=NULL) free(CNames[i]);
			if(NNames!=NULL) free(NNames);
			if(CNames!=NULL) free(CNames);
			if(NNums!=NULL) free(NNums);
			if(CNums!=NULL) free(CNums);
			if(fp!=NULL) fclose(fp);
			return -1;
		}
	}
	if(!error) {
		if(nbonds>0) printf("Bond info: \n");
		for(i=0; i<nbonds; i++) printf("%s %d--%s %d\n", NNames[i], NNums[i], CNames[i], CNums[i]);
		printf("\n");
	}
	if(fp!=NULL) fclose(fp);
	*p_NNames = NNames;
	*p_NNums = NNums;
	*p_CNames = CNames;
	*p_CNums = CNums; 
	return nbonds; 
}


int parse_shifts(char* Parname, int *p_nshift, int maxchar) {
	FILE* fp = fopen(Parname, "r");
	char line[maxchar];
	
	SHIFTNDX = NULL;
	SHIFT = NULL;

	*p_nshift = 0;
	int error = 0;
	if(fp==NULL) error = 1;
	while(!error && (fgets(line, maxchar, fp)!=NULL)) {
		if(strncmp(line, "SHIFT", 5)==0) {
			if(sscanf(line, "%*s %d", p_nshift)!=1) error = 1;
			else {
				SHIFT = (real*) malloc((*p_nshift)*sizeof(real));
				if(SHIFT==NULL) error = 2;
				SHIFTNDX = (int*) malloc((*p_nshift)*sizeof(int));
				if(SHIFT==NULL) error = 2;
				int i;
				for(i=0; i<(*p_nshift); i++) {
					if(fgets(line, maxchar, fp)==NULL) error = 1;
					else {
						if(sscanf(line, "%d %f", &SHIFTNDX[i], &SHIFT[i])!=2) error = 1;
					}
					if(error) break;
				}
			}
		}
	}

	if(error) return 0;
	else return 1;
}


/********************************************************************************
* Spectral calculations.							*
********************************************************************************/

int gen_ham_2Q(GNREAL *Ham1Q, int nosc, GNREAL *Ham2Q, int n2Q, real delta) {
	int m,n,p;
	int ndx;
	
	// States are numbered in order 00, 10, 11, 20, 21, 22,...,(nosc-1)(nosc-1)
	// The index of state mn is (m+1)*m/2 + n.
	for(m=0; m<n2Q; m++) for(n=0; n<n2Q; n++) Ham2Q[m*n2Q+n] = 0.0;
	for(m=0; m<nosc; m++) {
		// Double excitation at site m
		Ham2Q[NDX2Q(m,m)*n2Q + NDX2Q(m,m)] = 2*Ham1Q[m*nosc+m]-delta;
		// Single excitations at sites m and n
		for(n=0; n<m; n++) Ham2Q[ NDX2Q(m,n)*n2Q + NDX2Q(m,n) ] = Ham1Q[m*nosc+m]+Ham1Q[n*nosc+n];
		// Coupling between mm and mn states: <mm|H|mn> = sqrt(2)*<m|H|n>
		for(n=0; n<m; n++) Ham2Q[ NDX2Q(m,m)*n2Q + NDX2Q(m,n) ] = sqrt(2.0)*Ham1Q[m*nosc+n];
		// Coupling between mm and nm states: <mm|H|nm> = sqrt(2)*<m|H|n>
		for(n=m+1; n<nosc; n++) Ham2Q[ NDX2Q(m,m)*n2Q + NDX2Q(n,m) ] = sqrt(2.0)*Ham1Q[m*nosc+n];
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

int gen_dip_2Q(GNREAL *Dip1Q, GNREAL *Dip2Q, int nosc, int n2Q) {
	int i,m,n,M;
	real fac;
	real sqrt2 = sqrt(2.0);
	// cjfeng 04/05/2016
	// Interchange the loop and the matrix is now transposed to be (nbuffer*3)*n2Q*nosc. 
	int num = n2Q*nosc;
	for(M=0; M<num; M++) Dip2Q[M] = 0.0;
	
	// cjfeng 04/05/2016
	// Interchange the loop and the matrix is now transposed to be (nbuffer*3)*n2Q*nosc.
	for(n=0; n<nosc; n++) {
		for(m=0; m<nosc; m++) {
			if(m!=n) fac = 1.0;
			else fac = sqrt2;
			if(m<n) M = ((n*(n+1)/2)+m);
			else M = ((m*(m+1)/2)+n);
			
			Dip2Q[M*nosc+m] = fac*Dip1Q[n];
		}
	}
	return 1;
}

GNREAL orient(GNREAL **Dip1[3], GNREAL **Dip2[3], int tid, int ndxA, int ndxB, int ndxC, int ndxD, int pol) {
	real Val = 0.0;
	int a,b;
	// cjfeng 06/27/2016
	// Improving locality.
	GNREAL M0, M1, M2, M3;
	M0 =  M_ijkl_IJKL[pol][0];
	M1 =  M_ijkl_IJKL[pol][1];
	M2 =  M_ijkl_IJKL[pol][2];
	M3 =  M_ijkl_IJKL[pol][3];
	for(a=0; a<3; a++) {
		GNREAL D1aA, D1aB, D1bB, D2aC, D2bC, D2aD, D2bD;
		D1aA = Dip1[a][tid][ndxA];
		D1aB = Dip1[a][tid][ndxB];
		D2aC = Dip2[a][tid][ndxC];
		D2aD = Dip2[a][tid][ndxD];
		Val += M0*D1aA*D1aB*D2aC*D2aD;
		// Val += M_ijkl_IJKL[pol][0]*Dip1[a][tid][ndxA]*Dip1[a][tid][ndxB]*Dip2[a][tid][ndxC]*Dip2[a][tid][ndxD];
		for(b=0; b<3; b++) {
			if(b!=a) {
				D1bB = Dip1[b][tid][ndxB];
				D2bC = Dip2[b][tid][ndxC];
				D2bD = Dip2[b][tid][ndxD];
				Val += M1*D1aA*D1aB*D2bC*D2bD;
				Val += M2*D1aA*D1bB*D2bC*D2aD;
				Val += M3*D1aA*D1bB*D2aC*D2bD;
				// Val += M_ijkl_IJKL[pol][1]*Dip1[a][tid][ndxA]*Dip1[a][tid][ndxB]*Dip2[b][tid][ndxC]*Dip2[b][tid][ndxD];
				// Val += M_ijkl_IJKL[pol][2]*Dip1[a][tid][ndxA]*Dip1[b][tid][ndxB]*Dip2[b][tid][ndxC]*Dip2[a][tid][ndxD];
				// Val += M_ijkl_IJKL[pol][3]*Dip1[a][tid][ndxA]*Dip1[b][tid][ndxB]*Dip2[a][tid][ndxC]*Dip2[b][tid][ndxD];
			}
		}
	}
	return Val; 
}

int gen_pert_prop( GNCOMP *U1Q, GNCOMP *U2Q, int nosc, int n2Q, double expfac, double delta) {
	int i,j,ip,jp;
	GNCOMP fac;
	double sqrt2inv = 0.707106781186547;
	int n2Qsq = n2Q*n2Q;
	for(i=0; i<n2Qsq; i++) U2Q[i] = 0.0;
	for(i=0; i<nosc; i++) {
		for(j=0; j<=i; j++) {
			for(ip=0; ip<nosc; ip++) {
				for(jp=0; jp<=ip; jp++) {
					fac = 1.0;
					// negative sign on sin() term since delta is provided as a positive argument.
					if(i==j) fac *= sqrt2inv;
					if(ip==jp) fac *= sqrt2inv;
					if(i==j) fac *= ( cos(0.5*expfac*delta) - I*sin(0.5*expfac*delta) );
					if(ip==jp) fac *= ( cos(0.5*expfac*delta) - I*sin(0.5*expfac*delta) );
					U2Q[NDX2Q(i,j)*n2Q+NDX2Q(ip,jp)] = fac*( U1Q[i*nosc+ip]*U1Q[j*nosc+jp] + U1Q[j*nosc+ip]*U1Q[i*nosc+jp] );
				}
			}
		}
	}
	return 1;
}

int gen_perturb_2Q_energies(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evals2Q, int nosc, real delta ) {
	int m,n,i,j;
	real val, anh, c;
	for(m=0; m<nosc; m++) {
		for(n=0; n<=m; n++) {
			val = Evals1Q[m] + Evals1Q[n]; // harmonic 
			anh = 0.0;
			for(i=0; i<nosc; i++) anh += Evecs1Q[i*nosc+n]*Evecs1Q[i*nosc+m]*Evecs1Q[i*nosc+n]*Evecs1Q[i*nosc+m];
			if(m!=n) anh = anh * 2.0;
			Evals2Q[NDX2Q(m,n)] = val-delta*anh;
		}
	}
	return 1;
}

int gen_perturb_2Q_vectors(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evecs2Q, GNREAL *Evals2Q, int nosc, int n2Q, real delta ) {
	GNREAL *Temp2Q = malloc(n2Q*n2Q*sizeof(GNREAL));
	real cutoff = 1.0; // cm-1
	int m,n,i,j,mp,np;
	real sqrt2 = sqrt(2.0);
	real sqrt2inv = (1.0)/sqrt2;
	real val, anh, c, gap;

	for(i=0; i<n2Q; i++) for( j=0; j<n2Q; j++) Evecs2Q[i*n2Q+j] = 0.0;
	// Harmonic vectors
	for(m=0; m<nosc; m++) {
		for(n=0; n<=m; n++) {
			for(i=0; i<nosc; i++) {
				for(j=0; j<=i; j++) {
					Temp2Q[NDX2Q(i,j)*n2Q + NDX2Q(m,n)] = Evecs1Q[i*nosc+m]*Evecs1Q[j*nosc+n] + Evecs1Q[j*nosc+m]*Evecs1Q[i*nosc+n]*(1.0*(m!=n) + sqrt2inv*(m==n))*(1.0*(i!=j) + sqrt2inv*(i==j));
				}
			}
		}
	}

	for(i=0; i<n2Q; i++) for(j=0; j<n2Q; j++) Evecs2Q[i*n2Q+j] = Temp2Q[i*n2Q+j];
	// And the anharmonic correction
	for( m=0; m<nosc; m++) {
		for(n=0; n<=m; n++) {
			for(mp=0; mp<nosc; mp++) {
				for(np=0; np<=mp; np++) {
					anh = 0;
					for(i=0; i<nosc; i++) anh += Evecs1Q[i*nosc+m]*Evecs1Q[i*nosc+n]*Evecs1Q[i*nosc+mp]*Evecs1Q[i*nosc+np];
					anh *= delta*(1.0*(mp!=np) + sqrt2inv*(mp==np) )*(1.0*(m!=n) + sqrt2inv*(m==n));
					gap = Evals1Q[m] + Evals1Q[n] - Evals1Q[mp] - Evals1Q[np];
					if( ((gap>0) && (gap>cutoff)) || ((gap<0) && (gap<-cutoff)) ) {
						for(i=0; i<n2Q; i++) Evecs2Q[i*n2Q+NDX2Q(m,n)] -= (anh/gap)*Temp2Q[i*n2Q+NDX2Q(mp,np)];
					}
				}
			}
		}
	}
	for(m=0; m<n2Q; m++) {
		val = 0.0;
		for(i=0; i<n2Q; i++) val += Evecs2Q[i*n2Q+m]*Evecs2Q[i*n2Q+m];
		val = (1.0)/sqrt(val);
		for(i=0; i<n2Q; i++) Evecs2Q[i*n2Q+m] *= val;
	}
	free(Temp2Q);
	return 1;
}

int gen_perturb_2Q_matrix(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evecs2Q, GNREAL *Evals2Q, GNREAL *Temp2Q, int nosc, int n2Q, real delta ) {
	int m,n,i,j,mp,np;
	real sqrt2 = sqrt(2.0);
	real sqrt2inv = (1.0)/sqrt2;
	real val, anh, c, fac;

	for(i=0; i<n2Q; i++) for( j=0; j<n2Q; j++) Temp2Q[i*n2Q+j] = 0.0;
	for(m=0; m<nosc; m++) {
		for(n=0; n<=m; n++) {
			for(mp=0; mp<nosc; mp++) {
				for(np=0; np<=mp; np++) {
					for(i=0; i<nosc; i++) {
						fac = 2.0*delta;
						if(m==n) fac *= sqrt2inv;
						if(mp==np) fac *= sqrt2inv;
						Temp2Q[NDX2Q(m,n)*n2Q+NDX2Q(mp,np)] -= fac*Evecs1Q[i*nosc+m]*Evecs1Q[i*nosc+n]*Evecs1Q[i*nosc+mp]*Evecs1Q[i*nosc+np];
					}
				}
			}
		}
	}
	return 1;
}

int calc_2dir(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], GNREAL **ExDip2QAr[3], int tid, int nosc, int n2Q, int npts, GNREAL res, GNREAL start, GNREAL stop, GNCOMP *REPH[4], GNCOMP *NREPH[4], int POL[4], int npol, int reph, int nreph) {
	int i, a, b, c;
	int ndx1, ndx3;

	GNREAL *Evals1Q = Evals1QAr[tid];
	GNREAL *Evals2Q = Evals2QAr[tid];

	for(a=0; a<nosc; a++) {
		ndx1 = floor( (Evals1Q[a]-start)/res + 0.5 );
		if( (ndx1>0) && (ndx1<npts) ) {
			for(b=0; b<nosc; b++) {
				// Pathway 1 - rephasing & non-rephasing
				ndx3 = floor( (Evals1Q[b]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) {
					if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1*npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);
					if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1*npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);
				}

				// Pathway 2 - rephasing
				ndx3 = floor( (Evals1Q[b]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1*npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, a, b, POL[i]);

				// Pathway 2 - non-rephasing
				ndx3 = floor( (Evals1Q[a]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1*npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, b, a, POL[i]);
				for(c=0; c<n2Q; c++) {
					// Pathway 3 - rephasing
					ndx3 = floor( ( (Evals2Q[c]-Evals1Q[a])-start)/res + 0.5 );
					if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1*npts+ndx3] -= orient(ExDip1QAr, ExDip2QAr, tid, a, b, b*n2Q+c, a*n2Q+c, POL[i]);

					// Pathway 3 - non-rephasing
					ndx3 = floor( ( (Evals2Q[c]-Evals1Q[b])-start)/res + 0.5 );
					if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1*npts+ndx3] -= orient(ExDip1QAr, ExDip2QAr, tid, a, b, a*n2Q+c, b*n2Q+c, POL[i]);
				}
			}
		}
	}
	return 1;
}


int calc_2dir_pert(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], int tid, int nosc, int n2Q, int npts, GNREAL res, GNREAL start, GNREAL stop, GNCOMP *REPH[4], GNCOMP *NREPH[4], int POL[4], int npol, int reph, int nreph) {
	int i, a, b, c;
	int ndx, ndx1, ndx3;

	GNREAL *Evals1Q = Evals1QAr[tid];
	GNREAL *Evals2Q = Evals2QAr[tid];

	for(a=0; a<nosc; a++) {
		ndx1 = floor( (Evals1Q[a]-start)/res + 0.5 );
		if( (ndx1>0) && (ndx1<npts) ) {
			for(b=0; b<nosc; b++) {
				// Pathway 1 - rephasing & non-rephasing
				ndx3 = floor( (Evals1Q[b]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) {
					if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1*npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);
					if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1*npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);
				}

				// Pathway 2 - rephasing
				ndx3 = floor( (Evals1Q[b]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1*npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, a, b, POL[i]);

				// Pathway 2 - non-rephasing
				ndx3 = floor( (Evals1Q[a]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1*npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, b, a, POL[i]);

				if(a>b) ndx = ((a*(a+1)/2)+b);
				else ndx = ((b*(b+1)/2)+a);
				
				// Pathway 3 - rephasing (0,0) -> (0,a) -> (a,a) -> (ab,a)
				ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[a])-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1*npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);

				// Pathway 3 - rephasing (0,0) -> (0,a) -> (b,a) -> (ab,a)
				ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[a])-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1*npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, b, a, b, POL[i]);

				// Pathway 3 - non-rephasing (0,0) -> (0,a) -> (a,a) -> (a,ab)
				ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[a])-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1*npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);

				// Pathway 3 - non-rephasing (0,0) -> (0,a) -> (b,a) -> (b,ab)
				ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[b])-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1*npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, b, b, a, POL[i]);
			}	
		}
	}
	return 1;
}

int make_fnames( char* Infoname, char* Hamname, char* hamnm, char Dipnames[3][maxchar], char* dipxnm, char* dipynm, char* dipznm, char* Sitesname, char* sitesnm, char* axisnm, char* ftirnm, char* lognm, char polnm[4][16], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnm[10][maxchar], char* Parname, char* paramnm, char* outname, char* deffnm, int do_traj ) {
	int i;		// Index running over all possible polarizations.
	if(outname!=NULL) {
		printf("Output file name base: %s\n", outname);
		strcpy(ftirnm, outname);				// Copying output name base to ftir filename.
		strcpy(lognm, outname);					// Copying output name base to log filename.
		strcpy(axisnm, outname);				// Copying output name base to frequency axis filename.
		for(i=0; i<4; i++) {						// Running over polarizations
			strcpy(rephnm[i], outname);		// Copying output name base to rephasing spectra filenames
			strcpy(nrephnm[i], outname);	// Copying output name base to non-rephasing spectra filenames
		}
		// Adding "_" between output file name base and suffices 
		// when name base is not referring to a folder.
		if(outname[strlen(outname)-1]!='/') {
			strcat(ftirnm, "_");
			strcat(lognm, "_");
			strcat(axisnm, "_");
			for(i=0; i<4; i++) {
				strcat(rephnm[i], "_");
				strcat(nrephnm[i], "_");
			}
		}
	} else {
		ftirnm[0] = '\0';
		lognm[0] = '\0';
		axisnm[0] = '\0';
		for(i=0; i<4; i++) {
			rephnm[i][0] = '\0';
			nrephnm[i][0] = '\0';
		}
	}
	strcat(ftirnm, "ftir.txt");			// Adding suffix to ftir file
	strcat(lognm, "log.txt");				// Adding suffix to log file
	strcat(axisnm, "waxis.txt");		// Adding suffix to frequency axis file

	for(i=0; i<4; i++) {
		strcat(rephnm[i], "reph_");		// Adding suffix to rephasing spectra files
		strcat(rephnm[i], polnm[i]);
		strcat(nrephnm[i], "nreph_");	// Adding suffix to non-rephasing spectra file
		strcat(nrephnm[i], polnm[i]);
	}


	// Input file names
	// First info file
	if(deffnm!=NULL) {
		strcpy(Infoname, deffnm);			// Copying input file name base
		if(deffnm[strlen(deffnm)-1]!='/') strcat(Infoname, "_"); // Adding "_" between deffnm and suffix.
	} else {
		Infoname[0] = '\0';
	}
	strcat(Infoname, "info.txt");		// Adding suffix

	// Hamiltonian file
	if( (deffnm!=NULL) && (!strcmp(hamnm, "ham.txt")) ) {
		printf("Made it inside\n");
		strcpy(Hamname, deffnm);
		if(deffnm[strlen(deffnm)-1]!='/') strcat(Hamname, "_");
	} else Hamname[0] = '\0';
	strcat(Hamname, hamnm);
	
	// Dipole moment files	
	if( (deffnm!=NULL) && (!strcmp(dipxnm,"dipx.txt")) ) {
                strcpy(Dipnames[0], deffnm);
                if(deffnm[strlen(deffnm)-1]!='/') strcat(Dipnames[0], "_");
        } else Dipnames[0][0] = '\0';
        strcat(Dipnames[0], dipxnm);

        if( (deffnm!=NULL) && (!strcmp(dipynm,"dipy.txt")) ) {
                strcpy(Dipnames[1], deffnm);
                if(deffnm[strlen(deffnm)-1]!='/') strcat(Dipnames[1], "_");
        } else Dipnames[1][0] = '\0';
        strcat(Dipnames[1], dipynm);

        if( (deffnm!=NULL) && (!strcmp(dipznm,"dipz.txt")) ) {
                strcpy(Dipnames[2], deffnm);
                if(deffnm[strlen(deffnm)-1]!='/') strcat(Dipnames[2], "_");
        } else Dipnames[2][0] = '\0';
        strcat(Dipnames[2], dipznm);

	// Sitename and Parname are unused if not specified
	// The default file name is never appended.
	Sitesname[0] = '\0';
	strcat(Sitesname, sitesnm);
	
	Parname[0] = '\0';
	strcat(Parname, paramnm);

	if(do_traj) {
		if(outname!=NULL) {
			strcpy(Trajnm[0], outname);
			if(outname[strlen(outname)-1]!='/') strcat(Trajnm[0], "_");
		} else {
			Trajnm[0][0] = '\0';
		}
		strcat(Trajnm[0], "tstamp");
		strncpy(Trajnm[1], ftirnm, strlen(ftirnm)-4);
		for(i=2; i<6; i++) strncpy(Trajnm[i], rephnm[i-2], strlen(rephnm[i-2])-4);
		for(i=6; i<10; i++) strncpy(Trajnm[i], nrephnm[i-6], strlen(nrephnm[i-6])-4);
		for(i=0; i<10; i++) {
			strcat(Trajnm[i], "_traj.txt");
		}
	}
	return 0;
}

int open_all( char* Hamname, char Dipnames[3][maxchar], char* Sitesname, char* axisnm, char* ftirnm, char* lognm, int npol, int POL[4], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnames[10][maxchar], int do2d, int do_traj) {
	int error=0;
	int i,j;
	if(!error) {
		hfp = fopen(Hamname, "r");
		if(hfp==NULL) error = 2;
	} else hfp = NULL;
	for(i=0; i<3; i++) {
		if(!error) {
			Dfp[i] = fopen(Dipnames[i], "r");
			if(Dfp[i]==NULL) error = 2;
		} else Dfp[i] = NULL;
	}
	if(strlen(Sitesname)>0) {
		sfp = fopen(Sitesname, "r");
		if(sfp==NULL) error = 2;
	} else sfp = NULL;
	if(!error) {
		afp = fopen(axisnm, "w");
		if(afp==NULL) error = 2;
	} else afp = NULL;
	if(!error) {
		ffp = fopen(ftirnm, "w");
		if(ffp==NULL) error = 2;
	} else ffp = NULL;
	if(!error) {
		lfp = fopen(lognm, "w");
		if(lfp==NULL) error = 2;
	} else lfp = NULL;
	if( (!error) && (do_traj) ) {
		// Delete file names that don't get used
		for(i=0; i<4; i++) {
			int use_pol = 0;
			for(j=0; j<npol; j++) if(POL[j]==i) use_pol = 1;
			if(use_pol==0) {
				Trajnames[2+i][0] = '\0';
				Trajnames[6+i][0] = '\0';
			}
		}
		if(!do2d) for(i=2; i<10; i++) Trajnames[i][0] = '\0';
		for(i=0; i<10; i++) {
			if( (strlen(Trajnames[i])>0) && (!error) ) {
				printf("%s\n", Trajnames[i]);
				Trajfp[i] = fopen(Trajnames[i], "w");
				if(Trajfp[i]==NULL) error = 2;
			} else Trajfp[i] = NULL;
		}
	} else for(i=0; i<10; i++) Trajfp[i] = NULL;

	if(do2d) {
		for(i=0; i<npol; i++) {
			if(!error) {
				rfp[POL[i]] = fopen(rephnm[POL[i]], "w");
				if(rfp[POL[i]]==NULL) error = 2;
			} else rfp[POL[i]] = NULL;
			if(!error) {
				nrfp[POL[i]] = fopen(nrephnm[POL[i]], "w");
				if(nrfp[POL[i]]==NULL) error = 2;
			} else nrfp[POL[i]] = NULL;
		}
	} else {
		for(i=0; i<npol; i++) {
			rfp[POL[i]] = NULL; 
			nrfp[POL[i]] = NULL;
		}
	}
	if(!error) printf("Opened all files.\n");
	return error;
}

int allocate_all( int window, int win2d, int winzpad, int nosc, int nbuffer, int nthreads, int n2Q, int pert, int do2d, GNREAL tstep, GNREAL TauP, int npts, int nise ) {
	int error=0;
	int i,j,k;
	int npopdec = 2*(win2d-window);	// equivalent with (window + T2/tstep + 1)
	int windowsq = window*window;
	int noscsq = nosc*nosc;
	int n2Qsq = n2Q*n2Q;
	int nptssq = npts*npts;
	// int npopdec = 2*(window+win2d-2*window);
	const double pi = 3.14159265;

	int stick;
	if(nise) stick = 0;
	else stick = 1;

	if(!stick) {
		if(!error) {
			CorrFunc = (GNCOMP*) malloc(sizeof(GNCOMP)*window);
			if(CorrFunc==NULL) error = 1;
			else for(i=0; i<window; i++) {
				CorrFunc[i] = 0.0;
			}
		} else CorrFunc = NULL;
		if(!error) {
			NetCorrFunc = (GNCOMP*) malloc(sizeof(GNCOMP)*window);
			if(NetCorrFunc==NULL) error = 1;
			else for(i=0; i<window; i++) {
				NetCorrFunc[i] = 0.0;
			}
		} else NetCorrFunc = NULL;

		if(!error) {
			popdecay1Q = (GNREAL*) malloc(sizeof(GNREAL)*npopdec);
			if(popdecay1Q==NULL) error = 1;
			else for(i=0; i<npopdec; i++) popdecay1Q[i] = exp(-i*tstep/(2*TauP));
		} else popdecay1Q = NULL;
		if((!error) && (do2d)) {
			popdecay2Q = (GNREAL*) malloc(sizeof(GNREAL)*npopdec);
			if(popdecay2Q==NULL) error = 1;
			else for(i=0; i<npopdec; i++) popdecay2Q[i] = exp(-i*tstep/(2*TauP));
		} else popdecay2Q = NULL;
		for(i=0; i<4; i++) {
			if((!error) && (do2d)) {
				REPH[i] = (GNCOMP* ) malloc(windowsq*sizeof(GNCOMP));
				if(REPH[i]==NULL) error = 1;
				else for(j=0; j<windowsq; j++) REPH[i][j] = 0.0;
			} else REPH[i] = NULL;
			if((!error) && (do2d)) {
				NREPH[i] = (GNCOMP* ) malloc(windowsq*sizeof(GNCOMP));
				if(NREPH[i]==NULL) error = 1;
				else for(j=0; j<windowsq; j++) NREPH[i][j] = 0.0;
			} else NREPH[i] = NULL;
		}
		for(i=0; i<4; i++) {
			if((!error) && (do2d)) {
				NetREPH[i] = (GNCOMP* ) malloc(windowsq*sizeof(GNCOMP));
				if(NetREPH[i]==NULL) error = 1;
				else for(j=0; j<windowsq; j++) NetREPH[i][j] = 0.0;
			} else NetREPH[i] = NULL;
			if((!error) && (do2d)) {
				NetNREPH[i] = (GNCOMP* ) malloc(windowsq*sizeof(GNCOMP));
				if(NetNREPH[i]==NULL) error = 1;
				else for(j=0; j<windowsq; j++) NetNREPH[i][j] = 0.0;
			} else NetNREPH[i] = NULL;
		}

	
		// One-quantum arrays
		if(!error) {
			U1Q = (GNCOMP*) malloc(noscsq*sizeof(GNCOMP));
			if(U1Q==NULL) error = 1;
		} else U1Q = NULL;
		if(!error) {
			U1Qs = (GNCOMP*) malloc(noscsq*sizeof(GNCOMP));
			if(U1Qs==NULL) error = 1;
		} else U1Qs = NULL;
		for(i=0; i<3; i++) {
			if(!error) {
				cDip1Q[i] = (GNCOMP*) malloc(nosc*sizeof(GNCOMP));
				if(cDip1Q[i]==NULL) error = 1;
			} else cDip1Q[i] = NULL;
		}
		if(!error) {
			U1QMat = (GNCOMP**) malloc(nbuffer*sizeof(GNCOMP*));
			if(U1QMat==NULL) error = 1;
		} else U1QMat = NULL;
		if(U1QMat!=NULL) {
			for(i=0; i<nbuffer; i++) {
				if(!error) {
					U1QMat[i] = (GNCOMP*) malloc(noscsq*sizeof(GNCOMP));
					if(U1QMat[i]==NULL) error = 1;
				} else U1QMat[i] = NULL;
			}
		}
		for(i=0; i<3; i++) {
			if((!error) && (do2d) ) {
				psi_a[i] = (GNCOMP**) malloc(win2d*sizeof(GNCOMP*));
				if(psi_a[i]==NULL) error = 1;
			} else psi_a[i] = NULL;
			if((!error) && (do2d) ) {
				psi_b1[i] = (GNCOMP**) malloc(win2d*sizeof(GNCOMP*));
				if(psi_b1[i]==NULL) error = 1;
			} else psi_b1[i] = NULL;
			if((!error) && (do2d) ) {
				psi_b12[i] = (GNCOMP**) malloc(win2d*sizeof(GNCOMP*));
				if(psi_b12[i]==NULL) error = 1;
			} else psi_b12[i] = NULL;
			for(j=0; j<win2d; j++) {
				if((!error) && (do2d) ) {
					psi_a[i][j] = (GNCOMP*) malloc(nosc*sizeof(GNCOMP));
					if(psi_a[i][j]==NULL) error = 1;
				} else if(psi_a[i]!=NULL) psi_a[i][j] = NULL;
				if((!error) && (do2d) ) {
					psi_b1[i][j] = (GNCOMP*) malloc(nosc*sizeof(GNCOMP));
					if(psi_b1[i][j]==NULL) error = 1;
				} else if(psi_b1[i]!=NULL) psi_b1[i][j] = NULL;
				if((!error) && (do2d) ) {
					psi_b12[i][j] = (GNCOMP*) malloc(nosc*sizeof(GNCOMP));
					if(psi_b12[i][j]==NULL) error = 1;
				} else if(psi_b12[i]!=NULL) psi_b12[i][j] = NULL;
			}
			if( (!error) && (do2d) ) {
				for(j=0; j<win2d; j++) {	
					for(k=0; k<nosc; k++) {
						psi_a[i][j][k] = 0.0;
						psi_b1[i][j][k] = 0.0;
						psi_b12[i][j][k] = 0.0;
					}
				}
			}
		}
		for(i=0; i<9; i++) {
			if((!error) && (do2d) ) {
				psi_cb[i] = (GNCOMP*) malloc(n2Q*sizeof(GNCOMP));
				if(psi_cb[i]==NULL) error = 1;
			} else psi_cb[i] = NULL;
			if((!error) && (do2d) ) {
				psi_ca[i] = (GNCOMP*) malloc(n2Q*sizeof(GNCOMP));
				if(psi_ca[i]==NULL) error = 1;
			} else psi_ca[i] = NULL;
		}
		if((!error) && (do2d) ) {
			psi2Q = (GNCOMP*) malloc(n2Q*sizeof(GNCOMP));
			if(psi2Q==NULL) error = 1;
		} else psi2Q = NULL;
		
		// Two-quantum arrays
		if((!error) && (do2d) ) {
			U2Qs = (GNCOMP*) malloc(n2Qsq*sizeof(GNCOMP));
			if(U2Qs==NULL) error = 1;
		} else U2Qs = NULL;
	
		if((!error) && (do2d) ) {
			U2QMat = (GNCOMP**) malloc(nbuffer*sizeof(GNCOMP*));
			if(U2QMat==NULL) error = 1;
		} else U2QMat = NULL;
		if(U2QMat!=NULL) {
			for(i=0; i<nbuffer; i++) {
				if((!error) && (do2d) ) {
					U2QMat[i] = (GNCOMP*) malloc(n2Qsq*sizeof(GNCOMP));
					if(U2QMat[i]==NULL) error = 1;
				} else U2QMat[i] = NULL;
			}
		}
		ftir = NULL;
		netftir = NULL;
		for(i=0; i<3; i++) ExDip1QAr[i] = NULL;
		for(i=0; i<3; i++) ExDip2QAr[i] = NULL;
	} else {
		for(i=0; i<4; i++) {
			if((!error) && (do2d) ) {
				REPH[i] = (GNCOMP* ) malloc(nptssq*sizeof(GNCOMP));
				if(REPH[i]==NULL) error = 1;
				else for(j=0; j<nptssq; j++) REPH[i][j] = 0.0;
			} else REPH[i] = NULL;
			if((!error) && (do2d) ) {
				NREPH[i] = (GNCOMP* ) malloc(nptssq*sizeof(GNCOMP));
				if(NREPH[i]==NULL) error = 1;
				else for(j=0; j<nptssq; j++) NREPH[i][j] = 0.0;
			} else NREPH[i] = NULL;
		}
		for(i=0; i<4; i++) {
			if((!error) && (do2d) ) {
				NetREPH[i] = (GNCOMP* ) malloc(nptssq*sizeof(GNCOMP));
				if(NetREPH[i]==NULL) error = 1;
				else for(j=0; j<nptssq; j++) NetREPH[i][j] = 0.0;
			} else NetREPH[i] = NULL;
			if((!error) && (do2d) ) {
				NetNREPH[i] = (GNCOMP* ) malloc(nptssq*sizeof(GNCOMP));
				if(NetNREPH[i]==NULL) error = 1;
				else for(j=0; j<nptssq; j++) NetNREPH[i][j] = 0.0;
			} else NetNREPH[i] = NULL;
		}

		if(!error) {
			ftir = (GNREAL*) malloc(npts*sizeof(GNREAL));
			if(ftir==NULL) error = 1;
			else for(i=0; i<npts; i++) ftir[i] = 0.0;
		} else ftir = NULL;
		if(!error) {
			netftir = (GNREAL*) malloc(npts*sizeof(GNREAL));
			if(netftir==NULL) error = 1;
			else for(i=0; i<npts; i++) netftir[i] = 0.0;
		} else netftir = NULL;
		for(i=0; i<3; i++) {
			if(!error) {
				ExDip1QAr[i] = (GNREAL**) malloc(nthreads*sizeof(GNREAL*));
				if(ExDip1QAr[i]==NULL) error = 1;
			} else ExDip1QAr[i] = NULL;
			if((!error) && (do2d) && ( (!pert) || pertvec )) {
				ExDip2QAr[i] = (GNREAL**) malloc(nthreads*sizeof(GNREAL*));
				if(ExDip2QAr[i]==NULL) error = 1;
			} else ExDip2QAr[i] = NULL;
			for(j=0; j<nthreads; j++) {
				if(!error) {
					ExDip1QAr[i][j] = (GNREAL*) malloc(nosc*sizeof(GNREAL));
					if(ExDip1QAr[i][j]==NULL) error = 1;
				} else if(ExDip1QAr[i]!=NULL) ExDip1QAr[i][j] = NULL;
	
				if((!error) && (do2d) && ( (!pert) || pertvec ))  {
					ExDip2QAr[i][j] = (GNREAL*) malloc(nosc*n2Q*sizeof(GNREAL));
					if(ExDip2QAr[i][j]==NULL) error = 1;
				} else if(ExDip2QAr[i]!=NULL) ExDip2QAr[i][j] = NULL;
			}
		}

		CorrFunc = NULL;
		NetCorrFunc = NULL;
		popdecay1Q = NULL;
		popdecay2Q = NULL;
	
		// One-quantum arrays
		U1Q = NULL;
		U1Qs = NULL;
		for(i=0; i<3; i++) cDip1Q[i] = NULL;
		U1QMat = NULL;
		for(i=0; i<3; i++) {
			psi_a[i] = NULL;
			psi_b1[i] = NULL;
			psi_b12[i] = NULL;
		}
		for(i=0; i<9; i++) psi_cb[i] = NULL;
		psi2Q = NULL;
		U2Qs = NULL;
		U2QMat = NULL;
	}
	if(!error) {
		SitesBuffer = (GNREAL*) malloc(sizeof(GNREAL)*nosc);
		if(SitesBuffer==NULL) error = 1;
	} else SitesBuffer = NULL;

	if(!error) {
		FTin1D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*winzpad);
		if(FTin1D==NULL) error = 1;
	} else FTin1D = NULL;
	if(!error) {
		FTout1D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*winzpad);
		if(FTout1D==NULL) error = 1;
	} else FTout1D = NULL;
	if(!error) {
		FTplan1D = fftw_plan_dft_1d(winzpad, FTin1D, FTout1D, FFTW_BACKWARD, FFTW_ESTIMATE);
	} else FTplan1D = NULL;

	if((!error) && (do2d) ) {
		FTin2D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*winzpad*winzpad);
		if(FTin2D==NULL) error = 1;
	} else FTin2D = NULL;
	if((!error) && (do2d) ) {
		FTout2D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*winzpad*winzpad);
		if(FTout2D==NULL) error = 1;
	} else FTout2D = NULL;
	if((!error) && (do2d) ) {
		FTplan2D = fftw_plan_dft_2d(winzpad, winzpad, FTin2D, FTout2D, FFTW_BACKWARD, FFTW_ESTIMATE);
	} else FTplan2D = NULL;


	if(!error) {
		hann = (GNREAL*) malloc(sizeof(GNREAL)*window);
		if(hann==NULL) error = 1;
		else {
			// Initialize
			if(window==1) hann[0] = 1;
			else {
				for(i=0; i<window; i++) hann[i] = 0.5*(1+cos((pi*i)/(window-1)));
				// And normalize s.t. sum(hann) = 1.
				GNREAL val = 0.0;
				for(i=0; i<window; i++) val += hann[i];
				for(i=0; i<window; i++) hann[i] /= val;
			}
		}
	} else hann = NULL;

		
	if(!error) {
		Ham1QMat = (GNREAL**) malloc(nbuffer*sizeof(GNREAL*));
		if(Ham1QMat==NULL) error = 1;
	} else Ham1QMat = NULL;
	if(Ham1QMat!=NULL) {
		for(i=0; i<nbuffer; i++) {
			if(!error) {
				Ham1QMat[i] = (GNREAL*) malloc(noscsq*sizeof(GNREAL));
				if(Ham1QMat[i]==NULL) error = 1;
			} else Ham1QMat[i] = NULL;
		}
	}
	for(i=0; i<3; i++) {
		if(!error) {
			Dip1QMat[i] = (GNREAL**) malloc(nbuffer*sizeof(GNREAL*));
			if(Dip1QMat[i]==NULL) error = 1;
		} else Dip1QMat[i] = NULL;
		if(Dip1QMat[i]!=NULL) {
			for(j=0; j<nbuffer; j++) {
				if(!error) {
					Dip1QMat[i][j] = (GNREAL*) malloc(nosc*sizeof(GNREAL));
					if(Dip1QMat[i][j]==NULL) error = 1;
				} else Dip1QMat[i][j] = NULL;
			}
		}

	}

	if(!error) {
		Ham1QAr = (GNREAL**) malloc(nthreads*sizeof(GNREAL*));
		if(Ham1QAr==NULL) error = 1;
	} else Ham1QAr = NULL;
	if(!error) {
		Evals1QAr = (GNREAL**) malloc(nthreads*sizeof(GNREAL*));
		if(Evals1QAr==NULL) error = 1;
	} else Evals1QAr = NULL;

	for(i=0; i<3; i++) {
		if(!error) {
			Dip1QAr[i] = (GNREAL**) malloc(nthreads*sizeof(GNREAL*));
			if(Dip1QAr[i]==NULL) error = 1;
		} else Dip1QAr[i] = NULL;
	}

	for(i=0; i<nthreads; i++) {
		if(!error) {
			Ham1QAr[i] = (GNREAL*) malloc(noscsq*sizeof(GNREAL));
			if(Ham1QAr[i]==NULL) error = 1;
		} else if(Ham1QAr!=NULL) Ham1QAr[i] = NULL;
		if(!error) {
			Evals1QAr[i] = (GNREAL*) malloc(nosc*sizeof(GNREAL));
			if(Evals1QAr[i]==NULL) error = 1;
		} else if(Evals1QAr!=NULL) Evals1QAr[i] = NULL;

		for(j=0; j<3; j++) {
			if(!error) {
				Dip1QAr[j][i] = (GNREAL*) malloc(nosc*sizeof(GNREAL));
				if(Dip1QAr[j][i]==NULL) error = 1;
			} else if(Dip1QAr[j]!=NULL) Dip1QAr[j][i] = NULL;
		}
	}


	if((!error) && (do2d) && ((!pert) || (pertvec)) ) {
		Ham2QAr = (GNREAL**) malloc(nthreads*sizeof(GNREAL*));
		if(Ham2QAr==NULL) error = 1;
	} else Ham2QAr = NULL;
	if((!error) && (do2d) ) {
		Evals2QAr = (GNREAL**) malloc(nthreads*sizeof(GNREAL*));
		if(Evals2QAr==NULL) error = 1;
	} else Evals2QAr = NULL;

	for(i=0; i<3; i++) {
		if((!error) && (do2d) && ( (!pert) || pertvec )) {
			Dip2QAr[i] = (GNREAL**) malloc(nthreads*sizeof(GNREAL*));
			if(Dip2QAr[i]==NULL) error = 1;
		} else Dip2QAr[i] = NULL;
	}

	// cjfeng 04/05/2016
	// Interchange the order of the array to make
	// Dip2QMat a nbuffer*3*n2Q*nosc array
	if((!error) && (do2d) ) {
		Dip2QMat = (GNREAL***) malloc(nbuffer*sizeof(GNREAL**));
		if(Dip2QMat==NULL) error = 1;
		else {
			for(j=0; j<nbuffer; j++) {
				if(!error) {
					Dip2QMat[j] = (GNREAL**) malloc(3*sizeof(GNREAL*));
					if(Dip2QMat[j]==NULL) error = 1;
					else {
						if(!error) {
							for(i=0; i<3; i++) {
								Dip2QMat[j][i] = (GNREAL*) malloc(n2Q*nosc*sizeof(GNREAL));
								if(Dip2QMat[j][i]==NULL) error = 1;
							}
						} 
						else Dip2QMat[j][i] = NULL;
					}
				}
				else Dip2QMat[j] = NULL;
			}
		} 
	}
	else Dip2QMat = NULL;

	for(i=0; i<nthreads; i++) {
		if((!error) && (do2d) && ((!pert) || (pertvec)) ) {
			Ham2QAr[i] = (GNREAL*) malloc(n2Qsq*sizeof(GNREAL));
			if(Ham2QAr[i]==NULL) error = 1;
		} else if(Ham2QAr!=NULL) Ham2QAr[i] = NULL;
		if((!error) && (do2d) ) {
			Evals2QAr[i] = (GNREAL*) malloc(n2Q*sizeof(GNREAL));
			if(Evals2QAr[i]==NULL) error = 1;
		} else if(Evals2QAr!=NULL) Evals2QAr[i] = NULL;

		for(j=0; j<3; j++) {
			if((!error) && (do2d) && ((!pert) || pertvec )) {
				// cjfeng 04/05/2016
				// Tranpose the matrix here
				// Dip2QAr[j][i] = (GNREAL*) malloc(n2Q*nosc*sizeof(GNREAL));
				// The original line
				Dip2QAr[j][i] = (GNREAL*) malloc(nosc*n2Q*sizeof(GNREAL));
				if(Dip2QAr[j][i]==NULL) error = 1;
			} else if(Dip2QAr[j]!=NULL) Dip2QAr[j][i] = NULL;
		}
	}
	return error; 
}

// Close files, free memory and exit gracefully
int graceful_exit( int error, int nbuffer, int win2d, int nthreads, int npol, int nise, int nosc ) {
	int i,j;
	int stick;
	if(nise) stick = 0;
	else stick = 1;

	if(error==1) printf("Error allocating memory.\n");
	if(error==2) printf("Error opening files.\n");

	if(NNames!=NULL) for(i=0; i<nosc; i++) if(NNames[i]!=NULL) free(NNames[i]);
	if(CNames!=NULL) for(i=0; i<nosc; i++) if(CNames[i]!=NULL) free(CNames[i]);
	if(NNames!=NULL) free(NNames);
	if(CNames!=NULL) free(CNames);
	if(NNums!=NULL) free(NNums);
	if(CNums!=NULL) free(CNums);

	if(ifp!=NULL) fclose(ifp);
	if(hfp!=NULL) fclose(hfp);
	if(sfp!=NULL) fclose(sfp);
	for(i=0; i<3; i++) if(Dfp[i]!=NULL) fclose(Dfp[i]);
	for(i=0; i<4; i++) if(REPH[i]!=NULL) free(REPH[i]);
	for(i=0; i<4; i++) if(NetREPH[i]!=NULL) free(NetREPH[i]);
	for(i=0; i<4; i++) if(NREPH[i]!=NULL) free(NREPH[i]);
	for(i=0; i<4; i++) if(NetNREPH[i]!=NULL) free(NetNREPH[i]);
	if(ffp!=NULL) fclose(ffp);
	if(lfp!=NULL) fclose(lfp);
	if(afp!=NULL) fclose(afp);
	for(i=0; i<10; i++) if(Trajfp[i]!=NULL) fclose(Trajfp[i]);
	for(i=0; i<npol; i++) if(rfp[POL[i]]!=NULL) fclose(rfp[POL[i]]);
	for(i=0; i<npol; i++) if(nrfp[POL[i]]!=NULL) fclose(nrfp[POL[i]]);
	if(U1Q!=NULL) free(U1Q);
	if(U1Qs!=NULL) free(U1Qs);
	for(i=0; i<3; i++) if(cDip1Q[i]!=NULL) free(cDip1Q[i]);
	if(U1QMat!=NULL) for(i=0; i<nbuffer; i++) if(U1QMat[i]!=NULL) free(U1QMat[i]);
	if(U1QMat!=NULL) free(U1QMat);
	if(U2Qs!=NULL) free(U2Qs);
	if(U2QMat!=NULL) for(i=0; i<nbuffer; i++) if(U2QMat[i]!=NULL) free(U2QMat[i]);
	if(U2QMat!=NULL) free(U2QMat);
	if(CorrFunc!=NULL) free(CorrFunc);
	if(NetCorrFunc!=NULL) free(NetCorrFunc);
	if(hann!=NULL) free(hann);
	if(popdecay1Q!=NULL) free(popdecay1Q);
	if(popdecay2Q!=NULL) free(popdecay2Q);
	if(FTin1D!=NULL) free(FTin1D);
	if(FTout1D!=NULL) free(FTout1D);
	if(FTplan1D!=NULL) fftw_destroy_plan(FTplan1D);
	if(FTin2D!=NULL) fftw_free(FTin2D);
	if(FTout2D!=NULL) fftw_free(FTout2D);
	if(FTplan2D!=NULL) fftw_destroy_plan(FTplan2D);
	for(i=0; i<3; i++) if(psi_a[i]!=NULL) for(j=0; j<win2d; j++) if(psi_a[i][j]!=NULL) free(psi_a[i][j]);
	for(i=0; i<3; i++) if(psi_b1[i]!=NULL) for(j=0; j<win2d; j++) if(psi_b1[i][j]!=NULL) free(psi_b1[i][j]);
	for(i=0; i<3; i++) if(psi_b12[i]!=NULL) for(j=0; j<win2d; j++) if(psi_b12[i][j]!=NULL) free(psi_b12[i][j]);
	for(i=0; i<3; i++) if(psi_a[i]!=NULL) free(psi_a[i]);
	for(i=0; i<3; i++) if(psi_b1[i]!=NULL) free(psi_b1[i]);
	for(i=0; i<3; i++) if(psi_b12[i]!=NULL) free(psi_b12[i]);
	for(i=0; i<9; i++) if(psi_ca[i]!=NULL) free(psi_ca[i]);
	for(i=0; i<9; i++) if(psi_cb[i]!=NULL) free(psi_cb[i]);
	if(psi2Q!=NULL) free(psi2Q);
	if(ftir!=NULL) free(ftir);
	if(netftir!=NULL) free(netftir);
	for(j=0; j<3; j++) if(ExDip2QAr[j]!=NULL) for(i=0; i<nthreads; i++) if(ExDip2QAr[j][i]!=NULL) free(ExDip2QAr[j][i]);
	for(j=0; j<3; j++) if(ExDip2QAr[j]!=NULL) free(ExDip2QAr[j]);

	if(Ham1QMat!=NULL) for(i=0; i<nbuffer; i++) if(Ham1QMat[i]!=NULL) free(Ham1QMat[i]);
	if(Ham1QMat!=NULL) free(Ham1QMat);
	for(i=0; i<3; i++) {
		if(Dip1QMat[i]!=NULL) for(j=0; j<nbuffer; j++) if(Dip1QMat[i][j]!=NULL) free(Dip1QMat[i][j]);
		if(Dip1QMat[i]!=NULL) free(Dip1QMat[i]);
	}

	if(SHIFTNDX!=NULL) free(SHIFTNDX);
	if(SHIFT!=NULL) free(SHIFT);
	for(i=0; i<nthreads; i++) {
		if(Ham1QAr!=NULL) if(Ham1QAr[i]!=NULL) free(Ham1QAr[i]);
		if(Evals1QAr!=NULL) if(Evals1QAr[i]!=NULL) free(Evals1QAr[i]);
		for(j=0; j<3; j++) if(Dip1QAr[j]!=NULL) if(Dip1QAr[j][i]!=NULL) free(Dip1QAr[j][i]);
		for(j=0; j<3; j++) if(ExDip1QAr[j]!=NULL) if(ExDip1QAr[j][i]!=NULL) free(ExDip1QAr[j][i]);
		if(Ham2QAr!=NULL) if(Ham2QAr[i]!=NULL) free(Ham2QAr[i]);
		if(Evals2QAr!=NULL) if(Evals2QAr[i]!=NULL) free(Evals2QAr[i]);
	}
	if(Ham1QAr!=NULL) free(Ham1QAr);
	if(Evals1QAr!=NULL) free(Evals1QAr);
	for(j=0; j<3; j++) if(Dip1QAr[j]!=NULL) free(Dip1QAr[j]);
	for(j=0; j<3; j++) if(ExDip1QAr[j]!=NULL) free(ExDip1QAr[j]);

	if(Ham2QAr!=NULL) free(Ham2QAr);
	if(Evals2QAr!=NULL) free(Evals2QAr);
	for(j=0; j<3; j++) if(Dip2QAr[j]!=NULL) for(i=0; i<nthreads; i++) if(Dip2QAr[j][i]!=NULL) free(Dip2QAr[j][i]);
	for(j=0; j<3; j++) if(Dip2QAr[j]!=NULL) free(Dip2QAr[j]);

	// cjfeng 04/05/2016 	
	// The modified code due to changing the allocation of Dip2QMat
	if(Dip2QMat!=NULL) for(j=0; j<nbuffer; j++) if(Dip2QMat[j]!=NULL) for(i=0; i<3; i++) if(Dip2QMat[j][i]!=NULL) free(Dip2QMat[j][i]);
	if(Dip2QMat!=NULL) for(j=0; j<nbuffer; j++) if(Dip2QMat[j]!=NULL) free(Dip2QMat[j]);
	if(Dip2QMat!=NULL) free(Dip2QMat);

	if(SitesBuffer!=NULL) free(SitesBuffer);
	return 0;
}

/******************************\
 * Main program								*
\******************************/

int main ( int argc, char * argv[] ) {

/***************************************************************
 * Initial variable declaration and applying default settings  *
 ***************************************************************/
	double starttime, stoptime;		// Start and stop time during spectral simulation.
	char* deffnm = NULL;					// Input file name base
	char* outname = NULL;					// Output file name base
	char* hamnm = "ham.txt";			// Hamiltonian file suffix
	char* dipxnm = "dipx.txt";		// Dipole moment x-component file suffix
	char* dipynm = "dipy.txt";		// Dipole moment y-component file suffix
	char* dipznm = "dipz.txt";		// Dipole moment z-component file suffix
	char* sitesnm = "";						// Site energy file
	char* paramnm = "";						// Isotope shift file
	int nise = 0;									// Static averaging by default (0)
	int reph = 1;									// Calculate rephasing spectrum
	int nreph = 1;								// Calculate non-rephasing spectrum
	int zzzz = 1;									// ZZZZ polarization
	int zzyy = 1;									// ZZYY polarization
	int zyyz = 0;									// No ZYYZ polarization
	int zyzy = 0;									// No ZYZY poloarization
	int do2d = 0;									// No 2D IR spectral simulation.
	int pert = 0;									// No first-order perturbative correction on site energies
	int nthreads = 1;							// Number of threads
	int nread = 4*nthreads;				// nread for determining buffer length.
	int tscan = 0; 								// Scan time for NISE or Averaging window time for TAA in fs
	int window;										// Window size for TAA or linear NISE.
	int win2d;										// Window size for 2D NISE.
	int winzpad;									// Zero-padding length
	int whann = 1;								// Hann window applied on the response in NISE time domain or averaged Hamiltonian in TAA.
	int dump = -1;								// IR spectra dump time in ps
	int tstep = 0;								// Trajectory time step in fs
	int T2 = 0;										// Waiting time in fs, only valid in 2D NISE so far.
	int skip = 1;									// Number of skips between reference frames.
	real delta = 16.0;						// Anharmonicity in wavenumbers (Weak anharmonic approximation applied in the program.)
	int TauP = 1300;							// Amide I lifetime in fs

	int nosc = 0;									// Number of oscillators

	real wstart = 1500.0;					// Starting frequency of interest in wavenumber.
	real wstop  = 1800.0;					// Ending frequency of interest in wavenumber.
	real wres = 1;								// Frequency resolution in wavenumber.
	real winterp = 1;							// Output frequency spacing in wavenumber, but for NISE, it determines 
																// only zero-padding length, not the true resolution.

	double c = 2.9979e-5;					// Speed of light in cm/fsec
	double pi = 3.14159265;				// The ration of a circle's circumference

	int error = 0;								// Integer indicating error information.
	int nbuffer, nbuffer1d;				// nbuffer determines the size of Hamiltonian and dipole moment array sizes 
																// while nbuffer1d is used only for FTIR correlation function scanning.
	int npol;											// Number of polarizations to be simulated.

	// cjfeng 07/12/2016
	// Include switches on simulating pathways
	int reph1 = 1;
	int reph2 = 1;
	int reph3 = 1;
	int nreph1 = 1;
	int nreph2 = 1;
	int nreph3 = 1;
	// cjfeng 07/12/2016
	// Include switch on printing time-domain response
	int tresponse = 0;

	// cjfeng 06/27/2016
	// Added usage to track memory usage.
	struct rusage r_usage;				// resource usage, please refer to getrusage manual page.

	// A list of command line file flags
	t_filenm fnm[] = { };

	// A list of additional command line arguments
	t_pargs pa [] = {
		{"-deffnm", FALSE, etSTR, {&deffnm},
		 "Input file name base"},
		{"-outname", FALSE, etSTR, {&outname}, 
		 "Base for output file name." }, 
		{"-ham", FALSE, etSTR, {&hamnm},
		 "Hamiltonian trajectory file."},
		{"-dipx", FALSE, etSTR, {&dipxnm},
		 "Dipole moment x-component file."},
		{"-dipy", FALSE, etSTR, {&dipynm}, 
		 "Dipole moment y-component file."},
		{"-dipz", FALSE, etSTR, {&dipznm},
		 "Dipole moment z-component file."},
		{"-sites", FALSE, etSTR, {&sitesnm},
		 "Site energy file (replaces diagonal elements in Hamiltonian)."},
		{"-shift", FALSE, etSTR, {&paramnm},
		 "Isotope shift file (optional). Note that oscillator indexing begins at zero, not one."},
		{"-reph", FALSE, etBOOL, {&reph},
		 "Calculate rephasing spectrum"},
		{"-nreph", FALSE, etBOOL, {&nreph},
		 "Calculate non-rephasing spectrum"},
		{"-zzzz", FALSE, etBOOL, {&zzzz},
		 "Calculate ZZZZ polarization"},
		{"-zzyy", FALSE, etBOOL, {&zzyy},
		 "Calculate ZZYY polarization"},
		{"-zyyz", FALSE, etBOOL, {&zyyz},
		 "Calculate ZYYZ polarization"},
		{"-zyzy", FALSE, etBOOL, {&zyzy},
		 "Calculate ZYZY polarization"},
		{"-2dir", FALSE, etBOOL, {&do2d},
		 "Calculate 2DIR spectrum"},
		{"-nise", FALSE, etBOOL, {&nise},
		 "Static averaging (no dynamics)"},
		{"-pert", FALSE, etBOOL, {&pert},
		 "Use perturbative approximation for energies"},
		//{"-pertvec", FALSE, etBOOL, {&pertvec},
		// "Use perturbative approximation for vectors"},
		{"-skip", FALSE, etINT, {&skip}, 
		 "Number of skips between reference frames"},
		{"-nt", FALSE, etINT, {&nthreads},
		 "Number of threads"},
		{"-delta", FALSE, etREAL, {&delta},
		 "Anharmonicity (wavenumbers)"},
		{"-taup", FALSE, etINT, {&TauP},
		 "Population lifetime (fs)"},
		{"-tau2", FALSE, etINT, {&T2},
		 "Waiting time (fs). Only for NISE simulations."},
		{"-tstep", FALSE, etINT, {&tstep},
		 "Trajectory time step (fs)."},
		{"-tscan", FALSE, etINT, {&tscan},
		 "Scan time for NISE or averaging time for TAA (fs)"},
		{"-dump", FALSE, etINT, {&dump},
		 "FTIR dump time (ps)"},
		{"-wstart", FALSE, etREAL, {&wstart},
		 "Frequency range start"},
		{"-wstop", FALSE, etREAL, {&wstop},
		 "Frequency range stop"},
		{"-winterp", FALSE, etREAL, {&winterp},
		 "Output frequency spacing. NB: for NISE calculations This determines only zero-pad length, NOT the true resolution (which is determined by the -tscan flag). "},
		{"-whann", FALSE, etBOOL, {&whann},
		 "Use Hann window in time domain. In NISE simulations, this windows the response. In static calculations, this weights the averaged Hamiltonians."}, 
		// cjfeng 07/12/2016
		// Including switches on pathways
		// {"-reph1", FALSE, etBOOL, {&reph1},
		//  "Simulate rephasing pathway 1 in NISE."},
		// {"-reph2", FALSE, etBOOL, {&reph2},
		//  "Simulate rephasing pathway 2 in NISE."},
		// {"-reph3", FALSE, etBOOL, {&reph3},
		//  "Simulate rephasing pathway 3 in NISE."},
		// {"-nreph1", FALSE, etBOOL, {&nreph1},
		//  "Simulate non-rephasing pathway 1 in NISE."},
		// {"-nreph2", FALSE, etBOOL, {&nreph2},
		//  "Simulate non-rephasing pathway 2 in NISE."},
		// {"-nreph3", FALSE, etBOOL, {&nreph3},
		//  "Simulate non-rephasing pathway 3 in NISE."},
		{"-tresponse", FALSE, etBOOL, {&tresponse},
		 "Print time domain repsonse."},
		};	

	// The program description
	const char *desc[] = {"NB: 1. OpenMP parallel computation is not yet implemented in static averaging and time-averaging approximation scheme. "};

	// A description of known bugs
	const char *bugs[] = {""};

  output_env_t oenv;

	// parse_common_args() does exactly what it sounds like. The arg list is provided in the arrays t_filenm and t_pargs defined above. 
	parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_BE_NICE, asize(fnm), fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv);

	// Vector perturbation is slow.  Discontinue as of 1/31/2016
	//if(pertvec) pert = 1;

	// Variable assignments after passing arguments.
	nread = 4*nthreads;
	npol = zzzz + zzyy + zyyz + zyzy;	// Number of polarizations	

	// Checking reasonble simulation conditions are provided.
	if( (tstep<=0) && (tscan>0) ) {
		printf("Please supply a (positive) trajectory time step in fsec (-tstep flag).\n");
		graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
		return 0;
	} else if (tscan==0) tstep = 20;
	if( (tscan<0) ) {
		printf("Please supply a (non-negative) scan time step in fsec (-tscan flag).\n");
		printf("(Enter 0 for a single frame in TAA).\n");
		graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
		return 0;
	}
	if( (T2!=0) && (!nise) ) {
		printf("Error: Waiting time T2 must be zero for static spectral calculations.\n");
		graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
		return 0;
	}

	if( ((reph+nreph)==0) && (do2d) ) {
		printf("Nothing to calculate! Please request either rephasing or non-rephasing.\n");
		graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
		return 0;
	}
	if( (npol==0) && (do2d) ) {
		printf("Nothing to calculate! Please select at least one polarization condition.\n");
		printf("Note that the default values for ZYYZ and ZYZY are false.\n");
		graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
		return 0;
	} else {
		int count = 0;
		if(zzzz) {
			POL[count] = 0;
			count++;
		}
		if(zzyy) {
			POL[count] = 1;
			count++;
		}
		if(zyyz) {
			POL[count] = 2;
			count++;
		}
		if(zyzy) {
			POL[count] = 3;
			count++;
		}
	}
	
/***************************************************************************************
 Hamiltonian and Dipole files: check for consistency and open pointers 
***************************************************************************************/

	int i,j,k;
	char ftirnm[maxchar];			// FTIR file name
	char axisnm[maxchar];			// Frequency axis file name
	char lognm[maxchar];			// Log file name
	char rephnm[4][maxchar];	// Rephasing spectra file names
	char nrephnm[4][maxchar];	// Non-rephasing spectra file names
	char Trajnm[10][maxchar];	// Spectral trajectory file names
	char polnm[4][16] = { "zzzz.txt", "zzyy.txt", "zyyz.txt", "zyzy.txt" };
	int vals, n2Q;						// vals: number of entries, n2Q: number of 2Q states.
	int nframes;							// nframes: number of frames
	// Infoname: info file name, Hamname: Hamiltonian file name, Parname: isotope shift file name, Sitesname: site energy file name
	char Infoname[maxchar], Hamname[maxchar], Dipnames[3][maxchar], Parname[maxchar], Sitesname[maxchar];
	int nshift;								// Number of isotope shifts
	int npts;									// Number of frequency bins in static calculations.
	int do_traj = (dump>=1);	// Dumping spectral data as a spectral trajectory or not

	double wo = 0.0;
	double expfac = -2.0*pi*tstep*c;	// Exponential factor 

	// Define array dimensions. 
	wres = winterp;
	npts = ceil((wstop-wstart)/wres)+1;	// Used for static calculations.
	
	// For NISE FTIR and TAA 2DIR calculations, window is the number
	// of frames involved in each spectral calculation (the correlation
	// function window in FTIR or the time averaging window in TAA).
	// For NISE 2DIR calculations, this is the correlation function 
	// window for each dimension; the actual actual number of frames 
	// stored is win2d defined below. 
	window = ((int) tscan/tstep) + 1;
	
	// win2d and is used only for NISE 2D calculations
	win2d = 2*window + (int) T2/tstep;
	
	// For static calculations, winzpad is to prevent wraparound error 
	// when dressing with a lifetime-limited lorentzian. 
	if(nise) winzpad = (int) (1)/(winterp*c*tstep) + 1;
	else winzpad = 5*npts;
	
	// The buffer length must be nread-1 frames longer than the required time window
	// so that we can hold nread complete windows in memory at once. 
	// cjfeng 06/29/2016
	// nbuffer definition changes when 2D NISE simulation is performed, resulting in 
	// incorrect scanning of FTIR correlation function. nbuffer1d is declared to
	// correct the indicies when computing CorrFunc[n].
	if(do2d && nise) nbuffer = win2d + nread - 1;
	else nbuffer = window + nread - 1;
	nbuffer1d = window + nread -1;

	// Define frequency axis (only actually used for NISE calculations). 
	GNREAL dw = (1.0) / (winzpad*tstep*c); 	// Transform resolution in cm-1
	int maxit = 1e+9;
	int ndxstart = -1;
	GNREAL tol = dw/2.0;
	GNREAL wdif = 2*tol;
	while( (wdif>tol)  && (ndxstart<maxit) ) {
		ndxstart++;
		if( (ndxstart*dw-wstart)>=0 ) wdif = ndxstart*dw-wstart;
		else if( (ndxstart*dw-wstart)<0 ) wdif = wstart-ndxstart*dw;
	}	
	if(wdif>tol) printf("Error finding start frequency. Frequency axis should not be trusted!\n");
	int nprint = (int) ( wstop - wstart )/dw;		// Number of frequency bins to be printed.

	// Set file names for opening. 
	make_fnames( Infoname, Hamname, hamnm, Dipnames, dipxnm, dipynm, dipznm, Sitesname, sitesnm, axisnm, ftirnm, lognm, polnm, rephnm, nrephnm, Trajnm, Parname, paramnm, outname, deffnm, do_traj );

	// Read info file
	int info  = read_info( ifp, Infoname, &NNames, &NNums, &CNames, &CNums);
	if(info<0) { 
		graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
		return 0;
	} else if(info>0) {
		printf("Successfully read info file. Will be looking for %d oscillators in input files. \n", info);
	}

	// Read Param file
	if(strlen(paramnm)>0) {
		printf("Reading parameters from file %s.\n", Parname);
		if(!parse_shifts(Parname, &nshift, maxchar)) {
			printf("Error parsing parameter file %s. Please check input.\n", Parname);
			graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
			return 0;
		}
		printf("Loaded %d oscillator shifts: \n", nshift);
		if(info!=0) for(i=0; i<nshift; i++) printf("%s %d (index %d) :\t%6.2f cm-1.\n", NNames[SHIFTNDX[i]], NNums[SHIFTNDX[i]], SHIFTNDX[i], SHIFT[i]);
		else for(i=0; i<nshift; i++) printf("Oscillator %d: \t%6.2f cm-1.\n", SHIFTNDX[i], SHIFT[i]);
	} else {
		nshift = 0;
	}

	// Check validity of Hamiltonian file.
	nframes = count_lines(Hamname);
	printf("Located %d lines in input file %s\n", nframes, Hamname);

	if(nframes<nbuffer-nread+1) {
		printf("Error: Not enough (%d) frames for accurate averaging with requested window size (%d).\n", nframes, window);
		printf("At least %d frames have to be provided to compute spectra for %d thread(s).\n", nbuffer-nread+1, nthreads);
		graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
		return 0;
	}

	vals = count_entries(Hamname);
	nosc = floor(sqrt(vals)+0.5);
	n2Q = nosc*(nosc+1)/2;
	printf("Located %d oscillators in input file %s\n", nosc, Hamname);
	if( (nosc!=info) && (info!=0) ) {
		printf("Error! Info file specifies %d oscillators, but found %d in Hamiltonian file. \n", info, nosc);
		graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
		return 0;
	}
	for(i=0; i<nshift; i++) {
		if(SHIFTNDX[i]>=nosc) {
			printf("Error! Requested shift index (oscillator %d) is larger than total number of oscillators.\n", SHIFTNDX[i]);
			graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
			return 0;
		}
	}

	// Checking the consistency between Hamiltonian file and site energy file.
	if( (strlen(Sitesname)>0) ) {
		vals = count_entries(Sitesname);
		if(nosc!=vals) {			// Checking number of oscillators
			printf("Error! Different number of oscillators (%d vs. %d) located in Hamiltonian file %s and sites file %s\n", nosc, vals, Hamname, Sitesname); 
			graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
			return 0;
		}
		vals = count_lines(Sitesname);
		if((nframes!=vals)) {	// Checking number of frames
			printf("Error! Different number of lines (%d vs. %d) in Hamiltonian file %s and sites file %s\n", nframes, vals, Hamname, Sitesname);
			graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
			return 0;
		}
	}
	
	// Checking the consistency between Hamiltonian file and dipole moment files.
	for(i=0; i<3; i++) {
		vals = count_entries(Dipnames[i]);
		if(vals!=nosc) {			// Checking number of oscillators
			printf("Error! Different number of oscillators (%d vs. %d) in Hamiltonian file %s and dipole file %s\n", vals, nosc, Hamname, Dipnames[i]);
			graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
			return 0;
		}
		vals = count_lines(Dipnames[i]);
		if(vals!=nframes) {		// Checking number of frames
			printf("Error! Different number of lines (%d vs. %d) in Hamiltonian file %s and dipole file %s\n", vals, nframes, Hamname, Dipnames[i]);
			graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
			return 0;
		}
	}
	
	// Assigning dump frequency
	if(dump<1) dump = nframes;
	else dump = (int) ( 1000*dump / tstep );	// dump frequency in ps

	// Open files
	if(!error) error = open_all( Hamname, Dipnames, Sitesname, axisnm, ftirnm, lognm, npol, POL, rephnm, nrephnm, Trajnm, do2d, do_traj );
	if(!error) printf("Finished opening files.\n");

	// Allocate memory
	if(!pert && do2d) printf("Allocating memory for %d two-quantum states.\n", n2Q);
	error = allocate_all( window, win2d, winzpad, nosc, nbuffer, nthreads, n2Q, pert, do2d, tstep, TauP, npts, nise );
	if(!error) printf("Finished allocating memory.\n");

	// cjfeng 07/12/2016
	// Time-domain response function
	if(tresponse && do2d) {
		if (outname!=NULL) {
			for(i=0; i<8; i++) {
				strcpy(fidname[i],outname);
				strcat(fidname[i],"_");
				strcat(fidname[i],fidbase[i]);
			}
		}
		else for(i=0; i<8; i++) strcpy(fidname[i],fidbase[i]);
		for (i=0; i<npol; i++) {
			printf("%s\n",fidname[POL[i]]);
			printf("%s\n",fidname[POL[i]+4]);
		}
		for(i=0; i<npol; i++) {
			fid_reph[i] = fopen(fidname[POL[i]],"w");
			fid_nreph[i] = fopen(fidname[POL[i]+4],"w");
		}
	}

	// cjfeng 06/27/2016
	// Tracking memory usage.
	getrusage(RUSAGE_SELF, &r_usage);
	printf("Memory usage = %ld kB.\n", r_usage.ru_maxrss);

/***************************************************************************************\
 * Starting spectral simulation part
 ***************************************************************************************/
	// And go!
	int alldone = 0;	// Flag for indicating completing all of the calculations.
	int readframe = 0;// End point of the read frames
	int fr = 0;				// index for referencing reference point during spectral simulation
	int frame = 0;
	int nAvgd = 0; // Used only for static calculation
	// Setting number of threads for further parallel computation
	#if OMP_PARALLEL 
	omp_set_num_threads(nthreads);
	#endif

	// Measuring start time.
	starttime = omp_get_wtime();
	if(!error) printf("Beginning trajectory parsing...\n");
	// Step through trajectory until all frames have been processed
	while( (!error) & (!alldone) ) {
		int justread = 0;		// location of the last read frame
		// Read 1Q data for the next nread frames and store in Ham1QMat and Dip1QMat
		int frnum = MIN(nframes,readframe+nread);
		for(fr=readframe; fr<frnum; fr++) {
		// for(fr=readframe; fr<MIN(nframes,readframe+nread); fr++) {
			// cjfeng 06/27/2016
			// Reduce the number of modulus operations.
			// Modulus is required to fill in new frames when exceeding memory allocation size
			// of the arrays, determined by window or win2d.
			int xfr = fr%nbuffer;

			// First reading Hamiltonian
			if(!read_line(hfp, nosc, nosc, Ham1QMat[xfr])) {
				printf("Error reading from Hamiltonian file.\n");
				error = 1;
				break;
			}
			// Now dipole moments
			for(i=0; i<3; i++) if(!read_line(Dfp[i], nosc, 1, Dip1QMat[i][xfr])) {
				printf("Error reading from Dipole file.\n");
				error = 1;
				break;
			}
			// Site energies if specified
			if(sfp!=NULL) {
				if(!read_line(sfp, nosc, 1, SitesBuffer)) {
					printf("Error reading from site energy file.\n");
					error = 1; 
					break; 	
				} 
				else for(i=0; i<nosc; i++) Ham1QMat[xfr][i*nosc+i] = SitesBuffer[i];
			}
			// Add any specified shifts to the Hamiltonian matrix
			for(i=0; i<nshift; i++) {
				Ham1QMat[xfr][SHIFTNDX[i]*nosc+SHIFTNDX[i]] = Ham1QMat[xfr][SHIFTNDX[i]*nosc+SHIFTNDX[i]] + SHIFT[i]; 
			}
			// Generate two-quantum dipoles
			// 04/05/2016 cjfeng
			// Swapped the order of direction and frame due to changing the allocation of Dip2QMat.
			if( (do2d) ) {
				for(i=0; i<3; i++) gen_dip_2Q(Dip1QMat[i][xfr], Dip2QMat[xfr][i], nosc, n2Q);
			}
		}
		// Check if we've reached the end.
		if(fr==nframes) alldone = 1;

		// Note how many new frames have been read.
		justread = fr - readframe;

		// And update the total number of read frames. 
		readframe = fr;

/**********************************************
 * Numerical wavefunction propagation
 **********************************************/
		// Now process the new frames. 
		if( (!error) && (nise)) {
			// Again we use fr as our frame counter. It runs from readframe-justread 
			// (the first of the newly read frames) to readframe-1. 
			
			// cjfeng 03/31/2016
			// The parallel section won't be executed under single thread to avoid overhead.
			#if OMP_PARALLEL 
			#pragma omp parallel if(nthreads>1) shared(justread, Ham1QMat, Evals1QAr, U1QMat) private(fr, i, j) firstprivate(readframe, wo, nosc, expfac)
			#endif
			{
				// cjfeng 06/27/2016
				// Use temporary array instead for solving eigenvalue and eigenvectors.
				int noscsq = nosc*nosc;
				GNREAL *Evals1Qtmp, *Ham1Qtmp;
				Evals1Qtmp = (GNREAL*) malloc(nosc*sizeof(GNREAL));
				Ham1Qtmp = (GNREAL*) malloc(noscsq*sizeof(GNREAL));
				int *isuppz1Q = (int*) malloc(2*nosc*sizeof(int));
				#if OMP_PARALLEL
				#pragma omp for schedule(guided) nowait
				#endif
				for(fr=readframe-justread; fr<readframe; fr++) {
					// cjfeng 06/27/2016
					// Reduce modulus operations.
					int xfr = fr%nbuffer;
					int tid;
					double val;
					lapack_int info;
					tid = omp_get_thread_num();
					// Find one-quantum eigenvalues
					// Note that Ham1QMat[fr%nbuffer] now contains eigenvectors
					// cjfeng 06/27/2016
					// Use temporary array to solve eigenvalues and eigenvector,
					// and the Relatively Robust Representation to compute eigenvalues and eigenvectors.
					for (i=0; i<noscsq; i++) Ham1Qtmp[i] = Ham1QMat[xfr][i];
					info = SYEVR (LAPACK_ROW_MAJOR, 'V', 'A', 'U', nosc, Ham1Qtmp, nosc, wstart, wstop, 0, nosc-1, 0.00001, &nosc, Evals1Qtmp, Ham1Qtmp, nosc, isuppz1Q);

					// Note that our Hamiltonian is actually H/(h*c)
					// The exponent we calculate is -i*2*pi*tstep*c*Ham1Q
			
					// cjfeng 06/27/2016
					// Reduce for loop overhead.	
					for(i=0; i<noscsq; i++) U1QMat[xfr][i] = 0.0; 
					for(i=0; i<nosc; i++) {
						GNREAL cre, cim; 
						int ni = i*nosc;
						for(j=0; j<nosc; j++) {
							int k;
							int nj = j*nosc;
							cre = 0.0;
							cim = 0.0;
							for(k=0; k<nosc; k++) {
								// cjfeng 06/27/2016
								// Use temporary array instead.
								cre += Ham1Qtmp[ni+k]*Ham1Qtmp[nj+k]*cos(expfac*(Evals1Qtmp[k]-wo));
								cim += Ham1Qtmp[ni+k]*Ham1Qtmp[nj+k]*sin(expfac*(Evals1Qtmp[k]-wo));
							}
							U1QMat[xfr][ni+j] = cre + I*cim;
							// printf("%6.5f + (%6.5f)i\t", creal(U1QMat[fr%nbuffer][i*nosc+j]), cimag(U1QMat[fr%nbuffer][i*nosc+j]));
						}
						// printf("\n");
					}
					//printf("\n");
				}
				free(Evals1Qtmp);
				free(Ham1Qtmp);
				free(isuppz1Q);
			}
			// cjfeng 06/27/2016
			// Uncouple the 2D eigenvalue procedure from 1Q part.
			if ( do2d ) {	// Generating 2Q Propagator.
				#if OMP_PARALLEL 
				#pragma omp parallel if(nthreads>1) shared(justread, Ham1QMat, Evals2QAr, U1QMat, U2QMat) private(fr, i, j) firstprivate(readframe, wo, n2Q, do2d, pert, expfac, delta, nosc)
				#endif
				{
					// cjfeng 06/27/2016
					// Use temporary array.
					int n2Qsq = n2Q*n2Q;
					GNREAL *Evals2Qtmp, *Ham2Qtmp;
					Evals2Qtmp = (GNREAL*) malloc(n2Q*sizeof(GNREAL));
					Ham2Qtmp = (GNREAL*) malloc(n2Qsq*sizeof(GNREAL));
					int *isuppz2Q = (int*) malloc(2*n2Q*sizeof(int));
					#if OMP_PARALLEL
					#pragma omp for schedule(guided) nowait
					#endif
					for(fr=readframe-justread; fr<readframe; fr++) {
						// cjfeng 06/27/2016
						// Reduce modulus operations.
						int xfr = fr%nbuffer;
						int tid;
						double val;
						lapack_int info;
						tid = omp_get_thread_num();
						if( pert ) {
							// Skip eigenvalue calculation and use perturbative approximation to generate U2Q
							gen_pert_prop( U1QMat[xfr], U2QMat[xfr], nosc, n2Q, expfac, delta);
						}
						else {	// Two-quantum eigenvalues
					 		gen_ham_2Q(Ham1QMat[xfr], nosc, Ham2Qtmp, n2Q, delta);
							// printf("Starting eigenvalue calculation for thread %d\n", tid);

							info = SYEVR (LAPACK_ROW_MAJOR, 'V', 'A', 'U', n2Q, Ham2Qtmp, n2Q, wstart, wstop, 0, n2Q-1, 0.00001, &n2Q, Evals2Qtmp, Ham2Qtmp, n2Q, isuppz2Q);
							// printf("Finishing eigenvalue calculation for thread %d\n", tid);
							for(i=0; i<n2Qsq; i++) U2QMat[xfr][i] = 0.0; 
							for(i=0; i<n2Q; i++) {
								int ni = i*n2Q;
								GNREAL cre, cim; 
								for(j=0; j<n2Q; j++) {
									int nj = n2Q*nosc;
									int k;
									cre = 0.0;
									cim = 0.0;
									for(k=0; k<n2Q; k++) {
										cre += Ham2Qtmp[ni+k]*Ham2Qtmp[nj+k]*cos(expfac*(Evals2Qtmp[k]-wo));
										cim += Ham2Qtmp[ni+k]*Ham2Qtmp[nj+k]*sin(expfac*(Evals2Qtmp[k]-wo));
									}
									U2QMat[xfr][ni+j] = cre + I*cim;
									// cjfeng 07/05/2016
									// Absolute two-quantum frequency verified
									// printf("%6.5f + (%6.5f)i\t", creal(U2QMat[fr%nbuffer][i*n2Q+j]), cimag(U2QMat[fr%nbuffer][i*n2Q+j]));
								}
								// printf("\n");
							}
						} 
					}
					free(Evals2Qtmp);
					free(Ham2Qtmp);
					free(isuppz2Q);
				}
			}

			// FTIR correlation function computation.
			fr = -1;
			frame = readframe-nread;
			// Now do calculations for selected frames. 
			while( (fr!=-2) && (!error) ) {
				// Step through all recently read frames to see whether frame
				// is suitable as the end point of a calculation. 
				// The time window twin is either window or win2d, depending 
				// on whether do2d = 0 or 1. 
				int twin = nbuffer1d-nread+1;
				while(frame<readframe) {
					if((frame%100)==0) {
						printf("Frame: %d\n", frame);
						// cjfeng 06/27/2016
						// Tracking memory usage.
						getrusage(RUSAGE_SELF, &r_usage);
						printf("Memory usage = %ld kB.\n", r_usage.ru_maxrss);
						fprintf(lfp, "Frame: %d\n", frame);
						fflush(lfp);
					}
					// We check if at least window or win2d frames have been read 
					// and if frame is a multiple of skip. 
					if( ((frame%skip)==0) && (frame>=twin-1) ) {
						// frame will be the endpoint of the calculation. 
						// It will start at frame fr = frame-twin+1.
						fr = (frame-twin+1)%nbuffer;
						break;
					} else frame++;
				}
				// If we got all the way through the fr loop without finding anything
				// set fr = -2 to break out of the outside loop.
				if(frame==readframe) fr = -2;
				if( (fr>=0) && (nise) ) {
					// Do a dynamic calculation using fr as our starting point and augment frame.
					frame++;  // We don't reference frame further in the calculation. 

					int xfr = 0;
					int n,k;
					GNCOMP cval;
					int noscsq = nosc*nosc;
					
					// cjfeng 06/27/2016
					// Reduce for loop overhead.
					
					for(i=0; i<noscsq; i++) U1Q[i] = 0.0;
					for(i=0; i<nosc; i++) U1Q[i*nosc+i] = 1.0;
					
					for(n=0; n<window; n++) {
						xfr = (fr+n)%nbuffer;
						cval = 0.0;
						#if OMP_PARALLEL
						#pragma omp parallel if(nosc >= NBCOMP && nthreads>1) shared(nosc, cDip1Q, U1Q, Dip1QMat, xfr, fr) private(i,j,k)
						// cjfeng 04/06/2016
						// reduction used to avoid race condition
						#pragma omp for schedule(guided) collapse(2) reduction(+:cval)
						#endif
						for(i=0; i<3; i++) {
							for(j=0; j<nosc; j++) {
								int nj = j*nosc;
								cDip1Q[i][j] = 0.0;
								for(k=0; k<nosc; k++) cDip1Q[i][j] += U1Q[nj+k]*Dip1QMat[i][fr][k];
								cval += Dip1QMat[i][xfr][j]*cDip1Q[i][j];
							}
						}
						CorrFunc[n] += cval;
						// cjfeng 04/06/2016
						// U1Qs is transposed

						// And multiply by U1QMat[xfr] to extend propagation by one frame. The result goes in U1Q.
						// mmult_comp(U1QMat[xfr], U1Qs, U1Q, nosc, nthreads);

						trans_comp(U1Q,U1Qs,nosc,nosc,nthreads);
						if(nosc>=NBCOMP) {
							mmult_comp_block(U1QMat[xfr], U1Qs, U1Q, nosc, nthreads);
						}
						else {	
							mmult_comp_trans(U1QMat[xfr], U1Qs, U1Q, nosc, nthreads);
						}
					}
				}
			}
			// 2D part
			if (do2d) {
				fr = -1;
				frame = readframe-nread;
				// Now do calculations for selected frames. 
				while( (fr!=-2) && (!error) ) {
					// Step through all recently read frames to see whether frame
					// is suitable as the end point of a calculation. 
					// The time window twin is either window or win2d, depending 
					// on whether do2d = 0 or 1. 
					int twin = nbuffer-nread+1;
					while(frame<readframe) {
						// We check if at least window or win2d frames have been read 
						// and if frame is a multiple of skip. 
						if( ((frame%skip)==0) && (frame>=twin-1) ) {
							// frame will be the endpoint of the calculation. 
							// It will start at frame fr = frame-twin+1.
							fr = (frame-twin+1)%nbuffer;
							break;
						} else frame++;
					}
					// If we got all the way through the fr loop without finding anything
					// set fr = -2 to break out of the outside loop.
					if(frame==readframe) fr = -2;
					if( (fr>=0) && (nise) ) {
						// Do a dynamic calculation using fr as our starting point and augment frame.
						int n,k;
						frame++;  // We don't reference frame further in the calculation. 
					
						// For 2D spectra, we need 1Q propagators from 
						// 	t0-->t1
						// 	t0-->t1+t2
						// 	t0-->t1+t2+t3
						// 	t1-->t1+t2
						// 	t1-->t1+t2+t3
						// 	t1+t2-->t1+t2+t3
						// The 2Q propagator will be needed only from 
						// t1+t2 to t1+t2+t3.
					
						// xfr1, xfr2, and xfr3 will point to the location 
						// of the tau1, tau2, and tau3 frames, respectively. 
						int tau, tau1, tau2, tau3;
						int xfr0, xfr1, xfr2, xfr3;
						tau2 = T2/tstep;
						xfr0 = fr%nbuffer;

						// Initialize first frame of psi_a[i] array to be simply Dip1Qmat[i][0]
						// cjfeng 06/27/2016
						// Start using xfrs.
						for(i=0; i<3; i++) for(j=0; j<nosc; j++) psi_a[i][0][j] = Dip1QMat[i][xfr0][j];
						for(i=0; i<3; i++) {
							for(tau=0; tau<win2d-1; tau++) {
								// Fill in psi_a array. 
								// psi_a[i][tau+1] = U1QMat[(fr+tau)%nbuffer]*psi_a[i][tau];
								mvmult_comp(U1QMat[(fr+tau)%nbuffer], psi_a[i][tau], psi_a[i][tau+1], nosc, nthreads);
							}
						}

						// cjfeng 04/05/2016
						// Change the order of Dip2QMat to be a nbuffer*3*n2Q*nosc array.
						for(tau1=0; tau1<window; tau1++) {

							// cjfeng 06/27/2016
							// Start using xfrs
							xfr1 = (fr+tau1)%nbuffer;
							xfr2 = (fr+tau1+tau2)%nbuffer;
							// cjfeng 04/05/2016
							// Initialize psi_b[i][tau1].
							for(i=0; i<3; i++) for(j=0; j<nosc; j++) psi_b1[i][tau1][j] = Dip1QMat[i][xfr1][j];

							// Fill in psi_b1 array
							// psi_b1[i][tau+1] = U1QMat[(fr+tau)%nbuffer]*psi_b1[i][tau];
							for(i=0; i<3; i++) {
								for(tau=tau1; tau<tau1+tau2+window-1; tau++) {
									mvmult_comp(U1QMat[(fr+tau)%nbuffer], psi_b1[i][tau], psi_b1[i][tau+1], nosc, nthreads);
								}
							}

							// Initialize psi_b12[i][tau1+tau2].
							for(i=0; i<3; i++) for(j=0; j<nosc; j++) psi_b12[i][tau1+tau2][j] = Dip1QMat[i][xfr2][j];

							// Fill in psi_b12 array
							for(i=0; i<3; i++) {
								for(tau=tau1+tau2; tau<tau1+tau2+window-1; tau++) {
									mvmult_comp(U1QMat[(fr+tau)%nbuffer], psi_b12[i][tau], psi_b12[i][tau+1], nosc, nthreads);
								}
							}

							// Initialize the 2Q wavefunctions. 
							// psi_ca[i*3+j] = Dip2Q[j][tau1+tau2]*psi_a[i][tau1+tau2]
							// psi_cb[i*3+j] = Dip2Q[j][tau1+tau2]*psi_b1[i][tau1+tau2]
							
							for(i=0; i<3; i++) {
								for(j=0; j<3; j++) {
									mvmult_comp_trans_x(Dip2QMat[tau1+tau2][j], psi_a[i][tau1+tau2], psi_ca[i*3+j], nosc, n2Q, nthreads);
									mvmult_comp_trans_x(Dip2QMat[tau1+tau2][j], psi_b1[i][tau1+tau2], psi_cb[i*3+j], nosc, n2Q, nthreads);
								}
							}

							for(tau3=0; tau3<window; tau3++) {

								// cjfeng 06/27/2016
								// Start using xfr3.
								xfr3 = (fr+tau1+tau2+tau3)%nbuffer;

								GNREAL popfac1Q = popdecay1Q[tau1+2*tau2+tau3];
								GNREAL popfac2Q = popfac1Q*popdecay2Q[tau3];
								GNCOMP cval1 = 0.0;
								GNCOMP cval2 = 0.0;
								int P,p,l;
								// Computing nonlinear response function
								int tau1window = tau1*window;
								for(P=0; P<npol; P++) {
									int a, b;
									// cjfeng 04/27/2016
									// Adding local variable to reduce memory access
									GNCOMP REPH_private;
									GNCOMP NREPH_private;
									REPH_private = 0.0;
									NREPH_private = 0.0;
									// p sums over the forms iiii, iijj, ijji, and ijij
									for(p=0; p<4; p++) {
										for(a=0; a<3; a++) {
											for(b=0; b<3; b++) {
												if(p==0) { i=a; j=a; k=a; l=a; }
												else if(p==1) { i=a; j=a; k=b; l=b; }
												else if(p==2) { i=a; j=b; k=b; l=a; }
												else { i=a; j=b; k=a; l=b; };		// p==3
												// For p==0, we only add the signal when a==b. 
												if( (p!=0) || (a==b) ) {
													// cjfeng 07/07/2016
													// Increase locality of M_ijkl_IJKL
													GNREAL orient_fac = M_ijkl_IJKL[P][p];
													// First rephasing pathway:
													//
													// 	| b 0 |
													// 	| 0 0 |
													// 	| 0 a |
													// 	| 0 0 |
													//
													// The contribution to the rephasing spectrum is the orientionally averaged value
													// 	( Dip1Q[t3+t2+t1]*psi_b12[tau1+tau2+tau3] ) * conj( Dip1Q[tau1] * psi_a[tau1] )

													// cjfeng 07/12/2016
													// Include switches on simulating pathways
													// if (reph1) {
														cval1 = 0.0; cval2 = 0.0;
													
														for(n=0; n<nosc; n++) cval1 += Dip1QMat[l][xfr3][n]*psi_b12[k][tau1+tau2+tau3][n];
														for(n=0; n<nosc; n++) cval2 += Dip1QMat[j][xfr1][n]*psi_a[i][tau1][n];
														REPH_private += orient_fac*cval1*conj(cval2)*popfac1Q;
														// REPH[POL[P]][tau1*window+tau3] += M_ijkl_IJKL[P][p]*cval1*conj(cval2)*popfac1Q;
													// }

													// Second rephasing pathway:
													//
													// 	| b 0 |
													// 	| b a |
													// 	| 0 a |
													// 	| 0 0 |
													//
													// The contribution to the rephasing spectrum is the orientionally averaged value
													// 	( Dip1Q[tau3+tau2+tau1]*psi_b1[tau1+tau2+tau3] ) * conj( Dip1Q[tau2+tau1]*psi_a[tau2+tau1] )
													// cjfeng 07/12/2016
													// Include switches on simulating pathways
													// if (reph2) {
														cval1 = 0.0; cval2 = 0.0;
														for(n=0; n<nosc; n++) cval1 += Dip1QMat[l][xfr3][n]*psi_b1[j][tau1+tau2+tau3][n];
														for(n=0; n<nosc; n++) cval2 += Dip1QMat[k][xfr2][n]*psi_a[i][tau1+tau2][n];
														REPH_private += orient_fac*cval1*conj(cval2)*popfac1Q;
														// REPH[POL[P]][tau1*window+tau3] += M_ijkl_IJKL[P][p]*cval1*conj(cval2)*popfac1Q;
													// }

													// Final rephasing pathway:
													//
													// 	| c a |
													// 	| b a |
													// 	| 0 a |
													// 	| 0 0 |
													//
													// The contribution to the rephasing spectrum is the orientionally averaged value
													// 	( psi_cb[j*3+k] ) * ( Dip2Q[tau3+tau2+tau1]*psi_a[tau3+tau2+tau1] )^\dagger
													// cjfeng 07/12/2016
													// Include switches on simulating pathways
													// if (reph3) {
													
													// cjfeng 04/05/2016
													// Changed the order of Dip2QMat	
													// to be a nbuffer*3*n2Q*nosc array
														mvmult_comp_trans_x(Dip2QMat[xfr3][l], psi_a[i][tau1+tau2+tau3], psi2Q, nosc, n2Q, nthreads);
														cval1 = 0.0; 
														for(n=0; n<n2Q; n++) cval1 += psi_cb[j*3+k][n]*conj(psi2Q[n]);
														REPH_private -= orient_fac*cval1*popfac2Q;
														// REPH[POL[P]][tau1*window+tau3] -= M_ijkl_IJKL[P][p]*cval1*popfac2Q;
													// }

													// First non-rephasing pathway: 
													// 	
													// 	| b 0 |
													// 	| 0 0 |
													// 	| a 0 |
													// 	| 0 0 |
													//
													// The contribution to the non-rephasing spectrum is the orientationally averaged value
													// 	( Dip1Q[tau1+tau2+tau3]*psi_b12[tau1+tau2+tau3] ) * ( Dip1Q[tau1]*psi_a[tau1] ) 
													// cjfeng 07/12/2016
													// Include switches on simulating pathways
													// if (nreph1) {
														cval1 = 0.0; cval2 = 0.0;
														for(n=0; n<nosc; n++) cval1 += Dip1QMat[l][xfr3][n]*psi_b12[k][tau1+tau2+tau3][n];
														for(n=0; n<nosc; n++) cval2 += Dip1QMat[j][xfr1][n]*psi_a[i][tau1][n];
														NREPH_private += orient_fac*cval1*cval2*popfac1Q;
														// NREPH[POL[P]][tau1*window+tau3] += M_ijkl_IJKL[P][p]*cval1*cval2*popfac1Q;
													// } 

													// Second non-rephasing pathway:
													// 	| a 0 |
													// 	| a b |
													// 	| a 0 |
													// 	| 0 0 |
													//
													// The contribution to the non-rephasing spectrum is the orientionally averaged value
													// 	( Dip1Q[tau1+tau2+tau3]*psi_a[tau1+tau2+tau3] ) * conj( Dip1Q[tau1+tau2]*psi_b1[tau1+tau2] )
													// cjfeng 07/12/2016
													// Include switches on simulating pathways
													// if (nreph2) {
														cval1 = 0.0; cval2 = 0.0;
														for(n=0; n<nosc; n++) cval1 += Dip1QMat[l][xfr3][n]*psi_a[i][tau1+tau2+tau3][n];
														for(n=0; n<nosc; n++) cval2 += Dip1QMat[k][xfr2][n]*psi_b1[j][tau1+tau2][n];
														NREPH_private += orient_fac*cval1*conj(cval2)*popfac1Q;
														// NREPH[POL[P]][tau1*window+tau3] += M_ijkl_IJKL[P][p]*cval1*conj(cval2)*popfac1Q;
													// }

													// Final non-rephasing pathway:
													// 	| c b |
													// 	| a b |
													// 	| a 0 |
													// 	| 0 0 |
													//
													// The contribution to the non-rephasing spectrum is the orientationally averaged value
													// 	( psi_ca[k*3+i] ) * ( Dip2Q[tau1+tau2+tau3]*psi_b1[tau1+tau2+tau3] )
												
													// cjfeng 07/12/2016
													// Include switches on simulating pathways
													// if (nreph3) {
													// cjfeng 04/05/2016
													// Swapped the order of Dip2QMat.
														mvmult_comp_trans_x(Dip2QMat[xfr3][l], psi_b1[j][tau1+tau2+tau3], psi2Q, nosc, n2Q, nthreads);

														cval1 = 0.0;
													// Swapping i and k
													// cjfeng 06/08/2016
														for(n=0; n<n2Q; n++) cval1 += psi_ca[i*3+k][n]*conj(psi2Q[n]);
														NREPH_private -= orient_fac*cval1*popfac2Q;
														// NREPH[POL[P]][tau1*window+tau3] -= M_ijkl_IJKL[P][p]*cval1*popfac2Q;
													// }
												}
											}
										}
									}
									// cjfeng 04/27/2016
									REPH[POL[P]][tau1window+tau3] += REPH_private;
									// REPH[POL[P]][tau1*window+tau3] += REPH_private;
									NREPH[POL[P]][tau1window+tau3] += NREPH_private;
									// NREPH[POL[P]][tau1*window+tau3] += NREPH_private;
								}
								// cjfeng 07/13/2016
								// Propagate 2Q wavefunctions after computing correlation function. 
								// psi_ca[i*3+j] = U2Q[tau3]*psi_ca[i*3+j]
								// psi_cb[i*3+j] = U2Q[tau3]*psi_cb[i*3+j]
								if(pert) {
									gen_pert_prop( U1QMat[xfr3], U2Qs, nosc, n2Q, expfac, delta);
									for(i=0; i<9; i++) {
										for(j=0; j<n2Q; j++) psi2Q[j] = psi_ca[i][j];
										mvmult_comp(U2Qs, psi2Q, psi_ca[i], n2Q, nthreads);
										for(j=0; j<n2Q; j++) psi2Q[j] = psi_cb[i][j];
										mvmult_comp(U2Qs, psi2Q, psi_cb[i], n2Q, nthreads);
									}
								} else {
									for(i=0; i<9; i++) {
										for(j=0; j<n2Q; j++) psi2Q[j] = psi_ca[i][j];
										mvmult_comp(U2QMat[xfr3], psi2Q, psi_ca[i], n2Q, nthreads);
										for(j=0; j<n2Q; j++) psi2Q[j] = psi_cb[i][j];
										mvmult_comp(U2QMat[xfr3], psi2Q, psi_cb[i], n2Q, nthreads);
									}	
								}
							}
						}
					}
				}
			}
			 
		}

/*****************************************************
 * Static averaging and time-averaging approximation *
 *****************************************************/ 
		else if ((!error) && !nise) {
			fr = -1;
			frame = readframe-nread;
			// Now do calculations for selected frames. 
			while( (fr!=-2) && (!error) ) {
				// Step through all recently read frames to see whether frame
				// is suitable as the end point of a calculation. 
				// The time window twin is either window or win2d, depending 
				// on whether do2d = 0 or 1. 
				int twin = nbuffer-nread+1;
				while(frame<readframe) {
					if((frame%100)==0) {
						printf("Frame: %d\n", frame);
						// cjfeng 06/27/2016
						// Tracking memory usage.
						getrusage(RUSAGE_SELF, &r_usage);
						printf("Memory usage = %ld kB.\n", r_usage.ru_maxrss);
						fprintf(lfp, "Frame: %d\n", frame);
						fflush(lfp);
					}
					// We check if at least window or win2d frames have been read 
					// and if frame is a multiple of skip. 
					if( ((frame%skip)==0) && (frame>=twin-1) ) {
						// frame will be the endpoint of the calculation. 
						// It will start at frame fr = frame-twin+1.
						fr = (frame-twin+1)%nbuffer;
						break;
					} else frame++;
				}
				// If we got all the way through the fr loop without finding anything
				// set fr = -2 to break out of the outside loop.
				if(frame==readframe) fr = -2;
				else if( (fr>=0) ) {
					int k;
					int noscsq = nosc*nosc;
					int n2Qsq = n2Q*n2Q;
					int noscn2Q = nosc*n2Q;
					
					// Generate averaged Hamiltonian for a window starting at frame fr 
					// and augment frame. The average Hamiltonian is stored in Ham1QAr. 
					// If Ham1QAr gets filled up, we stop and do a spectral calculation 
					// for all stored frames before proceeding. 
					frame++;  // We don't reference frame further in the calculation. 
					// Generate averaged Hamiltonian. 
					for(j=0; j<noscsq; j++) Ham1QAr[nAvgd][j] = 0.0;
					for(i=0; i<3; i++) for(j=0; j<nosc; j++) Dip1QAr[i][nAvgd][j] = 0.0;
					if( (do2d) && ((!pert) || pertvec) ) {
						for(i=0; i<3; i++) for(j=0; j<noscn2Q; j++) Dip2QAr[i][nAvgd][j] = 0.0;
					}
					if(whann) {
						for(i=0; i<window; i++) for(j=0; j<noscsq; j++) Ham1QAr[nAvgd][j] += Ham1QMat[(fr+i)%nbuffer][j]*hann[i];
						for(i=0; i<3; i++) for(j=0; j<nosc; j++) Dip1QAr[i][nAvgd][j] = Dip1QMat[i][(fr+((int) (window/2)))%nbuffer][j];
						// cjfeng 04/05/2016
						// Dip2QAr is not transposed compared to the original form.
						// However, Dip2QMat is now a nbuffer*3*n2Q*nosc array.
						// It will be slower than before due to swapping the order of Dip2QMat.
						if( (do2d) && ((!pert) || pertvec) ) for(i=0; i<3; i++) for(j=0; j<nosc; j++) for(k=0; k<n2Q; k++) Dip2QAr[i][nAvgd][j*n2Q+k] = Dip2QMat[(fr+((int) (window/2)))%nbuffer][i][k*nosc+j];
					} else {
						for(i=0; i<window; i++) for(j=0; j<noscsq; j++) Ham1QAr[nAvgd][j] += Ham1QMat[(fr+i)%nbuffer][j]/window;
						for(i=0; i<3; i++) for(j=0; j<nosc; j++) Dip1QAr[i][nAvgd][j] = Dip1QMat[i][(fr+((int) (window/2)))%nbuffer][j];
						// cjfeng 04/05/2016
						// Dip2QAr is not transposed compared to the original form.
						if( (do2d) && ( (!pert) || pertvec )) for(i=0; i<3; i++) for(j=0; j<nosc; j++) for(k=0; k<n2Q; k++) Dip2QAr[i][nAvgd][j*n2Q+k] = Dip2QMat[(fr+((int) (window/2)))%nbuffer][i][k*nosc+j];
					}
					nAvgd++;

					// Declaring local index array for eigensolver
					int *isuppz1Q = malloc(2*nosc*sizeof(int));
					int *isuppz2Q = malloc(2*n2Q*sizeof(int));

					// If we've filled up the Ham1QAr array, do a calculation and re-start. 
					if(nAvgd==nthreads) {
						int thr;
						lapack_int info;
						for(thr=0; thr<nthreads; thr++) {
							int tid = omp_get_thread_num();
							// If needed, generate 2Q Hamiltonian. This MUST be done before eigenvalue calculation. 
							if( (do2d) && (!pert) ) gen_ham_2Q(Ham1QAr[tid], nosc, Ham2QAr[tid], n2Q, delta);
							// Find one-quantum eigenvalues
							// Note that Ham1QAr[tid] now contains eigenvectors
						
							// cjfeng 06/30/2016
							// RRR algorithm
							info = SYEVR( LAPACK_ROW_MAJOR, 'V', 'A', 'U', nosc, Ham1QAr[tid], nosc, wstart, wstop, 0, nosc-1, 0.00001, &nosc, Evals1QAr[tid], Ham1QAr[tid], nosc, isuppz1Q);
			
							// Calculate one-quantum dipole moments
							for(i=0; i<3; i++) mvmult_real_serial_trans(Ham1QAr[tid], Dip1QAr[i][tid], ExDip1QAr[i][tid], nosc, nthreads);
							// Calculate FTIR spectrum
							int n,N,d,ndx,k,K;
							GNREAL osc; 
							
							for(n=0; n<nosc; n++) {

								osc = 0.0;
								for(d=0; d<3; d++) osc += ExDip1QAr[d][tid][n]*ExDip1QAr[d][tid][n];
								ndx = floor( (Evals1QAr[tid][n]-wstart)/wres + 0.5 );
								if(ndx<npts && ndx>=0) ftir[ndx] += osc;
							}
							if( (do2d) && (!pert) ) {
								// Two-quantum eigenvalues
								// Note that Ham2QAr[tid] now contains eigenvectors
							
								// cjfeng 06/30/2016
								// RRR algorithm
								info = SYEVR( LAPACK_ROW_MAJOR, 'V', 'A', 'U', n2Q, Ham2QAr[tid], n2Q, wstart, wstop, 0, n2Q-1, 0.00001, &n2Q, Evals2QAr[tid], Ham2QAr[tid], n2Q, isuppz2Q);

								for(i=0; i<3; i++) {
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) ExDip2QAr[i][tid][n*n2Q+N] = 0.0;
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) for(k=0; k<nosc; k++) ExDip2QAr[i][tid][n*n2Q+N] += Ham1QAr[tid][k*nosc+n]*Dip2QAr[i][tid][k*n2Q+N];
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) Dip2QAr[i][tid][n*n2Q+N] = ExDip2QAr[i][tid][n*n2Q+N];
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) ExDip2QAr[i][tid][n*n2Q+N] = 0.0;
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) for(K=0; K<n2Q; K++) ExDip2QAr[i][tid][n*n2Q+N] += Dip2QAr[i][tid][n*n2Q+K]*Ham2QAr[tid][K*n2Q+N];
								}
								// Now calculate the 2D spectrum. 
								calc_2dir(Evals1QAr, Evals2QAr, ExDip1QAr, ExDip2QAr, tid, nosc, n2Q, npts, wres, wstart, wstop, REPH, NREPH, POL, npol, reph, nreph);
							} else if( (do2d) && (!pertvec) ) {
								// Generate first-order anharmonically perturbed 2Q energies.
								// Remember that Ham1QAr[tid] now contains eigenvectors. 
								gen_perturb_2Q_energies(Ham1QAr[tid], Evals1QAr[tid], Evals2QAr[tid], nosc, delta);
								// And calculate the 2D spectrum. 
								calc_2dir_pert(Evals1QAr, Evals2QAr, ExDip1QAr, tid, nosc, n2Q, npts, wres, wstart, wstop, REPH, NREPH, POL, npol, reph, nreph);
							} else if( (do2d) && (pertvec) ) {
								// Generate first-order anharmonically perturbed 2Q energies.
								// Remember that Ham1QAr[tid] now contains eigenvectors. 
								gen_perturb_2Q_energies(Ham1QAr[tid], Evals1QAr[tid], Evals2QAr[tid], nosc, delta);
								// Generate Eigenvectors
								//gen_perturb_2Q_matrix(Ham1QAr[tid], Evals1QAr[tid], Ham2QAr[tid], Evals2QAr[tid], Tau2QAr[tid], nosc, n2Q, delta);
								gen_perturb_2Q_vectors(Ham1QAr[tid], Evals1QAr[tid], Ham2QAr[tid], Evals2QAr[tid], nosc, n2Q, delta);
								// gen_perturb_2Q_vectors(Ham1QAr[tid], Evals1QAr[tid], Ham2QAr[tid], Evals2QAr[tid], Tau2QAr[tid], nosc, n2Q, delta);
								// And calculate the 2D spectrum. 
								for(i=0; i<3; i++) {
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) ExDip2QAr[i][tid][n*n2Q+N] = 0.0;
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) for(k=0; k<nosc; k++) ExDip2QAr[i][tid][n*n2Q+N] += Ham1QAr[tid][k*nosc+n]*Dip2QAr[i][tid][k*n2Q+N];
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) Dip2QAr[i][tid][n*n2Q+N] = ExDip2QAr[i][tid][n*n2Q+N];
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) ExDip2QAr[i][tid][n*n2Q+N] = 0.0;
									for(n=0; n<nosc; n++) for(N=0; N<n2Q; N++) for(K=0; K<n2Q; K++) ExDip2QAr[i][tid][n*n2Q+N] += Dip2QAr[i][tid][n*n2Q+K]*Ham2QAr[tid][K*n2Q+N];
								}
								// Now calculate the 2D spectrum. 
								calc_2dir(Evals1QAr, Evals2QAr, ExDip1QAr, ExDip2QAr, tid, nosc, n2Q, npts, wres, wstart, wstop, REPH, NREPH, POL, npol, reph, nreph);
							}
						}
						nAvgd = 0;
					}
					// Release memory
					free(isuppz1Q);
					free(isuppz2Q);
				}


			}

		}

		fr = -1;
		frame = readframe-nread;
		// Now do calculations for selected frames. 
		while( (fr!=-2) && (!error) ) {
			// Step through all recently read frames to see whether 
			// we should dump the spectrum. 
			// The time window twin is either window or win2d, depending 
			// on whether do2d = 0 or 1. 
			int twin = nbuffer-nread+1;
			while(frame<readframe) {
				// We check if at least window or win2d frames have been read 
				// and if frame is a multiple of dump. 
				if( ((frame%dump)==0) && (frame>=twin-1) ) {
					// frame will be the endpoint of the calculation. 
					// It will start at frame fr = frame-twin+1.
					fr = (frame-twin+1)%nbuffer;
					break;
				} else frame++;
			}
			// If we got all the way through the fr loop without finding anything
			// set fr = -2 to break out of the outside loop.
			if(frame==readframe) fr = -2;
			if( (fr>=0) && (frame!=0) ) printf("Dumping at frame %d.\n", frame);
			if( (fr>=0) && (frame!=0) ) fprintf(Trajfp[0], "%6.10f\n", (((GNREAL) tstep)*((GNREAL) frame))/1000.0);
			// Dump calculated spectra in trajectory files
			if( (fr>=0) && (nise) && (frame!=0)) {
				frame++; 
				
				// First do FTIR
				for(i=0; i<window; i++) NetCorrFunc[i] += CorrFunc[i];
				if(whann==1) for(i=0; i<window; i++) CorrFunc[i] = hann[i]*popdecay1Q[i]*CorrFunc[i];
				else for(i=0; i<window; i++) CorrFunc[i] = popdecay1Q[i]*CorrFunc[i];
				for(i=0; i<winzpad; i++) FTin1D[i][0] = 0.0;
				for(i=0; i<winzpad; i++) FTin1D[i][1] = 0.0;
				for(i=0; i<window; i++) FTin1D[i][0] = creal(CorrFunc[i]);
				for(i=0; i<window; i++) FTin1D[i][1] = cimag(CorrFunc[i]);
				// FTIR spectrum now stored in FTout1D[.][0]
				fftw_execute(FTplan1D);
				for(i=0; i<nprint; i++) { 
					if(i==nprint-1) fprintf(Trajfp[1], "%6.10f\n", FTout1D[(ndxstart+i)%winzpad][0]);
					else fprintf(Trajfp[1], "%6.10f, ", FTout1D[(ndxstart+i)%winzpad][0]);
				}
				for(i=0; i<window; i++) CorrFunc[i] = 0.0;

				// Now 2DIR
				if(do2d) {
					int winzpad2=winzpad*winzpad;
					GNCOMP cval = 0.0;
					int p;
					// FT and print rephasing spectra
					for(p=0; p<npol; p++) {
						for(i=0; i<window; i++) for(j=0; j<window; j++) NetREPH[POL[p]][i*window+j] += REPH[POL[p]][i*window+j];
						// cjfeng 06/27/2016
						// Reduce loop overhead.
						for(i=0; i<winzpad2; i++) {
							FTin2D[i][0] = 0.0;
							FTin2D[i][1] = 0.0;
						}
						// for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = 0.0;
						// for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = 0.0;
						// Apply Hann window, if requested. 
						if(whann) {
							for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][0] = creal(REPH[POL[p]][i*window+j])*hann[i];
							for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][1] = cimag(REPH[POL[p]][i*window+j])*hann[i];
						} else {
							for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][0] = creal(REPH[POL[p]][i*window+j]);
							for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][1] = cimag(REPH[POL[p]][i*window+j]);
						}
						fftw_execute(FTplan2D);
			
						// For rephasing only: flip spectrum along w1
						for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = FTout2D[(winzpad-1-i)*winzpad+j][0];
						for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = FTout2D[(winzpad-1-i)*winzpad+j][1];
						// cjfeng 06/27/2016
						// Reduce loop overhead.
						//
						for(i=0; i<winzpad2; i++) {
							FTout2D[i][0] = FTin2D[i][0];
							FTout2D[i][1] = FTin2D[i][1];
						}
						// for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTout2D[i*winzpad+j][0] = FTin2D[i*winzpad+j][0];
						// for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTout2D[i*winzpad+j][1] = FTin2D[i*winzpad+j][1];
						int ndx1, ndx3;
						for(i=0; i<nprint; i++) {
							ndx1 = (ndxstart+i)%winzpad;
							for(j=0; j<nprint; j++) {
								ndx3 = (ndxstart+j)%winzpad;
								cval = FTout2D[ndx1*winzpad+ndx3][0] + FTout2D[ndx1*winzpad+ndx3][1]*I;
								if(cimag(cval)<0) fprintf(Trajfp[2+POL[p]], "%6.10f%6.10fi\t", creal(cval), cimag(cval));
								else fprintf(Trajfp[2+POL[p]], "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
							}
						}
						fprintf(Trajfp[2+POL[p]], "\n");
						for(i=0; i<window; i++) for(j=0; j<window; j++) REPH[POL[p]][i*window+j] = 0.0;
					}
		
					// FT and print non-rephasing spectra
					for(p=0; p<npol; p++) {
						for(i=0; i<window; i++) for(j=0; j<window; j++) NetNREPH[POL[p]][i*window+j] += NREPH[POL[p]][i*window+j];
						// cjfeng 06/27/2016
						// Reduce loop overhead.
						for(i=0; i<winzpad2; i++) {
							FTin2D[i][0] = 0.0;
							FTin2D[i][1] = 0.0;
						}
						// for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = 0.0;
						// for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = 0.0;
						// Apply Hann window, if requested
						if(whann) {
							for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][0] = creal(NREPH[POL[p]][i*window+j])*hann[i];
							for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][1] = cimag(NREPH[POL[p]][i*window+j])*hann[i];
						} else {
							for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][0] = creal(NREPH[POL[p]][i*window+j]);
							for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][1] = cimag(NREPH[POL[p]][i*window+j]);
						}
						fftw_execute(FTplan2D);
						int ndx1, ndx3;
						for(i=0; i<nprint; i++) {
							ndx1 = (ndxstart+i)%winzpad;
							for(j=0; j<nprint; j++) {
								ndx3 = (ndxstart+j)%winzpad;
								cval = FTout2D[ndx1*winzpad+ndx3][0] + FTout2D[ndx1*winzpad+ndx3][1]*I;
								if(cimag(cval)<0) fprintf(Trajfp[6+POL[p]], "%6.10f%6.10fi\t", creal(cval), cimag(cval));
								else fprintf(Trajfp[6+POL[p]], "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
							}
						}
						fprintf(Trajfp[6+POL[p]], "\n");
						for(i=0; i<window; i++) for(j=0; j<window; j++) NREPH[POL[p]][i*window+j] = 0.0;
					}
				}
		
			} else if( (fr>=0) && (!nise) && (frame!=0) ) {
				frame++; 

				// First do FTIR
				for(i=0; i<npts; i++) netftir[i] += ftir[i];
				for(i=0; i<winzpad; i++) FTin1D[i][0] = 0.0;
				for(i=0; i<winzpad; i++) FTin1D[i][1] = 0.0;
				for(i=0; i<npts; i++) FTin1D[i][0] = ftir[i];
				// cjfeng 07/05/2016
				// Should it be FORWARD or -iwt by another plan?
				fftw_execute(FTplan1D);
				for(i=0; i<winzpad; i++) FTin1D[i][0] = FTout1D[i][0]*exp(-i/(2*TauP*c*wres*winzpad));
				for(i=0; i<winzpad; i++) FTin1D[i][1] = -FTout1D[i][1]*exp(-i/(2*TauP*c*wres*winzpad));
				fftw_execute(FTplan1D);
				// cjfeng 07/05/2016
				// Normalization?
				for(i=0; i<winzpad; i++) FTout1D[i][1] = -FTout1D[i][1];
				for(i=0; i<npts; i++) ftir[i] = FTout1D[i][0];
				for(i=0; i<npts; i++) {
					if(i==npts-1) fprintf(Trajfp[1], "%6.10f\n", ftir[i]);
					else fprintf(Trajfp[1], "%6.10f, ", ftir[i]);
				}
				for(i=0; i<npts; i++) ftir[i] = 0.0;

				// Now 2DIR
				if(do2d) {
					int p;
					for(p=0; p<npol; p++) {
						for(i=0; i<npts; i++) for(j=0; j<npts; j++) NetREPH[POL[p]][i*npts+j] += creal(REPH[POL[p]][i*npts+j]);
						for(i=0; i<npts; i++) for(j=0; j<npts; j++) NetNREPH[POL[p]][i*npts+j] += creal(NREPH[POL[p]][i*npts+j]);

						if(TauP>0) {
							// Dress rephasing spectrum with Lorentzian profile. Note that it is so far purely real. 
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = 0.0;
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = 0.0;
							for(i=0; i<npts; i++) for(j=0; j<npts; j++) FTin2D[i*winzpad+j][0] = creal(REPH[POL[p]][i*npts+j]);
							// cjfeng 07/05/2016
							// Should it be also FORWARD, or -iwt?
							fftw_execute(FTplan2D);
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = FTout2D[i*winzpad+j][0]*exp(-((winzpad-i)+j)/(2*TauP*c*wres*winzpad));
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = -FTout2D[i*winzpad+j][1]*exp(-((winzpad-i)+j)/(2*TauP*c*wres*winzpad));
							// cjfeng 07/05/2016
							// Normalization?
							fftw_execute(FTplan2D);
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTout2D[i*winzpad+j][1] = -FTout2D[i*winzpad+j][1];
							for(i=0; i<npts; i++) for(j=0; j<npts; j++) REPH[POL[p]][i*npts+j] = FTout2D[i*winzpad+j][0] + I*FTout2D[i*winzpad+j][1];
			
							// Dress non-rephasing spectrum with Lorentzian profile. Note that it is so far purely real. 
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = 0.0;
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = 0.0;
							for(i=0; i<npts; i++) for(j=0; j<npts; j++) FTin2D[i*winzpad+j][0] = creal(NREPH[POL[p]][i*npts+j]);
							// cjfeng 07/05/2016
							// Should it be also FORWARD, or -iwt?
							fftw_execute(FTplan2D);
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = FTout2D[i*winzpad+j][0]*exp(-(i+j)/(2*TauP*c*wres*winzpad));
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = -FTout2D[i*winzpad+j][1]*exp(-(i+j)/(2*TauP*c*wres*winzpad));
							// cjfeng 07/05/2016
							// Normalization?
							fftw_execute(FTplan2D);
							for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTout2D[i*winzpad+j][1] = -FTout2D[i*winzpad+j][1];
							for(i=0; i<npts; i++) for(j=0; j<npts; j++) NREPH[POL[p]][i*npts+j] = FTout2D[i*winzpad+j][0] + I*FTout2D[i*winzpad+j][1];	
						}
						GNCOMP cval;
						for(j=0; j<npts; j++) {
							for(k=0; k<npts; k++) {
								cval = REPH[POL[p]][j*npts+k];
								if(cimag(cval)<0) fprintf(Trajfp[2+POL[p]], "%6.10f%6.10fi\t", creal(cval), cimag(cval));
								else fprintf(Trajfp[2+POL[p]], "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
								// Changed print syntax to reflect complex data -- MER 7/21/2015
								cval = NREPH[POL[p]][j*npts+k];
								if(cimag(cval)<0) fprintf(Trajfp[6+POL[p]], "%6.10f%6.10fi\t", creal(cval), cimag(cval));
								else fprintf(Trajfp[6+POL[p]], "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
								// Changed print syntax to reflect complex data -- MER 7/21/2015
							}
						}
						fprintf(Trajfp[2+POL[p]], "\n");
						fprintf(Trajfp[6+POL[p]], "\n");
						for(i=0; i<npts; i++) for(j=0; j<npts; j++) REPH[POL[p]][i*npts+j] = 0.0;
						for(i=0; i<npts; i++) for(j=0; j<npts; j++) NREPH[POL[p]][i*npts+j] = 0.0;
					}
				}
			} else if(frame==0) frame++;
		}

	}
		
	if( (!error) && (nise) ) {
		for(i=0; i<window; i++) CorrFunc[i] += NetCorrFunc[i];
		if(whann==1) for(i=0; i<window; i++) CorrFunc[i] = hann[i]*popdecay1Q[i]*CorrFunc[i];
		else for(i=0; i<window; i++) CorrFunc[i] = popdecay1Q[i]*CorrFunc[i];
		for(i=0; i<winzpad; i++) FTin1D[i][0] = 0.0;
		for(i=0; i<winzpad; i++) FTin1D[i][1] = 0.0;
		for(i=0; i<window; i++) FTin1D[i][0] = creal(CorrFunc[i]);
		for(i=0; i<window; i++) FTin1D[i][1] = cimag(CorrFunc[i]);

		// FTIR spectrum now stored in FTout1D[.][0]
		fftw_execute(FTplan1D);
	
		// Print FTIR spectrum and frequency axis
		for(i=0; i<nprint; i++) {
			fprintf(afp, "%6.10f\n", (ndxstart+i)*dw);
			fprintf(ffp, "%6.10f\n", FTout1D[(ndxstart+i)%winzpad][0]);
		}

		if(do2d) {
			GNCOMP cval = 0.0;
			int p;
			// FT and print rephasing spectra
			for(p=0; p<npol; p++) {
				for(i=0; i<window; i++) for(j=0; j<window; j++) REPH[POL[p]][i*window+j] += NetREPH[POL[p]][i*window+j];
				for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = 0.0;
				for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = 0.0;
				// Apply Hann window, if requested. 
				if(whann) {
					for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][0] = creal(REPH[POL[p]][i*window+j])*hann[i];
					for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][1] = cimag(REPH[POL[p]][i*window+j])*hann[i];
				} else {
					for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][0] = creal(REPH[POL[p]][i*window+j]);
					for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][1] = cimag(REPH[POL[p]][i*window+j]);
				}
				// cjfeng 07/12/2016
				// Print REPH FID
				if (tresponse) {
					for(i=0; i<winzpad; i++) {
						for(j=0; j<winzpad; j++) {
							GNREAL re, im;
							re = FTin2D[i*winzpad+j][0];
							im = FTin2D[i*winzpad+j][1];
							fprintf(fid_reph[POL[p]], "%6.10f\t", re);
							// if(im<0) fprintf(fid_reph[POL[p]], "%6.10f%6.10fi\t", re, im);
							// else fprintf(fid_reph[POL[p]], "%6.10f+%6.10fi\t", re, im);
						}
						fprintf(fid_reph[POL[p]], "\n");
					}
				}

				fftw_execute(FTplan2D);
	
				// For rephasing only: flip spectrum along w1
				for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = FTout2D[(winzpad-1-i)*winzpad+j][0];
				for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = FTout2D[(winzpad-1-i)*winzpad+j][1];
				for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTout2D[i*winzpad+j][0] = FTin2D[i*winzpad+j][0];
				for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTout2D[i*winzpad+j][1] = FTin2D[i*winzpad+j][1];
				int ndx1, ndx3;
				for(i=0; i<nprint; i++) {
					ndx1 = (ndxstart+i)%winzpad;
					for(j=0; j<nprint; j++) {
						ndx3 = (ndxstart+j)%winzpad;
						cval = FTout2D[ndx1*winzpad+ndx3][0] + FTout2D[ndx1*winzpad+ndx3][1]*I;
						if(cimag(cval)<0) fprintf(rfp[POL[p]], "%6.10f%6.10fi\t", creal(cval), cimag(cval));
						else fprintf(rfp[POL[p]], "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
					}
					fprintf(rfp[POL[p]], "\n");
				}
			}

			// FT and print non-rephasing spectra
			for(p=0; p<npol; p++) {
				for(i=0; i<window; i++) for(j=0; j<window; j++) NREPH[POL[p]][i*window+j] += NetNREPH[POL[p]][i*window+j];
				for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = 0.0;
				for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = 0.0;
				// Apply Hann window, if requested
				if(whann) {
					for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][0] = creal(NREPH[POL[p]][i*window+j])*hann[i];
					for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][1] = cimag(NREPH[POL[p]][i*window+j])*hann[i];
				} else {
					for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][0] = creal(NREPH[POL[p]][i*window+j]);
					for(i=0; i<window; i++) for(j=0; j<window; j++) FTin2D[i*winzpad+j][1] = cimag(NREPH[POL[p]][i*window+j]);
				}	
				// cjfeng 07/12/2016
				// Print NREPH FID
				if (tresponse) {
					for(i=0; i<winzpad; i++) {
						for(j=0; j<winzpad; j++) {
							GNREAL re, im;
							re = FTin2D[i*winzpad+j][0];
							im = FTin2D[i*winzpad+j][1];
							fprintf(fid_nreph[POL[p]], "%6.10f\t", re);
							// if(im<0) fprintf(fid_nreph[POL[p]], "%6.10f%6.10fi\t", re, im);
							// else fprintf(fid_nreph[POL[p]], "%6.10f+%6.10fi\t", re, im);
						}
						fprintf(fid_nreph[POL[p]], "\n");
					}
				}

				fftw_execute(FTplan2D);
				int ndx1, ndx3;
				for(i=0; i<nprint; i++) {
					ndx1 = (ndxstart+i)%winzpad;
					for(j=0; j<nprint; j++) {
						ndx3 = (ndxstart+j)%winzpad;
						cval = FTout2D[ndx1*winzpad+ndx3][0] + FTout2D[ndx1*winzpad+ndx3][1]*I;
						if(cimag(cval)<0) fprintf(nrfp[POL[p]], "%6.10f%6.10fi\t", creal(cval), cimag(cval));
						else fprintf(nrfp[POL[p]], "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
					}
					fprintf(nrfp[POL[p]], "\n");
				}
			}
		}
	} else if( (!error) && (!nise) ) {
		int p;
		for(i=0; i<npts; i++) ftir[i] += netftir[i];
		for(i=0; i<winzpad; i++) FTin1D[i][0] = 0.0;
		for(i=0; i<winzpad; i++) FTin1D[i][1] = 0.0;
		for(i=0; i<npts; i++) FTin1D[i][0] = ftir[i];
		// cjfeng 07/05/2016
		// Should it be FORWARD or -iwt?
		fftw_execute(FTplan1D);
		for(i=0; i<winzpad; i++) FTin1D[i][0] = FTout1D[i][0]*exp(-i/(2*TauP*c*wres*winzpad));
		for(i=0; i<winzpad; i++) FTin1D[i][1] = -FTout1D[i][1]*exp(-i/(2*TauP*c*wres*winzpad));
		fftw_execute(FTplan1D);
		// cjfeng 07/05/2016
		// Normalization?
		for(i=0; i<winzpad; i++) FTout1D[i][1] = -FTout1D[i][1];
		for(i=0; i<npts; i++) ftir[i] = FTout1D[i][0];
		
		if(do2d) {
			for(p=0; p<npol; p++) {
				// Dress rephasing spectrum with Lorentzian profile. Note that it is so far purely real. 
				for(i=0; i<npts; i++) for(j=0; j<npts; j++) REPH[POL[p]][i*npts+j] += NetREPH[POL[p]][i*npts+j];
				if(TauP>0) {
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = 0.0;
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = 0.0;
					for(i=0; i<npts; i++) for(j=0; j<npts; j++) FTin2D[i*winzpad+j][0] = creal(REPH[POL[p]][i*npts+j]);
					// cjfeng 07/12/2016
					// Print REPH FID
					if (tresponse) {
						for(i=0; i<winzpad; i++) {
							for(j=0; j<winzpad; j++) {
								GNREAL re, im;
								re = FTin2D[i*winzpad+j][0];
								im = FTin2D[i*winzpad+j][1];
								fprintf(fid_reph[POL[p]], "%6.10f\t", re);
								// if(im<0) fprintf(fid_reph[POL[p]], "%6.10f%6.10fi\t", re, im);
								// else fprintf(fid_reph[POL[p]], "%6.10f+%6.10fi\t", re, im);
							}
							fprintf(fid_reph[POL[p]], "\n");
						}
					}
					
					// cjfeng 07/05/2016
					// Should it be FORWARD or -iwt?
					fftw_execute(FTplan2D);
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = FTout2D[i*winzpad+j][0]*exp(-((winzpad-i)+j)/(2*TauP*c*wres*winzpad));
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = -FTout2D[i*winzpad+j][1]*exp(-((winzpad-i)+j)/(2*TauP*c*wres*winzpad));
					// cjfeng 07/05/2016
					// Normalization?
					fftw_execute(FTplan2D);
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTout2D[i*winzpad+j][1] = -FTout2D[i*winzpad+j][1];
					for(i=0; i<npts; i++) for(j=0; j<npts; j++) REPH[POL[p]][i*npts+j] = FTout2D[i*winzpad+j][0] + I*FTout2D[i*winzpad+j][1];
				}
				// Dress non-rephasing spectrum with Lorentzian profile. Note that it is so far purely real. 
				for(i=0; i<npts; i++) for(j=0; j<npts; j++) NREPH[POL[p]][i*npts+j] += NetNREPH[POL[p]][i*npts+j];
				if(TauP>0) {
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = 0.0;
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = 0.0;
					for(i=0; i<npts; i++) for(j=0; j<npts; j++) FTin2D[i*winzpad+j][0] = creal(NREPH[POL[p]][i*npts+j]);
					// cjfeng 07/12/2016
					// Print NREPH FID
					if (tresponse) {
						for(i=0; i<winzpad; i++) {
							for(j=0; j<winzpad; j++) {
								GNREAL re, im;
								re = FTin2D[i*winzpad+j][0];
								im = FTin2D[i*winzpad+j][1];
								fprintf(fid_nreph[POL[p]], "%6.10f\t", re);
								// if(im<0) fprintf(fid_nreph[POL[p]], "%6.10f%6.10fi\t", re, im);
								// else fprintf(fid_nreph[POL[p]], "%6.10f+%6.10fi\t", re, im);
							}
							fprintf(fid_nreph[POL[p]], "\n");
						}
					}
					// cjfeng 07/05/2016
					// Should it be FORWARD or -iwt?
					fftw_execute(FTplan2D);
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][0] = FTout2D[i*winzpad+j][0]*exp(-(i+j)/(2*TauP*c*wres*winzpad));
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTin2D[i*winzpad+j][1] = -FTout2D[i*winzpad+j][1]*exp(-(i+j)/(2*TauP*c*wres*winzpad));
					// cjfeng 07/05/2016
					// Normalization?
					fftw_execute(FTplan2D);
					for(i=0; i<winzpad; i++) for(j=0; j<winzpad; j++) FTout2D[i*winzpad+j][1] = -FTout2D[i*winzpad+j][1];
					for(i=0; i<npts; i++) for(j=0; j<npts; j++) NREPH[POL[p]][i*npts+j] = FTout2D[i*winzpad+j][0] + I*FTout2D[i*winzpad+j][1];
				}
			}
		}

		// print the frequency axis and ftir spectrum
		for(i=0; i<npts; i++) fprintf(afp, "%6.10f\n", wstart+i*wres);
		for(i=0; i<npts; i++) fprintf(ffp, "%6.10f\n", ftir[i]);
	
		// Now print the 2D spectrum. 

		if(do2d) {
			GNCOMP cval = 0.0;
			for(i=0; i<npol; i++) {
				for(j=0; j<npts; j++) {
					for(k=0; k<npts; k++) {
						cval = REPH[POL[i]][j*npts+k];
						// Adding plus sign if the imaginary part is positive.
						if(cimag(cval)<0) fprintf(rfp[POL[i]], "%6.10f%6.10fi\t", creal(cval), cimag(cval));
						else fprintf(rfp[POL[i]], "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
						// Changed print syntax to reflect complex data -- MER 7/21/2015
						cval = NREPH[POL[i]][j*npts+k];
						// Adding plus sign if the imaginary part is positive.
						if(cimag(cval)<0) fprintf(nrfp[POL[i]], "%6.10f%6.10fi\t", creal(cval), cimag(cval));
						else fprintf(nrfp[POL[i]], "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
						// Changed print syntax to reflect complex data -- MER 7/21/2015
					}
					fprintf(rfp[POL[i]], "\n");
					fprintf(nrfp[POL[i]], "\n");
				}
			}
		}

	}
	// Measuring stop time and computing the total time for simulation.
	stoptime = omp_get_wtime();
	printf("Time: %2.4f seconds \n", stoptime-starttime);
	fprintf(lfp, "Time: %2.4f seconds \n", stoptime-starttime);
	fflush(lfp);

	// cjfeng 07/12/2016
	// Close 2D FID files
	if (tresponse && do2d) {
		for(i=0; i<npol; i++) fclose(fid_reph[i]);
		for(i=0; i<npol; i++) fclose(fid_nreph[i]);
	}
	graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
	return 0;
}

// cjfeng 07/05/2016
// Remained issues:
// 1. Is winzpad correct and reasonable?
// 2. Normalization of Fourier transform is not accounted.
// 3. Parallel scaling is not solved. What is the optimal automatic algorithm to determine the best nthreads?
// 4. Parallel computation on static averaging and Time-averaging approximation is not implemented.
