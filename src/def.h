#ifndef DEF_H
#define DEF_H

#include "fftw3.h"      // fftw for Fourier Transform
#include "macros.h"
#include "gstat.h"
#include "string2.h"
#include "lapacke.h"    // LAPACKE as C interface for solving eigen problems.
#include "complex.h"    // Complex utilities
// cjfeng 06/27/2016
// Tracking memory usage of g_spec
#include <sys/resource.h>

// cjfeng 07/14/2016
// Tracking execution time and parallel computation
#if OMP_PARALLEL 
  #include "omp.h"
#else
  #include "time.h"
#endif

#define THRESHOLD 16

#if DOUBLE_PRECISION
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

// Global variable declaration
// OpenBLAS function
// cjfeng 11/16/2016
extern void cgemv_(char*, int*, int*, GNCOMP*, GNCOMP*, int*, GNCOMP*, int*, GNCOMP*, GNCOMP*, int*);
extern void sgemv_(char*, int*, int*, GNCOMP*, GNCOMP*, int*, GNCOMP*, int*, GNCOMP*, GNCOMP*, int*);

extern int maxchar;

// Now spectra
extern GNCOMP *CorrFunc;
extern GNCOMP *NetCorrFunc;
extern GNREAL *popdecay;
extern GNREAL *popdecay2Q;
extern GNREAL *hann;
extern GNREAL *ftir;
extern GNREAL *netftir;
extern GNCOMP *REPH[4], *NREPH[4];
extern GNCOMP *NetREPH[4], *NetNREPH[4];

// Generic wavefunctions
extern GNCOMP ***psi_a;
extern GNCOMP ***psi_b1;
extern GNCOMP ***psi_b12;
extern GNCOMP *psi_cb[9];
extern GNCOMP *psi_ca[9];
extern GNCOMP *psi2Q;

extern GNREAL **Ham1QMat;
extern GNREAL ***Dip1QMat;
extern GNCOMP **U1QMat;

extern GNREAL **Ham1QAr;
extern GNREAL **Evals1QAr;
extern GNREAL **Dip1QAr[3];
extern GNREAL **ExDip1QAr[3];

extern GNREAL ***Dip2QMat;
// cjfeng 12/19/2018
// Using sparse matrix
extern GNREAL ***Dip2QSpMat;
extern int *Dip2Qrow;
extern int *Dip2Qcol;
extern GNCOMP **U2QMat;

// cjfeng 05/07/2018
// 2Q trajectories
extern GNREAL **Ham2QMat;

extern GNREAL **Ham2QAr;
extern GNREAL **Evals2QAr;
extern GNREAL **Dip2QAr[3];
extern GNREAL **ExDip2QAr[3];

#define f15 ( (1.0)/(5.0) )
#define f115 ( (1.0)/(15.0) )
#define f215 ( (2.0)/(15.0) )
#define f130 (-(1.0)/(30.0) )

#define getName(var) #var

#endif
