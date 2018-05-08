#ifndef DEF_H
#define DEF_H

#include "fftw3.h"			// fftw for Fourier Transform

// #include "physics.h"		// Physical constants used in Gromacs
#include "macros.h"
#include "gstat.h"
#include "string2.h"

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

// Global variable declaration
// OpenBLAS function
// cjfeng 11/16/2016
extern void cgemv_(char*, int*, int*, GNCOMP*, GNCOMP*, int*, GNCOMP*, int*, GNCOMP*, GNCOMP*, int*);
extern void sgemv_(char*, int*, int*, GNCOMP*, GNCOMP*, int*, GNCOMP*, int*, GNCOMP*, GNCOMP*, int*);

extern int pertvec; // first-order perturbative correction on both site energies and vectors

extern int maxchar;

extern int POL[4];
extern int* SHIFTNDX;
extern real* SHIFT;

// File pointers
extern FILE *hfp, *Dfp[3], *sfp, *afp, *ffp, *lfp, *rfp[4], *nrfp[4], *ifp;
// cjfeng 05/07/2018
// hfp2Q: Two-quantum Hamiltonian, ham2Q.txt
// Dfp2Q: Two-quantum dipole moments, dip[xyz]2Q.txt
extern FILE *hfp2Q, *Dfp2Q[3];
// Trajfp is an array of trajectory file pointers: time-stamp, FTIR, rzzzz, rzzyy, rzyyz, rzyzy, nrzzzz, nrzzyy, nrzyyz, nrzyzy
extern FILE *Trajfp[10];

extern int *NNums, *CNums;
extern char **NNames, **CNames;

// Now spectra
extern GNCOMP *CorrFunc;
extern GNCOMP *NetCorrFunc;
extern GNREAL *popdecay;
extern GNREAL *hann;
extern GNREAL *ftir;
extern GNREAL *netftir;
extern fftw_complex *FTin1D, *FTout1D;
extern fftw_plan FTplan1D;
extern GNCOMP *REPH[4], *NREPH[4];
extern GNCOMP *NetREPH[4], *NetNREPH[4];
extern fftw_complex *FTin2D, *FTout2D;
extern fftw_plan FTplan2D;

extern GNCOMP *U1Q;
extern GNCOMP *U1Qs;
extern GNCOMP *cDip1Q[3];

// Generic wavefunctions
extern GNCOMP **psi_a[3];
extern GNCOMP **psi_b1[3];
extern GNCOMP **psi_b12[3];
extern GNCOMP *psi_cb[9];
extern GNCOMP *psi_ca[9];
extern GNCOMP *psi2Q;

extern GNREAL **Ham1QMat;
extern GNREAL **Dip1QMat[3];
extern GNCOMP **U1QMat;

extern GNREAL **Ham1QAr;
extern GNREAL **Evals1QAr;
extern GNREAL **Dip1QAr[3];
extern GNREAL **ExDip1QAr[3];

extern GNCOMP *U2Qs;

extern GNREAL ***Dip2QMat;
extern GNCOMP **U2QMat;

// cjfeng 05/07/2018
// 2Q trajectories
extern GNREAL **Ham2QMat;

extern GNREAL **Ham2QAr;
extern GNREAL **Evals2QAr;
extern GNREAL **Dip2QAr[3];
extern GNREAL **ExDip2QAr[3];

extern GNREAL *SitesBuffer;

#define f15 ( (1.0)/(5.0) )
#define f115 ( (1.0)/(15.0) )
#define f215 ( (2.0)/(15.0) )
#define f130 (-(1.0)/(30.0) )

extern GNREAL M_ijkl_IJKL[4][4];

// cjfeng 07/05/2017
#define getName(var) #var

#endif
