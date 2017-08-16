#include "mem_helper.h"

int allocate_all( int window, int win2d, int winzpad, int nosc, int nbuffer, int nthreads, int n2Q, int pert, int do2d, GNREAL tstep, GNREAL TauP, int npts, int nise ) {
	int error = 0;
	int i,j,k;
	int npopdec = 2*(win2d-window);	// equivalent with (window + T2/tstep + 1)
	int windowsq = window*window;
	int noscsq = nosc*nosc;
	int n2Qsq = n2Q*n2Q;
	int nptssq = npts*npts;
	const double pi = 3.14159265;

	int stick;
	if(nise) stick = 0;
	else stick = 1;

	if(!stick) {
		if(!error) {
			error = set_gncomp_1d_array(&CorrFunc, 1, window);
			if(!error) for(i=0; i<window; i++) CorrFunc[i] = 0.0;
		} else CorrFunc = NULL;
		if(!error) {
			error = set_gncomp_1d_array(&NetCorrFunc, 1, window);
			if(!error) for(i=0; i<window; i++) NetCorrFunc[i] = 0.0;
		} else NetCorrFunc = NULL;

		if(!error) {
			error = set_gnreal_1d_array(&popdecay, 1, npopdec);
		} else popdecay = NULL;

		for(i=0; i<4; i++) {
			if((!error) && (do2d)) {
				error = set_gncomp_1d_array(&REPH[i], 1, windowsq);
				if(!error) for(j=0; j<windowsq; j++) REPH[i][j] = 0.0;
			} else REPH[i] = NULL;
		}
		for(i=0; i<4; i++) {
			if((!error) && (do2d)) {
				error = set_gncomp_1d_array(&NREPH[i], 1, windowsq);
				if(!error) for(j=0; j<windowsq; j++) NREPH[i][j] = 0.0;
			} else NREPH[i] = NULL;
		}
		for(i=0; i<4; i++) {
			if((!error) && (do2d)) {
				error = set_gncomp_1d_array(&NetREPH[i], 1, windowsq);
				if(!error) for(j=0; j<windowsq; j++) NetREPH[i][j] = 0.0;
			} else NetREPH[i] = NULL;
		}
		for(i=0; i<4; i++) {
			if((!error) && (do2d)) {
				error = set_gncomp_1d_array(&NetNREPH[i], 1, windowsq);
				if(!error) for(j=0; j<windowsq; j++) NetNREPH[i][j] = 0.0;
			} else NetNREPH[i] = NULL;
		}
	
		// One-quantum arrays
		if(!error) error = set_gncomp_1d_array(&U1Q, 1, noscsq);
		if(!error) error = set_gncomp_1d_array(&U1Qs, 1, noscsq);
		for(i=0; i<3; i++) if(!error) set_gncomp_1d_array(&cDip1Q[i], 1, nosc);
		if(!error) error = set_gncomp_2d_array(&U1QMat, nbuffer, noscsq);

		for(i=0; i<3; i++) if((!error) && (do2d) ) error = set_gncomp_2d_array(&psi_a[i], win2d, nosc);
		for(i=0; i<3; i++) if((!error) && (do2d) ) error = set_gncomp_2d_array(&psi_b1[i], win2d, nosc);
		for(i=0; i<3; i++) if((!error) && (do2d) ) error = set_gncomp_2d_array(&psi_b12[i], win2d, nosc);
		
		for(i=0; i<3; i++) {
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

		// Two-quantum arrays
		for(i=0; i<9; i++) if((!error) && (do2d) ) error = set_gncomp_1d_array(&psi_ca[i], 1, win2d);
		for(i=0; i<9; i++) if((!error) && (do2d) ) error = set_gncomp_1d_array(&psi_cb[i], 1, win2d);
		if((!error) && (do2d) ) error = set_gncomp_1d_array(&psi2Q, 1, n2Q);
		if((!error) && (do2d) ) error = set_gncomp_1d_array(&U2Qs, 1, n2Qsq);
		if((!error) && (do2d) ) error = set_gncomp_2d_array(&U2QMat, nbuffer, n2Qsq);
	
		// Only used for static averaging / time averaging approximation
		ftir = NULL;
		netftir = NULL;
		for(i=0; i<3; i++) ExDip1QAr[i] = NULL;
		for(i=0; i<3; i++) ExDip2QAr[i] = NULL;
	} else { // Static averaging / time averaging approximation
		for(i=0; i<4; i++) {
			if((!error) && (do2d) ) {
				error = set_gncomp_1d_array(&REPH[i], 1, nptssq);
				if(!error) for(j=0; j<nptssq; j++) REPH[i][j] = 0.0;
			} else REPH[i] = NULL;
		}
		for(i=0; i<4; i++) {
			if((!error) && (do2d) ) {
				error = set_gncomp_1d_array(&NREPH[i], 1, nptssq);
				if(!error) for(j=0; j<nptssq; j++) NREPH[i][j] = 0.0;
			} else NREPH[i] = NULL;
		}
		for(i=0; i<4; i++) {
			if((!error) && (do2d) ) {
				error = set_gncomp_1d_array(&NetREPH[i], 1, nptssq);
				if(!error) for(j=0; j<nptssq; j++) NetREPH[i][j] = 0.0;
			} else NetREPH[i] = NULL;
		}
		for(i=0; i<4; i++) {
			if((!error) && (do2d) ) {
				error = set_gncomp_1d_array(&NetNREPH[i], 1, nptssq);
				if(!error) for(j=0; j<nptssq; j++) NetNREPH[i][j] = 0.0;
			} else NetNREPH[i] = NULL;
		}

		if(!error) {
			error = set_gnreal_1d_array(&ftir, 1, npts);
			if(!error) for(i=0; i<npts; i++) ftir[i] = 0.0;
		} else ftir = NULL;
		if(!error) {
			error = set_gnreal_1d_array(&netftir, 1, npts);
			if(!error) for(i=0; i<npts; i++) netftir[i] = 0.0;
		} else netftir = NULL;
		for(i=0; i<3; i++) if(!error) error = set_gnreal_2d_array(&ExDip1QAr[i], nthreads, nosc);
		for(i=0; i<3; i++) if(!error) error = set_gnreal_2d_array(&ExDip2QAr[i], nthreads, nosc*n2Q);

		// Only used for NISE.
		CorrFunc = NULL;
		NetCorrFunc = NULL;
		popdecay = NULL;
	
		// One-quantum arrays
		U1Q = NULL;
		U1Qs = NULL;
		for(i=0; i<3; i++) cDip1Q[i] = NULL;
		U1QMat = NULL;
		// Two-quantum arrays
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
	if(!error) error = set_gnreal_1d_array(&SitesBuffer, 1, nosc);
	else SitesBuffer = NULL;

	// FFT
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
		error = set_gnreal_1d_array(&hann, 1, window);
		if(!error) {
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

	if(!error) error = set_gnreal_2d_array(&Ham1QMat, nbuffer, noscsq);

	for(i=0; i<3; i++) {
		if(!error) error = set_gnreal_2d_array(&Dip1QMat[i], nbuffer, nosc);
		else Dip1QMat[i] = NULL;
	}

	if(!error) error = set_gnreal_2d_array(&Ham1QAr, nthreads, noscsq);
	else Ham1QAr = NULL;

	if(!error) error = set_gnreal_2d_array(&Evals1QAr, nthreads, nosc);
	else Evals1QAr = NULL;

	for(i=0; i<3; i++) {
		if(!error) error = set_gnreal_2d_array(&Dip1QAr[i], nthreads, nosc);
		else Dip1QAr[i] = NULL;
	}
	
	if(!error) error = set_gnreal_2d_array(&Ham2QAr, nthreads, n2Qsq);
	else Ham2QAr = NULL;

	if(!error) error = set_gnreal_2d_array(&Evals2QAr, nthreads, n2Q);
	else Evals2QAr = NULL;

	for(i=0; i<3; i++) {
		if((!error) && (do2d) && ( (!pert) || pertvec )) set_gnreal_2d_array(&Dip2QAr[i], nthreads, nosc*n2Q);
		else Dip2QAr[i] = NULL;
	}

	if((!error) && (do2d) ) error = set_gnreal_3d_array(&Dip2QMat, nbuffer, 3, nosc*n2Q);
	else Dip2QMat = NULL;

	return error; 
}

// Close files, free memory and exit gracefully
int graceful_exit( int error, int nbuffer, int win2d, int nthreads, int npol, int nise, int nosc) {
	int i;
	int stick;
	if(nise) stick = 0;
	else stick = 1;

	if(error==1) printf("Error allocating memory.\n");
	if(error==2) printf("Error opening files.\n");
	// Char array
	if(NNames!=NULL) for(i=0; i<nosc; i++) if(NNames[i]!=NULL) free(NNames[i]);
	if(CNames!=NULL) for(i=0; i<nosc; i++) if(CNames[i]!=NULL) free(CNames[i]);
	if(NNames!=NULL) free(NNames);
	if(CNames!=NULL) free(CNames);
	if(NNums!=NULL) free(NNums);
	if(CNums!=NULL) free(CNums);
	// File pointers
	if(ifp!=NULL) fclose(ifp);
	if(hfp!=NULL) fclose(hfp);
	if(sfp!=NULL) fclose(sfp);
	if(ffp!=NULL) fclose(ffp);
	if(lfp!=NULL) fclose(lfp);
	if(afp!=NULL) fclose(afp);
	for(i=0; i<10; i++) if(Trajfp[i]!=NULL) fclose(Trajfp[i]);
	for(i=0; i<npol; i++) if(rfp[POL[i]]!=NULL) fclose(rfp[POL[i]]);
	for(i=0; i<npol; i++) if(nrfp[POL[i]]!=NULL) fclose(nrfp[POL[i]]);
	for(i=0; i<3; i++) if(Dfp[i]!=NULL) fclose(Dfp[i]);
	//GNCOMP arrays
	for(i=0; i<4; i++) unset_gncomp_1d_array(REPH[i]);
	for(i=0; i<4; i++) unset_gncomp_1d_array(NREPH[i]);
	for(i=0; i<4; i++) unset_gncomp_1d_array(NetREPH[i]);
	for(i=0; i<4; i++) unset_gncomp_1d_array(NetNREPH[i]);
	unset_gncomp_1d_array(U1Q);
	unset_gncomp_1d_array(U1Qs);
	unset_gncomp_2d_array(U1QMat, nbuffer);
	unset_gncomp_1d_array(U2Qs);
	unset_gncomp_2d_array(U2QMat, nbuffer);
	unset_gncomp_1d_array(CorrFunc);
	unset_gncomp_1d_array(NetCorrFunc);

	for(i=0; i<3; i++) unset_gncomp_1d_array(cDip1Q[i]);
	for(i=0; i<3; i++) unset_gncomp_2d_array(psi_a[i], win2d);
	for(i=0; i<3; i++) unset_gncomp_2d_array(psi_b1[i], win2d);
	for(i=0; i<3; i++) unset_gncomp_2d_array(psi_b12[i], win2d);

	for(i=0; i<9; i++) unset_gncomp_1d_array(psi_ca[i]);
	for(i=0; i<9; i++) unset_gncomp_1d_array(psi_cb[i]);
	unset_gncomp_1d_array(psi2Q);

	//GNREAL arrays
	unset_gnreal_1d_array(hann);
	unset_gnreal_1d_array(popdecay);
	
	unset_gnreal_1d_array(ftir);
	unset_gnreal_1d_array(netftir);
	for(i=0; i<3; i++) unset_gnreal_2d_array(Dip1QAr[i], nthreads);
	for(i=0; i<3; i++) unset_gnreal_2d_array(Dip2QAr[i], nthreads);
	for(i=0; i<3; i++) unset_gnreal_2d_array(ExDip1QAr[i], nthreads);
	for(i=0; i<3; i++) unset_gnreal_2d_array(ExDip2QAr[i], nthreads);
	for(i=0; i<3; i++) unset_gnreal_2d_array(Dip1QMat[i], nbuffer);
	unset_gnreal_2d_array(Ham1QMat, nbuffer);

	if(SHIFTNDX!=NULL) free(SHIFTNDX);
	if(SHIFT!=NULL) free(SHIFT);
	for(i=0; i<nthreads; i++) {
		if(Ham1QAr!=NULL) if(Ham1QAr[i]!=NULL) free(Ham1QAr[i]);
		if(Evals1QAr!=NULL) if(Evals1QAr[i]!=NULL) free(Evals1QAr[i]);
		if(Ham2QAr!=NULL) if(Ham2QAr[i]!=NULL) free(Ham2QAr[i]);
		if(Evals2QAr!=NULL) if(Evals2QAr[i]!=NULL) free(Evals2QAr[i]);
	}
	if(Ham1QAr!=NULL) free(Ham1QAr);
	if(Evals1QAr!=NULL) free(Evals1QAr);

	if(Ham2QAr!=NULL) free(Ham2QAr);
	if(Evals2QAr!=NULL) free(Evals2QAr);

	// cjfeng 04/05/2016 	
	// The modified code due to changing the allocation of Dip2QMat
	unset_gnreal_3d_array(Dip2QMat, nbuffer, 3);

	unset_gnreal_1d_array(SitesBuffer);

	// FFT arrays
	if(FTin1D!=NULL) free(FTin1D);
	if(FTout1D!=NULL) free(FTout1D);
	if(FTplan1D!=NULL) fftw_destroy_plan(FTplan1D);
	if(FTin2D!=NULL) fftw_free(FTin2D);
	if(FTout2D!=NULL) fftw_free(FTout2D);
	if(FTplan2D!=NULL) fftw_destroy_plan(FTplan2D);
	return 0;
}

// cjfeng 07/19/2017
int set_eig_array(eig_array *eig, int n) {
	eig->n = n;
	set_gnreal_1d_array(&(eig->Ham), eig->n, eig->n);
	set_gnreal_1d_array(&(eig->Evals), eig->n, 1);
	set_int_1d_array(&(eig->isuppz), eig->n, 2);
	return 0;
}

// cjfeng 07/19/2017
int unset_eig_array(eig_array *eig) {
	unset_int_1d_array(eig->isuppz);
	unset_gnreal_1d_array(eig->Ham);
	unset_gnreal_1d_array(eig->Evals);
	return 0;
}

// cjfeng 04/11/2017
// Subroutines
int set_gnreal_1d_array(GNREAL **mat, int nrows, int ncols) {
	if ( nrows*ncols == 0 ) {
		printf("The size of %s is zero.\n", getName(mat));
		return 1;
	}
	else {
		*mat = (GNREAL *) malloc(nrows*ncols*sizeof(GNREAL));
		if ( *mat == NULL) {
			printf("Out of memory for allocating %s.\n", getName(mat));
			return 1;
		}
	}
	return 0;
}

int set_gnreal_2d_array(GNREAL ***mat, int dim1, int dim2) {
	if ( dim1*dim2 == 0) {
		printf("The size of %s is zero.\n", getName(mat));
		return 1;
	}
	else {
		int i;
		*mat = (GNREAL **) malloc(dim1*sizeof(GNREAL *));
		if ( *mat == NULL) {
			printf("Out of memory for allocating %s.\n", getName(mat));
			return 1;
		}
		else {
			for(i=0; i<dim1; i++) {
				(*mat)[i] = (GNREAL *) malloc(dim2*sizeof(GNREAL));
				if ( (*mat)[i] == NULL) {
					int k;
        	for(k=0; k<i; k++) free(*mat[k]);
        	free(*mat);
					printf("Out of memory for allocating the second dimention of %s[%d].\n", getName(mat), i);
					return 1;
				}
			}
		}
	}
	return 0;
}

int set_gnreal_3d_array(GNREAL ****mat, int dim1, int dim2, int dim3) {
	if ( dim1*dim2*dim3 == 0) {
		printf("The size of %s is zero.\n", getName(mat));
		return 1;
	}
	else {
		int i, j;
		*mat = (GNREAL ***) malloc(dim1*sizeof(GNREAL **));
		if ( *mat == NULL) {
			printf("Out of memory for allocating %s.\n", getName(mat));
			return 1;
		}
		else {
			for(i=0; i<dim1; i++) {
				(*mat)[i] = (GNREAL **) malloc(dim2*sizeof(GNREAL *));
				if ( (*mat)[i] == NULL) {
					int k;
					for(k=0; k<i; k++) free(*mat[k]);
					free(*mat);
					printf("Out of memory for allocating %s[%d].\n", getName(mat), i);
					return 1;
				}
				else {
					for(j=0; j<dim2; j++) {
						(*mat)[i][j] = (GNREAL *) malloc(dim3*sizeof(GNREAL));
						if ( (*mat)[i][j] == NULL) {
							int k, l;
							for(k=0; k<i; k++) {
								for(l=0; l<j; l++) free(*mat[k][l]);
								free(*mat[k]);
							}
							free(*mat);
							printf("Out of memory for allocating %s[%d][%d].\n", getName(mat), i, j);
							return 1;
						}
					}
				}
			}
		}
	}
	return 0;
}

int unset_gnreal_1d_array(GNREAL *mat) {
	if ( mat !=NULL ) free(mat);
	return 0;
}

int unset_gnreal_2d_array(GNREAL **mat, int dim1) {
	if ( mat !=NULL ) {
		int i;
		for(i=0; i<dim1; i++) if( mat[i] !=NULL ) free(mat[i]);
		free(mat);
	}
	return 0;
}

int unset_gnreal_3d_array(GNREAL ***mat, int dim1, int dim2) {
	if ( mat !=NULL ) {
		int i,j;
		for(i=0; i<dim1; i++) {
			if( mat[i] !=NULL ) {
				for(j=0; j<dim2; j++) {
					if( mat[i][j] != NULL) free(mat[i][j]);
				}
				free(mat[i]);
			}
		}
		free(mat);
	}
	return 0;
}

int set_gncomp_1d_array(GNCOMP **mat, int nrows, int ncols) {
	if ( nrows*ncols == 0 ) {
		printf("The size of %s is zero.\n", getName(mat));
		return 1;
	}
	else {
		*mat = (GNCOMP *) malloc(nrows*ncols*sizeof(GNCOMP));
		if ( *mat == NULL) {
			printf("Out of memory for allocating %s.\n", getName(mat));
			return 1;
		}
	}
	return 0;
}

int set_gncomp_2d_array(GNCOMP ***mat, int dim1, int dim2) {
	if ( dim1*dim2 == 0) {
		printf("The size of %s is zero.\n", getName(mat));
		return 1;
	}
	else {
		int i;
		*mat = (GNCOMP **) malloc(dim1*sizeof(GNCOMP *));
		if ( *mat == NULL) {
			printf("Out of memory for allocating %s.\n", getName(mat));
			return 1;
		}
		else {
			for(i=0; i<dim1; i++) {
				(*mat)[i] = (GNCOMP *) malloc(dim2*sizeof(GNCOMP));
				if ( (*mat)[i] == NULL) {
					int k;
					for(k=0; k<i; k++) free(*mat[k]);
					free(*mat);
					printf("Out of memory for allocating %s[%d].\n", getName(mat), i);
					return 1;
				}
			}
		}
	}
	return 0;
}

int unset_gncomp_1d_array(GNCOMP *mat) {
	if ( mat !=NULL ) free(mat);
	return 0;
}

int unset_gncomp_2d_array(GNCOMP **mat, int dim1) {
	if ( mat !=NULL ) {
		int i;
		for(i=0; i<dim1; i++) if( mat[i] !=NULL ) free(mat[i]);
		free(mat);
	}
	return 0;
}

int set_int_1d_array(int **mat, int nrows, int ncols) {
	if ( nrows*ncols == 0 ) {
		printf("The size of %s is zero.\n", getName(mat));
		return 1;
	}
	else {
		*mat = (int *) malloc(nrows*ncols*sizeof(int));
		if ( *mat == NULL) {
			printf("Out of memory for allocating %s.\n", getName(mat));
			return 1;
		}
	}
	return 0;
}

int unset_int_1d_array(int *mat) {
	if ( mat !=NULL ) free(mat);
	return 0;
}

// cjfeng 07/24/2017
int set_char_1d_array(char **mat, int dim1) {
	if ( dim1 == 0) {
		printf("The size of %s is zero.\n", getName(mat));
		return 1;
	}
	else {
		*mat = (char *) malloc(dim1*sizeof(char));
		if ( *mat == NULL) {
			printf("Out of memory for allocating the %s.\n", getName(mat));
			return 1;
		}
	}
	return 0;
}

// cjfeng 04/25/2017
int set_char_2d_array(char ***mat, int dim1, int dim2) {
	if ( dim1*dim2 == 0) {
		printf("The size of %s is zero.\n", getName(mat));
		return 1;
	}
	else {
		int i;
		*mat = (char **) malloc(dim1*sizeof(char *));
		if ( *mat == NULL) {
			printf("Out of memory for allocating the %s.\n", getName(mat));
			return 1;
		}
		else {
			for(i=0; i<dim1; i++) {
				(*mat)[i] = (char *) malloc(dim2*sizeof(char));
				if ( (*mat)[i] == NULL) {
					int j;
					for(j=0; j<i; j++) free(*mat[j]);
					free(*mat);
					printf("Out of memory for allocating %s[%d].\n", getName(mat), i);
					return 1;
				}
			}
		}
	}
	return 0;
}

// cjfeng 07/24/2017
int unset_char_1d_array(char *mat) {
	if ( mat !=NULL ) free(mat);
	return 0;
}

int unset_char_2d_array(char **mat, int dim1) {
	if ( mat !=NULL ) {
		int i;
		for(i=0; i<dim1; i++) if( mat[i] !=NULL ) free(mat[i]);
		free(mat);
	}
	return 0;
}
