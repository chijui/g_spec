#include "spec_mod.h"

/*********************************
* Spectral parameter generations *
**********************************/
int gen_ham_2Q(GNREAL *Ham1Q, int nosc, GNREAL *Ham2Q, int n2Q, real delta) {
	int m, n, p;
	int n2Qsq = n2Q * n2Q;
	GNREAL sqrt2 = sqrt(2.0);
	// States are numbered in order 00, 10, 11, 20, 21, 22,...,(nosc-1)(nosc-1)
	// The index of state mn is (m+1)*m/2 + n.
	// cjfeng 04/14/2017
	for(m=0; m<n2Qsq; m++) Ham2Q[m] = 0.0;

	for(m=0; m<nosc; m++) {
		int mnosc = m*nosc;
		int NDXm2Q = NDX2Q(m,m)*n2Q;
		// Double excitation at site m
		Ham2Q[NDXm2Q + NDX2Q(m,m)] = 2*Ham1Q[mnosc+m]-delta;
		// Single excitations at sites m and n
		for(n=0; n<m; n++) Ham2Q[ NDX2Q(m,n)*n2Q + NDX2Q(m,n) ] = Ham1Q[mnosc+m]+Ham1Q[n*nosc+n];
		// Coupling between mm and mn states: <mm|H|mn> = sqrt(2)*<m|H|n>
		for(n=0; n<m; n++) Ham2Q[ NDXm2Q + NDX2Q(m,n) ] = sqrt2*Ham1Q[mnosc+n];
		// Coupling between mm and nm states: <mm|H|nm> = sqrt(2)*<m|H|n>
		for(n=m+1; n<nosc; n++) Ham2Q[ NDXm2Q + NDX2Q(n,m) ] = sqrt2*Ham1Q[mnosc+n];
		// Coupling between nm and mp states (n<p<m): <mn|H|mp> = <n|H|p>
		// cjfeng 07/03/2017
		// Reverse the order of n and p to reserve Cache access.
		for(n=0; n<m-1; n++) {
			int NDXn2Q = NDX2Q(m,n) * n2Q;
			int nnosc = n*nosc;
			for(p=n+1; p<m; p++) Ham2Q[NDXn2Q + NDX2Q(m,p) ] = Ham1Q[nnosc+p];
		}
		// Coupling between pm and mn states (n<m<p): <mn|H|pm> = <n|H|p>
		for(n=0; n<m; n++) {
			int NDXn2Q = NDX2Q(m,n)*n2Q;
			int nnosc = n*nosc;
			for(p=m+1; p<nosc; p++) Ham2Q[ NDXn2Q + NDX2Q(p,m) ] = Ham1Q[nnosc+p];
		}
		
		// Coupling between pm and nm states (m<n<p): <nm|H|pm> = <n|H|p>
		for(n=m+1; n<nosc; n++) {
			int NDXn2Q = NDX2Q(m,n)*n2Q;
			int nnosc = n*nosc;
			for(p=n+1; p<nosc; p++) Ham2Q[ NDXn2Q + NDX2Q(p,m) ] = Ham1Q[nnosc+p];
		}
	}
	for(m=0; m<n2Q; m++) {
		int m2Q = m * n2Q;
		for(n=0; n<m; n++) Ham2Q[m2Q+n] += Ham2Q[n*n2Q+m];
	}
	for(m=0; m<n2Q; m++) {
		int m2Q = m * n2Q;
		for(n=m+1; n<n2Q; n++) Ham2Q[m2Q+n] = Ham2Q[n*n2Q+m];
	}
	return 1;
}

int gen_dip_2Q(GNREAL *Dip1Q, GNREAL *Dip2Q, int nosc, int n2Q) {
	int m, n, M;
	int n2Qnosc = n2Q * nosc;
	GNREAL fac;
	GNREAL sqrt2 = sqrt(2.0);
	for(M=0; M<n2Qnosc; M++) Dip2Q[M] = 0.0;
	
	for(n=0; n<nosc; n++) {
		GNREAL Dipval = Dip1Q[n];
		int nn2Q = n*(n+1)/2;
		for(m=0; m<nosc; m++) {
			if(m!=n) fac = 1.0;
			else fac = sqrt2;
			if(m<n) M = (nn2Q+m);
			else M = ((m*(m+1)/2)+n);
			
			Dip2Q[M*nosc+m] = fac * Dipval;
		}
	}
	return 1;
}

// cjfeng 04/11/2017
int gen_ExDip2Q(GNREAL ***ExDip2QAr, GNREAL ***Dip2QAr, GNREAL **Ham1QAr, GNREAL **Ham2QAr, int nosc, int n2Q, int tid) {
	int i, n, N, k;
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

int gen_avg_ham(GNREAL **ham, GNREAL *avg_ham, GNREAL *hann, int nosc, int fr, int whann, int window, int nbuffer) {
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

int gen_avg_dip(GNREAL ***dip1Q, GNREAL ***dip2Q, GNREAL ***avg_dip1Q, GNREAL ***avg_dip2Q, int nosc, int n2Q, int fr, int nAvgd, int do2d, int pert, int pertvec, int whann, int window, int nbuffer) {
	int i, j, k;
	int nosc_n2Q = nosc * n2Q;
	int xfr_half_wind = (fr + ((int) (window/2))) % nbuffer;
	// Initialize averaged dipole moments.
	for(i=0; i<3; i++) for(j=0; j<nosc; j++) avg_dip1Q[i][nAvgd][j] = 0.0;
	if( (do2d) && ((!pert) || pertvec) ) for(i=0; i<3; i++) for(j=0; j<nosc_n2Q; j++) avg_dip2Q[i][nAvgd][j] = 0.0;
	if (whann) {
		// one-quantum dipole moment
		for(i=0; i<3; i++) for(j=0; j<nosc; j++) avg_dip1Q[i][nAvgd][j] = dip1Q[i][xfr_half_wind][j];
		// two-quantum dipole moment
		if( (do2d) && ( (!pert) || pertvec) ) {
			for(i=0; i<3; i++) for(j=0; j<nosc; j++) {
				int jj = j*n2Q;
				for(k=0; k<n2Q; k++) avg_dip2Q[i][nAvgd][jj+k] = dip2Q[xfr_half_wind][i][k*nosc+j];
			}
		}
	}
	else {
		// one-quantum dipole moment
		for(i=0; i<3; i++) for(j=0; j<nosc; j++) avg_dip1Q[i][nAvgd][j] = dip1Q[i][xfr_half_wind][j];
		// two-quantum dipole moment
		if( (do2d) && ( (!pert) || pertvec) ) {
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

// Modified on 04/14/2017
int gen_perturb_2Q_vectors(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evecs2Q, GNREAL *Evals2Q, int nosc, int n2Q, real delta ) {
	GNREAL *Temp2Q = malloc(n2Q*n2Q*sizeof(GNREAL));
	GNREAL cutoff = 1.0; // cm-1
	int m, n, i, j, mp, np;
	GNREAL sqrt2 = sqrt(2.0);
	GNREAL sqrt2inv = (1.0)/sqrt2;
	GNREAL val, anh, gap;
	int n2Qsq = n2Q * n2Q;
	// cjfeng 04/14/2017
	for(i=0; i<n2Qsq; i++) Evecs2Q[i] = 0.0;
	
	// Harmonic vectors
	for(i=0; i<nosc; i++) {
		int ni = i * nosc;
		for(j=0; j<=i; j++) {
			int nj = j * nosc;
			int Nij = NDX2Q(i,j)*n2Q;
			for(m=0; m<nosc; m++) {
				for(n=0; n<=m; n++) {
					Temp2Q[Nij + NDX2Q(m,n)] = Evecs1Q[ni+m]*Evecs1Q[nj+n] + Evecs1Q[nj+m]*Evecs1Q[ni+n]*(1.0*(m!=n) + sqrt2inv*(m==n))*(1.0*(i!=j) + sqrt2inv*(i==j));
				}
			}
		}
	}

	for(i=0; i<n2Qsq; i++) Evecs2Q[i] = Temp2Q[i];

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
						// cjfeng 04/14/2017
						// It may have a bug about accessing i or j in either term.
						for(j=0; j<n2Q; j++) Evecs2Q[j*n2Q+NDX2Q(m,n)] -= (anh/gap)*Temp2Q[j*n2Q+NDX2Q(mp,np)];
					}
				}
			}
		}
	}

	// Probability normalization
	// cjfeng 04/14/2017
	// Use Temp2Q as eigenvector matrix transpose
	trans_real(Evecs2Q, Temp2Q, n2Q, n2Q, 1);
	// Compute the normalization factor
	for(m=0; m<n2Q; m++) {
		val = 0.0;
		int nm = m * n2Q;
		for(i=0; i<n2Q; i++) val += Temp2Q[nm+i] * Temp2Q[nm+i];
		val = (1.0) / sqrt(val);
		for(i=0; i<n2Q; i++) Temp2Q[nm+i] *= val;
	}
	trans_real(Temp2Q, Evecs2Q, n2Q, n2Q, 1);
	
	// for(m=0; m<n2Q; m++) {
	// 	val = 0.0;
	// 	for(i=0; i<n2Q; i++) val += Evecs2Q[i*n2Q+m]*Evecs2Q[i*n2Q+m];
	// 	val = (1.0)/sqrt(val);
	// 	for(i=0; i<n2Q; i++) Evecs2Q[i*n2Q+m] *= val;
	// }
	free(Temp2Q);
	return 1;
}

int gen_perturb_2Q_matrix(GNREAL *Evecs1Q, GNREAL *Evals1Q, GNREAL *Evecs2Q, GNREAL *Evals2Q, GNREAL *Temp2Q, int nosc, int n2Q, real delta ) {
	int m,n,i,mp,np;
	GNREAL sqrt2 = sqrt(2.0);
	GNREAL sqrt2inv = (1.0)/sqrt2;
	GNREAL fac;
	int n2Qsq = n2Q * n2Q;

	for(i=0; i<n2Qsq; i++) Temp2Q[i] = 0.0;

	for(i=0; i<nosc; i++) {
		int ni = i * nosc;
		for(m=0; m<nosc; m++) {
			for(n=0; n<=m; n++) {
				int NDXn2Q = NDX2Q(m,n)*n2Q;
				for(mp=0; mp<nosc; mp++) {
					for(np=0; np<=mp; np++) {
						fac = 2.0*delta;
						if(m==n) fac *= sqrt2inv;
						if(mp==np) fac *= sqrt2inv;
						Temp2Q[NDXn2Q+NDX2Q(mp,np)] -= fac*Evecs1Q[ni+m]*Evecs1Q[ni+n]*Evecs1Q[ni+mp]*Evecs1Q[ni+np];
					}
				}
			}
		}
	}
	return 1;
}

int gen_popdecay(GNREAL *popdecay, GNREAL tstep, GNREAL TauP, int window, int win2d) {
	int i;
	int npopdec = 2*(win2d-window);
	GNREAL expfac = tstep/(2*TauP);
	for(i=0; i<npopdec; i++) popdecay[i] = exp(-i*expfac);
	return 0;
}

int gen_popdecay_elec(GNREAL **popdecay_elec, GNREAL **elec, GNREAL tstep, int nosc, int fr, int nelecpts, int window, int win2d) {
	int n, t; // n: oscillator index, t: time index 
	int npopdec = 2*(win2d-window);
	GNREAL k; // k: instantaneous rate in 1/fs
	// Initial assignment: Uniform excitation.
	for(n=0; n<nosc; n++) popdecay_elec[0][n] = 1;

	for(t=1; t<npopdec; t++) {
		int xfr = (fr + t) % nelecpts;
		for(n=0; n<nosc; n++) {
			// k = 0.001 * ( (elec[xfr][n]*elec[xfr][n]) + 0.0016) / 0.0023;	// Based on the equation fit from Ala-Ala MA absolute surface data.
			k = 0.0009579 * exp( - (elec[xfr][n] + 0.0206) / 0.002269) + 0.0007325;	// Based on the equation fit from Ala-Ala MA absolute surface data.
			popdecay_elec[t][n] = popdecay_elec[t-1][n] * exp(-k*tstep/2);	// Analogous to exp(-tstep/(2*TauP));
		}
	}
	return 0;
}

int gen_BLAS_gemv_opt(BLAS_gemv_opt *gemv_opt) {
	gemv_opt->ta = 'N';
	gemv_opt->inc = 1;
	gemv_opt->alpha = 1;
	gemv_opt->beta = 0.0;
	return 0;
}

/********************************************************************************
* Spectral calculations.							*
********************************************************************************/

// cjfeng 06/27/2017
// Convert the algorithm to matrix vector multiplication for speeding up correlation function computations.
// int calc_CorrFunc(GNCOMP *CorrFunc, GNCOMP **U1QMat, GNCOMP **psi_1d[3], GNREAL ***Dip1QMat, GNCOMP *cDip1Q[3], int nosc, int window, int fr, int nbuffer, double expfac, double wo, BLAS_gemv_opt *gemv_opt, int nthreads) {
// 	int noscsq = nosc * nosc;
// 	int i, j, n;
// 	int xfr;
// 	GNCOMP cval;
// 	initialize_psi(psi_1d, Dip1QMat, 0, fr, nosc, 3);
// 	
// 	for(n=0; n<window; n++) {
// 		xfr = (fr+n) % nbuffer;
// 		cval = 0.0;
// 		#if OMP_PARALLEL
// 		#pragma omp parallel if(nosc >= NBCOMP && nthreads>1) shared(nosc, cDip1Q, Dip1QMat, xfr, fr) private(i,j)
// 		#pragma omp for schedule(guided) collapse(2) reduction(+:cval)
// 		#endif
// 		for(i=0; i<3; i++) {
// 			for(j=0; j<nosc; j++) {
// 				cDip1Q[i][j] = psi_1d[i][0][j];
// 				cval += Dip1QMat[i][xfr][j] * cDip1Q[i][j];
// 			}
// 		}
// 		CorrFunc[n] += cval;
// 		
// 		// cjfeng 06/27/2017
// 		// propagate psi_1d
// 		#if BLAS_subroutines
// 		if(nosc > nosc_switch) for(i=0; i<3; i++) cgemv_(&(gemv_opt->ta), &nosc, &nosc, &(gemv_opt->alpha), U1QMat[xfr], &nosc, psi_1d[i][0], &(gemv_opt->inc), &(gemv_opt->beta), psi_1d[i][1], &(gemv_opt->inc));	
// 		else for(i=0; i<3; i++) mvmult_comp(U1QMat[xfr], psi_1d[i][0], psi_1d[i][1], nosc, nthreads);
// 		#else
// 		for(i=0; i<3; i++) mvmult_comp(U1QMat[xfr], psi_1d[i][0], psi_1d[i][1], nosc, nthreads);
// 		#endif
// 		for(i=0; i<3; i++) copy_gncomp(psi_1d[i][1], psi_1d[i][0], nosc, 1);
// 	}
// 	return 0;
// }

// int gen_UMat(GNCOMP **UMat, GNREAL *Ham, GNREAL *Evals, int N, int xfr, double expfac, double wo) {
// 	int i, j;
// 	int Nsq = N * N;
// 	// Initialize UMat
// 	for(i=0; i<Nsq; i++) UMat[xfr][i] = 0.0;
// 	// Compute U2Q elements
// 	for(i=0; i<N; i++) {
// 		int ni = i*N;
// 		GNREAL cre, cim;
// 		for(j=0; j<N; j++) {
// 			int k, nj = j*N;
// 			cre = 0.0;
// 			cim = 0.0;
// 			for(k=0; k<N; k++) {
// 				cre += Ham[ni+k] * Ham[nj+k] * cos( expfac * (Evals[k] - wo) );
// 				cim += Ham[ni+k] * Ham[nj+k] * sin( expfac * (Evals[k] - wo) );
// 			}
// 			UMat[xfr][ni+j] = cre + I*cim;
// 		}
// 	}
// 	return 0;
// }

GNREAL orient(GNREAL **Dip1[3], GNREAL **Dip2[3], int tid, int ndxA, int ndxB, int ndxC, int ndxD, int pol) {
	GNREAL Val = 0.0;
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

// int gen_pert_prop( GNCOMP *U1Q, GNCOMP *U2Q, int nosc, int n2Q, double expfac, double delta) {
// 	int i,j,ip,jp;
// 	GNCOMP fac;
// 	double sqrt2inv = 0.707106781186547;
// 	// negative sign on sin() term since delta is provided as a positive argument.
// 	GNCOMP facnew = sqrt2inv * ( cos(0.5*expfac*delta) - I*sin(0.5*expfac*delta) );
// 	// cjfeng 08/01/2016
// 	// No initialization required
// 	// for(i=0; i<n2Qsq; i++) U2Q[i] = 0.0;
// 	for(i=0; i<nosc; i++) {
// 		int ni = i*nosc;
// 		for(j=0; j<=i; j++) {
// 			int nj = j*nosc;
// 			int NDX2Qij = NDX2Q(i,j)*n2Q;
// 			for(ip=0; ip<nosc; ip++) {
// 				GNCOMP U1Qniip = U1Q[ni+ip];
// 				GNCOMP U1Qnjip = U1Q[nj+ip];
// 				for(jp=0; jp<=ip; jp++) {
// 					fac = 1.0;
// 					if(i==j) fac *= facnew;
// 					if(ip==jp) fac *= facnew;
// 					U2Q[NDX2Qij+NDX2Q(ip,jp)] = fac*( U1Qniip*U1Q[nj+jp] + U1Qnjip*U1Q[ni+jp] );
// 				}
// 			}
// 		}
// 	}
// 	return 1;
// }

int calc_2dir(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], GNREAL **ExDip2QAr[3], int tid, int nosc, int n2Q, int npts, GNREAL res, GNREAL start, GNREAL stop, GNCOMP *REPH[4], GNCOMP *NREPH[4], int POL[4], int npol, int reph, int nreph) {
	int i, a, b, c;
	int ndx1, ndx3;

	GNREAL *Evals1Q = Evals1QAr[tid];
	GNREAL *Evals2Q = Evals2QAr[tid];

	for(a=0; a<nosc; a++) {
		ndx1 = floor( (Evals1Q[a]-start)/res + 0.5 );
		if( (ndx1>0) && (ndx1<npts) ) {
			int ndx1npts = ndx1*npts;
			for(b=0; b<nosc; b++) {
				// Pathway 1 - rephasing & non-rephasing
				ndx3 = floor( (Evals1Q[b]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) {
					if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);
					if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);
				}

				// Pathway 2 - rephasing
				ndx3 = floor( (Evals1Q[b]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, a, b, POL[i]);

				// Pathway 2 - non-rephasing
				ndx3 = floor( (Evals1Q[a]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, b, a, POL[i]);
				for(c=0; c<n2Q; c++) {
					// Pathway 3 - rephasing
					ndx3 = floor( ( (Evals2Q[c]-Evals1Q[a])-start)/res + 0.5 );
					if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip2QAr, tid, a, b, b*n2Q+c, a*n2Q+c, POL[i]);

					// Pathway 3 - non-rephasing
					ndx3 = floor( ( (Evals2Q[c]-Evals1Q[b])-start)/res + 0.5 );
					if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip2QAr, tid, a, b, a*n2Q+c, b*n2Q+c, POL[i]);
				}
			}
		}
	}
	return 1;
}

int calc_2dir_pert(GNREAL **Evals1QAr, GNREAL **Evals2QAr, GNREAL **ExDip1QAr[3], int tid, int nosc, int n2Q, int npts, GNREAL res, GNREAL start, GNREAL stop, GNCOMP *REPH[4], GNCOMP *NREPH[4], int POL[4], int npol, int reph, int nreph) {
	int i, a, b;
	int ndx, ndx1, ndx3;

	GNREAL *Evals1Q = Evals1QAr[tid];
	GNREAL *Evals2Q = Evals2QAr[tid];

	for(a=0; a<nosc; a++) {
		ndx1 = floor( (Evals1Q[a]-start)/res + 0.5 );
		if( (ndx1>0) && (ndx1<npts) ) {
			int ndx1npts = ndx1*npts;
			for(b=0; b<nosc; b++) {
				// Pathway 1 - rephasing & non-rephasing
				ndx3 = floor( (Evals1Q[b]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) {
					if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);
					if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);
				}

				// Pathway 2 - rephasing
				ndx3 = floor( (Evals1Q[b]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, a, b, POL[i]);

				// Pathway 2 - non-rephasing
				ndx3 = floor( (Evals1Q[a]-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] += orient(ExDip1QAr, ExDip1QAr, tid, a, b, b, a, POL[i]);

				if(a>b) ndx = ((a*(a+1)/2)+b);
				else ndx = ((b*(b+1)/2)+a);
				
				// Pathway 3 - rephasing (0,0) -> (0,a) -> (a,a) -> (ab,a)
				ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[a])-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);

				// Pathway 3 - rephasing (0,0) -> (0,a) -> (b,a) -> (ab,a)
				ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[a])-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(reph) for(i=0; i<npol; i++) REPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, b, a, b, POL[i]);

				// Pathway 3 - non-rephasing (0,0) -> (0,a) -> (a,a) -> (a,ab)
				ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[a])-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, a, b, b, POL[i]);

				// Pathway 3 - non-rephasing (0,0) -> (0,a) -> (b,a) -> (b,ab)
				ndx3 = floor( ( (Evals2Q[ndx]-Evals1Q[b])-start)/res + 0.5 );
				if( (ndx3>0) && (ndx3<npts) ) if(nreph) for(i=0; i<npol; i++) NREPH[POL[i]][ndx1npts+ndx3] -= orient(ExDip1QAr, ExDip1QAr, tid, a, b, b, a, POL[i]);
			}	
		}
	}
	return 1;
}

// int initialize_psi(GNCOMP ***psi, GNREAL **Dip1QMat[3], int tinit, int xfr, int nosc, int nelem) {
// 	int i, j;
// 	for(i=0; i<nelem; i++) for(j=0; j<nosc; j++) psi[i][tinit][j] = Dip1QMat[i][xfr][j];
// 	return 0;
// }

// int initialize_psi2Q(GNCOMP ***psi, GNREAL ***Dip2QMat, GNCOMP **psi_2Q, int nosc, int n2Q, int nelem, int tau1, int tau2, int nthreads) {
// 	int i, j;
// 	for(i=0; i<nelem; i++) {
// 		int ii = i*3;
// 		for(j=0; j<nelem; j++) {
// 			mvmult_comp_trans_x(Dip2QMat[tau1+tau2][j], psi[i][tau1+tau2], psi_2Q[ii+j], nosc, n2Q, nthreads);
// 		}
// 	}
// 	return 0;
// }

// int propagate_psi1Q(GNCOMP **psi[3], GNREAL **Dip1QMat[3], GNCOMP **U1QMat, int nosc, int fr, int tau0, int xfr0, int nbuffer, int winsize, BLAS_gemv_opt *gemv_opt, int nthreads) {
// 	int i, tau;
// 	// Initial state
// 	initialize_psi(psi, Dip1QMat, tau0, xfr0, nosc, 3);
// 	#if BLAS_subroutines
// 	if(nosc > nosc_switch) {
// 		for(i=0; i<3; i++) for(tau=tau0; tau<tau0+winsize-1; tau++) cgemv_(&(gemv_opt->ta), &nosc, &nosc, &(gemv_opt->alpha), U1QMat[(fr+tau)%nbuffer], &nosc, psi[i][tau], &(gemv_opt->inc), &(gemv_opt->beta), psi[i][tau+1], &(gemv_opt->inc));
// 	}
// 	else {
// 		for(i=0; i<3; i++) for(tau=tau0; tau<tau0+winsize-1; tau++) mvmult_comp(U1QMat[(fr+tau)%nbuffer], psi[i][tau], psi[i][tau+1], nosc, nthreads );
// 	}
// 	#else
// 	for(i=0; i<3; i++) for(tau=tau0; tau<tau0+winsize-1; tau++) mvmult_comp(U1QMat[(fr+tau)%nbuffer], psi[i][tau], psi[i][tau+1], nosc, nthreads );
// 	#endif
// 	return 0;
// }

// int propagate_psi_2Q(GNCOMP **U2QMat, GNCOMP **psi_ca, GNCOMP **psi_cb, GNCOMP *psi2Q, int n2Q, int xfr3, BLAS_gemv_opt *gemv_opt, int nthreads) {
// 	int i, j;
// 	for(i=0; i<9; i++) {
// 		// First initialize psi_ca and propagate it.
// 		for(j=0; j<n2Q; j++) psi2Q[j] = psi_ca[i][j];
// 		#if BLAS_subroutines
// 		if (n2Q > nosc_switch) cgemv_(&(gemv_opt->ta), &n2Q, &n2Q, &(gemv_opt->alpha), U2QMat[xfr3], &n2Q, psi2Q, &(gemv_opt->inc), &(gemv_opt->beta), psi_ca[i], &(gemv_opt->inc));
// 		else mvmult_comp(U2QMat[xfr3], psi2Q, psi_ca[i], n2Q, nthreads);
// 		#else
// 		mvmult_comp(U2QMat[xfr3], psi2Q, psi_ca[i], n2Q, nthreads);
// 		#endif
// 		// Initialize psi_cb and propagate it.
// 		for(j=0; j<n2Q; j++) psi2Q[j] = psi_cb[i][j];
// 		#if BLAS_subroutines
// 		if (n2Q > nosc_switch) cgemv_(&(gemv_opt->ta), &n2Q, &n2Q, &(gemv_opt->alpha), U2QMat[xfr3], &n2Q, psi2Q, &(gemv_opt->inc), &(gemv_opt->beta), psi_cb[i], &(gemv_opt->inc));
// 		else mvmult_comp(U2QMat[xfr3], psi2Q, psi_cb[i], n2Q, nthreads);
// 		#else
// 		mvmult_comp(U2QMat[xfr3], psi2Q, psi_cb[i], n2Q, nthreads);
// 		#endif
// 	}
// 	return 0;
// }

// int propagate_psi_2Q_pert(GNCOMP **U1QMat, GNCOMP *U2Qs, GNCOMP **psi_ca, GNCOMP **psi_cb, GNCOMP *psi2Q, int nosc, int n2Q, real delta, double expfac, int xfr3, int nthreads) {
// 	int i, j;
// 	gen_pert_prop(U1QMat[xfr3], U2Qs, nosc, n2Q, expfac, delta);
// 	for(i=0; i<9; i++) {
// 		// First initialize psi_ca and propagate it.
// 		for(j=0; j<n2Q; j++) psi2Q[j] = psi_ca[i][j];
// 		mvmult_comp(U2Qs, psi2Q, psi_ca[i], n2Q, nthreads);
// 		// Initialize psi_cb and propagate it.
// 		for(j=0; j<n2Q; j++) psi2Q[j] = psi_cb[i][j];
// 		mvmult_comp(U2Qs, psi2Q, psi_cb[i], n2Q, nthreads);
// 	}
// 	return 0;
// }

// cjfeng 03/23/2017
// This routine aims for computing the nonlinear response function in nise algorithm.
// int calc_2dir_nise(GNCOMP ***psi_a, GNCOMP ***psi_b1, GNCOMP ***psi_b12, int *POL, int nosc, int n2Q, int npol, GNREAL popfac1Q, GNREAL popfac2Q, int tau1, int tau2, int tau3, int xfr1, int xfr2, int xfr3, int window, GNCOMP **REPH, GNCOMP **NREPH, int nthreads) {
// 	int tau1window = tau1 * window;
// 	int P, p, i, j, k, l, n, a, b;
// 	// cjfeng 07/02/2017
// 	GNCOMP *psi2Qa_aa, *psi2Qa_ba, *psi2Qb1_aa, *psi2Qb1_ab, *psi2Qb1_ba, *psi2Qb1_bb; 
// 	set_gncomp_1d_array(&psi2Qa_aa, n2Q, 1);
// 	set_gncomp_1d_array(&psi2Qa_ba, n2Q, 1);
// 	set_gncomp_1d_array(&psi2Qb1_aa, n2Q, 1);
// 	set_gncomp_1d_array(&psi2Qb1_ab, n2Q, 1);
// 	set_gncomp_1d_array(&psi2Qb1_ba, n2Q, 1);
// 	set_gncomp_1d_array(&psi2Qb1_bb, n2Q, 1);
// 	// cjfeng 07/03/2017
// 	// Change the order of loops
// 	// cjfeng 07/02/2017
// 	// Further optimzations to reduce number of 2Q computations.
// 	// Recognizing that i is always equal to a, and 2Q transitions are only related to i, and j,
// 	// we can pull out the i terms first and reduce operations.
// 	// p sums over the forms iiii, iijj, ijji, and ijij.
// 	for(a=0; a<3; a++) {
// 		// cjfeng 07/02/2017
// 		GNREAL *Dip2Qa = Dip2QMat[xfr3][a];
// 		mvmult_comp_trans_x(Dip2Qa, psi_a[a][tau1+tau2+tau3], psi2Qa_aa, nosc, n2Q, nthreads);
// 		mvmult_comp_trans_x(Dip2Qa, psi_b1[a][tau1+tau2+tau3], psi2Qb1_aa, nosc, n2Q, nthreads);
// 		for(b=0; b<3; b++) {
// 			// cjfeng 07/02/2017
// 			GNREAL *Dip2Qb = Dip2QMat[xfr3][b];
// 			if ( b != a ) {
// 				mvmult_comp_trans_x(Dip2Qb, psi_a[a][tau1+tau2+tau3], psi2Qa_ba, nosc, n2Q, nthreads);
// 				mvmult_comp_trans_x(Dip2Qa, psi_b1[b][tau1+tau2+tau3], psi2Qb1_ab, nosc, n2Q, nthreads);
// 				mvmult_comp_trans_x(Dip2Qb, psi_b1[a][tau1+tau2+tau3], psi2Qb1_ba, nosc, n2Q, nthreads);
// 				mvmult_comp_trans_x(Dip2Qb, psi_b1[b][tau1+tau2+tau3], psi2Qb1_bb, nosc, n2Q, nthreads);
// 			}
// 			else {
// 				for(n=0; n<n2Q; n++) psi2Qa_ba[n] = psi2Qa_aa[n];
// 				for(n=0; n<n2Q; n++) psi2Qb1_ab[n] = psi2Qb1_aa[n];
// 				for(n=0; n<n2Q; n++) psi2Qb1_ba[n] = psi2Qb1_aa[n];
// 				for(n=0; n<n2Q; n++) psi2Qb1_bb[n] = psi2Qb1_aa[n];
// 			}
// 			for(p=0; p<4; p++) {
// 				if(p==0) { i=a; j=a; k=a; l=a; }
// 				else if(p==1) { i=a; j=a; k=b; l=b; }
// 				else if(p==2) { i=a; j=b; k=b; l=a; }
// 				else { i=a; j=b; k=a; l=b; } // p==3
// 				// For p==0, we only add the signal when a==b.
// 				if( (p!=0) || (a==b) ) {
// 					GNCOMP REPH_private = 0.0;
// 					GNCOMP NREPH_private = 0.0;
// 					GNREAL *Dipj, *Dipk, *Dipl, *Dip2Q;
// 					Dipj = Dip1QMat[j][xfr1];
// 					Dipk = Dip1QMat[k][xfr2];
// 					Dipl = Dip1QMat[l][xfr3];
// 					Dip2Q = Dip2QMat[xfr3][l];
// 					int ik = i*3 + k;
// 					int jk = j*3 + k;
// 					// cjfeng 07/03/2017
// 					// Replaced the original code by sub-routines
// 					pathway_double_sides(Dipl, Dipj, psi_b12[k][tau1+tau2+tau3], psi_a[i][tau1], nosc, popfac1Q, &REPH_private);
// 					pathway_double_sides(Dipl, Dipk, psi_b1[j][tau1+tau2+tau3], psi_a[i][tau1+tau2], nosc, popfac1Q, &REPH_private);
// 					pathway_REPH_2Q(psi_cb[jk], psi2Qa_aa, psi2Qa_ba, n2Q, p, popfac2Q, &REPH_private);
// 					pathway_single_side(Dipl, Dipj, psi_b12[k][tau1+tau2+tau3], psi_a[i][tau1], nosc, popfac1Q, &NREPH_private);
// 					pathway_double_sides(Dipl, Dipk, psi_a[j][tau1+tau2+tau3], psi_b1[i][tau1+tau2], nosc, popfac1Q, &NREPH_private);
// 					pathway_NREPH_2Q(psi_ca[ik], psi2Qb1_aa, psi2Qb1_ba, psi2Qb1_ab, psi2Qb1_bb, n2Q, p, popfac2Q, &NREPH_private);
// 					// cjfeng 07/03/2017
// 					// Computing different polarized spectrum
// 					for(P=0; P<npol; P++) {
// 						int nP = POL[P];
// 						GNREAL orient_fac = M_ijkl_IJKL[P][p];
// 						REPH[nP][tau1window+tau3] += orient_fac * REPH_private;
// 						NREPH[nP][tau1window+tau3] += orient_fac * NREPH_private;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	unset_gncomp_1d_array(psi2Qa_aa);
// 	unset_gncomp_1d_array(psi2Qa_ba);
// 	unset_gncomp_1d_array(psi2Qb1_aa);
// 	unset_gncomp_1d_array(psi2Qb1_ab);
// 	unset_gncomp_1d_array(psi2Qb1_ba);
// 	unset_gncomp_1d_array(psi2Qb1_bb);
// 	return 0;
// }

// cjfeng 07/05/2017
// Grouping individual pathways into non-redundant functions
// int pathway_single_side(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val) {
// 	// Pathway involving only the ket side: 
// 	//
// 	// Non-rephasing pathway:
// 	//  | b 0 |       -- tau1+tau2+tau3
// 	//  | 0 0 |       -- tau1+tau2
// 	//  | a 0 |       -- tau1
// 	//  | 0 0 |       -- 0
// 	//  Dipl * psi_b1
// 	//  Dipj * psi_a    
// 	GNCOMP cval1, cval2;
// 	int n;
// 	cval1 = 0.0; cval2 = 0.0;
// 	for(n=0; n<nosc; n++) cval1 += Dip1[n] * psi1[n]; // Dipl and psi_a
// 	for(n=0; n<nosc; n++) cval2 += Dip2[n] * psi2[n]; // Dipk and psi_b1
// 	*val += cval1 * cval2 * popfac;
// 	return 0;
// }
// 
// int pathway_double_sides(GNREAL *Dip1, GNREAL *Dip2, GNCOMP *psi1, GNCOMP *psi2, int nosc, GNREAL popfac, GNCOMP *val) {
// 	// Pathways involving both sides:
// 	//
// 	// Rephasing pathways:
// 	//  | b 0 |         | b 0 |       -- tau1+tau2+tau3
// 	//  | 0 0 |         | b a |       -- tau1+tau2
// 	//  | 0 a |         | 0 a |       -- tau1
// 	//  | 0 0 |         | 0 0 |       -- 0
// 	//  Dipl * psi_b12  Dipl * psi_b1
// 	//  Dipj * psi_a    Dipk * psi_a 
// 	//
// 	// Non-rephasing pathway:
// 	//  | a 0 |       -- tau1+tau2+tau3
// 	//  | a b |       -- tau1+tau2
// 	//  | a 0 |       -- tau1
// 	//  | 0 0 |       --
// 	//  Dipl * psi_a
// 	//  Dipk * psi_b1
// 	// Note that the complex conjugate is required for the R.H.S.
// 	GNCOMP cval1, cval2;
// 	int n;
// 	cval1 = 0.0; cval2 = 0.0;
// 	for(n=0; n<nosc; n++) cval1 += Dip1[n] * psi1[n];
// 	for(n=0; n<nosc; n++) cval2 += Dip2[n] * psi2[n];
// 	*val += cval1 * conj(cval2) * popfac;
// 
// 	return 0;
// }
// 
// int pathway_REPH_2Q(GNCOMP *psi2Qcb, GNCOMP *psi2Qa_aa, GNCOMP *psi2Qa_ba, int n2Q, int p, GNREAL popfac2Q, GNCOMP *REPH_private) {
// 	// Final rephasing pathway:
// 	//
// 	//  | c a |
// 	//  | b a |
// 	//  | 0 a |
// 	//  | 0 0 |
// 	//
// 	// The contribution to the rephasing spectrum is the orientationally averaged value
// 	// ( psi_cb[j*3+k] ) * ( Dip2Q[tau1+tau2+tau3] * psi_a[tau1+tau2+tau3] )^\dagger
// 	GNCOMP cval1 = 0.0;
// 	int n;
// 	// cjfeng 07/02/2017
// 	if ( p == 0 || p == 2) for(n=0; n<n2Q; n++) cval1 += psi2Qcb[n] * conj(psi2Qa_aa[n]);
// 	else for(n=0; n<n2Q; n++) cval1 += psi2Qcb[n] * conj(psi2Qa_ba[n]);
// 	// cjfeng 07/02/2017
// 	// The original code
// 	// mvmult_comp_trans_x(Dip2Q, psi_a[i][tau1+tau2+tau3], psi2Q, nosc, n2Q, nthreads);
// 	// mvmult_comp_trans_x(Dip2QMat[xfr3][l], psi_a[i][tau1+tau2+tau3], psi2Q, nosc, n2Q, nthreads);
// 	// for(n=0; n<n2Q; n++) cval1 += psi_cb[jk][n] * conj(psi2Q[n]);
// 	*REPH_private -= cval1 * popfac2Q;
// 
// 	return 0;
// }
// 
// int pathway_NREPH_2Q(GNCOMP *psi2Qca, GNCOMP *psi2Qb1_aa, GNCOMP *psi2Qb1_ba, GNCOMP *psi2Qb1_ab, GNCOMP *psi2Qb1_bb, int n2Q, int p, GNREAL popfac2Q, GNCOMP *NREPH_private) {
// 	// Final non-rephasing pathway:
// 	//
// 	//  | c b |
// 	//  | a b |
// 	//  | a 0 |
// 	//  | 0 0 |
// 	//
// 	// The contribution to the non-rephasing spectrum is the orientationally averaged value
// 	// ( psi_ca[k*3+i] ) * ( Dip2QMat[tau1+tau2+tau3] * psi_b1[tau1+tau2+tau3] )
// 	GNCOMP cval1 = 0.0;
// 	int n;
// 	if ( p == 0 ) for(n=0; n<n2Q; n++) cval1 += psi2Qca[n] * conj(psi2Qb1_aa[n]);
// 	else if ( p == 1 ) for(n=0; n<n2Q; n++) cval1 += psi2Qca[n] * conj(psi2Qb1_ba[n]);
// 	else if ( p == 2 ) for(n=0; n<n2Q; n++) cval1 += psi2Qca[n] * conj(psi2Qb1_ab[n]);
// 	else for(n=0; n<n2Q; n++) cval1 += psi2Qca[n] * conj(psi2Qb1_bb[n]);
// 	*NREPH_private -= cval1 * popfac2Q;
// 
// 	return 0;
// }
