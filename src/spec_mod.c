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

int gen_avg_1Q_ham(GNREAL **ham, GNREAL *avg_ham, GNREAL *hann, int nosc, int fr, int whann, int window, int nbuffer) {
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
