#include "param_checker.h"

int check_tstep(int *tstep, int tscan) {
	if( (*tstep<=0) && (tscan>0) ) {
		printf("Please supply a (positive) trajectory time step in fs (-tstep flag).\n");
		return 3;
	}
	return 0;
}

int check_tscan(int tscan) {
	if( (tscan<0) ) {
		printf("Please supply a (non-negative) scan time step in fs (-tscan flag).\n");
		printf("(Enter 0 for a single frame in TAA).\n");
		return 3;
	}
	return 0;
}

int check_nise_tscan(int nise, int tscan) {
	if( (nise) && (tscan==0) ) {
		printf("Warning: the scan tiem is zero.\n");
		printf("There is no way to compute time correlation functions using NISE method.\n");
		printf("Please provide a (positive) scan time in fs (-tscan flag).\n");
		return 3;
	}
	return 0;
}

int check_nise_trotter(int nise, int trotter) {
	if ( (!nise) && (trotter) ) {
		printf("Trotter expansion only supports NISE method.\n");
		printf("Please specify -nise flag.\n");
		return 3;
	}
	return 0;
}

int check_T2(int T2, int nise) {
	if( (T2!=0) && (!nise) ) {
		printf("Error: Waiting time T2 must be zero for static spectral calculations.\n");
		return 3;
	}
	else if ( (T2<0) && (nise) ) {
		printf("Error: T2 should be greater or equal to zero.\n");
		return 3;
	}
	return 0;
}

int check_2d_spec_output(int *POL, int reph, int nreph, int npol, int zzzz, int zzyy, int zyyz, int zyzy,  int do2d) {
	if( ((reph+nreph)==0) && (do2d) ) {
		printf("Nothing to calculate! Please request either rephasing or non-rephasing.\n");
		return 3;
	}
	if( (npol==0) && (do2d) ) {
		printf("Nothing to calculate! Please select at least one polarization condition.\n");
		printf("Note that the default values for ZYYZ and ZYZY are false.\n");
		return 3;
	} 
	else {
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
	return 0;
}

int check_ham(char* Hamname, int *nframes, int *nosc, int *n2Q, int nbuffer, int nread, int info, int window, int nthreads) {
	int vals;
	*nframes = count_lines(Hamname);
	printf("Located %d lines in input file %s\n", *nframes, Hamname);
	vals = count_entries(Hamname);
	*nosc = floor(sqrt(vals)+0.5);
	*n2Q = (*nosc)*( (*nosc) + 1 )/2;
	printf("Located %d oscillators in input file %s\n", *nosc, Hamname);

	if(*nframes<nbuffer-nread+1) {
		printf("Error: Not enough (%d) frames for accurate averaging with requested window size (%d).\n", *nframes, window);
		printf("At least %d frames have to be provided to compute spectra for %d thread(s).\n", nbuffer-nread+1, nthreads);
		return 3;
	}
	if( (*nosc!=info) && (info!=0) ) {
		printf("Error! Info file specifies %d oscillators, but found %d in Hamiltonian file. \n", info, *nosc);
		return 3;
	}
	return 0;
}

int check_shift(int nshift, int *SHIFTNDX, int nosc) {
	int i;
	for(i=0; i<nshift; i++) {
		if(SHIFTNDX[i]>=nosc) {
			printf("Error! Requested shift index (oscillator %d) is larger than total number of oscillators.\n", SHIFTNDX[i]);
			return 3;
		}
	}
	return 0;
}

int check_sites(char* Sitesname, char* Hamname, int nosc, int nframes) {
	int vals;
	if( (strlen(Sitesname)>0) ) {
		vals = count_entries(Sitesname);
		if(nosc!=vals) {			// Checking number of oscillators
			printf("Error! Different number of oscillators (%d vs. %d) located in Hamiltonian file %s and sites file %s\n", nosc, vals, Hamname, Sitesname); 
			return 3;
		}
		vals = count_lines(Sitesname);
		if((nframes!=vals)) {	// Checking number of frames
			printf("Error! Different number of lines (%d vs. %d) in Hamiltonian file %s and sites file %s\n", nframes, vals, Hamname, Sitesname);
			return 3;
		}
	}

	return 0;
}

int check_dipoles(char Dipnames[3][maxchar], char *Hamname, int nosc, int nframes) {
	int i, vals;
	for(i=0; i<3; i++) {
		vals = count_entries(Dipnames[i]);
		if(vals!=nosc) {			// Checking number of oscillators
			printf("Error! Different number of oscillators (%d vs. %d) in Hamiltonian file %s and dipole file %s\n", vals, nosc, Hamname, Dipnames[i]);
			return 3;
		}
		vals = count_lines(Dipnames[i]);
		if(vals!=nframes) {		// Checking number of frames
			printf("Error! Different number of lines (%d vs. %d) in Hamiltonian file %s and dipole file %s\n", vals, nframes, Hamname, Dipnames[i]);
			return 3;
		}
	}
	return 0;
}
