#include "file_io.h"

/**************************
* File reading routines   *
***************************/
int read_line(FILE* fp, int nrows, int ncols, GNREAL *Mat) {
	int i, n = nrows * ncols;
	real val;
	for(i=0; i<n; i++) {
			if(fscanf(fp, "%f%*[ ,]", &val)==0) return 0;
			else Mat[i] = val;
	}
	return 1;
}

int read_info( FILE *fp, char* Infoname, char*** p_NNames, int** p_NNums, char*** p_CNames, int** p_CNums) {
	char line[1024];
	int i;
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
					// cjfeng 04/25/2017
					// Use mem_helper functions.
					if(!error) error = set_char_2d_array(&NNames, nbonds, nchar);
					if(!error) error = set_char_2d_array(&CNames, nbonds, nchar);
					if(!error) error = set_int_1d_array(&NNums, nbonds, 1);
					if(!error) error = set_int_1d_array(&CNums, nbonds, 1);
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
			// cjfeng 04/25/2017
			// Use mem_helper function
			if(NNames!=NULL) unset_char_2d_array(NNames, nbonds);
			if(CNames!=NULL) unset_char_2d_array(CNames, nbonds);
			if(NNums!=NULL) unset_int_1d_array(NNums);
			if(CNums!=NULL) unset_int_1d_array(CNums);
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

// cjfeng 03/23/2017
int read_ham(FILE *fp, GNREAL *ham, int nosc) {
	if(!read_line(fp, nosc, nosc, ham) ) {
		printf("Error reading from Hamiltonian file.\n");
		return 1;
	}
	else return 0;
}

// cjfeng 03/23/2017
int read_dip(FILE **fp, GNREAL ***Dip, int nosc, int xfr) {
	int i;
	for(i=0; i<3; i++) {
		if( !read_line(fp[i], nosc, 1, Dip[i][xfr]) ) {
			printf("Error reading from Dipole file.\n");
			return 1;
		}
	}
	return 0;
}

// cjfeng 03/23/2017
int read_sites(FILE *fp, GNREAL *sites, GNREAL **ham, int nosc, int xfr) {
	if(fp != NULL) {
		int i;
		if( !read_line(fp, nosc, 1, sites) ) {
			printf("Error reading from site energy file.\n");
			return 1;
		}
		else for(i=0; i<nosc; i++) ham[xfr][i*nosc+i] = sites[i];
	}
	return 0;
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

/***************************
* File name assignments    *
****************************/

int make_fnames( char* Infoname, char* Hamname, char* hamnm, char Dipnames[3][maxchar], char* dipxnm, char* dipynm, char* dipznm, char* Sitesname, char* sitesnm, char* axisnm, char* ftirnm, char* lognm, char polnm[4][16], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnm[10][maxchar], char* Parname, char* paramnm, char* outname, char* deffnm, int do_traj ) {
	// cjfeng 04/25/2017
	// Use subroutines to assign file names.
	// Input file names
	assign_input_files(Infoname, Hamname, hamnm, Dipnames, dipxnm, dipynm, dipznm, Sitesname, sitesnm, Parname, paramnm, deffnm);
	// Output file names
	copy_outname(axisnm, ftirnm, lognm, rephnm, nrephnm, Trajnm, outname, do_traj);
	add_suffix(axisnm, ftirnm, lognm, polnm, rephnm, nrephnm, Trajnm, outname, do_traj);

	return 0;
}

int assign_input_files( char* Infoname, char* Hamname, char* hamnm, char Dipnames[3][maxchar], char* dipxnm, char* dipynm, char* dipznm, char* Sitesname, char* sitesnm, char* Parname, char* paramnm, char* deffnm ) {
	char* hamsuffix = "ham.txt";
	char* dipxsuffix = "dipx.txt";
	char* dipysuffix = "dipy.txt";
	char* dipzsuffix = "dipz.txt";
	// Input file names
	// cjfeng 04/25/2017
	// Use subroutines to assign input file names.
	// First info file
	assign_prefix(Infoname, deffnm);
	strcat(Infoname, "info.txt");		// Adding suffix

	// Hamiltonian file
	if( (deffnm!=NULL) && (!strcmp(hamnm, hamsuffix)) ) printf("Made it inside\n");
	assign_input_name(Hamname, hamnm, deffnm, hamsuffix);
	
	// Dipole moment files	
	assign_input_name(Dipnames[0], dipxnm, deffnm, dipxsuffix);
	assign_input_name(Dipnames[1], dipynm, deffnm, dipysuffix);
	assign_input_name(Dipnames[2], dipznm, deffnm, dipzsuffix);

	// Sitename and Parname are unused if not specified
	// The default file name is never appended.
	Sitesname[0] = '\0';
	strcat(Sitesname, sitesnm);
	
	Parname[0] = '\0';
	strcat(Parname, paramnm);
	return 0;
}

int assign_input_name(char* fnm, char* new_suffix, char* prefix, char* default_suffix) {
	if(prefix!=NULL && (!strcmp(new_suffix, default_suffix)) ) {
		strcpy(fnm, prefix);		// Copying output name base to the filename.
		// Adding "_" between output file name base and suffices 
		// when name base is not referring to a folder.
		if(prefix[strlen(prefix)-1]!='/') strcat(fnm, "_");
	} 
	else fnm[0] = '\0';
	strcat(fnm, new_suffix);
	return 0; 
}

int assign_prefix(char* fnm, char* prefix) {
	if(prefix!=NULL) {
		strcpy(fnm, prefix);		// Copying output name base to the filename.
		// Adding "_" between output file name base and suffices 
		// when name base is not referring to a folder.
		if(prefix[strlen(prefix)-1]!='/') strcat(fnm, "_");
	} 
	else fnm[0] = '\0';
	return 0; 
}

int copy_outname(char* axisnm, char* ftirnm, char* lognm, char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnm[10][maxchar], char* outname, int do_traj ) {
	int i;
	if(outname!=NULL) printf("Output file name base: %s\n", outname);
	// cjfeng 04/25/2017
	// Use subroutines to assign outname for each file name.
	assign_prefix(ftirnm, outname);
	assign_prefix(lognm, outname);
	assign_prefix(axisnm, outname);
	for(i=0; i<4; i++) assign_prefix(rephnm[i], outname);
	for(i=0; i<4; i++) assign_prefix(nrephnm[i], outname);
	if(do_traj) assign_prefix(Trajnm[0], outname);
	return 0;
}

int add_suffix(char* axisnm, char* ftirnm, char* lognm, char polnm[4][16], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnm[10][maxchar], char* outname, int do_traj ) {
	int i;
	strcat(ftirnm, "ftir.txt");		// Adding suffix to ftir file
	strcat(lognm, "log.txt");		// Adding suffix to log file
	strcat(axisnm, "waxis.txt");		// Adding suffix to frequency axis file

	for(i=0; i<4; i++) {
		strcat(rephnm[i], "reph_");	// Adding suffix to rephasing spectra files
		strcat(rephnm[i], polnm[i]);
	}
	for(i=0; i<4; i++) {
		strcat(nrephnm[i], "nreph_");	// Adding suffix to non-rephasing spectra file
		strcat(nrephnm[i], polnm[i]);
	}
	if(do_traj) {
		// tstamp file
		strcat(Trajnm[0], "tstamp_traj.txt");
		// ftir trajectory file
		strncpy(Trajnm[1], ftirnm, strlen(ftirnm)-4);
		strcat(Trajnm[1], "_traj.txt");
		for(i=2; i<6; i++) {
			strncpy(Trajnm[i], rephnm[i-2], strlen(rephnm[i-2])-4);
			strcat(Trajnm[i], "_traj.txt");
		}
		for(i=6; i<10; i++) {
			strncpy(Trajnm[i], nrephnm[i-6], strlen(nrephnm[i-6])-4);
			strcat(Trajnm[i], "_traj.txt");
		}
	}
	return 0;
}

int file_opener(FILE **fp, char *fnm, char* ftype, int status) {
	int error = status;
	if(!error && strlen(fnm)>0 ) {
		*fp = fopen(fnm, ftype);
		if(*fp == NULL) error = 2;
	} else *fp = NULL;

	return error;
};

int open_all( char* Hamname, char Dipnames[3][maxchar], char* Sitesname, char* axisnm, char* ftirnm, char* lognm, int npol, int POL[4], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnames[10][maxchar], int do2d, int do_traj) {
	int error = 0;
	char* rtype = "r";
	char* wtype = "w";
	int i,j;
	// cjfeng 04/25/2017
	// Use subroutine to open files.
	error = file_opener(&hfp, Hamname, rtype, error);
	for(i=0; i<3; i++) error = file_opener(&Dfp[i], Dipnames[i], rtype, error);

	error = file_opener(&sfp, Sitesname, rtype, error); // Site energy file
	error = file_opener(&afp, axisnm, wtype, error);	// frequency axis file
	error = file_opener(&ffp, ftirnm, wtype, error);	// FTIR file
	error = file_opener(&lfp, lognm, wtype, error);	// log file
	
	if(do_traj) {
		if(!do2d) for(i=2; i<10; i++) Trajnames[i][0] = '\0';
		// Delete file names that don't get used
		else {
			for(i=0; i<4; i++) {
				int use_pol = 0;
				for(j=0; j<npol; j++) if(POL[j]==i) use_pol = 1;
				if(use_pol==0) {
					Trajnames[2+i][0] = '\0';
					Trajnames[6+i][0] = '\0';
				}
			}
		}
		for(i=0; i<10; i++) error = file_opener(&Trajfp[i], Trajnames[i], wtype, error); 
	}
	else for(i=0; i<10; i++) Trajfp[i] = NULL;

	if(do2d) {
		for(i=0; i<npol; i++) {
			int ni = POL[i];
			error = file_opener(&rfp[ni], rephnm[ni], wtype, error); // rephasing spectra
		}
		for(i=0; i<npol; i++) {
			int ni = POL[i];
			error = file_opener(&nrfp[ni], nrephnm[ni], wtype, error); // non-rephasing spectra
		}
	}	
	else {
		for(i=0; i<npol; i++) rfp[POL[i]] = NULL; 
		for(i=0; i<npol; i++) nrfp[POL[i]] = NULL;
	}

	if(!error) printf("Opened all files.\n");
	return error;
}

// cjfeng 03/23/2017
int print_frame(int frame, int frame_period, FILE *lfp) {
	if( !(frame % frame_period) ) {
		printf("Frame: %d\n", frame);
		fprintf(lfp, "Frame: %d\n", frame);
		fflush(lfp);
	}
	return 0;
}

// cjfeng 03/23/2017
int print_memory(int frame, int frame_period, struct rusage *r_usage, FILE *lfp) {
	if( !(frame % frame_period) ) {
		getrusage(RUSAGE_SELF, r_usage);
		printf("Memory usage = %ld kB.\n", r_usage->ru_maxrss);
		fprintf(lfp, "Memory usage = %ld kB.\n", r_usage->ru_maxrss);
	}
	return 0;
}

// cjfeng 03/24/2017
int print_nise_waxis(FILE *fp, GNREAL dw, int nprint, int ndxstart) {
	int i;
	for(i=0; i<nprint; i++) fprintf(fp, "%6.10f\n", (ndxstart+i) * dw);
	return 0;
}

int print_waxis(FILE *fp, real wres, int npts, real wstart) {
	int i;
	for(i=0; i<npts; i++) fprintf(fp, "%6.10f\n", wstart + i * wres);
	return 0;
}

// cjfeng 03/24/2017
int print_nise_ftir(FILE *fp, fftw_complex *FTout1D, int nprint, int ndxstart, int winzpad) {
	int i;
	// The FTIR spectrum is now stored in the real part, FTout1D[.][0]
	for(i=0; i<nprint; i++) fprintf(fp, "%6.10f\n", FTout1D[(ndxstart+i)%winzpad][0]);
	return 0;
}

int print_nise_ftir_traj(FILE *fp, fftw_complex *FTout1D, int nprint, int ndxstart, int winzpad) {
	int i;
	// The FTIR spectrum is now stored in the real part, FTout1D[.][0]
	for(i=0; i<nprint; i++) fprintf(fp, "%6.10f\t", FTout1D[(ndxstart+i)%winzpad][0]);
	fprintf(fp, "\n");
	return 0;
}

int print_ftir(FILE *fp, GNREAL *ftir, int npts) {
	int i;
	for(i=0; i<npts; i++) fprintf(fp, "%6.10f\n", ftir[i]);
	return 0;
}

int print_ftir_traj(FILE *fp, GNREAL *ftir, int npts) {
	int i;
	for(i=0; i<npts; i++) fprintf(fp, "%6.10f\t", ftir[i]);
	fprintf(fp, "\n");
	return 0;
}

int print_nise_2d(FILE *fp, fftw_complex *FTout2D, int nprint, int ndxstart, int winzpad) {
	int i, j, ndx1, ndx3;
	GNCOMP cval;
	for(i=0; i<nprint; i++) {
		ndx1 = ((ndxstart+i) % winzpad) * winzpad;
		for(j=0; j<nprint; j++) {
			ndx3 = (ndxstart+j) % winzpad;
			cval = FTout2D[ndx1 + ndx3][0] + FTout2D[ndx1 + ndx3][1]*I;
			if(cimag(cval) < 0) fprintf(fp, "%6.10f%6.10fi\t", creal(cval), cimag(cval));
			else fprintf(fp, "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
		}
		fprintf(fp, "\n");
	}
	return 0;
}

int print_nise_2d_traj(FILE *fp, fftw_complex *FTout2D, int nprint, int ndxstart, int winzpad) {
	int i, j, ndx1, ndx3;
	GNCOMP cval;
	for(i=0; i<nprint; i++) {
		ndx1 = ((ndxstart+i) % winzpad) * winzpad;
		for(j=0; j<nprint; j++) {
			ndx3 = (ndxstart+j) % winzpad;
			cval = FTout2D[ndx1 + ndx3][0] + FTout2D[ndx1 + ndx3][1]*I;
			if(cimag(cval) < 0) fprintf(fp, "%6.10f%6.10fi\t", creal(cval), cimag(cval));
			else fprintf(fp, "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
		}
	}
	fprintf(fp, "\n");
	return 0;
}

int print_2d(FILE *fp, GNCOMP **spec, int npts, int *POL, int p) {
	int j, k;
	GNCOMP cval;
	for(j=0; j<npts; j++) {
		int jj = j * npts;
		for(k=0; k<npts; k++) {
			cval = spec[POL[p]][jj+k];
			if(cimag(cval) < 0) fprintf(fp, "%6.10f%6.10fi\t", creal(cval), cimag(cval));
			else fprintf(fp, "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
		}
		fprintf(fp, "\n");
	}
	return 0;
}

int print_2d_traj(FILE *fp, GNCOMP **spec, int npts, int *POL, int p) {
	int j, k;
	GNCOMP cval;
	for(j=0; j<npts; j++) {
		int jj = j * npts;
		for(k=0; k<npts; k++) {
			cval = spec[POL[p]][jj+k];
			if(cimag(cval) < 0) fprintf(fp, "%6.10f%6.10fi\t", creal(cval), cimag(cval));
			else fprintf(fp, "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
		}
	}
	fprintf(fp, "\n");
	return 0;
}
