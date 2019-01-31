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

// cjfeng 12/09/2018
// Using ref_info, and fin
int read_info( fin *fin, res_info *res_info, int maxchar) {
  char line[maxchar];
  int i;
  int nbonds = 0;
  int nchar = 20;
  int error = 0; 
  fin->Info.fp = fopen(fin->Info.fnm, "r");
  if(fin->Info.fp==NULL) return 0;
  else {
    while( (!error) && (fgets(line, maxchar, fin->Info.fp)!=NULL)) {
      if(strncmp(line, "BONDS:", 6)==0) {
        if(sscanf(line, "%*s %d", &nbonds)!=1) {
          printf("Error reading number of bonds from info file. \n");
          error = 1;
        } else set_res_info(res_info, nbonds, nchar);

        for(i=0; i<nbonds; i++) {
          if( (fgets(line, maxchar, fin->Info.fp)==NULL) || (sscanf(line, "%s %d %s %d", res_info->NNames[i], &(res_info->NNums[i]), res_info->CNames[i], &(res_info->CNums[i]) )!=4) ) {
            error = 1;
            break;
          }
        }
      }
    }
    if(error) {
      printf("Error reading from info file.\n");
      unset_res_info(res_info);
      if(fin->Info.fp!=NULL) fclose(fin->Info.fp);
      return -1;
    }
  }
  if(!error) {
    if(nbonds>0) printf("Bond info: \n");
    for(i=0; i<nbonds; i++) printf("%s %d--%s %d\n", res_info->NNames[i], res_info->NNums[i], res_info->CNames[i], res_info->CNums[i]);
  }
  if(fin->Info.fp!=NULL) fclose(fin->Info.fp);
  return nbonds; 
}

// cjfeng 03/23/2017
int read_ham(FILE *fp, GNREAL *ham, int nosc) {
  if(!read_line(fp, nosc, nosc, ham) ) {
    printf("Error reading from Hamiltonian file.\n");
    return 1;
  }
  return 0;
}

// cjfeng 12/08/2018
// Using fin
// cjfeng 01/10/2019
// Unifying the read_dip function
int read_dip(f_pointer *fp, GNREAL **Dip, int nosc) {
  int i;
  for(i=0; i<3; i++) {
    if( !read_line(fp[i].fp, nosc, 1, Dip[i]) ) {
      printf("Error reading from Dipole file.\n");
      return 1;
    }
  }
  return 0;
}

// cjfeng 03/23/2017
int read_sites(FILE *fp, GNREAL *ham, int nosc) {
  int error = 0;
  GNREAL *sites;
  error = set_gnreal_1d_array(&sites, 1, nosc);
  if(error) return error;
  if(fp != NULL && sites != NULL ) {
    int i;
    if( !read_line(fp, nosc, 1, sites) ) {
      printf("Error reading from site energy file.\n");
      return 1;
    }
    else for(i=0; i<nosc; i++) ham[i*nosc+i] = sites[i];
  }
  unset_gnreal_1d_array(sites);
  return 0;
}

// cjfeng 12/08/2018
// Read param file using fnm_param
int read_param(fnm_param *fnm_param, char* Parname, res_info *res_info, shift_info *shift_info, int maxchar) {
  int i;
  int error = 0;
  if(strlen(fnm_param->paramnm)>0) {
    printf("Reading parameters from file %s.\n", Parname);
    if(!parse_shifts(Parname, shift_info, maxchar)) {
      printf("Error parsing parameter file %s. Please check input.\n", Parname);
      error = 3;
      return error;
    }
    int nshift = shift_info->nshift;
    printf("Loaded %d oscillator shifts: \n", nshift);
    if(res_info->nbonds!=0) {
      for(i=0; i<nshift; i++) printf("%s %d (index %d) :\t%6.2f cm-1.\n", res_info->NNames[shift_info->SHIFTNDX[i]], res_info->NNums[shift_info->SHIFTNDX[i]], shift_info->SHIFTNDX[i], shift_info->SHIFT[i]);
    } 
    else for(i=0; i<nshift; i++) printf("Oscillator %d: \t%6.2f cm-1.\n", shift_info->SHIFTNDX[i], shift_info->SHIFT[i]);
  } 
  else shift_info->nshift = 0;
  return error;
}

// cjfeng 12/09/2018
// Using res_info
int parse_shifts(char* Parname, shift_info *shift_info, int maxchar) {
  FILE* fp = fopen(Parname, "r");
  char line[maxchar];
  int i; 
  int error = 0;
  if(fp==NULL) error = 1;
  while(!error && (fgets(line, maxchar, fp)!=NULL)) {
    if(strncmp(line, "SHIFT", 5)==0) {
      if(sscanf(line, "%*s %d", &(shift_info->nshift) )!=1) error = 1;
      else {
        error = set_shift_info(shift_info);
        if (error) return 1;
        for(i=0; i<shift_info->nshift; i++) {
          if(fgets(line, maxchar, fp)==NULL) error = 1;
          else {
            if(sscanf(line, "%d %f", &(shift_info->SHIFTNDX[i]), &(shift_info->SHIFT[i]))!=2) error = 1;
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

// cjfeng 05/07/2018
// Added 2Q parts
// cjfeng 12/08/2018
// Using fnm_param, fin, and fout to replace individual file names
int make_fnames(fin *fin, fout *fout, fnm_param *fnm_param, char* Parname, int do_traj ) {
  printf("Making file names.\n");
  // Input file names
  assign_input_files(fin, fnm_param, Parname);
  // Output file names
  copy_outname(fout, fnm_param->outname, do_traj);
  printf("Finished copy_outname\n");
  add_suffix(fout, fnm_param->polnm, do_traj);
  printf("Finished making file names.\n");

  return 0;
}

// cjfeng 05/07/2018
// Added 2Q parts
// cjfeng 12/08/2018
// Using fnm_param, and fin
int assign_input_files( fin *fin, fnm_param *fnm_param, char* Parname ) {
  char* hamsuffix = "ham.txt";
  char* dipxsuffix = "dipx.txt";
  char* dipysuffix = "dipy.txt";
  char* dipzsuffix = "dipz.txt";
  // Input file names
  // First info file
  assign_prefix(fin->Info.fnm, fnm_param->deffnm);
  strcat(fin->Info.fnm, "info.txt"); // Adding suffix

  // Hamiltonian file
  if( (fnm_param->deffnm!=NULL) && (!strcmp(fnm_param->hamnm, hamsuffix)) ) printf("Made it inside\n");
  assign_input_name(fin->Ham1Q.fnm, fnm_param->hamnm, fnm_param->deffnm, hamsuffix);
  
  // Dipole moment files  
  assign_input_name(fin->Dip1Q[0].fnm, fnm_param->dipxnm, fnm_param->deffnm, dipxsuffix);
  assign_input_name(fin->Dip1Q[1].fnm, fnm_param->dipynm, fnm_param->deffnm, dipysuffix);
  assign_input_name(fin->Dip1Q[2].fnm, fnm_param->dipznm, fnm_param->deffnm, dipzsuffix);

  // Sitename and Parname are unused if not specified
  // The default file name is never appended.
  fin->Sites.fnm[0] = '\0';
  strcat(fin->Sites.fnm, fnm_param->sitesnm);
  
  Parname[0] = '\0';
  strcat(Parname, fnm_param->paramnm);

  // cjfeng 05/07/2018
  // 2Q parts, Ham2Qname and Dip2Qnames are unused if not specified.
  // The default file name is never appended.
  fin->Ham2Q.fnm[0] = '\0';
  strcat(fin->Ham2Q.fnm, fnm_param->ham2Qnm);
  
  // Dipole moment files  
  fin->Dip2Q[0].fnm[0] = '\0';
  fin->Dip2Q[1].fnm[0] = '\0';
  fin->Dip2Q[2].fnm[0] = '\0';
  strcat(fin->Dip2Q[0].fnm, fnm_param->dipx2Qnm);
  strcat(fin->Dip2Q[1].fnm, fnm_param->dipy2Qnm);
  strcat(fin->Dip2Q[2].fnm, fnm_param->dipz2Qnm);

  return 0;
}

int assign_input_name(char* fnm, char* new_suffix, char* prefix, char* default_suffix) {
  if(prefix!=NULL && (!strcmp(new_suffix, default_suffix)) ) {
    strcpy(fnm, prefix);    // Copying output name base to the filename.
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
    strcpy(fnm, prefix);    // Copying output name base to the filename.
    // Adding "_" between output file name base and suffices 
    // when name base is not referring to a folder.
    if(prefix[strlen(prefix)-1]!='/') strcat(fnm, "_");
  } 
  else fnm[0] = '\0';
  return 0; 
}

// cjfeng 12/09/2018
// Using fout
int copy_outname(fout *fout, char* outname, int do_traj ) {
  int i;
  if(outname!=NULL) printf("Output file name base: %s\n", outname);
  // cjfeng 04/25/2017
  // Use subroutines to assign outname for each file name.
  // cjfeng 12/09/2018
  assign_prefix(fout->waxis.fnm, outname);
  assign_prefix(fout->labs.fnm, outname);
  assign_prefix(fout->log.fnm, outname);
  for(i=0; i<4; i++) assign_prefix(fout->reph[i].fnm, outname);
  for(i=0; i<4; i++) assign_prefix(fout->nreph[i].fnm, outname);
  if(do_traj) for(i=0; i<10; i++) assign_prefix(fout->traj[i].fnm, outname);
  return 0;
}

// cjfeng 12/09/2018
// Using fout
int add_suffix(fout *fout, char *polnm[4], int do_traj ) {
  int i;
  strcat(fout->waxis.fnm, "waxis.txt");    // Adding suffix to frequency axis file
  strcat(fout->labs.fnm, "ftir.txt");    // Adding suffix to ftir file
  strcat(fout->log.fnm, "log.txt");    // Adding suffix to log file

  for(i=0; i<4; i++) {
    strcat(fout->reph[i].fnm, "reph_");  // Adding suffix to rephasing spectra files
    strcat(fout->reph[i].fnm, polnm[i]);
  }
  for(i=0; i<4; i++) {
    strcat(fout->nreph[i].fnm, "nreph_");  // Adding suffix to non-rephasing spectra file
    strcat(fout->nreph[i].fnm, polnm[i]);
  }
  if(do_traj) {
    // tstamp file
    strcat(fout->traj[0].fnm, "tstamp_traj.txt");
    // ftir trajectory file
    strcat(fout->traj[1].fnm, "ftir_traj.txt");
    for(i=2; i<6; i++) {
      strncpy(fout->traj[i].fnm, fout->reph[i-2].fnm, strlen(fout->reph[i-2].fnm)-4);
      strcat(fout->traj[i].fnm, "_traj.txt");
    }
    for(i=6; i<10; i++) {
      strncpy(fout->traj[i].fnm, fout->nreph[i-6].fnm, strlen(fout->nreph[i-6].fnm)-4);
      strcat(fout->traj[i].fnm, "_traj.txt");
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

// cjfeng 05/07/2018
// Added 2Q parts
// cjfeng 12/09/2018
// Using POL_info, spec_param, fin, and fout
int open_all(fin *fin, fout *fout, POL_info *pol_info, spec_param *spec_param) {
  int error = 0;
  char* rtype = "r";
  char* wtype = "w";
  int i,j;
  int npol = pol_info->npol;
  int POL[4];
  for(i=0; i<4; i++) POL[i] = pol_info->POL[i];

  int do2d = spec_param->do2d;
  int do_traj = spec_param->do_traj;
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;

  // cjfeng 04/25/2017
  // Use subroutine to open files.
  error = file_opener(&(fin->Ham1Q.fp), (fin->Ham1Q.fnm), rtype, error);
  for(i=0; i<3; i++) error = file_opener(&(fin->Dip1Q[i].fp), (fin->Dip1Q[i].fnm), rtype, error);

  error = file_opener(&(fin->Ham2Q.fp), (fin->Ham2Q.fnm), rtype, error);
  for(i=0; i<3; i++) error = file_opener(&(fin->Dip2Q[i].fp), (fin->Dip2Q[i].fnm), rtype, error);

  error = file_opener(&(fin->Sites.fp), (fin->Sites.fnm), rtype, error); // Site energy file

  error = file_opener(&(fout->waxis.fp), fout->waxis.fnm, wtype, error);  // frequency axis file
  error = file_opener(&(fout->labs.fp), fout->labs.fnm, wtype, error);  // FTIR file
  error = file_opener(&(fout->log.fp), fout->log.fnm, wtype, error);  // log file

  if(do2d) {
    if (reph) {
      for(i=0; i<npol; i++) {
        int ni = POL[i];
        error = file_opener(&(fout->reph[ni].fp), fout->reph[ni].fnm, wtype, error); // rephasing spectra
      }
    }
    else for(i=0; i<npol; i++) fout->reph[POL[i]].fp = NULL;
    if (nreph) {
      for(i=0; i<npol; i++) {
        int ni = POL[i];
        error = file_opener(&(fout->nreph[ni].fp), fout->nreph[ni].fnm, wtype, error); // non-rephasing spectra
      }
    } else for(i=0; i<npol; i++) fout->nreph[POL[i]].fp = NULL;
  }  
  else {
    for(i=0; i<npol; i++) fout->reph[POL[i]].fp = NULL; 
    for(i=0; i<npol; i++) fout->nreph[POL[i]].fp = NULL;
  }
  
  if(do_traj) {
    if(!do2d) for(i=2; i<10; i++) fout->traj[i].fnm[0] = '\0';
    // Delete file names that don't get used
    else {
      for(i=0; i<4; i++) {
        int use_pol = 0;
        for(j=0; j<npol; j++) if(POL[j]==i) use_pol = 1;
        if(use_pol==0 || !reph ) {
          fout->traj[2+i].fnm[0] = '\0';
        }
        if(use_pol==0 || !nreph ) {
          fout->traj[6+i].fnm[0] = '\0';
        }
      }
    }
    for(i=0; i<10; i++) error = file_opener(&(fout->traj[i].fp), fout->traj[i].fnm, wtype, error); 
  }
  else for(i=0; i<10; i++) fout->traj[i].fnm = NULL;


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
    fflush(lfp);
  }
  return 0;
}

// cjfeng 12/10/2018
int print_elapsed_time(double t, FILE *lfp) {
  printf("Time: %2.4f seconds \n", t);
  fprintf(lfp, "Time: %2.4f seconds \n", t);
  fflush(lfp);
  return 0;
}

void print_gmx_quote(FILE *lfp) {
  thanx(stdout);
  thanx(lfp);
  fflush(lfp);
}

// cjfeng 12/10/2018
int print_tstamp(int frame, int tstep, FILE *fp) {
  printf("Dumping at frame %d.\n", frame);
  GNREAL tstamp = ((GNREAL) tstep * (GNREAL) frame)/1000;
  fprintf(fp, "%6.10f\n", tstamp);
  fflush(fp);
  return 0;
}

// cjfeng 12/10/2018
// Using w_param
int print_nise_waxis(FILE *fp, w_param *w_param) {
  int i;
  GNREAL dw = w_param->dw;
  int nprint = w_param->nprint;
  int ndxstart = w_param->ndxstart;
  for(i=0; i<nprint; i++) fprintf(fp, "%6.10f\n", (ndxstart+i) * dw);
  fflush(fp);
  return 0;
}

// cjfeng 12/10/2018
// Using w_param
int print_waxis(FILE *fp, w_param *w_param) {
  int i;
  GNREAL wres = w_param->wres;
  int npts = w_param->npts;
  GNREAL wstart = w_param->wstart;
  for(i=0; i<npts; i++) fprintf(fp, "%6.10f\n", wstart + i * wres);
  fflush(fp);
  return 0;
}

// cjfeng 12/05/2018
// New subroutines using fftw_array
int print_nise_ftir(FILE *fp, fftw_array *FT1D, w_param *w_param) {
  int i;
  int nprint = w_param->nprint;
  int ndxstart = w_param->ndxstart;
  int nzpad = FT1D->nzpad;
  // The FTIR spectrum is now stored in the real part, FT1D->FTout[.][0]
  for(i=0; i<nprint; i++) fprintf(fp, "%6.10f\n", FT1D->FTout[(ndxstart+i)%nzpad][0]);
  fflush(fp);
  return 0;
}

int print_nise_ftir_traj(FILE *fp, fftw_array *FT1D, w_param *w_param) {
  int i;
  int nprint = w_param->nprint;
  int ndxstart = w_param->ndxstart;
  int nzpad = FT1D->nzpad;
  // The FTIR spectrum is now stored in the real part, FT1D->FTout1D[.][0]
  for(i=0; i<nprint; i++) fprintf(fp, "%6.10f\t", FT1D->FTout[(ndxstart+i)%nzpad][0]);
  fprintf(fp, "\n");
  fflush(fp);
  return 0;
}

int print_nise_2d(FILE *fp, fftw_array *FT2D, w_param *w_param) {
  int i, j, ndx1, ndx3;
  int nprint = w_param->nprint;
  int ndxstart = w_param->ndxstart;
  int nzpad = FT2D->nzpad;
  GNCOMP cval;
  for(i=0; i<nprint; i++) {
    ndx1 = ((ndxstart+i) % nzpad) * nzpad;
    for(j=0; j<nprint; j++) {
      ndx3 = (ndxstart+j) % nzpad;
      cval = FT2D->FTout[ndx1 + ndx3][0] + FT2D->FTout[ndx1 + ndx3][1]*I;
      if(cimag(cval) < 0) fprintf(fp, "%6.10f%6.10fi\t", creal(cval), cimag(cval));
      else fprintf(fp, "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
    }
    fprintf(fp, "\n");
  }
  fflush(fp);
  return 0;
}

int print_nise_2d_traj(FILE *fp, fftw_array *FT2D, w_param *w_param) {
  int i, j, ndx1, ndx3;
  int nzpad = FT2D->nzpad;
  int nprint = w_param->nprint;
  int ndxstart = w_param->ndxstart;
  GNCOMP cval;
  for(i=0; i<nprint; i++) {
    ndx1 = ((ndxstart+i) % nzpad) * nzpad;
    for(j=0; j<nprint; j++) {
      ndx3 = (ndxstart+j) % nzpad;
      cval = FT2D->FTout[ndx1 + ndx3][0] + FT2D->FTout[ndx1 + ndx3][1]*I;
      if(cimag(cval) < 0) fprintf(fp, "%6.10f%6.10fi\t", creal(cval), cimag(cval));
      else fprintf(fp, "%6.10f+%6.10fi\t", creal(cval), cimag(cval));
    }
  }
  fprintf(fp, "\n");
  fflush(fp);
  return 0;
}

int print_ftir(FILE *fp, GNREAL *ftir, int npts) {
  int i;
  for(i=0; i<npts; i++) fprintf(fp, "%6.10f\n", ftir[i]);
  fflush(fp);
  return 0;
}

int print_ftir_traj(FILE *fp, GNREAL *ftir, int npts) {
  int i;
  for(i=0; i<npts; i++) fprintf(fp, "%6.10f\t", ftir[i]);
  fprintf(fp, "\n");
  fflush(fp);
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
  fflush(fp);
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
  fflush(fp);
  return 0;
}
