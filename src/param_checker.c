#include "param_checker.h"

int check_tstep(int tstep) {
  if( tstep<=0 ) {
    printf("Please supply a positive trajectory time step in fs (-tstep flag).\n");
    return 3;
  }
  return 0;
}

// cjfeng 12/07/2018
// Merged tscan subroutines into a single subroutine
int check_tscan(int tscan, int nise) {
  if( (tscan<0) ) {
    printf("Please supply a (non-negative) scan time step in fs (-tscan flag).\n");
    printf("(Enter 0 for a single frame in TAA).\n");
    return 3;
  }
  else if( (tscan==0) && (nise) ) {
    printf("Warning: the scan tiem is zero.\n");
    printf("There is no way to compute time correlation functions using NISE method.\n");
    printf("Please provide a (positive) scan time in fs (-tscan flag).\n");
    return 3;
  }
  return 0;
}

int check_nise_flags(int nise, int trotter, int pert) {
  if ( (!nise) && (trotter) ) {
    printf("Trotter expansion only supports NISE method.\n");
    printf("Please specify -nise flag.\n");
    return 3;
  }
  else if ( trotter && pert ) {
    printf("Trotter expansion and perturbative approximation are not compatible to each other.\n");
    printf("Please provide either -trotter or -pert.\n");
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

// cjfeng 12/04/2018
// Using pol_info
int check_2d_spec_output(POL_info *pol_info, int reph, int nreph, int do2d) {
  if (do2d) {
    if( (reph+nreph)==0 ) {
      printf("Error: Nothing to calculate! Please request either rephasing or non-rephasing.\n");
      return 3;
    }
    if( pol_info->npol==0 ) {
      printf("Error: Nothing to calculate! Please select at least one polarization condition.\n");
      printf("Note that the default values for ZYYZ and ZYZY are false.\n");
      return 3;
    }
  } 
  return 0;
}

// cjfeng 01/03/2019
// Putting if2Qfiles into spec_param
// cjfeng 01/04/2019
// Added traj_param
int check_input_files(fin *fin, spec_param *spec_param, w_param *w_param, shift_info *shift_info, int *nframes, traj_param *traj_param, int nbonds) {
  int error = 0;
  int i;
  int if2Qfiles = spec_param->if2Qfiles;
  // Check validity of Hamiltonian file.
  error = check_ham(fin->Ham1Q.fnm, nframes, traj_param, nbonds, w_param->window, spec_param->nthreads);
  if(error) return error;
  // Check isotope shift file
  error = check_shift(shift_info, traj_param->nosc);
  if(error) return error;
  // Checking the consistency between Hamiltonian file and site energy file.
  error = check_sites(fin->Sites.fnm, fin->Ham1Q.fnm, traj_param->nosc, *nframes);
  if(error) return error;
  // Checking the consistency between Hamiltonian file and dipole moment files.
  for(i=0; i<3; i++) {
    error = check_dipole(fin->Dip1Q[i].fnm, fin->Ham1Q.fnm, traj_param->nosc, *nframes);
    if(error) return error;
  }
  // 2Q parts
  if(if2Qfiles) {
    error = check_ham2Q(fin->Ham2Q.fnm, nframes, traj_param, nbonds, w_param->win2d, spec_param->nthreads);
    if(error) return error;
    for(i=0; i<3; i++) {
      error = check_dip2Q(fin->Dip2Q[i].fnm, fin->Ham2Q.fnm, traj_param->n2Q, *nframes);
      if(error) return error;
    }
  }

  return error;
} 

// cjfeng 01/04/2019
// Added traj_param
int check_ham(char* Hamname, int *nframes, traj_param *traj_param, int info, int window, int nthreads) {
  int vals;
  int nbuffer = traj_param->nbuffer;
  int nread = traj_param->nread;
  *nframes = count_lines(Hamname);
  printf("Located %d lines in input file %s\n", *nframes, Hamname);
  vals = count_entries(Hamname);
  traj_param->nosc = floor(sqrt(vals)+0.5);
  traj_param->n2Q = (traj_param->nosc)*( (traj_param->nosc) + 1 )/2;
  printf("Located %d oscillators in input file %s\n", traj_param->nosc, Hamname);

  if(*nframes<nbuffer-nread+1) {
    printf("Error: Not enough (%d) frames for accurate averaging with requested window size (%d).\n", *nframes, window);
    printf("At least %d frames have to be provided to compute spectra for %d thread(s).\n", nbuffer, nthreads);
    return 3;
  }
  if( (traj_param->nosc!=info) && (info!=0) ) {
    printf("Error! Info file specifies %d oscillators, but found %d in Hamiltonian file. \n", info, traj_param->nosc);
    return 3;
  }
  return 0;
}

// cjfeng 01/04/2019
// Added traj_param
int check_ham2Q(char* Ham2Qname, int *nframes, traj_param *traj_param, int info, int win2d, int nthreads) {
  if( (strlen(Ham2Qname)>0) ) {
    int vals, n2Q_tmp;
    int nbuffer = traj_param->nbuffer;
    int nread = traj_param->nread;
    *nframes = count_lines(Ham2Qname);
    printf("Located %d lines in input file %s\n", *nframes, Ham2Qname);
    vals = count_entries(Ham2Qname);
    traj_param->n2Q = floor(sqrt(vals)+0.5);
    n2Q_tmp = (traj_param->nosc)*( (traj_param->nosc) + 1 )/2;
    printf("Located %d 2Q oscillators in input file %s\n", traj_param->n2Q, Ham2Qname);
    if ( traj_param->n2Q != n2Q_tmp ) {
      printf("n2Q (%d) from the input file %s is not equal to n2Q (%d) computed from the 1Q Ham.\n", traj_param->n2Q, Ham2Qname, n2Q_tmp);
      return 3;
    }

    if(*nframes<nbuffer-nread+1) {
      printf("Error: Not enough (%d) frames for accurate averaging with requested window size (%d).\n", *nframes, win2d);
      printf("At least %d frames have to be provided to compute spectra for %d thread(s).\n", nbuffer, nthreads);
      return 3;
    }
  }
  return 0;
}

// cjfeng 12/09/2018
// Using shift_info
int check_shift(shift_info *shift_info, int nosc) {
  int i;
  int nshift = shift_info->nshift;
  for(i=0; i<nshift; i++) {
    if(shift_info->SHIFTNDX[i]>=nosc) {
      printf("Error! Requested shift index (oscillator %d) is larger than total number of oscillators.\n", shift_info->SHIFTNDX[i]);
      return 3;
    }
  }
  return 0;
}

int check_sites(char* Sitesname, char* Hamname, int nosc, int nframes) {
  int vals;
  if( (strlen(Sitesname)>0) ) {
    vals = count_entries(Sitesname);
    if(nosc!=vals) {      // Checking number of oscillators
      printf("Error! Different number of oscillators (%d vs. %d) located in Hamiltonian file %s and sites file %s\n", nosc, vals, Hamname, Sitesname); 
      return 3;
    }
    vals = count_lines(Sitesname);
    if((nframes!=vals)) {  // Checking number of frames
      printf("Error! Different number of lines (%d vs. %d) in Hamiltonian file %s and sites file %s\n", nframes, vals, Hamname, Sitesname);
      return 3;
    }
  }

  return 0;
}

// cjfeng 12/08/2018
// Using fin
int check_dipole(char *Dipnames, char *Hamname, int nosc, int nframes) {
  int vals;
  vals = count_entries(Dipnames);
  if(vals!=nosc) {      // Checking number of oscillators
    printf("Error! Different number of oscillators (%d vs. %d) in Hamiltonian file %s and dipole file %s\n", vals, nosc, Hamname, Dipnames);
    return 3;
  }
  vals = count_lines(Dipnames);
  if(vals!=nframes) {    // Checking number of frames
    printf("Error! Different number of lines (%d vs. %d) in Hamiltonian file %s and dipole file %s\n", vals, nframes, Hamname, Dipnames);
    return 3;
  }
  return 0;
}

int check_dip2Q(char *Dip2Qnames, char *Ham2Qname, int n2Q, int nframes) {
  int vals;
  if (strlen(Dip2Qnames) == 0 ) {
    printf("Error! File names of 2Q dipole moments are not provided completely.\n");
    return 3;
  }
  vals = count_entries(Dip2Qnames);
  if(vals!=n2Q) {      // Checking number of oscillators
    printf("Error! Different number of 2Q oscillators (%d vs. %d) in Ham2Q file %s and dip2Q file %s\n", vals, n2Q, Ham2Qname, Dip2Qnames);
    return 3;
  }
  vals = count_lines(Dip2Qnames);
  if(vals!=nframes) {    // Checking number of frames
    printf("Error! Different number of lines (%d vs. %d) in Ham2Q file %s and dip2Q file %s\n", vals, nframes, Ham2Qname, Dip2Qnames);
    return 3;
  }
  
  return 0;
}

int check_nbonds(int nbonds) {
  int error = 0;
  if(nbonds>0) printf("Successfully read info file. Will be looking for %d oscillators in input files. \n", nbonds);
  else if (nbonds<0) error = 3;
  return error;
}

// cjfeng 12/07/2018
// Using a unified subroutine to check parameters instead of putting those into the main code.
int check_spec_param(spec_param *spec_param, POL_info *pol_info) {
  int error = 0;
  int tstep = spec_param->tstep;
  int tscan = spec_param->tscan;
  int nise = spec_param->nise;
  int trotter = spec_param->trotter;
  int pert = spec_param->pert;
  int T2 = spec_param->T2;
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;
  int do2d = spec_param->do2d;

  error = check_tstep(tstep);
  if(error) return error;

  error = check_tscan(tscan, nise);
  if(error) return error;
  
  error = check_nise_flags(nise, trotter, pert);
  if(error) return error;
  
  error = check_T2(T2, nise);
  if(error) return error;

  error = check_2d_spec_output(pol_info, reph, nreph, do2d);
  if(error) return error;

  return error;
}

int check_w_param(w_param *w_param) {
  int winzpad = w_param->winzpad;
  int window = w_param->window;
  if (winzpad < window) {
    printf("Error: the provided winterp gives shorter spectra array (%d) than requested window size (%d)\n", winzpad, window);
    printf("Please give smaller winterp or smaller tscan.\n");
    return 3;
  }
  return 0;
}

