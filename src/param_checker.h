#ifndef PARAM_CHECKER_H
#define PARAM_CHECKER_H

#ifndef DEF_H
#include "def.h"
#endif

int check_tstep(int *tstep, int tscan);
int check_tscan(int tscan);
int check_nise_tscan(int nise, int tscan);
int check_nise_trotter(int nise, int trotter);
int check_T2(int T2, int nise);
int check_2d_spec_output(int *POL, int reph, int nreph, int npol, int zzzz, int zzyy, int zyyz, int zyzy,  int do2d);

int check_ham(char* Hamname, int *nframes, int *nosc, int *n2Q, int nbuffer, int nread, int info, int window, int nthreads);
int check_shift(int nshift, int *SHIFTNDX, int nosc);
int check_sites(char* Sitesname, char* Hamname, int nosc, int nframes);
int check_dipoles(char Dipnames[3][maxchar], char *Hamname, int nosc, int nframes);

#endif
