#ifndef PARAM_CHECKER_H
#define PARAM_CHECKER_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

int check_tstep(int tstep);
int check_tscan(int tscan, int nise);
int check_nise_flags(int nise, int trotter, int pert);
int check_T2(int T2, int nise);
int check_2d_spec_output(POL_info *pol_info, int reph, int nreph, int do2d);

int check_input_files(fin *fin, spec_param *spec_param, w_param *w_param, shift_info *shift_info, int *nframes, traj_param *traj_param, int nbonds);
int check_ham(char* Hamname, int *nframes, traj_param *traj_param, int info, int window, int nthreads);
int check_shift(shift_info *shift_info, int nosc);
int check_sites(char* Sitesname, char* Hamname, int nosc, int nframes);
int check_dipole(char *Dipnames, char *Hamname, int nosc, int nframes);
int check_elecs(char* Elecname, char* Hamname, int nosc, int nframes);

int check_ham2Q(char* Ham2Qname, int *nframes, traj_param *traj_param, int info, int win2d, int nthreads);
int check_dip2Q(char *Dip2Qnames, char *Ham2Qname, int n2Q, int nframes);
int check_nbonds(int nbonds);

int check_spec_param(spec_param *spec_param, POL_info *pol_info);
int check_w_param(w_param *w_param);

#endif
