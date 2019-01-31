#ifndef FILE_IO_H
#define FILE_IO_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

int read_line(FILE* fp, int nrows, int ncols, GNREAL *Mat);
int read_info( fin *fin, res_info *res_info, int maxchar);

int read_ham(FILE *fp, GNREAL *ham, int nosc);
int read_dip(f_pointer *fp, GNREAL **Dip, int nosc);
int read_sites(FILE *fp, GNREAL *ham, int nosc);

int read_param(fnm_param *fnm_param, char* Parname, res_info *res_info, shift_info *shift_ino, int maxchar);
int parse_shifts(char* Parname, shift_info *shift_info, int maxchar);

int make_fnames(fin *fin, fout *fout, fnm_param *fnm_param, char* Parname, int do_traj );
int assign_input_files( fin *fin, fnm_param *fnm_param, char* Parname );
int assign_input_name(char* fnm, char* new_suffix, char* prefix, char* default_suffix);
int assign_prefix( char* fnm, char* prefix);

int copy_outname(fout *fout, char* outname, int do_traj );
int add_suffix(fout *fout, char *polnm[4], int do_traj );

int file_opener(FILE **fp, char *fnm, char* ftype, int status);

int open_all(fin *fin, fout *fout, POL_info *pol_info, spec_param *spec_param);

int print_frame(int frame, int frame_period, FILE *lfp);
int print_memory(int frame, int frame_period, struct rusage *r_usage, FILE *lfp);
int print_elapsed_time(double t, FILE *lfp);
void print_gmx_quote(FILE *lfp);
int print_tstamp(int frame, int tstep, FILE *fp);

int print_nise_waxis(FILE *fp, w_param *w_param);
int print_waxis(FILE *fp, w_param *w_param);

int print_nise_ftir(FILE *fp, fftw_array *FT1D, w_param *w_param);
int print_nise_ftir_traj(FILE *fp, fftw_array *FT1D, w_param *w_param);
int print_nise_2d(FILE *fp, fftw_array *FT2D, w_param *w_param);
int print_nise_2d_traj(FILE *fp, fftw_array *FT2D, w_param *w_param);

int print_ftir(FILE *fp, GNREAL *ftir, int npts);
int print_ftir_traj(FILE *fp, GNREAL *ftir, int npts);
int print_2d(FILE *fp, GNCOMP **spec, int npts, int *POL, int p);
int print_2d_traj(FILE *fp, GNCOMP **spec, int npts, int *POL, int p);

#endif
