#ifndef FILE_IO_H
#define FILE_IO_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

int read_line(FILE* fp, int nrows, int ncols, GNREAL *Mat);
int read_info( FILE *fp, char* Infoname, char*** p_NNames, int** p_NNums, char*** p_CNames, int** p_CNums);

int read_ham(FILE *fp, GNREAL *ham, int nosc);
int read_dip(FILE **fp, GNREAL ***Dip, int nosc, int xfr);
int read_sites(FILE *fp, GNREAL *sites, GNREAL **ham, int nosc, int xfr);

int parse_shifts(char* Parname, int *p_nshift, int maxchar);

int make_fnames( char* Infoname, char* Hamname, char* hamnm, char Dipnames[3][maxchar], char* dipxnm, char* dipynm, char* dipznm, char* Sitesname, char* sitesnm, char* axisnm, char* ftirnm, char* lognm, char polnm[4][16], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnm[10][maxchar], char* Parname, char* paramnm, char* outname, char* deffnm, int do_traj );

int assign_input_files( char* Infoname, char* Hamname, char* hamnm, char Dipnames[3][maxchar], char* dipxnm, char* dipynm, char* dipznm, char* Sitesname, char* sitesnm, char* Parname, char* paramnm, char* deffnm );
int assign_input_name(char* fnm, char* new_suffix, char* prefix, char* default_suffix);
int assign_prefix( char* fnm, char* prefix);
int copy_outname(char* axisnm, char* ftirnm, char* lognm, char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnm[10][maxchar], char* outname, int do_traj );
int add_suffix(char* axisnm, char* ftirnm, char* lognm, char polnm[4][16], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnm[10][maxchar], char* outname, int do_traj );

int file_opener(FILE **fp, char *fnm, char* ftype, int status);
int open_all( char* Hamname, char Dipnames[3][maxchar], char* Sitesname, char* axisnm, char* ftirnm, char* lognm, int npol, int POL[4], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnames[10][maxchar], int do2d, int do_traj);

int print_frame(int frame, int frame_period, FILE *lfp);
int print_memory(int frame, int frame_period, struct rusage *r_usage, FILE *lfp);

int print_nise_waxis(FILE *fp, GNREAL dw, int nprint, int ndxstart);
int print_waxis(FILE *fp, real wres, int npts, real wstart);

int print_nise_ftir(FILE *fp, fftw_complex *FTout1D, int nprint, int ndxstart, int winzpad);
int print_nise_ftir_traj(FILE *fp, fftw_complex *FTout1D, int nprint, int ndxstart, int winzpad);
int print_ftir(FILE *fp, GNREAL *ftir, int npts);
int print_ftir_traj(FILE *fp, GNREAL *ftir, int npts);

int print_nise_2d(FILE *fp, fftw_complex *FTout2D, int nprint, int ndxstart, int winzpad);
int print_nise_2d_traj(FILE *fp, fftw_complex *FTout2D, int nprint, int ndxstart, int winzpad);

int print_2d(FILE *fp, GNCOMP **spec, int npts, int *POL, int p);
int print_2d_traj(FILE *fp, GNCOMP **spec, int npts, int *POL, int p);

// cjfeng 05/07/2018
int read_dip2Q(FILE **fp, GNREAL **Dip, int n2Q);
int make_fnames_new( char* Infoname, char* Hamname, char* hamnm, char Dipnames[3][maxchar], char* dipxnm, char* dipynm, char* dipznm, char* Ham2Qname, char* ham2Qnm, char Dip2Qnames[3][maxchar], char* dipx2Qnm, char* dipy2Qnm, char* dipz2Qnm, char* Sitesname, char* sitesnm, char* axisnm, char* ftirnm, char* lognm, char polnm[4][16], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnm[10][maxchar], char* Parname, char* paramnm, char* outname, char* deffnm, int do_traj );
int assign_input_files_new( char* Infoname, char* Hamname, char* hamnm, char Dipnames[3][maxchar], char* dipxnm, char* dipynm, char* dipznm, char* Ham2Qname, char* ham2Qnm, char Dip2Qnames[3][maxchar], char* dipx2Qnm, char* dipy2Qnm, char* dipz2Qnm, char* Sitesname, char* sitesnm, char* Parname, char* paramnm, char* deffnm );
int open_all_new( char* Hamname, char Dipnames[3][maxchar], char* Ham2Qname, char Dip2Qnames[3][maxchar], char* Sitesname, char* axisnm, char* ftirnm, char* lognm, int npol, int POL[4], char rephnm[4][maxchar], char nrephnm[4][maxchar], char Trajnames[10][maxchar], int do2d, int do_traj);

#endif
