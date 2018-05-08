#ifndef MEM_HELPER_H
#define MEM_HELPER_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

int allocate_all( int window, int win2d, int winzpad, int nosc, int nbuffer, int nthreads, int n2Q, int pert, int do2d, GNREAL tstep, GNREAL TauP, int npts, int nise );
int graceful_exit( int error, int nbuffer, int win2d, int nthreads, int npol, int nise, int nosc); 
// cjfeng 05/07/2018
int allocate_all_new( int window, int win2d, int winzpad, int nosc, int nbuffer, int nthreads, int n2Q, int pert, int do2d, int if2Qfiles, GNREAL tstep, GNREAL TauP, int npts, int nise );
int graceful_exit_new( int error, int nbuffer, int win2d, int nthreads, int npol, int nise, int nosc); 

int set_eig_array(eig_array *eig, int n);
int unset_eig_array(eig_array *eig);

int set_gnreal_1d_array(GNREAL **mat, int nrows, int ncols);
int set_gnreal_2d_array(GNREAL ***mat, int dim1, int dim2);
int set_gnreal_3d_array(GNREAL ****mat, int dim1, int dim2, int dim3);
int unset_gnreal_1d_array(GNREAL *mat);
int unset_gnreal_2d_array(GNREAL **mat, int dim1);
int unset_gnreal_3d_array(GNREAL ***mat, int dim1, int dim2);

int set_gncomp_1d_array(GNCOMP **mat, int nrows, int ncols);
int set_gncomp_2d_array(GNCOMP ***mat, int dim1, int dim2);
int unset_gncomp_1d_array(GNCOMP *mat);
int unset_gncomp_2d_array(GNCOMP **mat, int dim1);

int set_int_1d_array(int **mat, int nrows, int ncols);
int unset_int_1d_array(int *mat);

int set_char_1d_array(char **mat, int dim1);
int set_char_2d_array(char ***mat, int dim1, int dim2);
int unset_char_1d_array(char *mat);
int unset_char_2d_array(char **mat, int dim1);

#endif
