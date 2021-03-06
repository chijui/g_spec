#ifndef MEM_HELPER_H
#define MEM_HELPER_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

#ifndef FILE_IO_H
#include "file_io.h"
#endif

int allocate_all(trotter_array *trotter1Q, trotter_array *trotter2Q, fftw_array *FT1D, fftw_array *FT2D, w_param *w_param, spec_param *spec_param, traj_param *traj_param );
int graceful_exit(int error, fin *fin, fout *fout, trotter_array *trotter1Q, trotter_array *trotter2Q, fftw_array *FT1D, fftw_array *FT2D, traj_param *traj_param, w_param *w_param, spec_param *spec_param, POL_info *pol_info, shift_info *shift_info);

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
int set_gncomp_3d_array(GNCOMP ****mat, int dim1, int dim2, int dim3);
int unset_gncomp_1d_array(GNCOMP *mat);
int unset_gncomp_2d_array(GNCOMP **mat, int dim1);
int unset_gncomp_3d_array(GNCOMP ***mat, int dim1, int dim2);

int set_int_1d_array(int **mat, int nrows, int ncols);
int unset_int_1d_array(int *mat);

int set_char_1d_array(char **mat, int dim1);
int set_char_2d_array(char ***mat, int dim1, int dim2);
int unset_char_1d_array(char *mat);
int unset_char_2d_array(char **mat, int dim1);

int initialize_POL_info(POL_info *pol_info);
int set_POL_info(POL_info *pol_info);

int set_fin(fin *fin, int maxchar);
int unset_fin(fin *fin);

int set_fout(fout *fout, int maxchar);
int unset_fout(fout *fout);

int initialize_res_info(res_info *res_info);
int set_res_info(res_info *res_info, int nbonds, int nchar);
int unset_res_info(res_info *res_info);
int initialize_shift_info(shift_info *shift_info);
int set_shift_info(shift_info *shift_info);
int unset_shift_info(shift_info *shift_info);

#endif
