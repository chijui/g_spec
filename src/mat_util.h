#ifndef MAT_UTIL_H
#define MAT_UTIL_H

#ifndef DEF_H
#include "def.h"
#endif

int copy_gncomp(GNCOMP *A, GNCOMP *B, int nrows, int ncols);
int zero_gncomp_mat(int nrows, int ncols, GNCOMP *Mat);
int gncomp_identity_mat(int N, GNCOMP *Mat);

int trans_comp(GNCOMP *A, GNCOMP *A_t, int nrows, int ncols, int nt);
int trans_real(GNREAL *A, GNREAL *A_t, int nrows, int ncols, int nt);
int mmult_comp_trans(GNCOMP *A, GNCOMP *B_t, GNCOMP *C, int dim, int nt);
int mvmult_real_serial_trans(GNREAL *A, GNREAL *vin, GNREAL *vout, int dim, int nt);
int mvmult_comp(GNCOMP *A, GNCOMP *vin, GNCOMP *vout, int dim, int nt);
int mvmult_comp_trans_x(GNREAL *A_t, GNCOMP *vin, GNCOMP *vout, int dimin, int dimout, int nt);

#endif
