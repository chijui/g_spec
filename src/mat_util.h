#ifndef MAT_UTIL_H
#define MAT_UTIL_H

#ifndef DEF_H
#include "def.h"
#endif

int copy_gncomp(GNCOMP *A, GNCOMP *B, int nrows, int ncols);

int mvmult_real_serial_trans(GNREAL *A, GNREAL *vin, GNREAL *vout, int dim, int nt);
int mvmult_comp(GNCOMP *A, GNCOMP *vin, GNCOMP *vout, int dim, int nt);
int mvmult_comp_trans_x(GNREAL *A_t, GNCOMP *vin, GNCOMP *vout, int dimin, int dimout, int nt);
int mvmult_comp_trans_Sp(int *A_row, int *A_col, GNREAL *A_Sp, GNCOMP *vin, GNCOMP *vout, int dimin, int dimout, int nt);

#endif
