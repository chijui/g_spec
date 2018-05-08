#ifndef TROTTER_H
#define TROTTER_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

#ifndef MEM_HELPER_H
#include "mem_helper.h"
#endif

int set_trotter_1Qarray(trotter_array *trotter, int nosc, int nbuffer);
int set_trotter_2Qarray(trotter_array *trotter, int nosc, int nbuffer);
int unset_trotter_array(trotter_array *trotter);
int construct_trotter_mat(trotter_array *trotter, GNREAL *Ham, int xfr, double expfac);
int construct_trotter_D(trotter_array *trotter, GNREAL *Ham, int xfr, GNREAL expfac);
int construct_trotter_index(trotter_array *trotter);
int construct_trotter_2Qindex(trotter_array *trotter, int nosc);
int construct_trotter_J(trotter_array *trotter, GNREAL *Ham, int xfr, GNREAL expfac);
int calc_CorrFunc_trotter(GNCOMP *CorrFunc, trotter_array *trotter, GNCOMP **psi[3], GNREAL ***Dip1QMat, int window, int fr, int nthreads);
int propagate_psi1Q_trotter(GNCOMP **psi[3], GNREAL **Dip1QMat[3], trotter_array *trotter, int tau0, int fr, int xfrtau, int winsize, int nthreads);
int propagate_psi2Q_trotter(trotter_array *trotter, GNCOMP **psi_c, GNCOMP *psi2Q, int xfr3, int nthreads);

#endif
