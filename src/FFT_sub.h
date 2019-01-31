#ifndef FFT_SUB_H
#define FFT_SUB_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

#ifndef MEM_HELPER_H
#include "mem_helper.h"
#endif

void initialize_fftw_array_pointer(fftw_array *fftw_array);
int allocate_fftw(fftw_complex **array, int dim, int nzpad);
int set_fftw_array(fftw_array *fftw_array, int dim, w_param *w_param);
int unset_fftw_array(fftw_array *fftw_array);

int initialize_fftw_array(fftw_array *fftw_array);
int windowing1D(GNCOMP *CorrFunc, GNREAL *hann, GNREAL *popdecay, int whann, int ndpts);
int fill_FT1D(fftw_array *FT1D, GNCOMP *CorrFunc );
int fill_FT2D(fftw_array *FT2D, GNCOMP *spec);
int windowing2D(fftw_array *FT2D, GNREAL *hann, int whann);
int zero_pad_FT1D(fftw_array *FT1D);
int FFT1D(GNCOMP *CorrFunc, GNREAL *hann, GNREAL *popdecay, fftw_array *FT1D, int whann);
int FFT2D(GNCOMP **spec, GNREAL *hann, fftw_array *FT2D, int *POL, int p, int whann);

int flip_w1(fftw_array *FT2D);

int dress_FT1D(GNREAL *ftir, fftw_array *FT1D, int npts, phys_const *Const);
int dress_FT2D(GNCOMP **spec, int *POL, fftw_array *FT2D, int p, int npts, phys_const *Const, int flip);

int FFT_linear(GNCOMP *CorrFunc, GNCOMP *NetCorrFunc, GNREAL *hann, GNREAL *popdecay, fftw_array *FT1D, w_param *w_param, fout *fout, int flag_traj);
int FFT_nonlinear(GNCOMP *REPH[4], GNCOMP *NetREPH[4], GNCOMP *NREPH[4], GNCOMP *NetNREPH[4], GNREAL *hann, fftw_array *FT2D, POL_info *pol_info, spec_param *spec_param, w_param *w_param, fout *fout, int flag_traj);
int FFT_static_linear(GNREAL *ftir, GNREAL *netftir, fftw_array *FT1D, spec_param *spec_param, int npts, phys_const *Const, fout *fout, int flag_traj);
int FFT_static_nonlinear(GNCOMP *REPH[4], GNCOMP *NetREPH[4], GNCOMP *NREPH[4], GNCOMP *NetNREPH[4], fftw_array *FT2D, POL_info *pol_info, spec_param *spec_param, int npts, phys_const *Const, fout *fout, int flag_traj);

#endif
