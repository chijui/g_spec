#ifndef SPEC_FFT_H
#define SPEC_FFT_H

#ifndef DEF_H
#include "def.h"
#endif

int FFT_ftir(GNCOMP *CorrFunc, GNREAL *hann, GNREAL *popdecay, fftw_plan FTplan1D, fftw_complex *FTin1D, fftw_complex *FTout1D, int whann, int window, int winzpad, double c, int TauP, real wres);
int FFT_2d(GNCOMP **spec, GNREAL *hann, int *POL, fftw_plan FTplan2D, fftw_complex *FTin2D, fftw_complex *FTout2D, int p, int whann, int window, int winzpad);
int FFT_flip_w1(fftw_complex *FTin2D, fftw_complex *FTout2D, int winzpad);

int dress_ftir_lorentzian(GNREAL *ftir,  fftw_plan FTplan1D, fftw_complex *FTin1D, fftw_complex *FTout1D, int npts, int winzpad, double c, int TauP, real wres);

int dress_reph_Lorentzian(GNCOMP **spec, int *POL, fftw_plan FTplan2D, fftw_complex *FTin2D, fftw_complex *FTout2D, int p, int npts, int winzpad, double c, int TauP, real wres);
int dress_nreph_Lorentzian(GNCOMP **spec, int *POL, fftw_plan FTplan2D, fftw_complex *FTin2D, fftw_complex *FTout2D, int p, int npts, int winzpad, double c, int TauP, real wres);

#endif
