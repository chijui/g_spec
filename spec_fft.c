#include "spec_fft.h"

int FFT_ftir(GNCOMP *CorrFunc, GNREAL *hann, GNREAL *popdecay, fftw_plan FTplan1D, fftw_complex *FTin1D, fftw_complex *FTout1D, int whann, int window, int winzpad, double c, int TauP, real wres) {
	int i;
	if(whann) for(i=0; i<window; i++) CorrFunc[i] = hann[i] * popdecay[i] * CorrFunc[i];
	else for(i=0; i<window; i++) CorrFunc[i] = popdecay[i] * CorrFunc[i];
	for(i=0; i<window; i++) {
		FTin1D[i][0] = creal(CorrFunc[i]);
		FTin1D[i][1] = cimag(CorrFunc[i]);
	}
	for(i=window; i<winzpad; i++) {
		FTin1D[i][0] = 0.0;
		FTin1D[i][1] = 0.0;
	}
	// FTIR spectrum now stored in FTout1D[.][0]
	fftw_execute(FTplan1D);
	// Reinitialize CorrFunc for next dump
	for(i=0; i<window; i++) CorrFunc[i] = 0.0;
	
	return 0;
}

int FFT_2d(GNCOMP **spec, GNREAL *hann, int *POL, fftw_plan FTplan2D, fftw_complex *FTin2D, fftw_complex *FTout2D, int p, int whann, int window, int winzpad) {
	int i, j;
	int np = POL[p];
	int winzpadsq = winzpad * winzpad;
	for(i=0; i<winzpadsq; i++) {
		FTin2D[i][0] = 0.0;
		FTin2D[i][1] = 0.0;
	}
	// Windowing and filling FTin2D
	if(whann) {
		for(i=0; i<window; i++) {
			GNREAL hann_private = hann[i];
			int iwin = i*window;
			int izpad = i*winzpad;
			for(j=0; j<window; j++) {
				FTin2D[izpad+j][0] = creal(spec[np][iwin+j]) * hann_private;
				FTin2D[izpad+j][1] = cimag(spec[np][iwin+j]) * hann_private;
			}
		}
	}
	else {
		for(i=0; i<window; i++) {
			int iwin = i*window;
			int izpad = i*winzpad;
			for(j=0; j<window; j++) {
				FTin2D[izpad+j][0] = creal(spec[np][iwin+j]);
				FTin2D[izpad+j][1] = cimag(spec[np][iwin+j]);
			}
		}
	}	
	fftw_execute(FTplan2D);

	return 0;
}

// For rephasing spectrum only, flip the spectrum along w1.
int FFT_flip_w1(fftw_complex *FTin2D, fftw_complex *FTout2D, int winzpad) {
	int i, j;
	int  winzpadsq = winzpad * winzpad;
	for(i=0; i<winzpad;i++) {
		int izpad = i*winzpad;
		int ii_r = (winzpad-1-i)*winzpad;
		for(j=0; j<winzpad;j++) {
			FTin2D[izpad+j][0] = FTout2D[ii_r+j][0];
			FTin2D[izpad+j][1] = FTout2D[ii_r+j][1];
		}
	}
	for(i=0; i<winzpadsq; i++) {
		FTout2D[i][0] = FTin2D[i][0];
		FTout2D[i][1] = FTin2D[i][1];
	}
	return 0;
}

int dress_ftir_lorentzian(GNREAL *ftir,  fftw_plan FTplan1D, fftw_complex *FTin1D, fftw_complex *FTout1D, int npts, int winzpad, double c, int TauP, real wres) {
	int i, j;
	GNREAL period = 2 * TauP * c * wres * winzpad;
	for(j=0; j<winzpad; j++) {
		FTin1D[j][0] = 0.0;
		FTin1D[j][1] = 0.0;	
	}
	for(i=0; i<npts; i++) FTin1D[i][0] = ftir[i];
	// cjfeng 07/05/2016
	// Should it be FORWARD or -iwt by another plan?
	fftw_execute(FTplan1D);
	for(j=0; j<winzpad; j++) {
		FTin1D[j][0] = FTout1D[j][0] * exp(-j / period);
		FTin1D[j][1] = -FTout1D[j][1] * exp(-j / period);
	}
	fftw_execute(FTplan1D);
	// cjfeng 07/05/2016
	// Normalization?
	for(j=0; j<winzpad; j++) FTout1D[j][1] = -FTout1D[j][1];
	for(i=0; i<npts; i++) ftir[i] = FTout1D[i][0];
	return 0;
} 

int dress_reph_Lorentzian(GNCOMP **spec, int *POL, fftw_plan FTplan2D, fftw_complex *FTin2D, fftw_complex *FTout2D, int p, int npts, int winzpad, double c, int TauP, real wres) {
	int i, j;
	int np = POL[p];
	int winzpadsq = winzpad * winzpad;
	GNREAL period = 2 * TauP * c * wres * winzpad;

	// Dress spectrum with Lorentzian profile only if vibrational lifetime is greater than zero.
	if(TauP>0) {
		for(i=0; i<winzpadsq; i++) {
			FTin2D[i][0] = 0.0;
			FTin2D[i][1] = 0.0;
		}
		// Dress spectrum with Lorentzian profile. Note the it is so far purely real.
		for(i=0; i<npts; i++) {
			int ii_pad = i * winzpad;
			int ii = i * npts;
			for(j=0; j<npts; j++) FTin2D[ii_pad+j][0] = creal(spec[np][ii+j]);
		}
		// cjfeng 07/05/2016
		// Should it be FORWARD or -iwt?
		fftw_execute(FTplan2D);
	
		for(i=0; i<winzpad; i++) {
			int ii_pad = i * winzpad;
			for(j=0; j<winzpad; j++) {
				int k = winzpad-i+j;
				FTin2D[ii_pad+j][0] = FTout2D[ii_pad+j][0] * exp( -k / period );
				FTin2D[ii_pad+j][1] = -FTout2D[ii_pad+j][1] * exp( -k / period );
			}
		}
		// cjfeng 07/05/2016
		// Normalization?
		fftw_execute(FTplan2D);

		for(i=0; i<winzpad; i++) {
			int ii_pad = i * winzpad;
			for(j=0; j<winzpad; j++) FTout2D[ii_pad+j][1] = -FTout2D[ii_pad+j][1];
		}
		for(i=0; i<npts; i++) {
			int ii = i * npts; 
			int ii_pad = i * winzpad;
			for(j=0; j<npts; j++) spec[np][ii+j] = FTout2D[ii_pad+j][0] + I*FTout2D[ii_pad+j][1];
		}
	}

	return 0;
}

int dress_nreph_Lorentzian(GNCOMP **spec, int *POL, fftw_plan FTplan2D, fftw_complex *FTin2D, fftw_complex *FTout2D, int p, int npts, int winzpad, double c, int TauP, real wres) {
	int i, j;
	int np = POL[p];
	int winzpadsq = winzpad * winzpad;
	GNREAL period = 2 * TauP * c * wres * winzpad;

	// Dress spectrum with Lorentzian profile only if vibrational lifetime is greater than zero.
	if(TauP>0) {
		for(i=0; i<winzpadsq; i++) {
			FTin2D[i][0] = 0.0;
			FTin2D[i][1] = 0.0;
		}
		// Dress spectrum with Lorentzian profile. Note the it is so far purely real.
		for(i=0; i<npts; i++) {
			int ii_pad = i * winzpad;
			int ii = i * npts;
			for(j=0; j<npts; j++) FTin2D[ii_pad+j][0] = creal(spec[np][ii+j]);
		}

		// cjfeng 07/05/2016
		// Should it be FORWARD or -iwt?
		fftw_execute(FTplan2D);
	
		for(i=0; i<winzpad; i++) {
			int ii_pad = i * winzpad;
			for(j=0; j<winzpad; j++) {
				int k = i+j;
				FTin2D[ii_pad+j][0] = FTout2D[ii_pad+j][0] * exp( -k / period );
				FTin2D[ii_pad+j][1] = -FTout2D[ii_pad+j][1] * exp( -k / period );
			}
		}
		// cjfeng 07/05/2016
		// Normalization?
		fftw_execute(FTplan2D);

		for(i=0; i<winzpad; i++) {
			int ii_pad = i * winzpad;
			for(j=0; j<winzpad; j++) FTout2D[ii_pad+j][1] = -FTout2D[ii_pad+j][1];
		}
		for(i=0; i<npts; i++) {
			int ii = i * npts; 
			int ii_pad = i * winzpad;
			for(j=0; j<npts; j++) spec[np][ii+j] = FTout2D[ii_pad+j][0] + I*FTout2D[ii_pad+j][1];
		}
	}

	return 0;
}
