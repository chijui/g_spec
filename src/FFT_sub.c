#include "FFT_sub.h"

// cjfeng 12/04/2018
// Using fftw_array to replace individual global fftw variables

// Initialize fftw_array
void initialize_fftw_array_pointer(fftw_array *fftw_array) {
  fftw_array->dim = 0;
  fftw_array->ndpts = 0;
  fftw_array->nzpad = 0;
  fftw_array->FTin = NULL;
  fftw_array->FTout = NULL;
  fftw_array->FTplan = NULL;
}

// cjfeng 12/05/2018
// Subroutine for memory allocation of fftw arrays
int allocate_fftw(fftw_complex **array, int dim, int nzpad) {
  if ( dim == 1 ) *array = (fftw_complex*) fftw_malloc( sizeof(fftw_complex)*nzpad );
  else if ( dim == 2 ) *array = (fftw_complex*) fftw_malloc ( sizeof(fftw_complex)*nzpad*nzpad );
  if ( *array == NULL ) {
    printf("Out of memory for allocating %s.\n", getName(array));
    return 1;
  }
  return 0;
}

// Allocation and de-allocation
int set_fftw_array(fftw_array *fftw_array, int dim, w_param *w_param) {
  int ndpts = w_param->window;
  int nzpad = w_param->winzpad;
  if ( ndpts == 0 || nzpad == 0 ) {
    printf("The array size is zero for %s.\n", getName(fftw_array));
    return 1;
  }
  else {
    int error = 0;
    fftw_array->dim = dim;
    fftw_array->ndpts = ndpts;
    fftw_array->nzpad = nzpad;
    if( !error ) error = allocate_fftw( &(fftw_array->FTin), fftw_array->dim, fftw_array->nzpad );
    else fftw_array->FTin = NULL;
    if( !error ) error = allocate_fftw( &(fftw_array->FTout), fftw_array->dim, fftw_array->nzpad );
    else fftw_array->FTout = NULL;
    
    if(!error) {
      if ( fftw_array->dim == 1) fftw_array->FTplan = fftw_plan_dft_1d(fftw_array->nzpad, fftw_array->FTin, fftw_array->FTout, FFTW_BACKWARD, FFTW_ESTIMATE);
      else if ( fftw_array->dim == 2 ) fftw_array->FTplan = fftw_plan_dft_2d(fftw_array->nzpad, fftw_array->nzpad, fftw_array->FTin, fftw_array->FTout, FFTW_BACKWARD, FFTW_ESTIMATE);
    } else fftw_array->FTplan = NULL;
  } 
  return 0;
}

int unset_fftw_array(fftw_array *fftw_array) {
  if( fftw_array->FTin != NULL ) fftw_free(fftw_array->FTin);
  if( fftw_array->FTout != NULL ) fftw_free(fftw_array->FTout);
  if( fftw_array->FTplan != NULL ) fftw_destroy_plan(fftw_array->FTplan);
  return 0;
}

// cjfeng 12/05/2018
// Initializing fftw_array
int initialize_fftw_array(fftw_array *fftw_array) {
  int i;
  int dim = fftw_array->dim;
  int nzpad;
  if (dim == 1) nzpad = fftw_array->nzpad;
  else if (dim == 2) nzpad = fftw_array->nzpad * fftw_array->nzpad;
  for(i=0; i<nzpad; i++) {
    fftw_array->FTin[i][0] = 0.0;
    fftw_array->FTin[i][1] = 0.0;
  }
  return 0;
}

// cjfeng 12/05/2018
// Applying windowing on 1D array
int windowing1D(GNCOMP *CorrFunc, GNREAL *hann, GNREAL *popdecay, int whann, int ndpts) {
  int i;
  if(whann) for(i=0; i<ndpts; i++) CorrFunc[i] = hann[i] * popdecay[i] * CorrFunc[i];
  else for(i=0; i<ndpts; i++) CorrFunc[i] = popdecay[i] * CorrFunc[i];
  return 0;
}

// cjfeng 12/05/2018
// Fill 1D fftw_array
int fill_FT1D(fftw_array *FT1D, GNCOMP *CorrFunc ) {
  int i;
  int ndpts = FT1D->ndpts;
  int nzpad = FT1D->nzpad;
  for(i=0; i<ndpts; i++) {
    FT1D->FTin[i][0] = creal(CorrFunc[i]);
    FT1D->FTin[i][1] = cimag(CorrFunc[i]);
  }
  return 0;
}

// cjfeng 12/05/2018
// Fill 2D fftw_array, note that *spec = spec[np];
int fill_FT2D(fftw_array *FT2D, GNCOMP *spec) {
  int i, j;
  int ndpts = FT2D->ndpts;
  int nzpad = FT2D->nzpad;
  
  // Windowing and filling FTin
  for(i=0; i<ndpts; i++) {
    int id = i*ndpts;
    int iz = i*nzpad;
    for(j=0; j<ndpts; j++) {
      FT2D->FTin[iz+j][0] = creal(spec[id+j]);
      FT2D->FTin[iz+j][1] = cimag(spec[id+j]);
    }
  }
  return 0;
}

// cjfeng 12/05/2018
// Applying windowing on 2D array
int windowing2D(fftw_array *FT2D, GNREAL *hann, int whann) {
  if(whann) {
    int i, j;
    int ndpts = FT2D->ndpts;
    int nzpad = FT2D->nzpad;
    for(i=0; i<ndpts; i++) {
      GNREAL hann_private = hann[i];
      int id = i*ndpts;
      int iz = i*nzpad;
      for(j=0; j<ndpts; j++) {
        FT2D->FTin[iz+j][0] = FT2D->FTin[iz+j][0] * hann_private;
        FT2D->FTin[iz+j][1] = FT2D->FTin[iz+j][1] * hann_private;
      }
    }
  }
  return 0;
}

// cjfeng 12/05/2018
// Zero padding FT1D
int zero_pad_FT1D(fftw_array *FT1D) {
  int i;
  int ndpts = FT1D->ndpts;
  int nzpad = FT1D->nzpad;
  for(i=ndpts; i<nzpad; i++) {
    FT1D->FTin[i][0] = 0.0;
    FT1D->FTin[i][1] = 0.0;
  }
  return 0;
}

// cjfeng 12/05/2018
// DFT subroutines
int FFT1D(GNCOMP *CorrFunc, GNREAL *hann, GNREAL *popdecay, fftw_array *FT1D, int whann) {
  int i;
  int ndpts = FT1D->ndpts;
  int nzpad = FT1D->nzpad;
  windowing1D(CorrFunc, hann, popdecay, whann, ndpts);
  fill_FT1D(FT1D, CorrFunc);
  zero_pad_FT1D(FT1D);
  // FTIR spectrum now stored in FT1D->FTout[.][0]
  fftw_execute(FT1D->FTplan);
  // Reinitialize CorrFunc for next dump
  for(i=0; i<ndpts; i++) CorrFunc[i] = 0.0;
  
  return 0;
}

// cjfeng 12/05/2018
// Currently we still pass through the whole spec, POL, and p, 
// but in the future it can be reduced to only have spec[POL[p]].
int FFT2D(GNCOMP **spec, GNREAL *hann, fftw_array *FT2D, int *POL, int p, int whann) {
  int i, j;
  int np = POL[p];
  int ndpts = FT2D->ndpts;
  int nzpad = FT2D->nzpad;

  initialize_fftw_array(FT2D);
  fill_FT2D(FT2D, spec[np]);
  windowing2D(FT2D, hann, whann);
  // Output stored in FT2D->FTout[.][0] and FT2D->FTout[.][1]
  fftw_execute(FT2D->FTplan);

  return 0;
}

// For rephasing spectrum only, flip the spectrum along w1.
int flip_w1(fftw_array *FT2D) {
  int i, j;
  int nzpad = FT2D->nzpad;
  int nzpadsq = nzpad * nzpad;
  for(i=0; i<nzpad; i++) {
    int iz = i * nzpad;
    int ii_r = (nzpad-1-i) * nzpad;
    for(j=0; j<nzpad;j++) {
      FT2D->FTin[iz+j][0] = FT2D->FTout[ii_r+j][0];
      FT2D->FTin[iz+j][1] = FT2D->FTout[ii_r+j][1];
    }
  }
  for(i=0; i<nzpadsq; i++) {
    FT2D->FTout[i][0] = FT2D->FTin[i][0];
    FT2D->FTout[i][1] = FT2D->FTin[i][1];
  }
  return 0;
}

// Subroutines for sum-over-state calculations
// cjfeng 01/04/2019
// Using phys_const
int dress_FT1D(GNREAL *ftir, fftw_array *FT1D, int npts, phys_const *Const) {
  int i, j;
  int nzpad = FT1D->nzpad;

  initialize_fftw_array(FT1D);
  for(i=0; i<npts; i++) FT1D->FTin[i][0] = ftir[i];

  fftw_execute(FT1D->FTplan);

  for(j=0; j<nzpad; j++) {
    GNREAL fac = exp( j * (GNREAL) Const->expfac );
    FT1D->FTin[j][0] = FT1D->FTout[j][0] * fac;
    FT1D->FTin[j][1] = -FT1D->FTout[j][1] * fac;
  }

  fftw_execute(FT1D->FTplan);

  for(j=0; j<nzpad; j++) FT1D->FTout[j][1] = -FT1D->FTout[j][1];
  for(i=0; i<npts; i++) ftir[i] = FT1D->FTout[i][0];

  return 0;
} 

// cjfeng 12/05/2018
// Merging the rephasing and non-rephasing into a single function
// Using flip to indicate rephasing (flip=1) or non-rephasing (flip=0);
// cjfeng 01/04/2019
// Using phys_const
int dress_FT2D(GNCOMP **spec, int *POL, fftw_array *FT2D, int p, int npts, phys_const *Const, int flip) {
  int i, j;
  int np = POL[p];
  int nzpad = FT2D->nzpad;
  int nzpadsq = nzpad * nzpad;

  initialize_fftw_array(FT2D);

  // Dress spectrum with Lorentzian profile. Note the it is so far purely real.
  for(i=0; i<npts; i++) {
    int iz = i * nzpad;
    int ii = i * npts;
    for(j=0; j<npts; j++) FT2D->FTin[iz+j][0] = creal(spec[np][ii+j]);
  }

  fftw_execute(FT2D->FTplan);
  
  for(i=0; i<nzpad; i++) {
    int iz = i * nzpad;
    for(j=0; j<nzpad; j++) {
      int k;
      if (flip==1) k = nzpad-i+j;
      else k = i+j;
      GNREAL fac = exp( k * (GNREAL) Const->expfac );
      FT2D->FTin[iz+j][0] = FT2D->FTout[iz+j][0] * fac;
      FT2D->FTin[iz+j][1] = -FT2D->FTout[iz+j][1] * fac;
    }
  }

  fftw_execute(FT2D->FTplan);

  for(i=0; i<nzpad; i++) {
    int iz = i * nzpad;
    for(j=0; j<nzpad; j++) FT2D->FTout[iz+j][1] = -FT2D->FTout[iz+j][1];
  }
  for(i=0; i<npts; i++) {
    int ii = i * npts; 
    int iz = i * nzpad;
    for(j=0; j<npts; j++) spec[np][ii+j] = FT2D->FTout[iz+j][0] + I*FT2D->FTout[iz+j][1];
  }
  return 0;
}

// cjfeng 12/09/2018
// Using a single sub-routine to perform FFT on linear response
int FFT_linear(GNCOMP *CorrFunc, GNCOMP *NetCorrFunc, GNREAL *hann, GNREAL *popdecay, fftw_array *FT1D, w_param *w_param, fout *fout, int flag_traj) {
  int i;
  int window = w_param->window;
  if (flag_traj) {
    for(i=0; i<window; i++) NetCorrFunc[i] += CorrFunc[i];
    FFT1D(CorrFunc, hann, popdecay, FT1D, w_param->whann);
    print_nise_ftir_traj(fout->traj[1].fp, FT1D, w_param);
  }
  else {
    for(i=0; i<window; i++) CorrFunc[i] += NetCorrFunc[i];
    FFT1D(CorrFunc, hann, popdecay, FT1D, w_param->whann);
    print_nise_ftir(fout->labs.fp, FT1D, w_param);
  }
  return 0;
}

// cjfeng 12/09/2018
// Using a single sub-routine to perform FFT on nonlinear response
int FFT_nonlinear(GNCOMP *REPH[4], GNCOMP *NetREPH[4], GNCOMP *NREPH[4], GNCOMP *NetNREPH[4], GNREAL *hann, fftw_array *FT2D, POL_info *pol_info, spec_param *spec_param, w_param *w_param, fout *fout, int flag_traj) {
  if(spec_param->do2d) {
    int i;
    int reph = spec_param->reph;
    int nreph = spec_param->nreph;
    int window = w_param->window;
    int whann = w_param->whann;
    int windowsq = window * window;
    int p;
    int npol = pol_info->npol;
    int POL[4];
    for(p=0; p<4; p++) POL[p] = pol_info->POL[p];
    // FT and print rephasing spectra
    if (reph) {
      for(p=0; p<npol; p++) {
        int np = POL[p];
        if(flag_traj) for(i=0; i<windowsq; i++) NetREPH[np][i] += REPH[np][i];
        else for(i=0; i<windowsq; i++) REPH[np][i] += NetREPH[np][i];
        FFT2D(REPH, hann, FT2D, POL, p, whann);
        // For rephasing only: flip spectrum along w1
        flip_w1(FT2D);
        if(flag_traj) {
          print_nise_2d_traj(fout->traj[2+np].fp, FT2D, w_param);
          for(i=0; i<windowsq; i++) REPH[np][i] = 0.0;
        }
        else print_nise_2d(fout->reph[np].fp, FT2D, w_param);
      }
    }
    // FT and print non-rephasing spectra
    if (nreph) {
      for(p=0; p<npol; p++) {
        int np = POL[p];
        if(flag_traj) for(i=0; i<windowsq; i++) NetNREPH[np][i] += NREPH[np][i];
        else for(i=0; i<windowsq; i++) NREPH[np][i] += NetNREPH[np][i];
        FFT2D(NREPH, hann, FT2D, POL, p, whann);
        if(flag_traj) {
          print_nise_2d_traj(fout->traj[6+np].fp, FT2D, w_param);
          for(i=0; i<windowsq; i++) NREPH[np][i] = 0.0;
        }
        else print_nise_2d(fout->nreph[np].fp, FT2D, w_param);
      }
    }
  }
  return 0;
}

// cjfeng 01/04/2019
// Using phys_const
int FFT_static_linear(GNREAL *ftir, GNREAL *netftir, fftw_array *FT1D, spec_param *spec_param, int npts, phys_const *Const, fout *fout, int flag_traj) {
  int i;
  if(flag_traj) {
    for(i=0; i<npts; i++) netftir[i] += ftir[i];
    dress_FT1D(ftir, FT1D, npts, Const);
    print_ftir_traj(fout->traj[1].fp, ftir, npts);
    for(i=0; i<npts; i++) ftir[i] = 0.0;
  }
  else {
    for(i=0; i<npts; i++) ftir[i] += netftir[i];
    dress_FT1D(ftir, FT1D, npts, Const);
    print_ftir(fout->labs.fp, ftir, npts);
  }
  return 0;
}
// cjfeng 01/04/2019
// Using phys_const
int FFT_static_nonlinear(GNCOMP *REPH[4], GNCOMP *NetREPH[4], GNCOMP *NREPH[4], GNCOMP *NetNREPH[4], fftw_array *FT2D, POL_info *pol_info, spec_param *spec_param, int npts, phys_const *Const, fout *fout, int flag_traj) {
  if(spec_param->do2d) {
    int nptssq = npts * npts;
    int reph = spec_param->reph;
    int nreph = spec_param->nreph;
    int i, p;
    int npol = pol_info->npol;
    int POL[4];
    for(p=0; p<4; p++) POL[p] = pol_info->POL[p];
    // First rephasing data
    if (reph) {
      for(p=0; p<npol; p++) {
        int np = POL[p];
        // cjfeng 04/11/2017
        if(flag_traj) for(i=0; i<nptssq; i++) NetREPH[np][i] += creal(REPH[np][i]);
        else for(i=0; i<nptssq; i++) REPH[np][i] += NetREPH[np][i];
        // cjfeng 12/05/2018
        // Using a unified function for both REPH and NREPH
        dress_FT2D(REPH, POL, FT2D, p, npts, Const, 1);
        if(flag_traj) {
          print_2d_traj(fout->traj[2+np].fp, REPH, npts, POL, p);            
          for(i=0; i<nptssq; i++) REPH[np][i] = 0.0;
        }
        else print_2d(fout->reph[np].fp, REPH, npts, POL, p);
      }
    }
    // Then non-rephasing data
    if (nreph) {
      for(p=0; p<npol; p++) {
        int np = POL[p];
        // cjfeng 04/11/2017
        if(flag_traj) for(i=0; i<nptssq; i++) NetNREPH[np][i] += creal(NREPH[np][i]);
        else for(i=0; i<nptssq; i++) NREPH[np][i] += NetNREPH[np][i];
        // cjfeng 12/05/2018
        // Using a unified function for both REPH and NREPH
        dress_FT2D(NREPH, POL, FT2D, p, npts, Const, 0);
        if(flag_traj) {
          print_2d_traj(fout->traj[6+np].fp, NREPH, npts, POL, p);            
          for(i=0; i<nptssq; i++) NREPH[np][i] = 0.0;
        }
        else print_2d(fout->nreph[np].fp, NREPH, npts, POL, p);
      }
    }
  }
  return 0;
}
