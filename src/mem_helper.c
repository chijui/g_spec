#include "mem_helper.h"

// cjfeng 05/07/2018
// Added Ham2QMat into allocations
// cjfeng 12/08/2018
// Using spec_param, and w_param
// cjfeng 01/03/2019
// Putting if2Qfiles into spec_param
// cjfeng 01/04/2019
// Incorporating fftw_array and trotter_array into allocations
// Added traj_param to replace nosc, nbuffer, and n2Q
int allocate_all(trotter_array *trotter1Q, trotter_array *trotter2Q, fftw_array *FT1D, fftw_array *FT2D, w_param *w_param, spec_param *spec_param, traj_param *traj_param ) {
  int error = 0;
  int i,j,k;
  int window = w_param->window;
  int win2d = w_param->win2d;
  int winzpad = w_param->winzpad;
  int npts = w_param->npts;

  int nosc = traj_param->nosc;
  int n2Q = traj_param->n2Q;
  int nbuffer = traj_param->nbuffer; 

  int npopdec = 2*(win2d-window);  // equivalent with (window + T2/tstep + 1)
  int nptssq;
  int noscsq = nosc*nosc;
  int n2Qsq = n2Q*n2Q;

  int nise = spec_param->nise;
  int do2d = spec_param->do2d;
  int reph = spec_param->reph;
  int nreph = spec_param->nreph;
  int pert = spec_param->pert;
  int trotter = spec_param->trotter;
  int nthreads = spec_param->nthreads;
  int if2Qfiles = spec_param->if2Qfiles;
  GNREAL tstep = spec_param->tstep;

  if(nise) nptssq = window * window;
  else nptssq = npts * npts;

  // cjfeng 01/03/2019
  // The output spectra are all common in various schemes
  if( do2d ) {
    for(i=0; i<4; i++) {
      if( (!error) && (do2d) && (reph) ) {
        error = set_gncomp_1d_array(&REPH[i], 1, nptssq);
        if(!error) for(j=0; j<nptssq; j++) REPH[i][j] = 0.0;
      } else REPH[i] = NULL;
    }
    for(i=0; i<4; i++) {
      if( (!error) && (nreph) ) {
        error = set_gncomp_1d_array(&NREPH[i], 1, nptssq);
        if(!error) for(j=0; j<nptssq; j++) NREPH[i][j] = 0.0;
      } else NREPH[i] = NULL;
    }
    for(i=0; i<4; i++) {
      if( (!error) && (reph) ) {
        error = set_gncomp_1d_array(&NetREPH[i], 1, nptssq);
        if(!error) for(j=0; j<nptssq; j++) NetREPH[i][j] = 0.0;
      } else NetREPH[i] = NULL;
    }
    for(i=0; i<4; i++) {
      if( (!error) && (nreph) ) {
        error = set_gncomp_1d_array(&NetNREPH[i], 1, nptssq);
        if(!error) for(j=0; j<nptssq; j++) NetNREPH[i][j] = 0.0;
      } else NetNREPH[i] = NULL;
    }
  }
  // FFT
  // cjfeng 01/04/2019
  // the fftw_array will be incorporated into this routine
  set_fftw_array(FT1D, 1, w_param);
  if(do2d) set_fftw_array(FT2D, 2, w_param);

  if(nise) {
    if(!error) {
      error = set_gncomp_1d_array(&CorrFunc, 1, window);
      if(!error) for(i=0; i<window; i++) CorrFunc[i] = 0.0;
    } else CorrFunc = NULL;
    if(!error) {
      error = set_gncomp_1d_array(&NetCorrFunc, 1, window);
      if(!error) for(i=0; i<window; i++) NetCorrFunc[i] = 0.0;
    } else NetCorrFunc = NULL;

    if(!error) {
      error = set_gnreal_1d_array(&popdecay, 1, npopdec);
    } else popdecay = NULL;
    if(!error) {
      error = set_gnreal_1d_array(&popdecay2Q, 1, npopdec);
    } else popdecay2Q = NULL;

    // One-quantum arrays
    // cjfeng 01/03/2019
    // Incorporating trotter_array into allocations
    // NISE Propagators
    // cjfeng 01/04/2019
    // Using U2QMat to replace U2Qs, and eliminate unnecessary allocation of U2QMat using -pert
    if (!trotter && !error) error = set_gncomp_2d_array(&U1QMat, nbuffer, noscsq);
    else U1QMat = NULL;
    if (do2d && !trotter && !pert && !error) error = set_gncomp_2d_array(&U2QMat, nbuffer, n2Qsq);
    else if (do2d && !trotter && pert && !error) error = set_gncomp_2d_array(&U2QMat, 1, n2Qsq);
    else U2QMat = NULL;
    // cjfeng 01/04/2019
    // Trotter propagators
    if (trotter && !error) {
      set_trotter_array(trotter1Q, nosc, nbuffer, 0);
      construct_trotter_index(trotter1Q);
      if (do2d) {
        set_trotter_array(trotter2Q, nosc, nbuffer, 1);
        construct_trotter_2Qindex(trotter2Q, nosc);
      }
    }
    if((!error) && do2d) error = set_gncomp_3d_array(&psi_a, win2d, 3, nosc);
    else psi_a = NULL;
    if((!error) && do2d) error = set_gncomp_3d_array(&psi_b1, win2d, 3, nosc);
    else psi_b1 = NULL;
    if( (!error) && do2d) error = set_gncomp_3d_array(&psi_b12, win2d, 3, nosc);
    else psi_b12 = NULL;
    if( !error && do2d ) for(i=0; i<win2d; i++) for(j=0; j<3; j++) for(k=0; k<nosc; k++) psi_a[i][j][k] = 0.0;
    if( !error && do2d ) for(i=0; i<win2d; i++) for(j=0; j<3; j++) for(k=0; k<nosc; k++) psi_b1[i][j][k] = 0.0;
    if( !error && do2d ) for(i=0; i<win2d; i++) for(j=0; j<3; j++) for(k=0; k<nosc; k++) psi_b12[i][j][k] = 0.0;

    // Two-quantum arrays
    for(i=0; i<9; i++) {
      if( (!error) && (do2d) ) error = set_gncomp_1d_array(&psi_ca[i], 1, n2Q);
      else psi_ca[i] = NULL;
    }
    for(i=0; i<9; i++) {
      if( (!error) && (do2d) ) error = set_gncomp_1d_array(&psi_cb[i], 1, n2Q);
      else psi_cb[i] = NULL;
    }
    if( (!error) && (do2d) ) error = set_gncomp_1d_array(&psi2Q, 1, n2Q);
    else psi2Q = NULL;
    
    // Only used for static averaging / time averaging approximation
    ftir = NULL;
    netftir = NULL;
    for(i=0; i<3; i++) ExDip1QAr[i] = NULL;
    for(i=0; i<3; i++) ExDip2QAr[i] = NULL;
  } 
  else { // Static averaging / time averaging approximation
    if(!error) error = set_gnreal_1d_array(&ftir, 1, npts);
    else ftir = NULL;
    if(!error) for(i=0; i<npts; i++) ftir[i] = 0.0;
    
    if(!error) error = set_gnreal_1d_array(&netftir, 1, npts);
    else netftir = NULL;
    if(!error) for(i=0; i<npts; i++) netftir[i] = 0.0;

    for(i=0; i<3; i++) {
      if(!error) error = set_gnreal_2d_array(&ExDip1QAr[i], nthreads, nosc);
      else ExDip1QAr[i] = NULL;
    }
    for(i=0; i<3; i++) {
      if(do2d && !error) error = set_gnreal_2d_array(&ExDip2QAr[i], nthreads, nosc*n2Q);
      else ExDip2QAr[i] = NULL;
    }

    // Only used for NISE.
    CorrFunc = NULL;
    NetCorrFunc = NULL;
    popdecay = NULL;
    popdecay2Q = NULL;
  
    // One-quantum arrays
    U1QMat = NULL;
    psi_a = NULL;
    psi_b1 = NULL;
    psi_b12 = NULL;
    // Two-quantum arrays
    for(i=0; i<9; i++) psi_ca[i] = NULL;
    for(i=0; i<9; i++) psi_cb[i] = NULL;
    psi2Q = NULL;
    U2QMat = NULL;
  }

  if(!error) error = set_gnreal_1d_array(&hann, 1, window);
  else hann = NULL;

  if(!error) error = set_gnreal_2d_array(&Ham1QMat, nbuffer, noscsq);
  else Ham1QMat = NULL;

  if(!error) error = set_gnreal_3d_array(&Dip1QMat, nbuffer, 3, nosc);
  else Dip1QMat = NULL;

  if(!error) error = set_gnreal_2d_array(&Ham1QAr, nthreads, noscsq);
  else Ham1QAr = NULL;

  if(!error) error = set_gnreal_2d_array(&Evals1QAr, nthreads, nosc);
  else Evals1QAr = NULL;

  for(i=0; i<3; i++) {
    if(!error) error = set_gnreal_2d_array(&Dip1QAr[i], nthreads, nosc);
    else Dip1QAr[i] = NULL;
  }
  
  if(!error) error = set_gnreal_2d_array(&Ham2QAr, nthreads, n2Qsq);
  else Ham2QAr = NULL;

  if(!error) error = set_gnreal_2d_array(&Evals2QAr, nthreads, n2Q);
  else Evals2QAr = NULL;

  for(i=0; i<3; i++) {
    if( (!error) && (do2d) && (!pert) ) set_gnreal_2d_array(&Dip2QAr[i], nthreads, nosc*n2Q);
    else Dip2QAr[i] = NULL;
  }

  if((!error) && (do2d) ) error = set_gnreal_3d_array(&Dip2QMat, nbuffer, 3, nosc*n2Q);
  else Dip2QMat = NULL;

  // cjfeng 12/19/2018
  // Using sparse matrix
  if((!error) && (do2d) ) error = set_gnreal_3d_array(&Dip2QSpMat, nbuffer, 3, nosc*nosc);
  else Dip2QSpMat = NULL;
  if((!error) && (do2d) ) error = set_int_1d_array(&Dip2Qrow, 1, nosc*nosc);
  else Dip2Qrow = NULL;
  if((!error) && (do2d) ) error = set_int_1d_array(&Dip2Qcol, 1, nosc*nosc);
  else Dip2Qcol = NULL;

  // cjfeng 05/07/2018
  if((!error) && (do2d) && (if2Qfiles) ) error = set_gnreal_2d_array(&Ham2QMat, nbuffer, n2Qsq);
  else if((!error) && (do2d) && !(if2Qfiles) ) error = set_gnreal_2d_array(&Ham2QMat, 1, n2Qsq);
  else Ham2QMat = NULL;

  return error; 
}

// Close files, free memory and exit gracefully
// cjfeng 05/07/2018
// Added 2Q file pointers
// cjfeng 12/09/2018
// Using pol_info, spec_param, w_param, fin, and fout
// cjfeng 01/04/2019
// Added trotter_array, and traj_param
int graceful_exit(int error, fin *fin, fout *fout, trotter_array *trotter1Q, trotter_array *trotter2Q, fftw_array *FT1D, fftw_array *FT2D, traj_param *traj_param, w_param *w_param, spec_param *spec_param, POL_info *pol_info, shift_info *shift_info) {
  int i;
  int npol = pol_info->npol;
  int POL[4];
  for(i=0; i<4; i++) POL[i] = pol_info->POL[i];
  int nbuffer = traj_param->nbuffer;
  int nise = spec_param->nise;
  int nthreads = spec_param->nthreads;
  int win2d = w_param->win2d;

  if(error==1) printf("Error allocating memory.\n");
  if(error==2) printf("Error opening files.\n");
  unset_fin(fin);
  unset_fout(fout);

  //GNCOMP arrays
  for(i=0; i<4; i++) unset_gncomp_1d_array(REPH[i]);
  for(i=0; i<4; i++) unset_gncomp_1d_array(NREPH[i]);
  for(i=0; i<4; i++) unset_gncomp_1d_array(NetREPH[i]);
  for(i=0; i<4; i++) unset_gncomp_1d_array(NetNREPH[i]);
  unset_gncomp_2d_array(U1QMat, nbuffer);
  // cjfeng 01/04/2019
  // Using U2QMat to replace U2Qs
  if( !(spec_param->pert) ) unset_gncomp_2d_array(U2QMat, nbuffer);
  else unset_gncomp_2d_array(U2QMat, 1);
  unset_gncomp_1d_array(CorrFunc);
  unset_gncomp_1d_array(NetCorrFunc);
  unset_gncomp_3d_array(psi_a, win2d, 3);
  unset_gncomp_3d_array(psi_b1, win2d, 3);
  unset_gncomp_3d_array(psi_b12, win2d, 3);

  for(i=0; i<9; i++) unset_gncomp_1d_array(psi_ca[i]);
  for(i=0; i<9; i++) unset_gncomp_1d_array(psi_cb[i]);
  unset_gncomp_1d_array(psi2Q);

  //GNREAL arrays
  unset_gnreal_1d_array(hann);
  unset_gnreal_1d_array(popdecay);
  unset_gnreal_1d_array(popdecay2Q);
  
  unset_gnreal_1d_array(ftir);
  unset_gnreal_1d_array(netftir);
  for(i=0; i<3; i++) unset_gnreal_2d_array(Dip1QAr[i], nthreads);
  for(i=0; i<3; i++) unset_gnreal_2d_array(Dip2QAr[i], nthreads);
  for(i=0; i<3; i++) unset_gnreal_2d_array(ExDip1QAr[i], nthreads);
  for(i=0; i<3; i++) unset_gnreal_2d_array(ExDip2QAr[i], nthreads);

  unset_gnreal_3d_array(Dip1QMat, nbuffer, 3);
  unset_gnreal_2d_array(Ham1QMat, nbuffer);

  unset_shift_info_shift(shift_info);

  for(i=0; i<nthreads; i++) {
    if(Ham1QAr!= NULL) if(Ham1QAr[i]!= NULL) free(Ham1QAr[i]);
    if(Evals1QAr!= NULL) if(Evals1QAr[i]!= NULL) free(Evals1QAr[i]);
    if(Ham2QAr!= NULL) if(Ham2QAr[i]!= NULL) free(Ham2QAr[i]);
    if(Evals2QAr!= NULL) if(Evals2QAr[i]!= NULL) free(Evals2QAr[i]);
  }
  if(Ham1QAr!= NULL) free(Ham1QAr);
  if(Evals1QAr!= NULL) free(Evals1QAr);

  if(Ham2QAr!= NULL) free(Ham2QAr);
  if(Evals2QAr!= NULL) free(Evals2QAr);

  // cjfeng 04/05/2016   
  // The modified code due to changing the allocation of Dip2QMat
  unset_gnreal_3d_array(Dip2QMat, nbuffer, 3);

  // cjfeng 12/19/2018
  // Using sparse matrix
  unset_gnreal_3d_array(Dip2QSpMat, nbuffer, 3);
  unset_int_1d_array(Dip2Qrow);
  unset_int_1d_array(Dip2Qcol);

  // cjfeng 05/07/2018
  // Ham2QMat
  if (spec_param->do2d && (spec_param->if2Qfiles) ) unset_gnreal_2d_array(Ham2QMat, nbuffer);
  if (spec_param->do2d && !(spec_param->if2Qfiles) ) unset_gnreal_2d_array(Ham2QMat, 1);

  // cjfeng 01/03/2019
  // FFT arrays
  unset_fftw_array( FT1D );
  if (spec_param->do2d) unset_fftw_array( FT2D );

  // cjfeng 01/04/2019
  // trotter_arrays
  if (spec_param->trotter) {
    unset_trotter_array(trotter1Q);
    if(spec_param->do2d) unset_trotter_array(trotter2Q);
  }
  
  return 0;
}

// cjfeng 07/20/2017
int set_Ham(Ham *ham, int n, int nbuffer) {
  ham->n = n;
  ham->nbuffer = nbuffer;
  set_gnreal_2d_array(&(ham->Mat), ham->nbuffer, ham->n*ham->n);
  return 0;
}

// cjfeng 07/20/2017
int unset_Ham(Ham *ham) {
  unset_gnreal_2d_array(ham->Mat, ham->nbuffer);
  return 0;
}

// cjfeng 07/20/2017
int set_Dip(Dip *dip, int n, int nbuffer) {
  dip->n = n;
  dip->nbuffer = nbuffer;
  set_gnreal_3d_array(&(dip->Mat), dip->nbuffer, 3, dip->n);
  return 0;
}

// cjfeng 07/20/2017
int unset_Dip(Dip *dip) {
  unset_gnreal_3d_array(dip->Mat, dip->nbuffer, 3);
  return 0;  
}

// cjfeng 07/19/2017
int set_eig_array(eig_array *eig, int n) {
  eig->n = n;
  set_gnreal_1d_array(&(eig->Ham), eig->n, eig->n);
  set_gnreal_1d_array(&(eig->Evals), eig->n, 1);
  set_int_1d_array(&(eig->isuppz), eig->n, 2);
  return 0;
}

// cjfeng 07/19/2017
int unset_eig_array(eig_array *eig) {
  unset_int_1d_array(eig->isuppz);
  unset_gnreal_1d_array(eig->Ham);
  unset_gnreal_1d_array(eig->Evals);
  return 0;
}

// cjfeng 04/11/2017
// Subroutines
int set_gnreal_1d_array(GNREAL **mat, int nrows, int ncols) {
  if ( nrows*ncols == 0 ) {
    printf("The size of %s is zero.\n", getName(mat));
    *mat = NULL;
    return 1;
  }
  else {
    *mat = (GNREAL *) malloc(nrows*ncols*sizeof(GNREAL));
    if ( *mat == NULL ) {
      printf("Out of memory for allocating %s.\n", getName(mat));
      return 1;
    }
  }
  return 0;
}

int set_gnreal_2d_array(GNREAL ***mat, int dim1, int dim2) {
  if ( dim1*dim2 == 0) {
    printf("The size of %s is zero.\n", getName(mat));
    *mat = NULL;
    return 1;
  }
  else {
    int i;
    *mat = (GNREAL **) malloc(dim1*sizeof(GNREAL *));
    if ( *mat == NULL ) {
      printf("Out of memory for allocating %s.\n", getName(mat));
      return 1;
    }
    else {
      for(i=0; i<dim1; i++) {
        (*mat)[i] = (GNREAL *) malloc(dim2*sizeof(GNREAL));
        if ( (*mat)[i] == NULL ) {
          int k;
          for(k=0; k<i; k++) free(*mat[k]);
          free(*mat);
          printf("Out of memory for allocating the second dimention of %s[%d].\n", getName(mat), i);
          *mat = NULL;
          return 1;
        }
      }
    }
  }
  return 0;
}

int set_gnreal_3d_array(GNREAL ****mat, int dim1, int dim2, int dim3) {
  if ( dim1*dim2*dim3 == 0) {
    printf("The size of %s is zero.\n", getName(mat));
    *mat = NULL;
    return 1;
  }
  else {
    int i, j;
    *mat = (GNREAL ***) malloc(dim1*sizeof(GNREAL **));
    if ( *mat == NULL) {
      printf("Out of memory for allocating %s.\n", getName(mat));
      return 1;
    }
    else {
      for(i=0; i<dim1; i++) {
        (*mat)[i] = (GNREAL **) malloc(dim2*sizeof(GNREAL *));
        if ( (*mat)[i] == NULL) {
          int k;
          for(k=0; k<i; k++) free(*mat[k]);
          free(*mat);
          printf("Out of memory for allocating %s[%d].\n", getName(mat), i);
          *mat = NULL;
          return 1;
        }
        else {
          for(j=0; j<dim2; j++) {
            (*mat)[i][j] = (GNREAL *) malloc(dim3*sizeof(GNREAL));
            if ( (*mat)[i][j] == NULL) {
              int k, l;
              for(k=0; k<i; k++) {
                for(l=0; l<j; l++) free(*mat[k][l]);
                free(*mat[k]);
              }
              free(*mat);
              printf("Out of memory for allocating %s[%d][%d].\n", getName(mat), i, j);
              *mat = NULL;
              return 1;
            }
          }
        }
      }
    }
  }
  return 0;
}

int unset_gnreal_1d_array(GNREAL *mat) {
  if ( mat != NULL ) free(mat);
  return 0;
}

int unset_gnreal_2d_array(GNREAL **mat, int dim1) {
  if ( mat != NULL ) {
    int i;
    for(i=0; i<dim1; i++) if( mat[i] != NULL ) free(mat[i]);
    free(mat);
  }
  return 0;
}

int unset_gnreal_3d_array(GNREAL ***mat, int dim1, int dim2) {
  if ( mat != NULL ) {
    int i,j;
    for(i=0; i<dim1; i++) {
      if( mat[i] != NULL ) {
        for(j=0; j<dim2; j++) {
          if( mat[i][j] != NULL ) free(mat[i][j]);
        }
        free(mat[i]);
      }
    }
    free(mat);
  }
  return 0;
}

int set_gncomp_1d_array(GNCOMP **mat, int nrows, int ncols) {
  if ( nrows*ncols == 0 ) {
    printf("The size of %s is zero.\n", getName(mat));
    *mat = NULL;
    return 1;
  }
  else {
    *mat = (GNCOMP *) malloc(nrows*ncols*sizeof(GNCOMP));
    if ( *mat == NULL ) {
      printf("Out of memory for allocating %s.\n", getName(mat));
      return 1;
    }
  }
  return 0;
}

int set_gncomp_2d_array(GNCOMP ***mat, int dim1, int dim2) {
  if ( dim1*dim2 == 0) {
    printf("The size of %s is zero.\n", getName(mat));
    *mat = NULL;
    return 1;
  }
  else {
    int i;
    *mat = (GNCOMP **) malloc(dim1*sizeof(GNCOMP *));
    if ( *mat == NULL ) {
      printf("Out of memory for allocating %s.\n", getName(mat));
      return 1;
    }
    else {
      for(i=0; i<dim1; i++) {
        (*mat)[i] = (GNCOMP *) malloc(dim2*sizeof(GNCOMP));
        if ( (*mat)[i] == NULL ) {
          int k;
          for(k=0; k<i; k++) free(*mat[k]);
          free(*mat);
          printf("Out of memory for allocating %s[%d].\n", getName(mat), i);
          *mat = NULL;
          return 1;
        }
      }
    }
  }
  return 0;
}

int set_gncomp_3d_array(GNCOMP ****mat, int dim1, int dim2, int dim3) {
  if ( dim1*dim2*dim3 == 0) {
    printf("The size of %s is zero.\n", getName(mat));
    *mat = NULL;
    return 1;
  }
  else {
    int i, j;
    *mat = (GNCOMP ***) malloc(dim1*sizeof(GNCOMP **));
    if ( *mat == NULL) {
      printf("Out of memory for allocating %s.\n", getName(mat));
      return 1;
    }
    else {
      for(i=0; i<dim1; i++) {
        (*mat)[i] = (GNCOMP **) malloc(dim2*sizeof(GNCOMP *));
        if ( (*mat)[i] == NULL) {
          int k;
          for(k=0; k<i; k++) free(*mat[k]);
          free(*mat);
          printf("Out of memory for allocating %s[%d].\n", getName(mat), i);
          *mat = NULL;
          return 1;
        }
        else {
          for(j=0; j<dim2; j++) {
            (*mat)[i][j] = (GNCOMP *) malloc(dim3*sizeof(GNCOMP));
            if ( (*mat)[i][j] == NULL) {
              int k, l;
              for(k=0; k<i; k++) {
                for(l=0; l<j; l++) free(*mat[k][l]);
                free(*mat[k]);
              }
              free(*mat);
              printf("Out of memory for allocating %s[%d][%d].\n", getName(mat), i, j);
              *mat = NULL;
              return 1;
            }
          }
        }
      }
    }
  }
  return 0;
}

int unset_gncomp_1d_array(GNCOMP *mat) {
  if ( mat != NULL ) free(mat);
  return 0;
}

int unset_gncomp_2d_array(GNCOMP **mat, int dim1) {
  if ( mat != NULL ) {
    int i;
    for(i=0; i<dim1; i++) if( mat[i] != NULL ) free(mat[i]);
    free(mat);
  }
  return 0;
}

int unset_gncomp_3d_array(GNCOMP ***mat, int dim1, int dim2) {
  if ( mat != NULL ) {
    int i,j;
    for(i=0; i<dim1; i++) {
      if( mat[i] != NULL ) {
        for(j=0; j<dim2; j++) {
          if( mat[i][j] != NULL ) free(mat[i][j]);
        }
        free(mat[i]);
      }
    }
    free(mat);
  }
  return 0;
}

int set_int_1d_array(int **mat, int nrows, int ncols) {
  if ( nrows*ncols == 0 ) {
    printf("The size of %s is zero.\n", getName(mat));
    *mat = NULL;
    return 1;
  }
  else {
    *mat = (int *) malloc(nrows*ncols*sizeof(int));
    if ( *mat == NULL ) {
      printf("Out of memory for allocating %s.\n", getName(mat));
      return 1;
    }
  }
  return 0;
}

int unset_int_1d_array(int *mat) {
  if ( mat != NULL ) free(mat);
  return 0;
}

// cjfeng 07/24/2017
int set_char_1d_array(char **mat, int dim1) {
  if ( dim1 == 0 ) {
    printf("The size of %s is zero.\n", getName(mat));
    *mat = NULL;
    return 1;
  }
  else {
    *mat = (char *) malloc(dim1*sizeof(char));
    if ( *mat == NULL ) {
      printf("Out of memory for allocating the %s.\n", getName(mat));
      return 1;
    }
  }
  return 0;
}

// cjfeng 04/25/2017
int set_char_2d_array(char ***mat, int dim1, int dim2) {
  if ( dim1*dim2 == 0 ) {
    printf("The size of %s is zero.\n", getName(mat));
    *mat = NULL;
    return 1;
  }
  else {
    int i;
    *mat = (char **) malloc(dim1*sizeof(char *));
    if ( *mat == NULL ) {
      printf("Out of memory for allocating the %s.\n", getName(mat));
      return 1;
    }
    else {
      for(i=0; i<dim1; i++) {
        (*mat)[i] = (char *) malloc(dim2*sizeof(char));
        if ( (*mat)[i] == NULL) {
          int j;
          for(j=0; j<i; j++) free(*mat[j]);
          free(*mat);
          printf("Out of memory for allocating %s[%d].\n", getName(mat), i);
          *mat = NULL;
          return 1;
        }
      }
    }
  }
  return 0;
}

// cjfeng 07/24/2017
int unset_char_1d_array(char *mat) {
  if ( mat != NULL ) free(mat);
  return 0;
}

int unset_char_2d_array(char **mat, int dim1) {
  if ( mat != NULL ) {
    int i;
    for(i=0; i<dim1; i++) if( mat[i] != NULL ) free(mat[i]);
    free(mat);
  }
  return 0;
}

// cjfeng 12/04/2018
// Using pol_info to replace zzzz, zzyy, zyyz, zyzy, npol, POL[4] and M_ijkl_IJKL.
int initialize_POL_info(POL_info *pol_info) {
  pol_info->zzzz = 1; // ZZZZ polarization
  pol_info->zzyy = 1; // ZZYY polarization
  pol_info->zyyz = 0; // No ZYYZ polarization
  pol_info->zyzy = 0; // No ZYZY polarization
  pol_info->npol = 2; // Number of polarizations to be simulated
  //                             zzzz  zzyy  zyyz  zyzy
  // pol_info->M_ijkl_IJKL = { {  f15, f115, f115, f115 },   // ZZZZ
  //                           { f115, f215, f130, f130 },   // ZZYY
  //                           { f115, f130, f215, f130 },   // ZYYZ
  //                           { f115, f130, f130, f215 } }; // ZYZY
  pol_info->M_ijkl_IJKL[0][0] =  f15; // zzzz->ZZZZ
  pol_info->M_ijkl_IJKL[0][1] = f115; // zzyy->ZZZZ
  pol_info->M_ijkl_IJKL[0][2] = f115; // zyyz->ZZZZ
  pol_info->M_ijkl_IJKL[0][3] = f115; // zyzy->ZZZZ
  pol_info->M_ijkl_IJKL[1][0] = f115; // zzzz->ZZYY
  pol_info->M_ijkl_IJKL[1][1] = f215; // zzyy->ZZYY
  pol_info->M_ijkl_IJKL[1][2] = f130; // zyyz->ZZYY
  pol_info->M_ijkl_IJKL[1][3] = f130; // zyzy->ZZYY
  pol_info->M_ijkl_IJKL[2][0] = f115; // zzzz->ZYYZ
  pol_info->M_ijkl_IJKL[2][1] = f130; // zzyy->ZYYZ
  pol_info->M_ijkl_IJKL[2][2] = f215; // zyyz->ZYYZ
  pol_info->M_ijkl_IJKL[2][3] = f130; // zyzy->ZYYZ
  pol_info->M_ijkl_IJKL[3][0] = f115; // zzzz->ZYZY
  pol_info->M_ijkl_IJKL[3][1] = f130; // zzyy->ZYZY
  pol_info->M_ijkl_IJKL[3][2] = f130; // zyyz->ZYZY
  pol_info->M_ijkl_IJKL[3][3] = f215; // zyzy->ZYZY
  return 0;
}

int set_POL_info(POL_info *pol_info) {
  pol_info->npol = pol_info->zzzz + pol_info->zzyy + pol_info->zyyz + pol_info-> zyzy; // Total number of polarizations to be simualted
  int count = 0;
  if(pol_info->zzzz) {
    pol_info->POL[count] = 0;
    count++;
  }
  if(pol_info->zzyy) {
    pol_info->POL[count] = 1;
    count++;
  }
  if(pol_info->zyyz) {
    pol_info->POL[count] = 2;
    count++;
  }
  if(pol_info->zyzy) {
    pol_info->POL[count] = 3;
    count++;
  }
  return 0;
}

void set_f_pointer(f_pointer *f_pointer, int maxchar, int error) {
  if(!error) {
    f_pointer->fp = NULL;
    error = set_char_1d_array( &(f_pointer->fnm), maxchar );
  }
}

void unset_f_pointer(f_pointer *f_pointer) {
  if( f_pointer->fp != NULL ) fclose(f_pointer->fp);
  unset_char_1d_array(f_pointer->fnm);
}

int set_fin(fin *fin, int maxchar) {
  int error = 0;
  int i;
  fin->maxchar = maxchar;
  set_f_pointer( &(fin->Ham1Q), fin->maxchar, error);
  if(error) return error;
  for(i=0; i<3; i++) {
    set_f_pointer( &(fin->Dip1Q[i]), fin->maxchar, error);
    if(error) return error;
  }
  set_f_pointer( &(fin->Info), fin->maxchar, error);
  if(error) return error;
  set_f_pointer( &(fin->Sites), fin->maxchar, error);
  if(error) return error;
  set_f_pointer( &(fin->Ham2Q), fin->maxchar, error);
  if(error) return error;
  for(i=0; i<3; i++) {
    set_f_pointer( &(fin->Dip2Q[i]), fin->maxchar, error);
    if(error) return error;
  }
  return error;
}

int unset_fin(fin *fin) {
  int i; 
  unset_f_pointer( &(fin->Ham1Q) );
  for(i=0; i<3; i++) unset_f_pointer( &(fin->Dip1Q[i]) );
  unset_f_pointer( &(fin->Sites) );
  unset_f_pointer( &(fin->Ham2Q) );
  for(i=0; i<3; i++) unset_f_pointer( &(fin->Dip2Q[i]) );
  return 0;
}

int set_fout(fout *fout, int maxchar) {
  int error = 0;
  int i;
  fout->maxchar = maxchar;
  set_f_pointer( &(fout->waxis), fout->maxchar, error);
  if(error) return error;
  set_f_pointer( &(fout->labs), fout->maxchar, error);
  if(error) return error;
  set_f_pointer( &(fout->log), fout->maxchar, error);
  if(error) return error;
  for(i=0; i<4; i++) {
    set_f_pointer( &(fout->reph[i]), fout->maxchar, error);
    if(error) return error;
  }
  for(i=0; i<4; i++) {
    set_f_pointer( &(fout->nreph[i]), fout->maxchar, error);
    if(error) return error;
  }
  for(i=0; i<10; i++) {
    set_f_pointer( &(fout->traj[i]), fout->maxchar, error);
    if(error) return error;
  }
  return error;
}

int unset_fout(fout *fout) {
  int i; 
  unset_f_pointer( &(fout->waxis) );
  unset_f_pointer( &(fout->labs) );
  unset_f_pointer( &(fout->log) );
  for(i=0; i<4; i++) unset_f_pointer( &(fout->reph[i]) );
  for(i=0; i<4; i++) unset_f_pointer( &(fout->nreph[i]) );
  for(i=0; i<10; i++) unset_f_pointer( &(fout->traj[i]) );
  return 0;
}

int initialize_res_info(res_info *res_info) {
  res_info->nbonds = 0;
  res_info->nchar = 20;
  res_info->NNums = NULL;
  res_info->CNums = NULL;
  res_info->NNames = NULL;
  res_info->CNames = NULL;
  return 0;
}

int set_res_info(res_info *res_info, int nbonds, int nchar) {
  int error = 0;
  res_info->nbonds = nbonds;
  res_info->nchar = nchar;
  if(!error) error = set_int_1d_array(&(res_info->NNums), res_info->nbonds, 1);
  else return error;
  if(!error) error = set_int_1d_array(&(res_info->CNums), res_info->nbonds, 1);
  else return error;
  if(!error) error = set_char_2d_array(&(res_info->NNames), res_info->nbonds, res_info->nchar);
  else return error;
  if(!error) error = set_char_2d_array(&(res_info->CNames), res_info->nbonds, res_info->nchar);
  else return error;
  return error;
}

int unset_res_info(res_info *res_info) {
  if(res_info->NNums!=NULL) unset_int_1d_array(res_info->NNums);
  if(res_info->CNums!=NULL) unset_int_1d_array(res_info->CNums);
  if(res_info->NNames!=NULL) unset_char_2d_array(res_info->NNames, res_info->nbonds);
  if(res_info->CNames!=NULL) unset_char_2d_array(res_info->CNames, res_info->nbonds);
  return 0;
}

int initialize_shift_info(shift_info *shift_info) {
  shift_info->nshift = 0;
  shift_info->SHIFTNDX = NULL;
  shift_info->SHIFT = NULL;
  return 0;
}

int set_shift_info(shift_info *shift_info) {
  int nshift = shift_info->nshift;
  int error = 0;
  if(!error) error = set_int_1d_array(&(shift_info->SHIFTNDX), nshift, 1);
  else return error;
  if(!error) error = set_gnreal_1d_array(&(shift_info->SHIFT), nshift, 1);
  else return error;

  return error;
}

int unset_shift_info_shift(shift_info *shift_info) {
  if(shift_info->SHIFTNDX!=NULL) unset_int_1d_array(shift_info->SHIFTNDX);
  if(shift_info->SHIFT!=NULL) unset_gnreal_1d_array(shift_info->SHIFT);
  return 0;
}
