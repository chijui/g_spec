#include "param_sub.h"

// cjfeng 05/07/2018
// Added two-quantum Hamiltonian and dipole moments
// Added more descriptions when 2Q Hamiltonian and dipole moments are supplied.
// cjfeng 06/22/2017
// Added Trotter flag.
// cjfeng 12/07/2018
// Using pol_info, spec_param, w_param, and fnm_param
// cjfeng 01/03/2019
// Using gmx_parse_args to warp all the details of argument parsing
int gmx_parse_args(int *argc, char *argv[], fnm_param *fnm_param, spec_param *spec_param, w_param *w_param, POL_info *pol_info, real *delta, output_env_t *oenv) {
  // A list of command line file flags
  t_filenm fnm[] = { };
  // A list of additional command line arguments
  t_pargs pa [] = {
    {"-deffnm", FALSE, etSTR, {&(fnm_param->deffnm)},
     "Input file name base"},
    {"-outname", FALSE, etSTR, {&(fnm_param->outname)}, 
     "Base for output file name." }, 
    {"-ham", FALSE, etSTR, {&(fnm_param->hamnm)},
     "Hamiltonian trajectory file."},
    {"-dipx", FALSE, etSTR, {&(fnm_param->dipxnm)},
     "Dipole moment x-component file."},
    {"-dipy", FALSE, etSTR, {&(fnm_param->dipynm)}, 
     "Dipole moment y-component file."},
    {"-dipz", FALSE, etSTR, {&(fnm_param->dipznm)},
     "Dipole moment z-component file."},
    {"-ham2Q", FALSE, etSTR, {&(fnm_param->ham2Qnm)},
     "Two-quantum Hamiltonian trajectory file."},
    {"-dipx2Q", FALSE, etSTR, {&(fnm_param->dipx2Qnm)},
     "Two-quantum dipole moment x-component file."},
    {"-dipy2Q", FALSE, etSTR, {&(fnm_param->dipy2Qnm)}, 
     "Two-quantum dipole moment y-component file."},
    {"-dipz2Q", FALSE, etSTR, {&(fnm_param->dipz2Qnm)},
     "Two-qunatum dipole moment z-component file."},
    {"-sites", FALSE, etSTR, {&(fnm_param->sitesnm)},
     "Site energy file (replaces diagonal elements in Hamiltonian)."},
    {"-shift", FALSE, etSTR, {&(fnm_param->paramnm)},
     "Isotope shift file (optional). Note that oscillator indexing begins at zero, not one."},
    {"-reph", FALSE, etBOOL, {&(spec_param->reph)},
     "Calculate rephasing spectrum"},
    {"-nreph", FALSE, etBOOL, {&(spec_param->nreph)},
     "Calculate non-rephasing spectrum"},
    {"-zzzz", FALSE, etBOOL, {&(pol_info->zzzz)},
     "Calculate ZZZZ polarization"},
    {"-zzyy", FALSE, etBOOL, {&(pol_info->zzyy)},
     "Calculate ZZYY polarization"},
    {"-zyyz", FALSE, etBOOL, {&(pol_info->zyyz)},
     "Calculate ZYYZ polarization"},
    {"-zyzy", FALSE, etBOOL, {&(pol_info->zyzy)},
     "Calculate ZYZY polarization"},
    {"-2dir", FALSE, etBOOL, {&(spec_param->do2d)},
     "Calculate 2DIR spectrum"},
    {"-nise", FALSE, etBOOL, {&(spec_param->nise)},
     "Static averaging (no dynamics)"},
    {"-pert", FALSE, etBOOL, {&(spec_param->pert)},
     "Use perturbative approximation for energies"},
    {"-skip", FALSE, etINT, {&(spec_param->skip)}, 
     "Number of skips between reference frames"},
    {"-nt", FALSE, etINT, {&(spec_param->nthreads)},
     "Number of threads"},
    {"-delta", FALSE, etREAL, {delta},
     "Anharmonicity (wavenumbers), not used when supplying 2Q trajectories."},
    {"-taup", FALSE, etINT, {&(spec_param->TauP)},
     "Population lifetime (fs)"},
    {"-taup2Q", FALSE, etINT, {&(spec_param->TauP2Q)},
     "Population lifetime of the 2Q states(fs)"},
    {"-tau2", FALSE, etINT, {&(spec_param->T2)},
     "Waiting time (fs). Only for NISE simulations."},
    {"-tstep", FALSE, etINT, {&(spec_param->tstep)},
     "Trajectory time step (fs)."},
    {"-tscan", FALSE, etINT, {&(spec_param->tscan)},
     "Scan time for NISE or averaging time for TAA (fs)"},
    {"-dump", FALSE, etINT, {&(spec_param->dump)},
     "FTIR dump time (ps)"},
    {"-wstart", FALSE, etREAL, {&(w_param->wstart)},
     "Frequency range start"},
    {"-wstop", FALSE, etREAL, {&(w_param->wstop)},
     "Frequency range stop"},
    {"-winterp", FALSE, etREAL, {&(w_param->winterp)},
     "Output frequency spacing. NB: for NISE calculations This determines only zero-pad length, NOT the true resolution (which is determined by the -tscan flag). "},
    {"-whann", FALSE, etBOOL, {&(w_param->whann)},
     "Don't Use Hann window in time domain. In NISE simulations, this windows the response. In static calculations, this weights the averaged Hamiltonians."}, 
    {"-trotter", FALSE, etBOOL, {&(spec_param->trotter)},
     "No trotter expansion on the propagator, only used in NISE method."},
  };

  // The program description
  const char *desc[] = {"NB: 1. OpenMP parallel computation is not yet implemented in static averaging, time-averaging approximation scheme and Trotter expansion."};
  
  // A description of known bugs
  const char *bugs[] = {"1. -dump may have unpredictable behavior.\n"};

  CopyRight(stderr, argv[0]);
  // parse_common_args() does exactly what it sounds like. The arg list is provided in the arrays t_filenm and t_pargs defined above. 
  parse_common_args(argc, argv, PCA_CAN_TIME | PCA_BE_NICE, asize(fnm), fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, oenv);

  return 0;
}

// cjfeng 12/07/2018
// Using spec_param to check all of the parameters
int initialize_spec_param(spec_param *spec_param) {
  spec_param->nise = 0;
  spec_param->trotter = 0;
  spec_param->do2d = 0;
  spec_param->reph = 1;
  spec_param->nreph = 1;
  spec_param->pert = 0;
  spec_param->tstep = 20;
  spec_param->tscan = 0;
  spec_param->T2 = 0;
  spec_param->TauP = 1300;
  spec_param->TauP2Q = 1300;
  spec_param->nthreads = 1;
  spec_param->dump = -1;
  spec_param->do_traj = (spec_param->dump>=1);
  spec_param->skip = 1; 
  spec_param->if2Qfiles = 0;
  return 0;
}

// cjfeng 12/08/2018
// Setting w_param
int initialize_w_param(w_param *w_param) {
  w_param->whann = 0;        // Applying Hann window on time domain response in NISE or average Hamiltonian in TAA.
  w_param->wstart = 1500.0;  // Starting frequency of interest in wavenumber
  w_param->wstop = 1800.0;   // Ending frequency of interest in wavenumber
  w_param->wres = 1;         // Frequency resolution in wavenumber
  w_param->winterp = 1;      // Output frequency spacing in wavenumber, but it determines only zero-padding length for NISE, not the true resolution
  return 0;
}

int set_w_param(w_param *w_param, spec_param *spec_param) {
  const double c = 2.9979e-5; // Speed of light in cm/fsec
  // For NISE FTIR and TAA 2DIR calculations, window is the number
  // of frames involved in each spectral calculation (the correlation
  // function window in FTIR or the time averaging window in TAA).
  // For NISE 2DIR calculations, this is the correlation function
  // window for each dimension; the actual actual number of frames
  // stored is win2d defined below.
  w_param->wres = w_param->winterp;
  w_param->window = ((int) spec_param->tscan/spec_param->tstep ) + 1; // equivalent with floor();
  w_param->win2d = 2*w_param->window + (int) spec_param->T2/spec_param->tstep;
  // Define array dimensions.
  w_param->npts = ceil( (w_param->wstop-w_param->wstart)/w_param->wres) + 1;
  if (spec_param->nise) {
    w_param->winzpad = (int) 1/(w_param->winterp * c * spec_param->tstep) + 1;
  } else w_param->winzpad = 5 * w_param->npts;
  // Define frequency resolution in NISE calculations
  w_param->dw = 1 / (w_param->winzpad*spec_param->tstep*c);
  set_waxis_nise(w_param);
  
  w_param->nprint = (int) ( w_param->wstop - w_param->wstart ) / w_param->dw;

  return 0;
}

int set_waxis_nise(w_param *w_param) {
  GNREAL dw = w_param->dw;
  int ndxstart = -1;
  int maxit = 1e+9;
  GNREAL tol = dw/2;
  GNREAL wdif = dw;
  
  while( (wdif>tol)  && (ndxstart<maxit) ) {
    ndxstart++;
    if( (ndxstart*dw-w_param->wstart)>=0 ) wdif = ndxstart*dw-w_param->wstart;
    else wdif = w_param->wstart-ndxstart*dw;
  }
  if(wdif>tol) printf("Error finding start frequency. Frequency axis should not be trusted!\n");
  w_param->ndxstart = ndxstart;

  return 0;
}

// cjfeng 12/08/2018
// Using fnm_param to construct file names
int initialize_fnm_param(fnm_param *fnm_param) {
  fnm_param->deffnm = NULL;         // Input file name base
  fnm_param->outname = NULL;        // Outname file name base
  fnm_param->hamnm = "ham.txt";     // Hamiltonian file suffix
  fnm_param->dipxnm = "dipx.txt";   // x-component file suffix of the dipole moment
  fnm_param->dipynm = "dipy.txt";   // y-component file suffix of the dipole moment
  fnm_param->dipznm = "dipz.txt";   // z-component file suffix of the dipole moment
  fnm_param->ham2Qnm = "";          // 2Q Hamiltonian
  fnm_param->dipx2Qnm = "";         // x-componet file suffix of the 2Q dipole moment
  fnm_param->dipy2Qnm = "";         // y-componet file suffix of the 2Q dipole moment
  fnm_param->dipz2Qnm = "";         // z-componet file suffix of the 2Q dipole moment
  fnm_param->sitesnm = "";          // Site energy file
  fnm_param->paramnm = "";          // Isotope shift file
  fnm_param->sites2Qnm = "";        // Site energy file for 2Q parts
  fnm_param->polnm[0] = "zzzz.txt"; // zzzz file suffix
  fnm_param->polnm[1] = "zzyy.txt"; // zzyy file suffix
  fnm_param->polnm[2] = "zyyz.txt"; // zyyz file suffix
  fnm_param->polnm[3] = "zyzy.txt"; // zyzy file suffix

  return 0;
}

// cjfeng 01/04/2019
int initialize_traj_param(traj_param *traj_param) {
  traj_param->nosc = 0;
  traj_param->n2Q = 0;
  traj_param->nread = 0;
  traj_param->nbuffer = 0;
  traj_param->nbuffer1d = 0;
  return 0;
}

int set_traj_param(traj_param *traj_param, spec_param *spec_param, w_param *w_param) {
  // The buffer length must be nread-1 frames longer than the required time window
  // so that we can hold nread complete windows in memory at once. 
  // cjfeng 06/29/2016
  // nbuffer definition changes when 2D NISE simulation is performed, resulting in 
  // incorrect scanning of FTIR correlation function. nbuffer1d is declared to
  // correct the indicies when computing CorrFunc[n].
  traj_param->nread = spec_param->nthreads;
  if(spec_param->do2d && spec_param->nise) traj_param->nbuffer = w_param->win2d + traj_param->nread - 1;
  else traj_param->nbuffer = w_param->window + traj_param->nread -1;
  traj_param->nbuffer1d = w_param->window + traj_param->nread -1;
  return 0;
}

int set_BLAS_gemv_opt(BLAS_gemv_opt *gemv_opt) {
  gemv_opt->ta = 'N';
  gemv_opt->inc = 1;
  gemv_opt->alpha = 1;
  gemv_opt->beta = 0.0;
  return 0;
}
