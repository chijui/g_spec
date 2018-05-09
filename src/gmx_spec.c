#include "def.h"
#include "types.h"
#include "count.h"
#include "mem_helper.h"
#include "file_io.h"
#include "mat_util.h"
#include "spec_mod.h"
#include "spec_fft.h"
#include "nise.h"
#include "trotter.h"

// Global variable definition. We need only one file for defining all the variables in the entire program
int pertvec = 0; // first-order perturbative correction on both site energies and vectors

int maxchar = 1024;

int POL[4];
int* SHIFTNDX = NULL;
real* SHIFT = NULL;

// File pointers
// cjfeng 07/20/2017
// Notes about file pointers, planning to use local pointers instead in the future.
// Input files:
// hfp: ham.txt
// Dfp: dip[xyz].txt
// ifp: info.txt
// sfp: sites.txt
// Output files
// afp: waxis.txt
// ffp: ftir.txt
// lfp: log.txt
// rfp: reph_${POL}.txt
// nrfp: nreph_${POL}.txt
FILE *hfp, *Dfp[3], *sfp, *afp, *ffp, *lfp, *rfp[4], *nrfp[4], *ifp;
// File pointers 
// cjfeng 05/07/2018
// hfp2Q: Two-quantum Hamiltonian, ham2Q.txt
// Dfp2Q: Two-quantum dipole moments, dip[xyz]2Q.txt
FILE *hfp2Q, *Dfp2Q[3]; 
// Trajfp is an array of trajectory file pointers: time-stamp, FTIR, rzzzz, rzzyy, rzyyz, rzyzy, nrzzzz, nrzzyy, nrzyyz, nrzyzy
FILE *Trajfp[10];

int *NNums, *CNums;
char **NNames, **CNames;

// Now spectra
GNCOMP *CorrFunc;
GNCOMP *NetCorrFunc;
GNREAL *popdecay;
GNREAL *hann;
GNREAL *ftir;
GNREAL *netftir;
fftw_complex *FTin1D, *FTout1D;
fftw_plan FTplan1D;
GNCOMP *REPH[4], *NREPH[4];
GNCOMP *NetREPH[4], *NetNREPH[4];
fftw_complex *FTin2D, *FTout2D;
fftw_plan FTplan2D;

GNCOMP *U1Q;
GNCOMP *U1Qs;
GNCOMP *cDip1Q[3];

// Generic wavefunctions
GNCOMP **psi_a[3];
GNCOMP **psi_b1[3];
GNCOMP **psi_b12[3];
GNCOMP *psi_cb[9];
GNCOMP *psi_ca[9];
GNCOMP *psi2Q;

GNREAL **Ham1QMat;
GNREAL **Dip1QMat[3];
GNCOMP **U1QMat;

GNREAL **Ham1QAr;
GNREAL **Evals1QAr;
GNREAL **Dip1QAr[3];
GNREAL **ExDip1QAr[3];

GNCOMP *U2Qs;

GNREAL ***Dip2QMat;
GNCOMP **U2QMat;
// cjfeng 05/07/2018
// 2Q trajectories
GNREAL **Ham2QMat;

GNREAL **Ham2QAr;
GNREAL **Evals2QAr;
GNREAL **Dip2QAr[3];
GNREAL **ExDip2QAr[3];


GNREAL *SitesBuffer;

//            zzzz   zzyy   zyyz   zyzy
GNREAL M_ijkl_IJKL[4][4] = { {  f15,  f115,  f115,  f115 },   // ZZZZ
           {  f115,  f215,  f130,  f130 },    // ZZYY
           {  f115,  f130,  f215,  f130 },    // ZYYZ
           {  f115,  f130,  f130,  f215 } };  // ZYZY

/******************************\
 * Main program               *
\******************************/

int main ( int argc, char * argv[] ) {

/***************************************************************
 * Initial variable declaration and applying default settings  *
 ***************************************************************/
  #if OMP_PARALLEL
    double starttime, stoptime;  // Start and stop time during spectral simulation.
  #else
    clock_t starttime, stoptime;
  #endif 
  char* deffnm = NULL;    // Input file name base
  char* outname = NULL;    // Output file name base
  char* hamnm = "ham.txt";  // Hamiltonian file suffix
  char* dipxnm = "dipx.txt";  // x-component file suffix of the dipole moment
  char* dipynm = "dipy.txt";  // y-component file suffix of the dipole moment
  char* dipznm = "dipz.txt";  // z-component file suffix of the dipole moment

  // cjfeng 05/07/2018
  char* ham2Qnm = ""; // 2Q Hamiltonian
  char* dipx2Qnm = ""; // x-componet file suffix of the 2Q dipole moment
  char* dipy2Qnm = ""; // y-componet file suffix of the 2Q dipole moment
  char* dipz2Qnm = ""; // z-componet file suffix of the 2Q dipole moment
  
  char* sitesnm = "";    // Site energy file
  char* paramnm = "";    // Isotope shift file
  // cjfeng 05/07/2018
  char* sites2Qnm = "";  // Site energy file for 2Q parts
  int nise = 0;      // Static averaging by default (0)
  int reph = 1;      // Calculate rephasing spectrum
  int nreph = 1;      // Calculate non-rephasing spectrum
  int zzzz = 1;      // ZZZZ polarization
  int zzyy = 1;      // ZZYY polarization
  int zyyz = 0;      // No ZYYZ polarization
  int zyzy = 0;      // No ZYZY poloarization
  int do2d = 0;      // No 2D IR spectral simulation.
  int pert = 0;      // No first-order perturbative correction on site energies
  // cjfeng 06/22/2017
  int trotter = 0;    // Not using Trotter expansion.
  int nthreads = 1;    // Number of threads
  int nread = 4*nthreads;    // nread for determining buffer length.
  int tscan = 0;       // Scan time for NISE or Averaging window time for TAA in fs
  int window;      // Window size for TAA or linear NISE.
  int win2d;      // Window size for 2D NISE.
  int winzpad;      // Zero-padding length
  int whann = 1;      // Hann window applied on the response in NISE time domain or averaged Hamiltonian in TAA.
  int dump = -1;      // IR spectra dump time in ps
  int tstep = 0;      // Trajectory time step in fs
  int T2 = 0;      // Waiting time in fs, only valid in 2D NISE so far.
  int skip = 1;      // Number of skips between reference frames.
  real delta = 16.0;    // Anharmonicity in wavenumbers (Weak anharmonic approximation applied in the program.)
  int TauP = 1300;    // Amide I lifetime in fs
  
  int nosc = 0;      // Number of oscillators
  // cjfeng 05/07/2018
  int if2Qfiles = 0; // Not loading 2Q trajectories from files

  real wstart = 1500.0;    // Starting frequency of interest in wavenumber.
  real wstop  = 1800.0;    // Ending frequency of interest in wavenumber.
  real wres = 1;      // Frequency resolution in wavenumber.
  real winterp = 1;    // Output frequency spacing in wavenumber, but for NISE, it determines 
                                // only zero-padding length, not the true resolution.
  
  double c = 2.9979e-5;    // Speed of light in cm/fsec
  double pi = 3.14159265;    // The ratio of a circle's circumference
  
  int error = 0;      // Integer indicating error information.
  int nbuffer, nbuffer1d;    // nbuffer determines the size of Hamiltonian and dipole moment array sizes 
          // while nbuffer1d is used only for FTIR correlation function scanning.
  int npol;      // Number of polarizations to be simulated.
  
  // cjfeng 06/27/2016
  // Added usage to track memory usage.
  struct rusage r_usage;    // resource usage, please refer to getrusage manual page.
  
  // cjfeng 11/16/2016
  // Added variables for OpenBLAS usage.
  BLAS_gemv_opt gemv_opt;
  gen_BLAS_gemv_opt(&gemv_opt);
  
  // cjfeng 07/05/2017
  // Trotter array
  trotter_array trotter1Q;
  trotter_array trotter2Q;
  
  // A list of command line file flags
  t_filenm fnm[] = { };

  // A list of additional command line arguments
  t_pargs pa [] = {
    {"-deffnm", FALSE, etSTR, {&deffnm},
     "Input file name base"},
    {"-outname", FALSE, etSTR, {&outname}, 
     "Base for output file name." }, 
    {"-ham", FALSE, etSTR, {&hamnm},
     "Hamiltonian trajectory file."},
    {"-dipx", FALSE, etSTR, {&dipxnm},
     "Dipole moment x-component file."},
    {"-dipy", FALSE, etSTR, {&dipynm}, 
     "Dipole moment y-component file."},
    {"-dipz", FALSE, etSTR, {&dipznm},
     "Dipole moment z-component file."},
    // cjfeng 05/07/2018
    // Added two-quantum Hamiltonian and dipole moments
    {"-ham2Q", FALSE, etSTR, {&ham2Qnm},
     "Two-quantum Hamiltonian trajectory file."},
    {"-dipx2Q", FALSE, etSTR, {&dipx2Qnm},
     "Two-quantum dipole moment x-component file."},
    {"-dipy2Q", FALSE, etSTR, {&dipy2Qnm}, 
     "Two-quantum dipole moment y-component file."},
    {"-dipz2Q", FALSE, etSTR, {&dipz2Qnm},
     "Two-qunatum dipole moment z-component file."},
    {"-sites", FALSE, etSTR, {&sitesnm},
     "Site energy file (replaces diagonal elements in Hamiltonian)."},
    {"-shift", FALSE, etSTR, {&paramnm},
     "Isotope shift file (optional). Note that oscillator indexing begins at zero, not one."},
    {"-reph", FALSE, etBOOL, {&reph},
     "Calculate rephasing spectrum"},
    {"-nreph", FALSE, etBOOL, {&nreph},
     "Calculate non-rephasing spectrum"},
    {"-zzzz", FALSE, etBOOL, {&zzzz},
     "Calculate ZZZZ polarization"},
    {"-zzyy", FALSE, etBOOL, {&zzyy},
     "Calculate ZZYY polarization"},
    {"-zyyz", FALSE, etBOOL, {&zyyz},
     "Calculate ZYYZ polarization"},
    {"-zyzy", FALSE, etBOOL, {&zyzy},
     "Calculate ZYZY polarization"},
    {"-2dir", FALSE, etBOOL, {&do2d},
     "Calculate 2DIR spectrum"},
    {"-nise", FALSE, etBOOL, {&nise},
     "Static averaging (no dynamics)"},
    {"-pert", FALSE, etBOOL, {&pert},
     "Use perturbative approximation for energies"},
    // {"-pertvec", FALSE, etBOOL, {&pertvec},
    // "Use perturbative approximation for vectors"},
    {"-skip", FALSE, etINT, {&skip}, 
     "Number of skips between reference frames"},
    {"-nt", FALSE, etINT, {&nthreads},
     "Number of threads"},
    {"-delta", FALSE, etREAL, {&delta},
    // cjfeng 05/07/2018
    // Added more descriptions when 2Q Hamiltonian and dipole moments are supplied.
     "Anharmonicity (wavenumbers), not used when supplying 2Q trajectories."},
    // "Anharmonicity (wavenumbers)"},
    {"-taup", FALSE, etINT, {&TauP},
     "Population lifetime (fs)"},
    {"-tau2", FALSE, etINT, {&T2},
     "Waiting time (fs). Only for NISE simulations."},
    {"-tstep", FALSE, etINT, {&tstep},
     "Trajectory time step (fs)."},
    {"-tscan", FALSE, etINT, {&tscan},
     "Scan time for NISE or averaging time for TAA (fs)"},
    {"-dump", FALSE, etINT, {&dump},
     "FTIR dump time (ps)"},
    {"-wstart", FALSE, etREAL, {&wstart},
     "Frequency range start"},
    {"-wstop", FALSE, etREAL, {&wstop},
     "Frequency range stop"},
    {"-winterp", FALSE, etREAL, {&winterp},
     "Output frequency spacing. NB: for NISE calculations This determines only zero-pad length, NOT the true resolution (which is determined by the -tscan flag). "},
    {"-whann", FALSE, etBOOL, {&whann},
     "Use Hann window in time domain. In NISE simulations, this windows the response. In static calculations, this weights the averaged Hamiltonians."}, 
    // cjfeng 06/22/2017
    // Added Trotter flag.
    {"-trotter", FALSE, etBOOL, {&trotter},
     "No trotter expansion on the propagator, only used in NISE method."},
  };

  // The program description
  const char *desc[] = {"NB: 1. OpenMP parallel computation is not yet implemented in static averaging, time-averaging approximation scheme and Trotter expansion."};
  
  // A description of known bugs
  const char *bugs[] = {""};
  
  output_env_t oenv;
  
  CopyRight(stderr, argv[0]);
  // parse_common_args() does exactly what it sounds like. The arg list is provided in the arrays t_filenm and t_pargs defined above. 
  parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_BE_NICE, asize(fnm), fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv);
  
  // Vector perturbation is slow.  Discontinue as of 1/31/2016
  //if(pertvec) pert = 1;
  
  // Variable assignments after passing arguments.
  nread = 4 * nthreads;
  npol = zzzz + zzyy + zyyz + zyzy;  // Number of polarizations  
  
  // OpenBLAS environment setup 
  // cjfeng 11/16/2016
  #if USE_OpenBLAS
    #if OMP_PARALLEL 
      openblas_set_num_threads(1);
    #else 
      #if BLAS_subroutines
        openblas_set_num_threads(nthreads);
      #endif
    #endif
  #endif
  // Checking if reasonble simulation conditions are provided.
  error = check_tstep(&tstep, tscan);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    return 0;
  }
  error = check_tscan(tscan);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    return 0;
  }
  error = check_nise_trotter(nise, trotter);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    return 0;
  }
  error = check_T2(T2, nise);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    return 0;
  }
  error = check_2d_spec_output(POL, reph, nreph, npol, zzzz, zzyy, zyyz, zyzy, do2d);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    return 0;
  }
  
/***************************************************************************************
 Hamiltonian and Dipole files: check for consistency and open pointers 
***************************************************************************************/

  int i;
  char ftirnm[maxchar];    // FTIR file name
  char axisnm[maxchar];    // Frequency axis file name
  char lognm[maxchar];    // Log file name
  // cjfeng 10/23/2017
  // Added number of realization file suffix;
  char nrealnm[maxchar]; 
  char rephnm[4][maxchar];  // Rephasing spectra file names
  char nrephnm[4][maxchar];  // Non-rephasing spectra file names
  char Trajnm[10][maxchar];  // Spectral trajectory file names
  char polnm[4][16] = { "zzzz.txt", "zzyy.txt", "zyyz.txt", "zyzy.txt" };
  int n2Q;      // n2Q: number of 2Q states.
  int nframes;      // nframes: number of frames
  // Infoname: info file name, Hamname: Hamiltonian file name, Parname: isotope shift file name, Sitesname: site energy file name
  char Infoname[maxchar], Hamname[maxchar], Dipnames[3][maxchar], Parname[maxchar], Sitesname[maxchar];
  // cjfeng 05/07/2018
  // Added file names for the 2Q parts
  char Ham2Qname[maxchar], Dip2Qnames[3][maxchar];
  // cjfeng 10/23/2017
  // Added number of realization file name string
  char nrealname[maxchar];
  int nshift;      // Number of isotope shifts
  int npts;      // Number of frequency bins in static calculations.
  int do_traj = (dump>=1);  // Dumping spectral data as a spectral trajectory or not
  
  double wo = 0.0;
  double expfac = -2.0*pi*tstep*c;  // Exponential factor 
  
  // Define array dimensions. 
  wres = winterp;
  npts = ceil((wstop-wstart)/wres)+1;  // Used for static calculations.
  
  // For NISE FTIR and TAA 2DIR calculations, window is the number
  // of frames involved in each spectral calculation (the correlation
  // function window in FTIR or the time averaging window in TAA).
  // For NISE 2DIR calculations, this is the correlation function 
  // window for each dimension; the actual actual number of frames 
  // stored is win2d defined below. 
  window = ((int) tscan/tstep) + 1;
  
  // win2d and is used only for NISE 2D calculations
  win2d = 2*window + (int) T2/tstep;
  
  // For static calculations, winzpad is to prevent wraparound error 
  // when dressing with a lifetime-limited lorentzian. 

  // cjfeng 05/08/2018
  if(nise) {
    winzpad = (int) (1)/(winterp*c*tstep) + 1;
    if (winzpad < window) {
      printf("Error: the provided winterp gives shorter spectra array (%d) than requested window size (%d)\n", winzpad, window);
      printf("Please give smaller winterp or smaller tscan.\n");
      graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
      return 0;
    }
  }
  else winzpad = 5*npts;
  
  // The buffer length must be nread-1 frames longer than the required time window
  // so that we can hold nread complete windows in memory at once. 
  // cjfeng 06/29/2016
  // nbuffer definition changes when 2D NISE simulation is performed, resulting in 
  // incorrect scanning of FTIR correlation function. nbuffer1d is declared to
  // correct the indicies when computing CorrFunc[n].
  if(do2d && nise) nbuffer = win2d + nread - 1;
  else nbuffer = window + nread - 1;
  nbuffer1d = window + nread -1;
  
  // Define frequency axis (only actually used for NISE calculations). 
  GNREAL dw = (1.0) / (winzpad*tstep*c);   // Transform resolution in cm-1
  int maxit = 1e+9;
  int ndxstart = -1;
  GNREAL tol = dw/2.0;
  GNREAL wdif = 2*tol;
  while( (wdif>tol)  && (ndxstart<maxit) ) {
    ndxstart++;
    if( (ndxstart*dw-wstart)>=0 ) wdif = ndxstart*dw-wstart;
    else wdif = wstart-ndxstart*dw;
  }  
  if(wdif>tol) printf("Error finding start frequency. Frequency axis should not be trusted!\n");
  int nprint = (int) ( wstop - wstart )/dw;    // Number of frequency bins to be printed.
  
  // Set file names for opening. 
  // cjfeng 05/07/2018
  // Added 2Q parts
  make_fnames_new( Infoname, Hamname, hamnm, Dipnames, dipxnm, dipynm, dipznm, Ham2Qname, ham2Qnm, Dip2Qnames, dipx2Qnm, dipy2Qnm, dipz2Qnm, Sitesname, sitesnm, axisnm, ftirnm, lognm, polnm, rephnm, nrephnm, Trajnm, Parname, paramnm, outname, deffnm, do_traj );
  // make_fnames( Infoname, Hamname, hamnm, Dipnames, dipxnm, dipynm, dipznm, Sitesname, sitesnm, axisnm, ftirnm, lognm, polnm, rephnm, nrephnm, Trajnm, Parname, paramnm, outname, deffnm, do_traj );

  if ( strlen(Ham2Qname)>0 ) if2Qfiles = 1;
  
  // Read info file
  int info  = read_info( ifp, Infoname, &NNames, &NNums, &CNames, &CNums);
  if(info<0) { 
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
    return 0;
  } else if(info>0) {
    printf("Successfully read info file. Will be looking for %d oscillators in input files. \n", info);
  }

  // Read Param file
  if(strlen(paramnm)>0) {
    printf("Reading parameters from file %s.\n", Parname);
    if(!parse_shifts(Parname, &nshift, maxchar)) {
      printf("Error parsing parameter file %s. Please check input.\n", Parname);
      // cjfeng 05/07/2018
      graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
      // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
      return 0;
    }
    printf("Loaded %d oscillator shifts: \n", nshift);
    if(info!=0) for(i=0; i<nshift; i++) printf("%s %d (index %d) :\t%6.2f cm-1.\n", NNames[SHIFTNDX[i]], NNums[SHIFTNDX[i]], SHIFTNDX[i], SHIFT[i]);
    else for(i=0; i<nshift; i++) printf("Oscillator %d: \t%6.2f cm-1.\n", SHIFTNDX[i], SHIFT[i]);
  } 
  else nshift = 0;

  // Check validity of Hamiltonian file.
  error = check_ham(Hamname, &nframes, &nosc, &n2Q, nbuffer, nread, info, window, nthreads);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
    return 0;
  }
  // cjfeng 05/07/2018
  // 2Q Ham
  error = check_ham2Q(Ham2Qname, &nframes, &nosc, &n2Q, nbuffer, nread, info, win2d, nthreads);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
    return 0;
  }

  // Check isotope shift file
  error = check_shift(nshift, SHIFTNDX, nosc);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
    return 0;
  }
  // Checking the consistency between Hamiltonian file and site energy file.
  error = check_sites(Sitesname, Hamname, nosc, nframes);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
    return 0;
  }
  // Checking the consistency between Hamiltonian file and dipole moment files.
  error = check_dipoles(Dipnames, Hamname, nosc, nframes);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
    return 0;
  }
  // cjfeng 05/07/2018
  // Added 2Q parts
  error = check_dip2Q(Dip2Qnames, Ham2Qname, n2Q, nframes);
  if(error) {
    // cjfeng 05/07/2018
    graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
    // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
    return 0;
  }

  // Assigning dump frequency
  if(dump<1) dump = nframes;
  else dump = (int) ( 1000*dump / tstep );  // dump frequency in ps

  // Open files
  // cjfeng 05/07/2018  
  if(!error) error = open_all_new( Hamname, Dipnames, Ham2Qname, Dip2Qnames, Sitesname, axisnm, ftirnm, lognm, npol, POL, rephnm, nrephnm, Trajnm, do2d, do_traj );
  // if(!error) error = open_all( Hamname, Dipnames, Sitesname, axisnm, ftirnm, lognm, npol, POL, rephnm, nrephnm, Trajnm, do2d, do_traj );
  if(!error) printf("Finished opening files.\n");

  // Allocate memory
  if(!pert && do2d) printf("Allocating memory for %d two-quantum states.\n", n2Q);
  // cjfeng 05/07/2018
  // Added Ham2QMat into allocations
  error = allocate_all_new( window, win2d, winzpad, nosc, nbuffer, nthreads, n2Q, pert, do2d, if2Qfiles, tstep, TauP, npts, nise );
  // error = allocate_all( window, win2d, winzpad, nosc, nbuffer, nthreads, n2Q, pert, do2d, tstep, TauP, npts, nise );
  
  // cjfeng 07/05/2017
  // Trotter array
  if(trotter) {
    set_trotter_1Qarray(&trotter1Q, nosc, nbuffer);
    construct_trotter_index(&trotter1Q);
  }
  // cjfeng 07/06/2017
  // Reduce coupling operations
  // if(trotter && do2d) set_trotter_1Qarray(&trotter2Q, n2Q, nbuffer);
  if(trotter && do2d) {
    set_trotter_2Qarray(&trotter2Q, nosc, nbuffer);
    construct_trotter_2Qindex(&trotter2Q, nosc);
  }
  if(!error) printf("Finished allocating memory.\n");

  // Assigning population decay
  // cjfeng 07/05/2017
  if(popdecay != NULL) gen_popdecay(popdecay, tstep, TauP, window, win2d);

  // cjfeng 06/27/2016
  // Tracking memory usage.
  getrusage(RUSAGE_SELF, &r_usage);
  printf("Memory usage = %ld kB.\n", r_usage.ru_maxrss);
  fprintf(lfp, "Memory usage = %ld kB.\n", r_usage.ru_maxrss);

/***************************************************************************************\
 * Starting spectral simulation part
 ***************************************************************************************/
  // And go!
  int alldone = 0;  // Flag for indicating completing all of the calculations.
  int readframe = 0;  // End point of the read frames
  int fr = 0;    // index for referencing reference point during spectral simulation
  int frame = 0;
  int nAvgd = 0;     // Used only for static calculation
  // Setting number of threads for further parallel computation
  #if OMP_PARALLEL 
    omp_set_num_threads(nthreads);
  #endif

  // Measuring start time.
  #if OMP_PARALLEL
    starttime = omp_get_wtime();
  #else
    starttime = clock();
  #endif
  // cjfeng 05/07/2018
  if(!error) printf("Beginning trajectory parsing...\n");
  // Step through trajectory until all frames have been processed
  while( (!error) & (!alldone) ) {
    int justread = 0;    // location of the last read frame
    // Read 1Q data for the next nread frames and store in Ham1QMat and Dip1QMat
    int fr_upper;
    if (nframes >= readframe+nread) fr_upper = readframe+nread;
    else fr_upper = nframes;
    for(fr=readframe; fr<fr_upper; fr++) {
      // cjfeng 06/27/2016
      // Reduce the number of modulus operations.
      // Modulus is required to fill in new frames when exceeding memory allocation size
      // of the arrays, determined by window or win2d.
      int xfr = fr % nbuffer;
      // First reading Hamiltonian
      error = read_ham(hfp, Ham1QMat[xfr], nosc);
      if(error) break;
      // Now dipole moments
      error = read_dip(Dfp, Dip1QMat, nosc, xfr);
      if(error) break;
      // Site energies if specified
      error = read_sites(sfp, SitesBuffer, Ham1QMat, nosc, xfr);
      // cjfeng 05/07/2018
      // 2Q Hamiltonian if specified
      if( do2d && if2Qfiles ) error = read_ham(hfp2Q, Ham2QMat[xfr], n2Q);
      if(error) break;
      // 2Q Dipole moments if specified
      if( do2d && if2Qfiles ) error = read_dip2Q(Dfp2Q, Dip2QMat[xfr], n2Q);
      if(error) break;

      // Add any specified shifts to the Hamiltonian matrix
      for(i=0; i<nshift; i++) Ham1QMat[xfr][SHIFTNDX[i]*nosc+SHIFTNDX[i]] += SHIFT[i]; 
      // Generate two-quantum dipoles
      if( (do2d) ) { 
        // cjfeng 05/07/2018
        // Skip the construction of 2Q dipole moments if loading from files.
        if ( (!if2Qfiles) ) for(i=0; i<3; i++) gen_dip_2Q(Dip1QMat[i][xfr], Dip2QMat[xfr][i], nosc, n2Q);
        // for(i=0; i<3; i++) gen_dip_2Q(Dip1QMat[i][xfr], Dip2QMat[xfr][i], nosc, n2Q);
      }
    }
    // Check if we've reached the end.
    if(fr==nframes) alldone = 1;

    // Note how many new frames have been read.
    justread = fr - readframe;

    // And update the total number of read frames. 
    readframe = fr;
    // cjfeng 07/05/2017
    // Printing frame and memory after reading instead of before response function calculations.
    print_frame(readframe, 100, lfp);
    print_memory(readframe, 100, &r_usage, lfp);

/**************************************
 * Numerical wavefunction propagation *
 **************************************/
    // Now process the new frames. 
    if( (!error) && nise) {
      // Again we use fr as our frame counter. It runs from readframe-justread 
      // (the first of the newly read frames) to readframe-1. 
      
      // cjfeng 06/22/2017
      // Conditioning for Trotter expansion
      if (!trotter) {
        // cjfeng 03/31/2016
        // The parallel section won't be executed under single thread to avoid overhead.
        #if OMP_PARALLEL 
        #pragma omp parallel if(nthreads>1) shared(justread, Ham1QMat, Evals1QAr, U1QMat) private(fr, i) firstprivate(readframe, wo, nosc, expfac)
        #endif
        {
          // cjfeng 07/19/2017
          // Use temporary struct for solving eigenvalue and eigenvectors.
          int noscsq = nosc*nosc;
          eig_array eig;
          set_eig_array(&eig, nosc);
          #if OMP_PARALLEL
          #pragma omp for schedule(guided) nowait
          #endif
          for(fr=readframe-justread; fr<readframe; fr++) {
            // cjfeng 06/27/2016
            // Reduce modulus operations.
            int xfr = fr%nbuffer;
            int tid = 0;
            lapack_int info;
            #if OMP_PARALLEL
              tid = omp_get_thread_num();
            #endif
            // Find one-quantum eigenvalues
            // Note that Ham1QMat[xfr] now contains eigenvectors
            // cjfeng 06/27/2016
            // Relatively Robust Representation to compute eigenvalues and eigenvectors.
            for (i=0; i<noscsq; i++) eig.Ham[i] = Ham1QMat[xfr][i];
            info = SYEVR (LAPACK_ROW_MAJOR, 'V', 'A', 'U', nosc, eig.Ham, nosc, wstart, wstop, 0, nosc-1, 0.00001, &nosc, eig.Evals, eig.Ham, nosc, eig.isuppz);

            // Note that our Hamiltonian is actually H/(h*c)
            // The exponent we calculate is -i*2*pi*tstep*c*Ham1Q
            // cjfeng 04/07/2016
            // Wrapping operations into subroutines
            gen_UMat(U1QMat, eig.Ham, eig.Evals, eig.n, xfr, expfac, wo);
          }
          unset_eig_array(&eig);
        }
      }
      // cjfeng 07/05/2017
      else {
        #if OMP_PARALLEL
        #pragma omp for schedule(guided) nowait
        #endif
        for(fr=readframe-justread; fr<readframe; fr++) {
          int xfr = fr%nbuffer;
          // cjfeng 07/19/2017
          construct_trotter_mat(&trotter1Q, Ham1QMat[xfr], xfr, expfac);
        }
      }
      // cjfeng 06/27/2016
      // Uncouple the 2D eigenvalue procedure from 1Q part.
      if ( do2d && (!pert) ) {  // Generating 2Q Propagator.
        // cjfeng 06/27/2017
        // Adding Trotter expansion
        if(!trotter) {
          #if OMP_PARALLEL 
          #pragma omp parallel if(nthreads>1) shared(justread, Ham1QMat, Evals2QAr, U1QMat, U2QMat) private(fr, i) firstprivate(readframe, wo, n2Q, do2d, pert, expfac, delta, nosc)
          #endif
          {
            // cjfeng 07/19/2017
            // Use temporary struct.
            int n2Qsq = n2Q*n2Q;
            eig_array eig;
            set_eig_array(&eig, n2Q);
            #if OMP_PARALLEL
            #pragma omp for schedule(guided) nowait
            #endif
            for(fr=readframe-justread; fr<readframe; fr++) {
              // cjfeng 06/27/2016
              // Reduce modulus operations.
              int xfr = fr%nbuffer;
              int tid = 0;
              lapack_int info;
              #if OMP_PARALLEL
                tid = omp_get_thread_num();
              #endif
              // cjfeng 05/07/2018
              // Skip the construction of 2Q Ham if loading from external files
              if ( !(if2Qfiles) ) gen_ham_2Q(Ham1QMat[xfr], nosc, eig.Ham, n2Q, delta);
              else { // Copy the hamiltonian from the 2Q trajectory
                for (i=0; i<n2Qsq; i++) eig.Ham[i] = Ham2QMat[xfr][i];  
              }
              // gen_ham_2Q(Ham1QMat[xfr], nosc, eig.Ham, n2Q, delta);
              // printf("Starting eigenvalue calculation for thread %d\n", tid);

              info = SYEVR(LAPACK_ROW_MAJOR, 'V', 'A', 'U', n2Q, eig.Ham, n2Q, wstart, wstop, 0, n2Q-1, 0.00001, &n2Q, eig.Evals, eig.Ham, n2Q, eig.isuppz);
              // printf("Finishing eigenvalue calculation for thread %d\n", tid);
              // cjfeng 04/07/2016
              // Wrapping operations into subroutines
              gen_UMat(U2QMat, eig.Ham, eig.Evals, eig.n, xfr, expfac, wo);
            }
            unset_eig_array(&eig);
          }
        }
        // cjfeng 07/06/2017
        else { // Trotter expansion
          int n2Qsq = n2Q*n2Q;
          GNREAL *Ham2Qtmp;
          set_gnreal_1d_array(&Ham2Qtmp, n2Qsq, 1);
          #if OMP_PARALLEL
          #pragma omp for schedule(guided) nowait
          #endif
          for(fr=readframe-justread; fr<readframe; fr++) {
            int xfr = fr%nbuffer;
            // cjfeng 05/07/2017
            // Skip the construction of 2Q Ham if loading from external files
            if ( !(if2Qfiles) ) gen_ham_2Q(Ham1QMat[xfr], nosc, Ham2Qtmp, n2Q, delta);
            else { // Copy the Hamiltonian from the 2Q trajectory
              for (i=0; i<n2Qsq; i++) Ham2Qtmp[i] = Ham2QMat[xfr][i];
            }
            // cjfeng 07/06/2017
            // Construct for only non-zero couplings
            construct_trotter_mat(&trotter2Q, Ham2Qtmp, xfr, expfac);
          }  
          unset_gnreal_1d_array(Ham2Qtmp);
        }
      }

      // FTIR correlation function computation.
      fr = -1;
      frame = readframe-nread;
      // Now do calculations for selected frames. 
      while( (fr!=-2) && (!error) ) {
        // cjfeng 07/05/2017
        // Use subroutine to count fr.
        int twin = nbuffer1d-nread+1;
        count_fr(&fr, frame, readframe, nread, twin, nbuffer, skip);
        if( (fr>=0) ) {
          // Do a dynamic calculation using fr as our starting point and augment frame.
          frame++;  // We don't reference frame further in the calculation. 
          // cjfeng 04/07/2017
          // Wrapping operations into subroutines.
          // cjfeng 06/27/2017
          // Set up wavefunction for propagation purpose
          // Eventually it will reduce the scaling from O(N^3) to O(N^2).
          GNCOMP **psi_1d[3];
          for(i=0; i<3; i++) set_gncomp_2d_array(&psi_1d[i], 2, nosc);
          // cjfeng 06/27/2017
          // Changed matrix multiplication to matrix vector multiplication.
          // Added Trotter expansion
          if(!trotter) calc_CorrFunc(CorrFunc, U1QMat, psi_1d, Dip1QMat, cDip1Q, nosc, window, fr, nbuffer, expfac, wo, &gemv_opt, nthreads);
          else calc_CorrFunc_trotter(CorrFunc, &trotter1Q, psi_1d, Dip1QMat, window, fr, nthreads);
          for(i=0; i<3; i++) unset_gncomp_2d_array(psi_1d[i], 2);
        }
      }
      // 2D part
      if (do2d) {
        fr = -1;
        frame = readframe-nread;
        // Now do calculations for selected frames. 
        while( (fr!=-2) && (!error) ) {
          // cjfeng 07/05/2017
          // Use subroutine to count fr.
          int twin = nbuffer-nread+1;
          count_fr(&fr, frame, readframe, nread, twin, nbuffer, skip);
          if( (fr>=0) ) {
      
            // Do a dynamic calculation using fr as our starting point and augment frame.
            frame++;  // We don't reference frame further in the calculation. 
          
            // For 2D spectra, we need 1Q propagators from 
            //   t0-->t1
            //   t0-->t1+t2
            //   t0-->t1+t2+t3
            //   t1-->t1+t2
            //   t1-->t1+t2+t3
            //   t1+t2-->t1+t2+t3
            // The 2Q propagator will be needed only from 
            // t1+t2 to t1+t2+t3.
          
            // xfr1, xfr2, and xfr3 will point to the location 
            // of the tau1, tau2, and tau3 frames, respectively. 
            int tau1, tau2, tau3;
            int xfr0, xfr1, xfr2, xfr3;
            tau2 = T2 / tstep;
            xfr0 = fr % nbuffer;
            // Initialize first frame of psi_a[i] array to be simply Dip1Qmat[i][0]
            // cjfeng 06/27/2016
            // Start using xfrs.
            // cjfeng 07/05/2017
            // Initialize and propagate psi_a by a generic propagator
            if(!trotter) propagate_psi1Q(psi_a, Dip1QMat, U1QMat, nosc, fr, 0, xfr0, nbuffer, win2d, &gemv_opt, nthreads);
            else propagate_psi1Q_trotter(psi_a, Dip1QMat, &trotter1Q, 0, fr, xfr0, win2d, nthreads);
            for(tau1=0; tau1<window; tau1++) {
              // cjfeng 06/27/2016
              // Start using xfrs
              xfr1 = (fr+tau1)%nbuffer;
              xfr2 = (fr+tau1+tau2)%nbuffer;
              // cjfeng 07/05/2017
              // Initialize and propagate psi_b1 (psi_b1[i][tau1]) by a generic propagator
              if(!trotter) propagate_psi1Q(psi_b1, Dip1QMat, U1QMat, nosc, fr, tau1, xfr1, nbuffer, window+tau2, &gemv_opt, nthreads);
              else propagate_psi1Q_trotter(psi_b1, Dip1QMat, &trotter1Q, tau1, fr, xfr1, window+tau2, nthreads);
              // cjfeng 07/05/2017
              // Initialize and propagate psi_b12 (psi_b12[i][tau1+tau2]) by generic propagator.
              if(!trotter) propagate_psi1Q(psi_b12, Dip1QMat, U1QMat, nosc, fr, tau1+tau2, xfr2, nbuffer, window, &gemv_opt, nthreads);
              else propagate_psi1Q_trotter(psi_b12, Dip1QMat, &trotter1Q, tau1+tau2, fr, xfr2, window, nthreads);
              // cjfeng 03/23/2017
              // Initialize the 2Q wavefunctions. 
              // psi_ca[i*3+j] = Dip2Q[j][tau1+tau2]*psi_a[i][tau1+tau2]
              // psi_cb[i*3+j] = Dip2Q[j][tau1+tau2]*psi_b1[i][tau1+tau2]
              initialize_psi2Q(psi_a, Dip2QMat, psi_ca, nosc, n2Q, 3, tau1, tau2, nthreads);
              initialize_psi2Q(psi_b1, Dip2QMat, psi_cb, nosc, n2Q, 3, tau1, tau2, nthreads);

              for(tau3=0; tau3<window; tau3++) {
                // cjfeng 06/27/2016
                // Start using xfr3.
                xfr3 = (fr+tau1+tau2+tau3)%nbuffer;

                GNREAL popfac1Q, popfac2Q;
                popfac1Q = popdecay[tau1+2*tau2+tau3];
                popfac2Q = popfac1Q*popdecay[tau3];
                // Computing nonlinear response function
                calc_2dir_nise(psi_a, psi_b1, psi_b12, POL, nosc, n2Q, npol, popfac1Q, popfac2Q, tau1, tau2, tau3, xfr1, xfr2, xfr3, window, REPH, NREPH, nthreads);
                // cjfeng 07/13/2016
                // Propagate 2Q wavefunctions after computing correlation function. 
                // psi_ca[i*3+j] = U2Q[tau3]*psi_ca[i*3+j]
                // psi_cb[i*3+j] = U2Q[tau3]*psi_cb[i*3+j]
                // cjfeng 03/23/2017
                // Propagate the 2Q wavefunctions after computing the response function.
                if(pert) propagate_psi_2Q_pert(U1QMat, U2Qs, psi_ca, psi_cb, psi2Q, nosc, n2Q, delta, expfac, xfr3, nthreads);
                else if (trotter) {
                  propagate_psi2Q_trotter(&trotter2Q, psi_ca, psi2Q, xfr3, nthreads);
                  propagate_psi2Q_trotter(&trotter2Q, psi_cb, psi2Q, xfr3, nthreads);
                } 
                else  propagate_psi_2Q(U2QMat, psi_ca, psi_cb, psi2Q, n2Q, xfr3, &gemv_opt, nthreads);
              }
            }
          }
        }
      }
    }

/*****************************************************
 * Static averaging and time-averaging approximation *
 *****************************************************/ 
    else if ((!error) && !nise) {
      fr = -1;
      frame = readframe-nread;
      // Now do calculations for selected frames. 
      while( (fr!=-2) && (!error) ) {
        int twin = nbuffer-nread+1;
        // cjfeng 07/05/2017
        // Use subroutine to count fr.
        count_fr(&fr, frame, readframe, nread, twin, nbuffer, skip);
        if( (fr>=0) ) {
          // Generate averaged Hamiltonian for a window starting at frame fr 
          // and augment frame. The average Hamiltonian is stored in Ham1QAr. 
          // If Ham1QAr gets filled up, we stop and do a spectral calculation 
          // for all stored frames before proceeding. 
          frame++;  // We don't reference frame further in the calculation. 
          // Generate averaged Hamiltonian. 
          gen_avg_ham(Ham1QMat, Ham1QAr[nAvgd], hann, nosc, fr, whann, window, nbuffer);
          // cjfeng 05/07/2018
          // Added 2Q part.
          gen_avg_ham(Ham2QMat, Ham2QAr[nAvgd], hann, n2Q, fr, whann, window, nbuffer);
          // Generate averaged dipole moments.
          gen_avg_dip(Dip1QMat, Dip2QMat, Dip1QAr, Dip2QAr, nosc, n2Q, fr, nAvgd, do2d, pert, pertvec, whann, window, nbuffer);
          nAvgd++;

          // Declaring local index array for eigensolver
          int *isuppz1Q, *isuppz2Q;
          set_int_1d_array(&isuppz1Q, nosc, 2);
          set_int_1d_array(&isuppz2Q, n2Q, 2);

          // If we've filled up the Ham1QAr array, do a calculation and re-start. 
          if(nAvgd==nthreads) {
            int thr;
            lapack_int info;
            for(thr=0; thr<nthreads; thr++) {
              int tid = 0;
              #if OMP_PARALLEL 
                tid = omp_get_thread_num();
              #endif 
              // If needed, generate 2Q Hamiltonian. This MUST be done before eigenvalue calculation. 
              if( (do2d) && (!pert) ) {
                // cjfeng 05/07/2018
                // Skip the construction of 2Q Hamiltonian when supplying external trajectory.
                if ( !(if2Qfiles) ) gen_ham_2Q(Ham1QAr[tid], nosc, Ham2QAr[tid], n2Q, delta);
              }
              // Find one-quantum eigenvalues
              // Note that Ham1QAr[tid] now contains eigenvectors
              info = SYEVR( LAPACK_ROW_MAJOR, 'V', 'A', 'U', nosc, Ham1QAr[tid], nosc, wstart, wstop, 0, nosc-1, 0.00001, &nosc, Evals1QAr[tid], Ham1QAr[tid], nosc, isuppz1Q);
      
              // Calculate one-quantum dipole moments
              for(i=0; i<3; i++) mvmult_real_serial_trans(Ham1QAr[tid], Dip1QAr[i][tid], ExDip1QAr[i][tid], nosc, nthreads);
              // Calculate FTIR spectrum
              int n, d, ndx;
              GNREAL osc; 
              
              for(n=0; n<nosc; n++) {
                osc = 0.0;
                for(d=0; d<3; d++) osc += ExDip1QAr[d][tid][n]*ExDip1QAr[d][tid][n];
                ndx = floor( (Evals1QAr[tid][n]-wstart)/wres + 0.5 );
                if(ndx<npts && ndx>=0) ftir[ndx] += osc;
              }
              if( (do2d) && (!pert) ) {
                // Two-quantum eigenvalues
                // Note that Ham2QAr[tid] now contains eigenvectors
                info = SYEVR( LAPACK_ROW_MAJOR, 'V', 'A', 'U', n2Q, Ham2QAr[tid], n2Q, wstart, wstop, 0, n2Q-1, 0.00001, &n2Q, Evals2QAr[tid], Ham2QAr[tid], n2Q, isuppz2Q);
                // cjfene 04/11/2017
                gen_ExDip2Q(ExDip2QAr, Dip2QAr, Ham1QAr, Ham2QAr, nosc, n2Q, tid);
                // Now calculate the 2D spectrum. 
                calc_2dir(Evals1QAr, Evals2QAr, ExDip1QAr, ExDip2QAr, tid, nosc, n2Q, npts, wres, wstart, wstop, REPH, NREPH, POL, npol, reph, nreph);
              } else if( (do2d) && (!pertvec) ) {
                // Generate first-order anharmonically perturbed 2Q energies.
                // Remember that Ham1QAr[tid] now contains eigenvectors. 
                gen_perturb_2Q_energies(Ham1QAr[tid], Evals1QAr[tid], Evals2QAr[tid], nosc, delta);
                // And calculate the 2D spectrum. 
                calc_2dir_pert(Evals1QAr, Evals2QAr, ExDip1QAr, tid, nosc, n2Q, npts, wres, wstart, wstop, REPH, NREPH, POL, npol, reph, nreph);
              } else if( (do2d) && (pertvec) ) {
                // Generate first-order anharmonically perturbed 2Q energies.
                // Remember that Ham1QAr[tid] now contains eigenvectors. 
                gen_perturb_2Q_energies(Ham1QAr[tid], Evals1QAr[tid], Evals2QAr[tid], nosc, delta);
                // Generate Eigenvectors
                //gen_perturb_2Q_matrix(Ham1QAr[tid], Evals1QAr[tid], Ham2QAr[tid], Evals2QAr[tid], Tau2QAr[tid], nosc, n2Q, delta);
                gen_perturb_2Q_vectors(Ham1QAr[tid], Evals1QAr[tid], Ham2QAr[tid], Evals2QAr[tid], nosc, n2Q, delta);
                // gen_perturb_2Q_vectors(Ham1QAr[tid], Evals1QAr[tid], Ham2QAr[tid], Evals2QAr[tid], Tau2QAr[tid], nosc, n2Q, delta);
                // And calculate the 2D spectrum. 
                // cjfeng 04/11/2017
                gen_ExDip2Q(ExDip2QAr, Dip2QAr, Ham1QAr, Ham2QAr, nosc, n2Q, tid);
                // Now calculate the 2D spectrum. 
                calc_2dir(Evals1QAr, Evals2QAr, ExDip1QAr, ExDip2QAr, tid, nosc, n2Q, npts, wres, wstart, wstop, REPH, NREPH, POL, npol, reph, nreph);
              }
            }
            nAvgd = 0;
          }
          // Release memory
          unset_int_1d_array(isuppz1Q);
          unset_int_1d_array(isuppz2Q);
        }
      }
    }

    fr = -1;
    frame = readframe-nread;

    // Now do calculations for selected frames. 
    while( (fr!=-2) && (!error) ) {
      // cjfeng 07/05/2017
      // Use subroutine to count fr.
      int twin = nbuffer-nread+1;
      count_fr(&fr, frame, readframe, nread, twin, nbuffer, dump);
      if( (fr>=0) && (frame!=0) ) printf("Dumping at frame %d.\n", frame);
      if( (fr>=0) && (frame!=0) ) fprintf(Trajfp[0], "%6.10f\n", (((GNREAL) tstep)*((GNREAL) frame))/1000.0);
      // Dump calculated spectra in trajectory files
      if( (fr>=0) && (nise) && (frame!=0)) {
        frame++; 
        // First do FTIR
        // cjfeng 03/23/2017
        for(i=0; i<window; i++) NetCorrFunc[i] += CorrFunc[i];
        FFT_ftir(CorrFunc, hann, popdecay, FTplan1D, FTin1D, FTout1D, whann, window, winzpad, c, TauP, wres);
        // cjfeng 03/24/2017
        print_nise_ftir_traj(Trajfp[1], FTout1D, nprint, ndxstart, winzpad);

        // Now 2DIR
        if(do2d) {
          int windowsq = window * window;
          int p;
          // FT and print rephasing spectra
          for(p=0; p<npol; p++) {
            int np = POL[p];
            for(i=0; i<windowsq; i++) NetREPH[np][i] += REPH[np][i];
            // cjfeng 03/24/2017
            FFT_2d(REPH, hann, POL, FTplan2D, FTin2D, FTout2D, p, whann, window, winzpad);
            FFT_flip_w1(FTin2D, FTout2D, winzpad);
            // cjfeng 03/24/2017
            print_nise_2d_traj(Trajfp[2+np], FTout2D, nprint, ndxstart, winzpad);
            for(i=0; i<windowsq; i++) REPH[np][i] = 0.0;
          }
    
          // FT and print non-rephasing spectra
          for(p=0; p<npol; p++) {
            int np = POL[p];
            for(i=0; i<windowsq; i++) NetNREPH[np][i] += NREPH[np][i];
            // cjfeng 03/24/2017
            FFT_2d(NREPH, hann, POL, FTplan2D, FTin2D, FTout2D, p, whann, window, winzpad);
            print_nise_2d_traj(Trajfp[6+np], FTout2D, nprint, ndxstart, winzpad);
            for(i=0; i<windowsq; i++) NREPH[np][i] = 0.0;
          }
        }
      }
      else if( (fr>=0) && (!nise) && (frame!=0) ) { // Static averaging or Time-averaging approximation scheme.
        frame++; 
        // First do FTIR
        // cjfeng 03/24/2017
        for(i=0; i<npts; i++) netftir[i] += ftir[i];
        dress_ftir_lorentzian(ftir, FTplan1D, FTin1D, FTout1D, npts, winzpad, c, TauP, wres);
        print_ftir_traj(Trajfp[1], ftir, npts);
        for(i=0; i<npts; i++) ftir[i] = 0.0;

        // Now 2DIR
        if(do2d) {
          int nptssq = npts * npts;
          int p;
          for(p=0; p<npol; p++) {
            int np = POL[p];
            // First rephasing data
            // cjfeng 04/11/2017
            for(i=0; i<nptssq; i++) NetREPH[np][i] += creal(REPH[np][i]);
            dress_reph_Lorentzian(REPH, POL, FTplan2D, FTin2D, FTout2D, p, npts, winzpad, c, TauP, wres);
            print_2d_traj(Trajfp[2+np], REPH, npts, POL, p);            
            // cjfeng 04/11/2017
            for(i=0; i<nptssq; i++) REPH[np][i] = 0.0;
            // Then non-rephasing data
            // cjfeng 04/11/2017
            for(i=0; i<nptssq; i++) NetNREPH[np][i] += creal(NREPH[np][i]);
            dress_nreph_Lorentzian(NREPH, POL, FTplan2D, FTin2D, FTout2D, p, npts, winzpad, c, TauP, wres);
            print_2d_traj(Trajfp[6+np], NREPH, npts, POL, p);            
            // cjfeng 04/11/2017
            for(i=0; i<nptssq; i++) NREPH[np][i] = 0.0;
          }
        }
      } else if(frame==0) frame++;
    }
  }

// Print the overall spectra  
  if( (!error) && (nise) ) {
    for(i=0; i<window; i++) CorrFunc[i] += NetCorrFunc[i];
    // cjfeng 03/24/2017
    // First, do FTIR.
    FFT_ftir(CorrFunc, hann, popdecay, FTplan1D, FTin1D, FTout1D, whann, window, winzpad, c, TauP, wres);
    // Print FTIR spectrum and frequency axis
    print_nise_ftir(ffp, FTout1D, nprint, ndxstart, winzpad);
    print_nise_waxis(afp, dw, nprint, ndxstart);

    if(do2d) {
      int p;
      int windowsq = window*window;
      // FT and print rephasing spectra
      for(p=0; p<npol; p++) {
        int np = POL[p];
        for(i=0; i<windowsq; i++) REPH[np][i] += NetREPH[np][i];
        FFT_2d(REPH, hann, POL, FTplan2D, FTin2D, FTout2D, p, whann, window, winzpad);
        // For rephasing only: flip spectrum along w1
        FFT_flip_w1(FTin2D, FTout2D, winzpad);
        // Print the rephasing spectrum.
        print_nise_2d(rfp[np], FTout2D, nprint, ndxstart, winzpad);
      }

      // FT and print non-rephasing spectra
      for(p=0; p<npol; p++) {
        int np = POL[p];
        // cjfeng 04/11/2017
        for(i=0; i<windowsq; i++) NREPH[np][i] += NetNREPH[np][i];
        FFT_2d(NREPH, hann, POL, FTplan2D, FTin2D, FTout2D, p, whann, window, winzpad);
        // Print the non-rephasing spectrum.
        print_nise_2d(nrfp[np], FTout2D, nprint, ndxstart, winzpad);
      }
    }
  } else if( (!error) && (!nise) ) { // Static averaging or Time averaging approximation scheme
    int p;
    for(i=0; i<npts; i++) ftir[i] += netftir[i];
    dress_ftir_lorentzian(ftir, FTplan1D, FTin1D, FTout1D, npts, winzpad, c, TauP, wres);

    if(do2d) {
      int nptssq = npts * npts;
      for(p=0; p<npol; p++) {
        int np = POL[p];
        // cjfeng 04/11/2017
        for(i=0; i<nptssq; i++) REPH[np][i] += NetREPH[np][i];
        // Dress rephasing spectrum with Lorentzian profile. Note that it is so far purely real. 
        dress_reph_Lorentzian(REPH, POL, FTplan2D, FTin2D, FTout2D, p, npts, winzpad, c, TauP, wres);
        // cjfeng 04/11/2017
        for(i=0; i<nptssq; i++) NREPH[np][i] += NetNREPH[np][i];
        // Dress non-rephasing spectrum with Lorentzian profile. Note that it is so far purely real. 
        dress_nreph_Lorentzian(NREPH, POL, FTplan2D, FTin2D, FTout2D, p, npts, winzpad, c, TauP, wres);
      }
    }

    // print the frequency axis and ftir spectrum
    print_waxis(afp, wres, npts, wstart);
    print_ftir(ffp, ftir, npts);
  
    // Now print the 2D spectrum. 
    if(do2d) {
      for(i=0; i<npol; i++) {
        int ni = POL[i];
        print_2d(rfp[ni], REPH, npts, POL, i);
        print_2d(nrfp[ni], NREPH, npts, POL, i);
      }
    }
  }
  // Measuring stop time and computing the total time for simulation.
  #if OMP_PARALLEL
    stoptime = omp_get_wtime();
    printf("Time: %2.4f seconds \n", stoptime-starttime);
    fprintf(lfp, "Time: %2.4f seconds \n", stoptime-starttime);
  #else
    stoptime = clock();
    printf("Time: %2.4f seconds \n", (double)(stoptime-starttime) / CLOCKS_PER_SEC );
    fprintf(lfp, "Time: %2.4f seconds \n", (double)(stoptime-starttime) / CLOCKS_PER_SEC );
  #endif
  thanx(stdout);
  thanx(lfp);
  fflush(lfp);

  // cjfeng 05/07/2018
  graceful_exit_new( error, nbuffer, win2d, nthreads, npol, nise, nosc );
  // graceful_exit( error, nbuffer, win2d, nthreads, npol, nise, nosc);
  if(trotter) {
    unset_trotter_array(&trotter1Q);
    if(do2d) unset_trotter_array(&trotter2Q);
  }
  return 0;
}

// cjfeng 05/07/2018
// Remained issues:
// 3. Parallel scaling is not solved. What is the optimal automatic algorithm to determine the best nthreads?
// 4. Parallel computation on static averaging, time-averaging approximation, and Trotter expanson is not implemented.
