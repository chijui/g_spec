#ifndef TYPES_H
#define TYPES_H

typedef struct {
  char ta;
  int inc;
  GNCOMP alpha;
  GNCOMP beta;
} BLAS_gemv_opt;

typedef struct {
  int  nosc;  // nosc*(nosc+1)/2 for 2Q
  int ncoup; // nosc * (nosc-1) /2 for 1Q, nosc^2 * (nosc-1)/2 for nonzero 2Q
  int nbuffer; 
  GNCOMP **D; // nbuffer * nosc
  GNREAL **J; // nbuffer * (2 * ncoup) 
  int *row; // row index, length: ncoup
  int *col; // column index, length: ncoup
} trotter_array;

typedef struct {
  int n;  // nosc or n2Q
  int nbuffer;
  GNREAL **Mat; // nbuffer * (nosc * nosc)
} Ham;

typedef struct {
  int n; // nosc or n2Q
  int nbuffer;
  GNREAL ***Mat; // nbuffer * 3 * nosc
} Dip;

typedef struct {
  int n;
  GNREAL *Ham; // nosc * nosc
  GNREAL *Evals;  // nosc
  int *isuppz; // 2 * nosc
} eig_array;

// cjfeng 12/04/2018
typedef struct {
  int zzzz; // ZZZZ polarization
  int zzyy; // ZZYY polarization
  int zyyz; // ZYYZ polarization
  int zyzy; // ZYZY polarization
  int npol; // Number of polarizations to be simualted
  int POL[4]; // Index table for finding the correct polarization
  GNREAL M_ijkl_IJKL[4][4]; // Orientational tensor
} POL_info;

// cjfeng 12/05/2018
typedef struct {
  int dim;             // Dimension of the Fourier transform (FT)
  int ndpts;           // Number of simulated data points along each dimension
  int nzpad;           // Number of data points used for zero-padding, also used for memory allocation.
  fftw_complex *FTin;  // Input for FT
  fftw_complex *FTout; // Output after FT
  fftw_plan FTplan;    // Plan of FT
} fftw_array;

// cjfeng 12/08/2018
// Using a generic file struct, which will be implemented later.
typedef struct {
  char* fnm; // File name
  FILE *fp; // File pointer
} f_pointer;

// cjfeng 12/08/2018
typedef struct {
  int maxchar;
  f_pointer Ham1Q;    // 1Q Hamiltonian
  f_pointer Dip1Q[3]; // 1Q Dipole moment, xyz
  f_pointer Info;     // Info
  f_pointer Sites;    // Site energy
  f_pointer Ham2Q;    // 2Q Hamiltonian, 
  f_pointer Dip2Q[3]; // 2Q Dipole moment, xyz
} fin;

// cjfeng 12/08/2018
typedef struct {
  int maxchar;
  f_pointer waxis;    // Frequency axis
  f_pointer labs;     // Linear absorption spectrum, e.g. FTIR 
  f_pointer log;      // Log file 
  f_pointer reph[4];  // Rephasing spectra, zzzz, zzyy, zyyz, zyzy. Detailed syntax depends on input flags, e.g. -nozzzz
  f_pointer nreph[4]; // Non-repshaing spectra, zzzz, zzyy, zyyz, and zyzy. Detailed syntax depends on input flags, e.g. -nozzzz
  f_pointer traj[10]; // Spectral trajectory, time-stamp, linear absorption, rzzzz, rzzyy, rzyyz, rzyzy, nrzzzz, nrzzyy, nrzyyz, nrzyzy 
} fout;

// cjfeng 12/07/2018
typedef struct {
  char* deffnm;      // Input file name base
  char* outname;     // Output file name base
  char* hamnm;       // 1Q Hamiltonian file
  char* dipxnm;      // x-component file suffix of the 1Q dipole moment
  char* dipynm;      // y-component file suffix of the 1Q dipole moment
  char* dipznm;      // z-component file suffix of the 1Q dipole moment
  char* ham2Qnm;     // 2Q Hamiltonian file
  char* dipx2Qnm;    // x-component file suffix of the 2Q dipole moment
  char* dipy2Qnm;    // y-component file suffix of the 2Q dipole moment
  char* dipz2Qnm;    // z-component file suffix of the 2Q dipole moment
  char* sitesnm;     // Site energy file
  char* paramnm;     // Isotope frequency shift file
  char* sites2Qnm;   // Site energy file for 2Q manifold.
  char* polnm[4];    // file suffix for various polarizations
} fnm_param;

// cjfeng 01/03/2019
typedef struct {
  int nise;      // Flag for static averaging (0) or NISE (1)
  int trotter;   // Flag for Suzuki-Trotter expansion
  int do2d;      // flag for 2D Spectral simulation.
  int reph;      // Flag for calculating rephasing spectrum
  int nreph;     // Flag for calculating non-rephasing spectrum
  int pert;      // Flag for first-order perturbative correction on site energies in 2Q manifold
  int tstep;     // Time step (fs) of the given trajectory
  int tscan;     // Scan time (fs) for NISE or average window time in Time-averaging approximation (TAA)
  int T2;        // Waiting time (fs), only valid in 2D NISE by far
  int TauP;      // Lifetime modeled by single exponential
  int TauP2Q;    // Lifetime of the 2Q states modeled by single exponential
  int nthreads;  // Number of threads
  int dump;      // Dump time (ps) for printing spectral trajectory
  int do_traj;   // Flog for dumping spectra data as a spectral trajectory
  int skip;      // Number of skips between starting time point for calculating spectra
  int if2Qfiles; // Flag of loading 2Q trajectories from files
} spec_param;

// cjfeng 12/08/2018
typedef struct {
  int whann;    // Flag for applying Hann window on time domain response in NISE or average Hamiltonian in TAA.
  int window;   // Window size for TAA or linear NISE.
  int win2d;    // Window size for TAA or 2D NISE.
  int winzpad;  // Zero-padding length
  int npts;     // Number of data points used for static calculations
  real wstart;  // Starting frequency of interest in wavenumber
  real wstop;   // Ending frequency of interest in wavenumber
  real wres;    // Frequency resolution in wavenumber
  real winterp; // Output frequency spacing in wavenumber, but it determines only zero-padding length for NISE, not the true resolution
  GNREAL dw;    // Transformed frequency resolution in wavenumber, used in NISE.
  int ndxstart; // Index of starting frequency
  int nprint;   // Number of frequency bins to be printed
} w_param;

// cjfeng 12/09/2018
// Parameters used for reading info 
typedef struct {
  int nbonds;    // Number of amide bonds, equivalent with nosc
  int nchar;     // Number of charactors
  int *NNums;    // Residue number of the N side on a given amide bond
  int *CNums;    // Residue number of the C side on a given amide bond
  char **NNames; // Residue name of the N side on a given amide bond
  char **CNames; // Residue name of the C side on a given amide bond
} res_info;

// cjfeng 12/09/2018
// Parameters used for reading site-specific frequency shift
typedef struct {
  int nshift;    // Number of site-specific frequency shift
  int *SHIFTNDX; // Index array for site-specific frequency shift
  GNREAL *SHIFT; // Site-specific frequency shift in wavenumber
} shift_info;

// cjfeng 01/04/2019
// Collecting a set of physcial constants into a struct
typedef struct {
  const double c;        // Speed of light in cm/fsec
  const double pi;       // The ratio of a circle's circumference
  const double sqrt2;    // root 2
  const double sqrt2inv; // 1/(root 2)
  double expfac;         // exponential factor used in NISE propagation and sum-over-state approaches
} phys_const;

// cjfeng 01/04/2019
// Parameters used for describing the information of Hamiltonian and transition dipole moments
typedef struct {
  int nosc;      // Number of oscillators in 1Q manifold
  int n2Q;       // Number of oscillators in 2Q manifold
  int nread;     // Number of frame to be read, used for determining buffer length
  int nbuffer;   // Buffer length, which determines the size of Hamiltonian and dipole moment array sizes
  int nbuffer1d; // Buffer length for two-point correlation function
} traj_param;

// cjfeng 01/08/2019
// Stopwatch for profiling different sub-routines of g_spec
typedef struct {
  double start;   // Start time in second
  double stop;    // Stop time in second
  double elapsed; // Elapsed time in second
  int n;       // number of elapsed measurements, used to track how many time the elapsed time has been added.
} stopwatch;

#endif
