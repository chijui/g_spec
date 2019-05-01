#include "def.h"
#include "types.h"
#include "count.h"
#include "profiler.h"
#include "file_io.h"
#include "param_checker.h"
#include "param_sub.h"
#include "mem_helper.h"
#include "mat_util.h"
#include "spec_mod.h"
#include "FFT_sub.h"
#include "nise.h"
#include "trotter.h"

// Global variable definition. We need only one file for defining all the variables in the entire program
int maxchar = 1024;

// Now spectra
GNCOMP *CorrFunc;
GNCOMP *NetCorrFunc;
GNREAL *popdecay;
GNREAL *popdecay2Q;
GNREAL *hann;
GNREAL *ftir;
GNREAL *netftir;
GNCOMP *REPH[4], *NREPH[4];
GNCOMP *NetREPH[4], *NetNREPH[4];

// Generic wavefunctions
GNCOMP ***psi_a;
GNCOMP ***psi_b1;
GNCOMP ***psi_b12;
GNCOMP *psi_cb[9];
GNCOMP *psi_ca[9];
GNCOMP *psi2Q;

GNREAL **Ham1QMat;
GNREAL ***Dip1QMat;
GNCOMP **U1QMat;

GNREAL **Ham1QAr;
GNREAL **Evals1QAr;
GNREAL **Dip1QAr[3];
GNREAL **ExDip1QAr[3];

GNREAL ***Dip2QMat;
// cjfeng 12/19/2018
// Using sparse matrix instead
GNREAL ***Dip2QSpMat;
int *Dip2Qrow, *Dip2Qcol;
GNCOMP **U2QMat;
// cjfeng 05/07/2018
// 2Q trajectories
GNREAL **Ham2QMat;

GNREAL **Ham2QAr;
GNREAL **Evals2QAr;
GNREAL **Dip2QAr[3];
GNREAL **ExDip2QAr[3];

/******************************\
 * Main program               *
\******************************/

int main ( int argc, char * argv[] ) {
/***************************************************************
 * Initial variable declaration and applying default settings  *
 ***************************************************************/
  // cjfeng 01/08/2019
  // Added clock_master for profiling
  clock_master clock_master;    // Master clock for profiling
  initialize_clock_master(&clock_master);
  count_time(&(clock_master.init.start));
  // cjfeng 12/10/2018
  // Using struct to assemble parameters
  fnm_param fnm_param;          // File name parameters
  spec_param spec_param;        // struct containing all the spectral parameters 
  w_param w_param;              // Frequency/windowing parameters
  traj_param traj_param;        // Trajectory parameters
  POL_info pol_info;            // Polarizations
  res_info res_info;            // Struct used to read info.txt
  shift_info shift_info;        // Struct used to read site-specific frequency shift
  fin fin;                      // Input file names and pointers
  fout fout;                    // Output file names and pointers
  real delta = 16.0;            // Anharmonicity in wavenumbers (Weak anharmonic approximation applied in the program.)
  BLAS_gemv_opt gemv_opt;       // Options for OpenBLAS matrix-vector multiplication
  // Initialize parameters
  initialize_fnm_param(&fnm_param);
  initialize_spec_param(&spec_param);
  initialize_w_param(&w_param);
  initialize_POL_info(&pol_info);
  initialize_res_info(&res_info);
  initialize_shift_info(&shift_info);
  initialize_traj_param(&traj_param);
  set_fin(&fin, maxchar);
  set_fout(&fout, maxchar);
  set_BLAS_gemv_opt(&gemv_opt);

  // Variable declarations
  int error = 0;                // Integer indicating error information.
  struct rusage r_usage;        // resource usage, used to track memory usage, please refer to getrusage manual page.
  
  // cjfeng 01/03/2019
  // Lumping parsing arguments into a single sub-routine
  output_env_t oenv;
  gmx_parse_args(&argc, argv, &fnm_param, &spec_param, &w_param, &pol_info, &delta, &oenv);

  // Variable assignments after passing arguments.
  set_POL_info(&pol_info);
  
  #if OMP_PARALLEL
    fftw_plan_with_nthreads(spec_param.nthreads);
  #endif
  fftw_array FT1D;
  fftw_array FT2D;

  // cjfeng 01/03/2019
  initialize_fftw_array_pointer(&FT1D);
  initialize_fftw_array_pointer(&FT2D);

  // cjfeng 07/05/2017
  // Trotter array
  trotter_array trotter1Q;      // 1Q propagator
  trotter_array trotter2Q;      // 2Q propagator
  initialize_trotter_array(&trotter1Q);
  initialize_trotter_array(&trotter2Q);
  
  // OpenBLAS environment setup 
  // cjfeng 11/16/2016
  #if USE_OpenBLAS
    #if OMP_PARALLEL 
      openblas_set_num_threads(1);
    #else 
      #if BLAS_subroutines
        openblas_set_num_threads(spec_param.nthreads);
      #endif
    #endif
  #endif

  error = check_spec_param(&spec_param, &pol_info);

  // cjfeng 01/04/2019
  // Added trotter_array, and traj_param
  if (error) {
    graceful_exit(error, &fin, &fout, &trotter1Q, &trotter2Q, &FT1D, &FT2D, &traj_param, &w_param, &spec_param, &pol_info, &shift_info );
    return 0;
  }
  
/***************************************************************************************
 Hamiltonian and Dipole files: check for consistency and open pointers 
***************************************************************************************/
  int nframes;  // nframes: number of frames
  char Parname[maxchar];

  spec_param.do_traj = (spec_param.dump>=1);  // Dumping spectral data as a spectral trajectory or not
  set_w_param(&w_param, &spec_param); 
  error = check_w_param(&w_param);
  if (error) {
    graceful_exit(error, &fin, &fout, &trotter1Q, &trotter2Q, &FT1D, &FT2D, &traj_param, &w_param, &spec_param, &pol_info, &shift_info );
    return 0;
  }

  set_traj_param(&traj_param, &spec_param, &w_param);
  
  // Set file names for opening. 
  make_fnames(&fin, &fout, &fnm_param, Parname, spec_param.do_traj );
  if ( strlen(fin.Ham2Q.fnm)>0 ) spec_param.if2Qfiles = 1;
  
  // Read info file
  int nbonds = read_info( &fin, &res_info, maxchar);
  error = check_nbonds(nbonds);
  if(error) {
    graceful_exit(error, &fin, &fout, &trotter1Q, &trotter2Q, &FT1D, &FT2D, &traj_param, &w_param, &spec_param, &pol_info, &shift_info );
    return 0;
  }

  // Read Param file
  error = read_param(&fnm_param, Parname, &res_info, &shift_info, maxchar);
  if(error) {
    graceful_exit(error, &fin, &fout, &trotter1Q, &trotter2Q, &FT1D, &FT2D, &traj_param, &w_param, &spec_param, &pol_info, &shift_info );
    return 0;
  }

  error = check_input_files(&fin, &spec_param, &w_param, &shift_info, &nframes, &traj_param, nbonds);
  if(error) {
    graceful_exit(error, &fin, &fout, &trotter1Q, &trotter2Q, &FT1D, &FT2D, &traj_param, &w_param, &spec_param, &pol_info, &shift_info );
    return 0;
  }

  // Assigning dump frequency
  if(spec_param.dump<1) spec_param.dump = nframes;
  else spec_param.dump = (int) ( 1000*spec_param.dump / spec_param.tstep );  // dump frequency in ps
  
  // Open files
  if(!error) error = open_all( &fin, &fout, &pol_info, &spec_param );
  if(!error) printf("Finished opening files.\n");

  // Allocate memory
  if(!spec_param.pert && spec_param.do2d) printf("Allocating memory for %d two-quantum states.\n", traj_param.n2Q);
  error = allocate_all(&trotter1Q, &trotter2Q, &FT1D, &FT2D, &w_param, &spec_param, &traj_param );
  if(error) {
    graceful_exit(error, &fin, &fout, &trotter1Q, &trotter2Q, &FT1D, &FT2D, &traj_param, &w_param, &spec_param, &pol_info, &shift_info );
    return 0;
  }
  else printf("Finished allocating memory.\n");

  // cjfeng 01/04/2019
  // Assign physical constant used in this program
  // Use phys_const, in order of c, pi, sqrt(2), 1/sqrt(2), and expfac
  // In NISE, expfac = -2 * Const.pi * Const.c * spec_param.tstep, while in sum-over-state schemes,
  // expfac = -1/(2 * spec_param.TauP * Const.c * w_param.wres * w_param.winzpad)
  phys_const Const = {2.9979e-5, 3.14159265, 1.41421356237, 0.707106781186547}; 
  if (spec_param.nise) Const.expfac = -2.0 * Const.pi * Const.c * (double) spec_param.tstep; 
  else Const.expfac = -1/(2.0 * (double) spec_param.TauP * Const.c * (double) w_param.wres * (double) w_param.winzpad); 

  // Assigning population decay and hann window
  if ( popdecay != NULL ) gen_popdecay(popdecay, &spec_param, &w_param, 0);
  if ( (spec_param.do2d) && (popdecay != NULL) ) gen_popdecay(popdecay2Q, &spec_param, &w_param, 1);
  if(hann != NULL) gen_hann(hann, w_param.window, &Const);
  count_time(&(clock_master.init.stop));
  calc_elapsed_time(&(clock_master.init));

/***************************************************************************************\
 * Starting spectral simulation part
 ***************************************************************************************/
  // And go!
  int alldone = 0;   // Flag for indicating completing all of the calculations.
  int readframe = 0; // End point of the read frames
  int fr = 0;        // index for referencing reference point during spectral simulation
  int frame = 0;
  int nAvgd = 0;     // Used only for static calculation
  // Setting number of threads for further parallel computation
  #if OMP_PARALLEL 
    omp_set_num_threads(spec_param.nthreads);
  #endif

  // Measuring start time.
  count_time(&(clock_master.run.start));
  // Print frequency axis first instead
  if(spec_param.nise) print_nise_waxis(fout.waxis.fp, &w_param);
  else print_waxis(fout.waxis.fp, &w_param);
  if(!error) printf("Beginning trajectory parsing...\n");
  // cjfeng 06/27/2016
  // Tracking memory usage.
  print_frame(readframe, 100, fout.log.fp);
  print_memory(readframe, 100, &r_usage, fout.log.fp);
  // Step through trajectory until all frames have been processed
  // cjfeng 12/08/2018
  // Using fin
  while( (!error) & (!alldone) ) {
    int justread = 0;    // location of the last read frame
    // Read 1Q data for the next nread frames and store in Ham1QMat and Dip1QMat
    int fr_upper;
    if (nframes >= readframe+traj_param.nread) fr_upper = readframe+traj_param.nread;
    else fr_upper = nframes;
    for(fr=readframe; fr<fr_upper; fr++) {
      // cjfeng 01/09/2019
      // Count read time
      count_time(&(clock_master.read.start));
      // Modulus is required to fill in new frames when exceeding memory allocation size
      // of the arrays, determined by window or win2d.
      int xfr = fr % traj_param.nbuffer;

      // First reading Hamiltonian
      error = read_ham(fin.Ham1Q.fp, Ham1QMat[xfr], traj_param.nosc);
      if(error) break;

      // Now dipole moments
      error = read_dip(fin.Dip1Q, Dip1QMat[xfr], traj_param.nosc);
      if(error) break;

      // Site energies if specified
      error = read_sites(fin.Sites.fp, Ham1QMat[xfr], traj_param.nosc);

      // cjfeng 05/07/2018
      // 2Q Hamiltonian if specified
      if( spec_param.do2d && spec_param.if2Qfiles ) error = read_ham(fin.Ham2Q.fp, Ham2QMat[xfr], traj_param.n2Q);
      if(error) break;

      // 2Q Dipole moments if specified
      if( spec_param.do2d && spec_param.if2Qfiles ) error = read_dip(fin.Dip2Q, Dip2QMat[xfr], traj_param.n2Q);
      if(error) break;
      count_time(&(clock_master.read.stop));
      calc_elapsed_time(&(clock_master.read));

      // Add any specified shifts to the Hamiltonian matrix
      // cjfeng 12/09/2018
      // Using subroutine to apply site-specific frequency shift
      apply_shift(Ham1QMat[xfr], &shift_info, traj_param.nosc);
      // Generate two-quantum dipoles
      // cjfeng 05/07/2018
      // Skip the construction of 2Q dipole moments if loading from files.
      // cjfeng 12/19/2018
      // Using sparse matrix for NISE calculations
      // cjfeng 01/09/2019
      count_time(&(clock_master.gen_Dip2Q.start));
      if( (spec_param.do2d) ) {
        if (spec_param.nise && !spec_param.if2Qfiles) {
          gen_dip_Sp2Q_ind(Dip2Qrow, Dip2Qcol, &traj_param);
          gen_dip_Sp2Q(Dip1QMat[xfr], Dip2QSpMat[xfr], &traj_param);
        }
        else gen_dip_2Q(Dip1QMat[xfr], Dip2QMat[xfr], &traj_param);
      }
      count_time(&(clock_master.gen_Dip2Q.stop));
      calc_elapsed_time(&(clock_master.gen_Dip2Q));
    }
    // Check if we've reached the end.
    if(fr==nframes) alldone = 1;

    // Note how many new frames have been read.
    justread = fr - readframe;

    // And update the total number of read frames. 
    readframe = fr;
    // cjfeng 07/05/2017
    // Printing frame and memory after reading instead of before response function calculations.
    // cjfeng 12/09/2018
    print_frame(readframe, 100, fout.log.fp);
    print_memory(readframe, 100, &r_usage, fout.log.fp);

/**************************************
 * Numerical wavefunction propagation *
 **************************************/
    // Now process the new frames. 
    if( (!error) && spec_param.nise) {
      // Again we use fr as our frame counter. It runs from readframe-justread 
      // (the first of the newly read frames) to readframe-1. 
      
      // cjfeng 06/22/2017
      // Implemented Trotter expansion
      // cjfeng 12/12/2018
      // Using gen_U1Q to wrap the code
      // cjfeng 01/04/2019
      // Added traj_param
      // cjfeng 01/09/2019
      count_time(&(clock_master.gen_U1Q.start));
      if (!spec_param.trotter) gen_U1Q(Ham1QMat, Evals1QAr, U1QMat, &spec_param, &w_param, &traj_param, justread, readframe, Const.expfac);
      else {
        #if OMP_PARALLEL
        #pragma omp for schedule(guided) nowait
        #endif
        for(fr=readframe-justread; fr<readframe; fr++) {
          int xfr = fr%traj_param.nbuffer;
          construct_trotter_mat(&trotter1Q, Ham1QMat[xfr], xfr, Const.expfac);
        }
      }
      count_time(&(clock_master.gen_U1Q.stop));
      calc_elapsed_time(&(clock_master.gen_U1Q));
      // cjfeng 06/27/2016
      // Uncouple the 2D eigenvalue procedure from 1Q part.
      // cjfeng 06/22/2017
      // Implemented Trotter expansion
      if ( spec_param.do2d && !(spec_param.pert) ) {  // Generating 2Q Propagator.
        count_time(&(clock_master.gen_U2Q.start));
        if(!spec_param.trotter) gen_U2Q(Ham1QMat, Ham2QMat, Evals2QAr, U2QMat, &spec_param, &w_param, &traj_param, justread, readframe, delta, Const.expfac);
        else {
          int n2Qsq = traj_param.n2Q*traj_param.n2Q;
          #if OMP_PARALLEL
          #pragma omp for schedule(guided) nowait
          #endif
          for(fr=readframe-justread; fr<readframe; fr++) {
            int xfr = fr%traj_param.nbuffer;
            // cjfeng 05/07/2017
            // Skip the construction of 2Q Ham if loading from external files
            if ( !(spec_param.if2Qfiles) ) {
              gen_ham_2Q(Ham1QMat[xfr], Ham2QMat[0], &traj_param, delta);
              construct_trotter_mat(&trotter2Q, Ham2QMat[0], xfr, Const.expfac);
            }
            // Copy the Hamiltonian from the 2Q trajectory
            else {
              construct_trotter_mat(&trotter2Q, Ham2QMat[xfr], xfr, Const.expfac);
            }
          }  
        }
        count_time(&(clock_master.gen_U2Q.stop));
        calc_elapsed_time(&(clock_master.gen_U2Q));
      }

      // FTIR correlation function computation.
      fr = -1;
      frame = readframe-traj_param.nread;
      // Now do calculations for selected frames. 
      while( (fr!=-2) && (!error) ) {
        // cjfeng 07/05/2017
        // Use subroutine to count fr.
        int twin = w_param.window;
        // cjfeng 12/07/2018
        // Using spec_param
        count_fr(&fr, frame, readframe, traj_param.nread, twin, traj_param.nbuffer, spec_param.skip);
        if( (fr>=0) ) {
          // cjfeng 01/09/2019
          count_time(&(clock_master.calc_R1.start));
          // Do a dynamic calculation using fr as our starting point and augment frame.
          frame++;  // We don't reference frame further in the calculation. 
          // cjfeng 04/07/2017
          // Wrapping operations into subroutines.
          // cjfeng 06/27/2017
          // Set up wavefunction for propagation purpose
          // Eventually it will reduce the scaling from O(N^3) to O(N^2).
          GNCOMP ***psi_1d;
          set_gncomp_3d_array(&psi_1d, 2, 3, traj_param.nosc);
          // cjfeng 06/27/2017
          // Changed matrix multiplication to matrix vector multiplication.
          // Added Trotter expansion
          // cjfeng 12/08/2018
          // Using spec_param and w_param
          if(!spec_param.trotter) calc_CorrFunc(CorrFunc, U1QMat, psi_1d, Dip1QMat, &traj_param, w_param.window, fr, Const.expfac, &gemv_opt, spec_param.nthreads);
          else calc_CorrFunc_trotter(CorrFunc, &trotter1Q, psi_1d, Dip1QMat, w_param.window, fr, spec_param.nthreads);
          unset_gncomp_3d_array(psi_1d, 2, 3);
          count_time(&(clock_master.calc_R1.stop));
          calc_elapsed_time(&(clock_master.calc_R1));
        }
      }
      // 2D part
      if (spec_param.do2d) {
        fr = -1;
        frame = readframe-traj_param.nread;
        // Now do calculations for selected frames. 
        while( (fr!=-2) && (!error) ) {
          // cjfeng 07/05/2017
          // Use subroutine to count fr.
          int twin = w_param.win2d;
          // cjfeng 12/07/2018
          // Using spec_param
          count_fr(&fr, frame, readframe, traj_param.nread, twin, traj_param.nbuffer, spec_param.skip);
          // cjfeng 06/27/2016
          // Start using xfrs
          // cjfeng 12/09/2018
          // Using pol_info, spec_param, and w_param
          if( (fr>=0) ) {
            // cjfeng 01/05/2019
            count_time(&(clock_master.calc_R3.start));
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
            //   t1+t2 to t1+t2+t3.
          
            // xfr1, xfr2, and xfr3 will point to the location 
            // of the tau1, tau2, and tau3 frames, respectively. 
            int tau1, tau2, tau3;
            int xfr0, xfr1, xfr2, xfr3;
            tau2 = spec_param.T2 / spec_param.tstep;
            xfr0 = fr % traj_param.nbuffer;
            // Initialize first frame of psi_a[0][i] array to be simply Dip1Qmat[0][i]
            // cjfeng 07/05/2017
            // Initialize and propagate psi_a by a generic propagator
            // cjfeng 01/04/2019
            // Added traj_param
            if(!spec_param.trotter) propagate_psi1Q(psi_a, Dip1QMat, U1QMat, &traj_param, fr, 0, xfr0, w_param.win2d, &gemv_opt, spec_param.nthreads);
            else propagate_psi1Q_trotter(psi_a, Dip1QMat, &trotter1Q, 0, fr, xfr0, w_param.win2d, spec_param.nthreads);
            int j, win2d = w_param.win2d;
            for(tau1=0; tau1<w_param.window; tau1++) {
              xfr1 = (fr+tau1)%traj_param.nbuffer;
              xfr2 = (fr+tau1+tau2)%traj_param.nbuffer;
              // cjfeng 07/05/2017
              // Initialize and propagate psi_b1 (psi_b1[i][tau1]) by a generic propagator
              if(!spec_param.trotter) propagate_psi1Q(psi_b1, Dip1QMat, U1QMat, &traj_param, fr, tau1, xfr1, w_param.window+tau2, &gemv_opt, spec_param.nthreads);
              else propagate_psi1Q_trotter(psi_b1, Dip1QMat, &trotter1Q, tau1, fr, xfr1, w_param.window+tau2, spec_param.nthreads);

              // cjfeng 07/05/2017
              // Initialize and propagate psi_b12 (psi_b12[i][tau1+tau2]) by generic propagator.
              if(!spec_param.trotter) propagate_psi1Q(psi_b12, Dip1QMat, U1QMat, &traj_param, fr, tau1+tau2, xfr2, w_param.window, &gemv_opt, spec_param.nthreads);
              else propagate_psi1Q_trotter(psi_b12, Dip1QMat, &trotter1Q, tau1+tau2, fr, xfr2, w_param.window, spec_param.nthreads);

              // cjfeng 03/23/2017
              // psi_ca[i*3+j] = Dip2Q[j][tau1+tau2]*psi_a[i][tau1+tau2]
              // psi_cb[i*3+j] = Dip2Q[j][tau1+tau2]*psi_b1[i][tau1+tau2]
              int i;
              // cjfeng 12/19/2018
              // Using sparse matrix
              // cjfeng 05/01/2018
              // Fixing the index of Dip2QMat and Dip2QSpMat from tau1+tau2 to xfr2
              if(!spec_param.if2Qfiles) {
                initialize_psi2Q_Sp(psi_a[tau1+tau2], Dip2QSpMat[xfr2], Dip2Qrow, Dip2Qcol, psi_ca, &traj_param, spec_param.nthreads);
                initialize_psi2Q_Sp(psi_b1[tau1+tau2], Dip2QSpMat[xfr2], Dip2Qrow, Dip2Qcol, psi_cb, &traj_param, spec_param.nthreads);
              }
              else {
                initialize_psi2Q(psi_a[tau1+tau2], Dip2QMat[xfr2], psi_ca, &traj_param, spec_param.nthreads);
                initialize_psi2Q(psi_b1[tau1+tau2], Dip2QMat[xfr2], psi_cb, &traj_param, spec_param.nthreads);
              }
              for(tau3=0; tau3<w_param.window; tau3++) {
                xfr3 = (fr+tau1+tau2+tau3)%traj_param.nbuffer;

                GNREAL popfac1Q, popfac2Q;
                popfac1Q = popdecay[tau1+2*tau2+tau3];
                // cjfeng 01/10/2019
                // Liang, C.; Jansen, T. L. C. J. Chem. Theory Comput. 2012, 8 (5), 1706–1713.
                // popfac2Q = popfac1Q;
                // Jansen, T. L. C.; Knoester, J. J. Phys. Chem. B 2006, 110, 22910–22916. 
                // popfac2Q = popfac1Q*popdecay[tau3];
                popfac2Q = popdecay[tau1+2*tau2] * popdecay2Q[tau3];
                // Computing third-order response function
                if (spec_param.if2Qfiles) calc_2dir_nise(psi_a, psi_b1, psi_b12, &pol_info, &traj_param, popfac1Q, popfac2Q, tau1, tau2, tau3, xfr1, xfr2, xfr3, w_param.window, REPH, NREPH, &spec_param, spec_param.nthreads);
                else calc_2dir_nise_Sp(psi_a, psi_b1, psi_b12, &pol_info, &traj_param, popfac1Q, popfac2Q, tau1, tau2, tau3, xfr1, xfr2, xfr3, w_param.window, REPH, NREPH, &spec_param, spec_param.nthreads);
                // cjfeng 07/13/2016
                // Propagate 2Q wavefunctions after computing correlation function. 
                // psi_ca[i*3+j] = U2Q[tau3]*psi_ca[i*3+j] // psi_cb[i*3+j] = U2Q[tau3]*psi_cb[i*3+j]
                // cjfeng 03/23/2017
                // Propagate the 2Q wavefunctions after computing the response function.
                // cjfeng 01/23/2019
                // Using single propagation function for -nise and -pert
                if (spec_param.trotter) {
                  propagate_psi2Q_trotter(&trotter2Q, psi_ca, xfr3, spec_param.nthreads);
                  propagate_psi2Q_trotter(&trotter2Q, psi_cb, xfr3, spec_param.nthreads);
                } 
                else {
                  if (spec_param.pert) {
                    gen_pert_prop(U1QMat[xfr3], U2QMat[0], &traj_param, Const.expfac, delta);
                    propagate_psi2Q(U2QMat, psi_ca, psi2Q, traj_param.n2Q, 0, &gemv_opt, spec_param.nthreads);
                    propagate_psi2Q(U2QMat, psi_cb, psi2Q, traj_param.n2Q, 0, &gemv_opt, spec_param.nthreads);
                  } 
                  else {
                    propagate_psi2Q(U2QMat, psi_ca, psi2Q, traj_param.n2Q, xfr3, &gemv_opt, spec_param.nthreads);
                    propagate_psi2Q(U2QMat, psi_cb, psi2Q, traj_param.n2Q, xfr3, &gemv_opt, spec_param.nthreads);
                  }
                }
              }
            }
            count_time(&(clock_master.calc_R3.stop));
            calc_elapsed_time(&(clock_master.calc_R3));
          }
        }
      }
    }

/*****************************************************
 * Static averaging and time-averaging approximation *
 *****************************************************/ 
    // cjfeng 05/07/2018
    // Added 2Q parts.
    // cjfeng 12/08/2018
    // Using pol_info, spec_param, and w_param
    else if ((!error) && !spec_param.nise) {
      fr = -1;
      frame = readframe-traj_param.nread;
      // Now do calculations for selected frames. 
      while( (fr!=-2) && (!error) ) {
        // cjfeng 01/05/2019
        int twin = traj_param.nbuffer-traj_param.nread+1;
        // cjfeng 07/05/2017
        // Use subroutine to count fr.
        count_fr(&fr, frame, readframe, traj_param.nread, twin, traj_param.nbuffer, spec_param.skip);
        if( (fr>=0) ) {
          // Generate averaged Hamiltonian for a window starting at frame fr 
          // and augment frame. The average Hamiltonian is stored in Ham1QAr. 
          // If Ham1QAr gets filled up, we stop and do a spectral calculation 
          // for all stored frames before proceeding. 
          frame++;  // We don't reference frame further in the calculation. 
          // cjfeng 01/09/2019
          count_time(&(clock_master.gen_avg.start));
          // Generate averaged Hamiltonian. 
          gen_avg_ham(Ham1QMat, Ham1QAr[nAvgd], hann, fr, &w_param, &traj_param);
          // Generate averaged dipole moments.
          gen_avg_dip(Dip1QMat, Dip2QMat, Dip1QAr, Dip2QAr, fr, nAvgd, &traj_param, &spec_param, &w_param);
          nAvgd++;
          count_time(&(clock_master.gen_avg.stop));
          calc_elapsed_time(&(clock_master.gen_avg));

          // Declaring local index array for eigensolver
          int *isuppz1Q, *isuppz2Q;
          set_int_1d_array(&isuppz1Q, traj_param.nosc, 2);
          if(spec_param.do2d) set_int_1d_array(&isuppz2Q, traj_param.n2Q, 2);

          // If we've filled up the Ham1QAr array, do a calculation and re-start. 
          if(nAvgd==spec_param.nthreads) {
            int thr;
            lapack_int info;
            for(thr=0; thr<spec_param.nthreads; thr++) {
              int tid = 0;
              #if OMP_PARALLEL 
                tid = omp_get_thread_num();
              #endif 
              // If needed, generate 2Q Hamiltonian. This MUST be done before eigenvalue calculation. 
              if( (spec_param.do2d) && !(spec_param.pert) ) {
                // Skip the construction of 2Q Hamiltonian when supplying external trajectory.
                // cjfeng 01/04/2019
                // Added traj_param
                count_time(&(clock_master.gen_avg.start));
                if ( !(spec_param.if2Qfiles) ) gen_ham_2Q(Ham1QAr[tid], Ham2QAr[tid], &traj_param, delta);
                else gen_avg_ham(Ham2QMat, Ham2QAr[tid], hann, fr, &w_param, &traj_param);
                count_time(&(clock_master.gen_avg.stop));
                calc_elapsed_time(&(clock_master.gen_avg));
              }
              // Find one-quantum eigenvalues
              // Note that Ham1QAr[tid] now contains eigenvectors
              count_time(&(clock_master.calc_R1.start));
              info = SYEVR(LAPACK_ROW_MAJOR, 'V', 'A', 'U', traj_param.nosc, Ham1QAr[tid], traj_param.nosc, w_param.wstart, w_param.wstop, 0, traj_param.nosc-1, 0.00001, &(traj_param.nosc), Evals1QAr[tid], Ham1QAr[tid], traj_param.nosc, isuppz1Q);
      
              // Calculate one-quantum dipole moments
              int i;
              for(i=0; i<3; i++) mvmult_real_serial_trans(Ham1QAr[tid], Dip1QAr[i][tid], ExDip1QAr[i][tid], traj_param.nosc, spec_param.nthreads);
              // Calculate FTIR spectrum
              // cjfeng 01/03/2019
              // Using calc_ftir subroutine
              calc_ftir(ftir, Evals1QAr, ExDip1QAr, &w_param, traj_param.nosc, tid);
              count_time(&(clock_master.calc_R1.stop));
              calc_elapsed_time(&(clock_master.calc_R1));
              if( spec_param.do2d ) {
                count_time(&(clock_master.calc_R3.start));
                eig_avg_Ham2Q(Ham1QAr, Ham2QAr, Evals1QAr, Evals2QAr, ExDip2QAr, isuppz2Q, &spec_param, &w_param, &traj_param, tid, delta);
                // Now calculate the 2D spectrum. 
                calc_2dir(Evals1QAr, Evals2QAr, ExDip1QAr, ExDip2QAr, tid, &traj_param, &w_param, REPH, NREPH, &pol_info, &spec_param);
                count_time(&(clock_master.calc_R3.stop));
                calc_elapsed_time(&(clock_master.calc_R3));
              }
            }
            nAvgd = 0;
          }
          // Release memory
          unset_int_1d_array(isuppz1Q);
          if(spec_param.do2d) unset_int_1d_array(isuppz2Q);
        }
      }
    }

    fr = -1;
    frame = readframe-traj_param.nread;

    // cjfeng 01/09/2019
    // Now do calculations for selected frames. 
    // cjfeng 12/09/2018
    // Using pol_info, fftw_array, spec_param, w_param, and fout
    while( (fr!=-2) && (!error) ) {
      int flag_traj = 1;
      // cjfeng 07/05/2017
      // Use subroutine to count fr.
      int twin = traj_param.nbuffer-traj_param.nread+1;
      count_fr(&fr, frame, readframe, traj_param.nread, twin, traj_param.nbuffer, spec_param.dump);
      if( (fr>=0) && (frame!=0) ) print_tstamp(frame, spec_param.tstep, fout.traj[0].fp);
      // Dump calculated spectra in trajectory files
      if( (fr>=0) && (spec_param.nise) && (frame!=0)) {
        count_time(&(clock_master.FFT.start));
        frame++; 
        FFT_linear(CorrFunc, NetCorrFunc, hann, popdecay, &FT1D, &w_param, &fout, flag_traj);
        FFT_nonlinear(REPH, NetREPH, NREPH, NetNREPH, hann, &FT2D, &pol_info, &spec_param, &w_param, &fout, flag_traj);
        count_time(&(clock_master.FFT.stop));
        calc_elapsed_time(&(clock_master.FFT));
      }
      else if( (fr>=0) && !(spec_param.nise) && (frame!=0) ) { // Static averaging or Time-averaging approximation scheme.
        count_time(&(clock_master.FFT.start));
        frame++; 
        FFT_static_linear(ftir, netftir, &FT1D, &spec_param, w_param.npts, &Const, &fout, flag_traj);
        FFT_static_nonlinear(REPH, NetREPH, NREPH, NetNREPH, &FT2D, &pol_info, &spec_param, w_param.npts, &Const, &fout, flag_traj); 
        count_time(&(clock_master.FFT.stop));
        calc_elapsed_time(&(clock_master.FFT));
      } else if(frame==0) frame++;
    }
  }

  // Print the overall spectra  
  // cjfeng 12/09/2018
  // Using pol_info, fftw_array, spec_param, w_param, and fout
  count_time(&(clock_master.FFT.start));
  int flag_traj = 0;
  if( (!error) && (spec_param.nise) ) {
    FFT_linear(CorrFunc, NetCorrFunc, hann, popdecay, &FT1D, &w_param, &fout, flag_traj);
    FFT_nonlinear(REPH, NetREPH, NREPH, NetNREPH, hann, &FT2D, &pol_info, &spec_param, &w_param, &fout, flag_traj);
  } else if( (!error) && !(spec_param.nise) ) { // Static averaging or Time averaging approximation scheme
    FFT_static_linear(ftir, netftir, &FT1D, &spec_param, w_param.npts, &Const, &fout, flag_traj);
    FFT_static_nonlinear(REPH, NetREPH, NREPH, NetNREPH, &FT2D, &pol_info, &spec_param, w_param.npts, &Const, &fout, flag_traj); 
  }
  count_time(&(clock_master.FFT.stop));
  calc_elapsed_time(&(clock_master.FFT));
  // cjfeng 01/08/2019
  // Measuring stop time and computing the total time for simulation.
  count_time(&(clock_master.run.stop));
  calc_elapsed_time(&(clock_master.run));
  print_profile_time(&clock_master, fout.log.fp);
  print_gmx_quote(fout.log.fp);

  graceful_exit(error, &fin, &fout, &trotter1Q, &trotter2Q, &FT1D, &FT2D, &traj_param, &w_param, &spec_param, &pol_info, &shift_info );
  return 0;
}
