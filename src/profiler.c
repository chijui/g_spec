#include "profiler.h"

int initialize_clock_master(clock_master *clock_master) {
  clock_master->init.elapsed = 0;       // Program initialization
  clock_master->init.n = 0; 
  clock_master->run.elapsed = 0;        // Program runtime
  clock_master->run.n = 0;
  clock_master->read.elapsed = 0;       // Reading a frame
  clock_master->read.n = 0;       
  clock_master->gen_Dip2Q.elapsed = 0;  // Generating 2Q dipole moments
  clock_master->gen_Dip2Q.n = 0;  
  clock_master->gen_U1Q.elapsed = 0;    // Generating 1Q Propagator
  clock_master->gen_U1Q.n = 0;    
  clock_master->gen_U2Q.elapsed = 0;    // Generating 2Q Propagator
  clock_master->gen_U2Q.n = 0;    
  clock_master->calc_R1.elapsed = 0;    // Calculating R1
  clock_master->calc_R1.n = 0;    
  clock_master->calc_R3.elapsed = 0;    // Calculating R3
  clock_master->calc_R3.n = 0;    
  clock_master->gen_avg.elapsed = 0;    // Generating averaged trajectory
  clock_master->gen_avg.n = 0;    
  clock_master->FFT.elapsed = 0;        // FFT
  clock_master->FFT.n = 0; 
  return 0;
}

int calc_elapsed_time(stopwatch *stopwatch) {
  stopwatch->elapsed += stopwatch->stop - stopwatch->start;
  stopwatch->n++;
  return 0;
}

// cjfeng 01/08/2019
// Use clock master to profile the program
int print_profile_time(clock_master *clock_master, FILE *lfp) {
  fprintf(lfp, "\n");
  fprintf(lfp, "Initialization Time: %2.4f seconds, n=%d\n", clock_master->init.elapsed, clock_master->init.n);      // Program initialization
  fprintf(lfp, "Read input Time: %2.4f seconds, n=%d\n", clock_master->read.elapsed, clock_master->read.n);          // Reading
  fprintf(lfp, "Gen_Dip2Q Time: %2.4f seconds, n=%d\n", clock_master->gen_Dip2Q.elapsed, clock_master->gen_Dip2Q.n); // Generating 2Q dipole moments
  fprintf(lfp, "Gen_U1Q Time: %2.4f seconds, n=%d\n", clock_master->gen_U1Q.elapsed, clock_master->gen_U1Q.n);       // Generating 1Q Propagator
  fprintf(lfp, "Gen_U2Q Time: %2.4f seconds, n=%d\n", clock_master->gen_U2Q.elapsed, clock_master->gen_U2Q.n);       // Generating 2Q Propagator
  fprintf(lfp, "Calc_R1 Time: %2.4f seconds, n=%d\n", clock_master->calc_R1.elapsed, clock_master->calc_R1.n);       // Calculating R1
  fprintf(lfp, "Calc_R3 Time: %2.4f seconds, n=%d\n", clock_master->calc_R3.elapsed, clock_master->calc_R3.n);       // Calculating R3
  fprintf(lfp, "Gen_Avg Time: %2.4f seconds, n=%d\n", clock_master->gen_avg.elapsed, clock_master->gen_avg.n);       // Generating averaged trajectory
  fprintf(lfp, "FFT Time: %2.4f seconds, n=%d\n", clock_master->FFT.elapsed, clock_master->FFT.n);                   // FFT
  fprintf(lfp, "Time: %2.4f seconds, n=%d\n", clock_master->run.elapsed, clock_master->run.n);                       // Program runtime
  fflush(lfp);
  printf("Time: %2.4f seconds \n", clock_master->run.elapsed);                                                        // Program runtime
  return 0;
}


