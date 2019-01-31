#ifndef PROFILER_H
#define PROFILER_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

typedef struct {
  stopwatch init;       // Program initialization
  stopwatch run;        // Program runtime
  stopwatch read;       // Reading a frame
  stopwatch gen_Dip2Q;  // Generating 2Q dipole moments
  stopwatch gen_U1Q;    // Generating 1Q Propagator
  stopwatch gen_U2Q;    // Generating 2Q Propagator
  stopwatch calc_R1;    // Calculating R1
  stopwatch calc_R3;    // Calculating R3
  stopwatch gen_avg;    // Generating averaged trajectory
  stopwatch FFT;        // FFT
} clock_master;

int initialize_clock_master(clock_master *clock_master);
int print_profile_time(clock_master *clock_master, FILE *lfp);

#endif

