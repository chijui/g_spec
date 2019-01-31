CC=gcc
OMP_PARALLEL=TRUE
DOUBLE_PREC=FALSE
BLAS_subroutines=FALSE
USE_OpenBLAS=FALSE

IDIR=<g_spec DIRECTORY>
LDIR=<g_spec DIRECTORY>

GROMACS_DIR=<GROMACS HOME DIRECTORY>
BLAS_DIR=<BLAS HOME DIRECTORY>
LAPACK_DIR=<LAPACK HOME DIRECTORY>
OpenBLAS_DIR=<OpenBLAS HOME DIRECTORY> # Optional

CFLAGS=-Wall -O3
LIBFLAGS=-g -lgmx -ldl -lfftw3 -lm

INC_DIRS=$(IDIR) -I $(GROMACS_DIR)/include/gromacs/
LIB_DIRS=$(LDIR) -L $(GROMACS_DIR)/lib

ifeq ($(USE_OpenBLAS), TRUE)
INC_DIRS += -I $(OpenBLAS_DIR)/include 
LIB_DIRS += -L $(OpenBLAS_DIR)/lib
LIBFLAGS += -lopenblas
CFLAGS += -DUSE_OpenBLAS
else
INC_DIRS += -I $(LAPACK_DIR)/lapacke/include 
LIB_DIRS += -L $(LAPACK_DIR)
LIBFLAGS += -lblas -llapack -llapacke
endif

ifeq ($(OMP_PARALLEL), TRUE)
CFLAGS += -fopenmp -DOMP_PARALLEL=TRUE
endif
ifeq ($(DOUBLE_PREC), TRUE)
CFLAGS += -DDOUBLE_PRECISION
endif
ifeq ($(BLAS_subroutines), TRUE)
CFLAGS += -DBLAS_subroutines
endif

OBJECTS=count.o profiler.o file_io.o param_checker.o param_sub.o mem_helper.o mat_util.o spec_mod.o FFT_sub.o nise.o trotter.o gmx_spec.o 

all: $(OBJECTS)
	$(CC) $(OBJECTS) -o g_spec -L $(LIB_DIRS) $(LIBFLAGS) $(CFLAGS)

count.o: count.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) count.c

profiler.o: profiler.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) profiler.c

file_io.o: file_io.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) file_io.c

param_checker.o: param_checker.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) param_checker.c

param_sub.o: param_sub.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) param_sub.c

mem_helper.o: mem_helper.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) mem_helper.c

mat_util.o: mat_util.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) mat_util.c

spec_mod.o: spec_mod.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) spec_mod.c

FFT_sub.o: FFT_sub.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) FFT_sub.c

nise.o: nise.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) nise.c

trotter.o: trotter.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) trotter.c

gmx_spec.o: gmx_spec.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) gmx_spec.c

clean:
	rm *.o
