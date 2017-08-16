CC=gcc
OMP_PARALLEL=TRUE
DOUBLE_PREC=FALSE
BLAS_subroutines=FALSE
USE_OpenBLAS=FALSE

GROMACS_DIR=<GROMACS HOME DIRECTORY>
BLAS_DIR=<BLAS HOME DIRECTORY>
LAPACK_DIR=<LAPACK HOME DIRECTORY>
OpenBLAS_DIR=<OpenBLAS HOME DIRECTORY> # Optional

CFLAGS=-Wall -O3

LIB_DIRS=$(GROMACS_DIR)/lib
INC_DIRS=$(GROMACS_DIR)/include/gromacs/
LIBFLAGS=-g -lgmx -ldl -lfftw3 -lm

ifeq ($(USE_OpenBLAS), TRUE)
LIB_DIRS += -L $(OpenBLAS_DIR)/lib
INC_DIRS += -I $(OpenBLAS_DIR)/include 
LIBFLAGS += -lopenblas
CFLAGS += -DUSE_OpenBLAS
else
LIB_DIRS += -L $(BLAS_DIR) -L $(LAPACK_DIR)
INC_DIRS += -I $(LAPACK_DIR)/lapacke/include 
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

OBJECTS=count.o file_io.o param_checker.o mem_helper.o mat_util.o spec_mod.o nise.o spec_fft.o trotter.o gmx_spec.o 

all: $(OBJECTS)
	$(CC) $(OBJECTS) -o g_spec -L $(LIB_DIRS) $(LIBFLAGS) $(CFLAGS)

count.o: count.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) count.c

file_io.o: file_io.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) file_io.c

param_checker.o: param_checker.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) param_checker.c

mem_helper.o: mem_helper.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) mem_helper.c

mat_util.o: mat_util.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) mat_util.c

spec_mod.o: spec_mod.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) spec_mod.c

nise.o: nise.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) nise.c

spec_fft.o: spec_fft.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) spec_fft.c

trotter.o: trotter.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) trotter.c

gmx_spec.o: gmx_spec.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) gmx_spec.c

clean:
	rm *.o
