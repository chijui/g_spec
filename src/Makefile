CC=gcc
OMP_PARALLEL=TRUE
GROMACS_DIR=<GROMACS HOME DIRECTORY>
LAPACK_DIR=<LAPACK HOME DIRECTORY>

LIB_DIRS=$(GROMACS_DIR)/lib -L $(LAPACK_DIR) 
INC_DIRS=$(GROMACS_DIR)/include/gromacs/ -I $(LAPACK_DIR)/lapacke/include 

ifeq ($(OMP_PARALLEL), TRUE)
CFLAGS=-Wall -O3 -fopenmp -DOMP_PARALLEL=TRUE 
LIBFLAGS=-lgmx -ldl -lfftw3 -lm -llapacke -llapack -lblas -lgfortran -lgomp
else
CFLAGS=-Wall -O3
LIBFLAGS=-lgmx -ldl -lfftw3 -lm -llapacke -llapack -lblas
endif


OBJECTS=gmx_spec.o

all: $(OBJECTS)
	$(CC) $(OBJECTS) -o g_spec -L $(LIB_DIRS) $(LIBFLAGS)

gmx_spec.o: gmx_spec.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) gmx_spec.c

clean:
	rm *.o
