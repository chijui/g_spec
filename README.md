# g_spec

g_spec is a spectral simulation program for converting molecular energy and
dipole trajectories into simulated linear and nonlinear optical spectra. In
conjunction with the open-source program g_amide, it enables the user to
translate  structural trajectories from a molecular dynamics (MD) simulation
in the GROMACS simulation package (www.gromacs.org) into Amide I (backbone
carbonyl stretch) vibrational spectra. 

The basic input to the g_spec program is a set of text files providing a
trajectory of exciton-like Hamiltonians and transition-dipole moment vectors.
In the context of protein Amide I calculations, such files can be generated
using the g_amide program referenced above. Using this input, g_spec
simulates linear absorption and nonlinear (two-dimensional) spectroscopic
signals using a variety of methods. 



DEPENDENCIES

g_spec calls several external libraries in its function. These libraries MUST
be available during installation or the code will not compile: 

(1) A working GROMACS installation for MD simulations. (See www.gromacs.org).

(2) A functional installation of FFTW3 for fast Fourier transforms (see
www.fftw.org).

(3) Working LAPACK and LAPACKE libraries for linear algebra processing (see
www.netlib.org/lapack). 

(4) Working BLAS libraries for linear algebra processing (see
www.netlib.org/blas). 



Although not required, g_spec can also make use of OpenMP libraries to
implement a certain level of parallelization into the calculation. To make use
of these features, a working OpenMP installation must also be available. 



INSTALLATION

(0) Choose an installation directory.

This may be any directory where you have read/write permissions, but should be
a permanent installation location (e.g. not your user download folder). For
example, if your username is mike, you might wish to install the program in an
"apps" directory such as /home/mike/apps/. Move the g_spec directory (where
this README.txt file is located) into the chosen installation directory and at
the command prompt cd to that location, e.g. in our example

        cd /home/mike/apps/g_spec/


(1) Prepare the Makefile.

The Makefile specifies machine parameters for compiling the code. Two
alternative makefiles are included in this distribution. The default file
named simple "Makefile" makes use of the gcc C-language compiler, installed
onmost linux-based operating systems. The second, Makefile_icc is configured
for compilation with the intel icc compiler. This compiler is able to
optimized the g_spec installation in a fashion that make make it run faster
on some systems. If the icc compiler is available on your system, we recommend
using it. To do so, enter at the command line

	cp Makefile_icc Makefile

Note that a backup copy of the gcc makefile is stored as Makefile_gcc if you
wish to rever to gcc installation. 

Next, modify the Makefile. E.g. to open the file using vi, type at the console

        vi src/Makefile

Modifications should be required onlyin the first four lines of the file.
Three lines should be checked. 

First, the flag OMP_PARALLEL should be set to either "TRUE" or (by default)
"FALSE". If set to true, this flag enables parallel processing using the
OpenMP library. Parallel processing will enable to code to execute
calculations faster if multiple processors are available, but installation
requires an OpenMP package accessible during compilation. (This is the -lgomp
flag under the LIBFLAGS variable in the Makefile). If you do not have an
OpenMP library available (or do not know how to access it), leave the flag set
to FALSE. 

Second, the GROMACS_DIR flag should be set to the absolute path of your
GROMACS installation directory, e.g. /home/mike/apps/gromacs/. 

Finally, the LAPACK_DIR flag must be set to the absolute path of the LAPACK
installation directory, e.g. /home/mike/apps/lapack. Note that this package
must include the lapacke header files for linking with C code. 


(2) Compile the code. In the same directory as your Makefile and source code
enter at the command line

	make 

The make program should locate and read your Makefile and the selected
compiler will compile the code. If you encounter errors, check that the
various libraries (-lxxx flags) and included header files (-I flags) are
available to the system. If you encounter difficulties, try using the simplest
settings possible first, using the gcc compiler (assuming it is available on
your system) with OMP_PARALLEL set to FALSE. 

(3) Run a test calculation! A set of example input files is included in the
directory g_spec/test/input, along with a brief tutorial describing basic
syntax for the program. 

