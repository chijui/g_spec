#ifndef TYPES_H
#define TYPES_H

typedef struct {
	char ta;
	int inc;
	GNCOMP alpha;
	GNCOMP beta;
} BLAS_gemv_opt;

typedef struct {
	int	nosc;	// nosc*(nosc+1)/2 for 2Q
	int ncoup; // nosc * (nosc-1) /2 for 1Q, nosc^2 * (nosc-1)/2 for nonzero 2Q
	int nbuffer; 
	GNCOMP **D; // nbuffer * nosc
	GNREAL **J; // nbuffer * (2 * ncoup) 
	int *row; // row index, length: ncoup
	int *col; // column index, length: ncoup
} trotter_array;

typedef struct {
	int n;
	GNREAL *Ham;	// nosc * nosc
	GNREAL *Evals;	// nosc
	int *isuppz;	// 2 * nosc
} eig_array;

#endif
