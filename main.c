#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "poisson.h"
#include "jacobi.h"
#include "utils.h"

#define GRID_SIZE 8
#define MAX_ITERS 2000
#define ERR 1.0E-07

int main(int argc, char* argv[]){
	double *uold_vals, *unew_vals, *f_vals;
	double **uold, **unew, **f;
	double ldiff, glob_diff;
	int myid, nprocs;
	int *dims, ndims = 2;
	int u, b, l, r;
	int nbrup, nbrdwn, nbrleft, nbrright;
	int *coords;
	MPI_Comm cart_comm;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	nx = nx = GRID_SIZE+1;
	uold_vals = (double*)calloc(nx*ny,sizeof(double));
	unew_vals = (double*)calloc(nx*ny,sizeof(double));
	f_vals = (double*)calloc(nx*ny,sizeof(double));
	
	init_arr(nx+1, ny+1, uold_vals, uold);
	init_arr(nx+1, ny+1, unew_vals, unew);
	init_arr(nx+1, ny+1, f_vals, f);

	dims[0] = size/2;
	dims[1] = 2;

	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, 0, 0, cart_comm);
	
	MPI_Cart_coords(cart_comm, myid, 2, coords);
	decomp2d(nx, ny, size, coords, &u, &b, &r, &l);

	init_range(unew, uold, f, u, b, r, l, fone, fone, fone, fone);

	MPI_Cart_shift(cart_comm, 0, 1, nbrleft, nbrright);
	MPI_Cart_shift(cart_comm, 0, 1, nbrdown, nbrup);

	MPI_Finalize();
	return 0;
}
