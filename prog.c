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
	int ys, ye, xs, xe;
	int nbrup, nbrdwn, nbrleft, nbrright;
	int *coords;
	int count;
	MPI_Comm cart_comm;
	MPI_Datatype coltype;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	// use getopt to assign grid size once running properly //
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
	decomp2d(nx, ny, size, coords, &ys, &ye, &xs, &xe);

	init_range(unew, uold, f, ys, ye, xs, xe, fone, fone, fone, fone);

	MPI_Cart_shift(cart_comm, 0, 1, nbrleft, nbrright);
	MPI_Cart_shift(cart_comm, 0, 1, nbrdown, nbrup);

	MPI_Type_vector((xe-xs+1),1,ny+1,MPI_DOUBLE,&coltype);
	MPI_Type_commit(coltype);

	count = 0;
	ldiff = 1000;
	while(count < MAX_ITERS){
		exchange(uold, ys, ye, xs, xe, cart_comm, nbrleft, brright, nbrup, nbrdow, coltype);
		solve(uold, f, unew, ys, ye, xs, xe, nx);

		exchange(unew, ys, ye, xs, xe, cart_comm, nbrleft, brright, nbrup, nbrdow, coltype);
		solve(unew, f, uold, ys, ye, xs, xe, nx);
	
		ldiff = griddiff(uold,unew,xs,xe,ys,ye);
		MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
    	if(myid==0 && count%2==0){
			printf("(myid %d) locdiff: %14.10lf; glob_diff: %14.10lf\n",myid, ldiff, glob_diff);
		}
		if(glob_diff < tol){
			if(myid==0){
				printf("iterative solve converged (tol: %E)\n",tol);
			}
	      	break;
    	}

    	printf("========================= Iteration %d complete ===========================\n",it);
    	fflush(stdout);
    	usleep(500);
    	MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}
