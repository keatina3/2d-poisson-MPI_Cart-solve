#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "poisson.h"
#include "jacobi.h"
#include "utils.h"

#define GRID_SIZE 9
#define MAX_ITERS 2000
#define ERR 1.0E-07

int main(int argc, char* argv[]){
	double *uold_vals, *unew_vals, *f_vals;
	double **uold, **unew, **f;
	double ldiff, glob_diff;
	int myid, nprocs;
	int dims[2], ndims = 2;
	int periods[2] = {0,0};
	int ys, ye, xs, xe;
	int nbrup, nbrdown, nbrleft, nbrright;
	int coords[2];
	int count;
	int nx, ny;
	MPI_Comm cart_comm;
	MPI_Datatype coltype;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	// use getopt to assign grid size once running properly //
	nx = ny = GRID_SIZE;

	uold_vals = (double*)calloc((nx+2)*(ny+2),sizeof(double));
	unew_vals = (double*)calloc((nx+2)*(ny+2),sizeof(double));
	f_vals = (double*)calloc((nx+2)*(ny+2),sizeof(double));
	uold = (double**)malloc((nx+2)*sizeof(double*));
	unew = (double**)malloc((nx+2)*sizeof(double*));
	f = (double**)malloc((nx+2)*sizeof(double*));

	init_arr(nx+2, ny+2, uold_vals, uold);
	init_arr(nx+2, ny+2, unew_vals, unew);
	init_arr(nx+2, ny+2, f_vals, f);

	//dims[0] = 2;
	//dims[1] = nprocs/2;
	calc_dims(nprocs, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 0, &cart_comm);
	
	MPI_Cart_coords(cart_comm, myid, 2, coords);
	decomp2d(nx, ny, dims[0], dims[1], coords, &xs, &xe, &ys, &ye);
	//printf("coords = (%d,%d), rank = %d, xs = %d, xe = %d, ys = %d, ye = %d\n",coords[0],coords[1], myid, xs, xe, ys, ye);
	MPI_Barrier(cart_comm);

	init_range(unew, uold, f, xs, xe, ys, ye, nx, ny, &fone, &fone, &fone, &fone);

	MPI_Cart_shift(cart_comm, 0, 1, &nbrleft, &nbrright);
	MPI_Cart_shift(cart_comm, 1, 1, &nbrdown, &nbrup);
	//printf("coords = (%d,%d), nbrleft = %d, nbrright = %d, nbrup = %d, nbrdown = %d\n", coords[0],coords[1], nbrleft, nbrright, nbrup, nbrdown);

	MPI_Type_vector((xe-xs+1), 1, ny+2, MPI_DOUBLE, &coltype);
	MPI_Type_commit(&coltype);

	count = 0;
	ldiff = 1000;
	while(count < MAX_ITERS){
		exchange(uold, xs, xe, ys, ye, cart_comm, nbrleft, nbrright, nbrup, nbrdown, coltype);
		solve(uold, f, unew, xs, xe, ys, ye, nx);

		exchange(unew, xs, xe, ys, ye, cart_comm, nbrleft, nbrright, nbrup, nbrdown, coltype);
		solve(unew, f, uold, xs, xe, ys, ye, nx);
	
		ldiff = griddiff(uold,unew,xs,xe,ys,ye);
		MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
    	if(myid==0 && count%2==0){
			printf("(myid %d) locdiff: %14.10lf; glob_diff: %14.10lf\n",myid, ldiff, glob_diff);
		}
		if(glob_diff < ERR){
			if(myid==0){
				printf("iterative solve converged (tol: %E, count: %d)\n",ERR, count);
			}
	      	break;
    	}

    	//printf("========================= Iteration %d complete ===========================\n",count);
    	fflush(stdout);
    	//usleep(500);
		count++;
    	MPI_Barrier(MPI_COMM_WORLD);
	}
	print_grid(unew, nx, ny, xs, xe, ys, ye, coords, dims, cart_comm);

	MPI_Finalize();
	return 0;
}
