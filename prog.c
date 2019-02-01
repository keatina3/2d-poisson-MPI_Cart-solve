#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <mpi.h>

#include "poisson.h"
#include "jacobi.h"
#include "utils.h"

#define GRID_SIZE 9

int main(int argc, char* argv[]){
	double *uold_vals, *unew_vals, *f_vals;
	double **uold, **unew, **f;
	int myid, nprocs;
	int dims[2], ndims = 2;
	int periods[2] = {0,0};
	int ys, ye, xs, xe;
	int nbrup, nbrdown, nbrleft, nbrright;
	int coords[2];
	int nx, ny;
	int opt;
    MPI_Comm cart_comm;
	MPI_Datatype coltype;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
    // getting grid size //
	nx = GRID_SIZE;
    if(myid==0){
        while((opt=getopt(argc,argv,"g:"))!=-1){
            switch(opt){
                case 'g':
                    nx = atoi(optarg);
                    printf("\nGrid size set to %d+2\n",nx);
                    break;
            }
        }
    }
    MPI_Bcast(&nx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    ny = nx;

    // allocating arrays to store grids //
	uold_vals = (double*)calloc((nx+2)*(ny+2),sizeof(double));
	unew_vals = (double*)calloc((nx+2)*(ny+2),sizeof(double));
	f_vals = (double*)calloc((nx+2)*(ny+2),sizeof(double));
	uold = (double**)malloc((nx+2)*sizeof(double*));
	unew = (double**)malloc((nx+2)*sizeof(double*));
	f = (double**)malloc((nx+2)*sizeof(double*));

	init_arr(nx+2, ny+2, uold_vals, uold);
	init_arr(nx+2, ny+2, unew_vals, unew);
	init_arr(nx+2, ny+2, f_vals, f);
    
    // getting dimensions of cartesian grid and setting up cart comminicator //
	calc_dims(nprocs, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 0, &cart_comm);
	
	MPI_Cart_coords(cart_comm, myid, 2, coords);
	decomp2d(nx, ny, dims[0], dims[1], coords, &xs, &xe, &ys, &ye);


    ///////////////// TESTING FOR ALL BOUNDARIES = 1.0 //////////////////////////
    
	if(myid==0){
        printf("=========================================================================================\n");
        printf("\n\nTESTING WITH BOUNDARY CONDITIONS = 1.0\n\n");
        printf("Initial grid with boundary conditions applied:\n\n");
    }
    init_range(unew, uold, f, xs, xe, ys, ye, nx, ny, &fone, &fone, &fone, &fone);
    print_grid(unew, nx, ny, xs, xe, ys, ye, coords, dims, cart_comm);
	
    // getting neighbours //
    MPI_Cart_shift(cart_comm, 0, 1, &nbrleft, &nbrright);
	MPI_Cart_shift(cart_comm, 1, 1, &nbrdown, &nbrup);
    
    // committing vector datatype to send columns //
	MPI_Type_vector((xe-xs+1), 1, ny+2, MPI_DOUBLE, &coltype);
	MPI_Type_commit(&coltype);
    
    // iterative Jacobi solver //
    itersolve(unew, uold, f, xs, xe, ys, ye, nx, nbrleft, nbrright, nbrup, nbrdown, myid, coltype, cart_comm);

	if(myid==0)
        printf("\nSolved PDE:\n\n");
	print_grid(unew, nx, ny, xs, xe, ys, ye, coords, dims, cart_comm);
    usleep(1000);
    MPI_Barrier(MPI_COMM_WORLD);

    ////////////////////////////////////////////////////////////////////////////////

    if(myid==0){
        printf("\n\nPress Enter to continue and test other function:\n");
        while( getchar() != '\n' );
    }
    MPI_Barrier(cart_comm);
    
    clear_arr(nx+2,ny+2,unew);
    clear_arr(nx+2,ny+2,uold);
    clear_arr(nx+2,ny+2,f);
    
    ///////////////// TESTING FOR BOUNDARY CONDITIONS IN NOTES /////////////////////
    
	if(myid==0){
        printf("=========================================================================================\n");
        printf("\n\nTESTING WITH BOUNDARY CONDITIONS FROM NOTES\n\n");
        printf("Initial grid with boundary conditions applied:\n\n");
    }
    MPI_Barrier(cart_comm);
    usleep(1000);
    fflush(stdout);
    init_range(unew, uold, f, xs, xe, ys, ye, nx, ny, &fiserles2, &fiserles3, &fiserles1, &fzero);
    print_grid(unew, nx, ny, xs, xe, ys, ye, coords, dims, cart_comm);
    
    itersolve(unew, uold, f, xs, xe, ys, ye, nx, nbrleft, nbrright, nbrup, nbrdown, myid, coltype, cart_comm);
    
	if(myid==0)
        printf("\nSolved PDE:\n\n");
	print_grid(unew, nx, ny, xs, xe, ys, ye, coords, dims, cart_comm);
    usleep(1000);
    MPI_Barrier(MPI_COMM_WORLD);
    //////////////////////////////////////////////////////////////////////////////
    
    free(f_vals);
    free(unew_vals);
    free(uold_vals);
    free(f);
    free(unew);
    free(uold);

	MPI_Finalize();
    	
    return 0;
}
