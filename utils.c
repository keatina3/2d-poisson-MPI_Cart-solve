#include <stdlib.h>
#include <mpi.h>
#include "utils.h"

int decomp1d(int n, int size, int rank, int *s, int *e){
    int nlocal, deficit;

    nlocal  = n / size;
    *s  = rank * nlocal + 1;
    deficit = n % size;
    *s  = *s + ((rank < deficit) ? rank : deficit);
    if (rank < deficit) nlocal++;
    *e      = *s + nlocal - 1;
    if (*e > n || rank == size-1) *e = n;
    
    return MPI_SUCCESS;
}

int decomp2d(int nx, int ny, int xprocs, int yprocs, int* coords, int *ys, int *ye, int *xs, int *xe){
    decomp1d(nx, xprocs, coords[0], xs, xe);
    decomp1d(ny, yprocs, coords[1], ys, ye);

    return MPI_SUCCESS;
}

void init_arr(int n, int m, double *x, double **x_ptr){
	int i;
	for(i=0;i<n;i++)
		x_ptr[i] = &x[i*m];
}

void init_range(double **unew, double **uold, double **f, int xs, int xe, int ys, int ye, int nx, int ny,
	double (*lbound)(int, int, int, int),
	double (*rbound)(int, int, int, int),
	double (*ubound)(int, int, int, int),
	double (*bbound)(int, int, int, int))
	{	
	int i;

	// lower boundary //
	if(ys == 1){
		for(i=xs;i<=xe;i++){
			uold[i][0] = bbound(i,0,nx,ny);
			unew[i][0] = bbound(i,0,nx,ny);
		}
	}
	// upper boundary //
	if(ye == ny){
		for(i=xs;i<=xe;i++){
			uold[i][ny+1] = ubound(i,ny+1,nx,ny);
			unew[i][ny+1] = ubound(i,ny+1,nx,ny);
		}
	}
	// left boundary //
	if(xs == 1){
		for(i=ys;i<=ye;i++){
			uold[0][i] = lbound(0,i,nx,ny);
			unew[0][i] = lbound(0,i,nx,ny);
		}
	}
	// right boundary //
	if(xe == nx){
		for(i=ys;i<=ye;i++){
			uold[nx+1][i] = lbound(nx+1,i,nx,ny);
			unew[nx+1][i] = lbound(nx+1,i,nx,ny);
		}
	}
}
