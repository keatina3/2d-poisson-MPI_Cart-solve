#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h>
#include "utils.h"

void calc_dims(int nprocs, int* dims){
	int root = (int)sqrt(nprocs);
	while(nprocs % root != 0)
		root--;
	dims[0] = nprocs/root;
	dims[1] = root;
}

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

int decomp2d(int nx, int ny, int xprocs, int yprocs, int* coords, int *xs, int *xe, int *ys, int *ye){
    decomp1d(nx, xprocs, coords[0], xs, xe);
    decomp1d(ny, yprocs, coords[1], ys, ye);

    return MPI_SUCCESS;
}

void init_arr(int n, int m, double *x, double **x_ptr){
	int i;
	for(i=0;i<n;i++)
		x_ptr[i] = &x[i*m];
}

void clear_arr(int n, int m, double **x){
    int i, j;
    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
            x[i][j] = 0.0;
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
		//printf("dbound, xs = %d, xe = %d, ys = %d, ye = %d\n",xs, xe, ys, ye);
		for(i=(xs-1);i<=(xe+1);i++){
			uold[i][0] = bbound(i,0,nx,ny);
			unew[i][0] = bbound(i,0,nx,ny);
		}
	}
	// upper boundary //
	if(ye == ny){
		//printf("ubound, xs = %d, xe = %d, ys = %d, ye = %d\n", xs, xe, ys, ye);
		for(i=(xs-1);i<=(xe+1);i++){
			uold[i][ny+1] = ubound(i,ny+1,nx,ny);
			unew[i][ny+1] = ubound(i,ny+1,nx,ny);
		}
	}
	// left boundary //
	if(xs == 1){
		//printf("lbound, xs = %d, xe = %d, ys = %d, ye = %d\n", xs, xe, ys, ye);
		for(i=ys;i<=ye;i++){
			uold[0][i] = lbound(0,i,nx,ny);
			unew[0][i] = lbound(0,i,nx,ny);
		}
	}
	// right boundary //
	if(xe == nx){
		//printf("rbound, xs = %d, xe = %d, ys = %d, ye = %d\n", xs, xe, ys, ye);
		for(i=ys;i<=ye;i++){
			uold[nx+1][i] = rbound(nx+1,i,nx,ny);
			unew[nx+1][i] = rbound(nx+1,i,nx,ny);
		}
	}
}

void print_grid(double **u, int nx, int ny, int xs, int xe, int ys, int ye, int* coords, int *dims, MPI_Comm comm){
	int i,j;
	
    // printing top boundary //
    for(i=0;i<dims[1];i++){
		if(xe==nx){
			if(coords[1] == i){
				if(ys==1)
					printf("%lf ", u[nx+1][0]); 
				for(j=ys;j<=ye;j++)
					printf("%lf ", u[nx+1][j]);
				if(ye==ny)
					printf("%lf\n", u[nx+1][ny+1]);
            }
		}
		fflush(stdout);
		usleep(500);
		MPI_Barrier(comm);
	}
    
    // printing main part of grid //
    for(i=nx;i>0;i--){
        print_row(u, nx, ny, i, xs, xe, ys, ye, coords, dims, comm);
        MPI_Barrier(comm);
    }
    
    // printing bottom boundary //	
    for(i=0;i<dims[1];i++){
		if(xs==1){
			if(coords[1] == i){
				if(ys==1)
					printf("%lf ", u[0][0]);
				for(j=ys;j<=ye;j++)
					printf("%lf ", u[0][j]);
				if(ye==ny)
					printf("%lf\n", u[0][ny+1]);
			}
        }
        fflush(stdout);
        usleep(500);
        MPI_Barrier(comm);
    }
}

void print_row(double **u, int nx, int ny, int row, int xs, int xe, int ys, int ye, int* coords, int *dims, MPI_Comm comm){
    int i,j;

    for(i=0;i<dims[1];i++){
		if(row>=xs && row <=xe){
			if(coords[1] == i){
				if(ys==1)
					printf("%lf ", u[row][0]);
				for(j=ys;j<=ye;j++)
					printf("%lf ", u[row][j]);
				if(ye==ny)
					printf("%lf\n", u[row][ny+1]);
		    }
        }
		fflush(stdout);
		usleep(500);
		MPI_Barrier(comm);
	}
}
