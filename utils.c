#include <stdlib.h>
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

int decomp2d(int nx, int ny, int xprocs, int ypocs, int* coords int *u, int *b, int *r, int *l){
    decomp1d(nx, xprocs, coords[0], l, r);
    decomp1d(ny, yprocs, coords[1], b, u);

    return MPI_SUCCESS;
}

void init_arr(int n, int m, double *x, double **x_ptr){
	int i;
	for(i=0;i<n;i++)
		x_ptr[i] = &x[i*m];
}

void init_range(double **unew, double **uold, double **f, int u, int b, int r, int l,
	double (*lbound)(int, int, int, int),
	double (*rbound)(int, int, int, int),
	double (*ubound)(int, int, int, int),
	double (*bbound)(int, int, int, int))
{	

}
