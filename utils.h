#ifndef _UTILS_H
#define _UTILS_H

int decomp1d( int n, int size, int rank, int *s, int *e );
int decomp2d(int nx, int ny, int xprocs, int yprocs, int* coord, int *ys, int *ye, int *xs, int *xe);
void init_arr(int n, int m, double *x, double **x_ptr);
void init_range(double **unew, double **uold, double **f, int ys, int ys, int xs, int xe,
	double (*lbound)(int, int, int, int), 
    double (*rbound)(int, int, int, int),
	double (*ubound)(int, int, int, int),
	double (*bbound)(int, int, int, int));

#endif
