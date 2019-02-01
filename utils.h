#ifndef _UTILS_H
#define _UTILS_H

void calc_dims(int nprocs, int* dims);
int decomp1d( int n, int size, int rank, int *s, int *e );
int decomp2d(int nx, int ny, int xprocs, int yprocs, int* coord, int *xs, int *xe, int *ys, int *ye);
void init_arr(int n, int m, double *x, double **x_ptr);
void clear_arr(int n, int m, double **x);
void init_range(double **unew, double **uold, double **f, int xs, int xe, int ys, int ye, int nx, int ny,
	double (*lbound)(int, int, int, int), 
    double (*rbound)(int, int, int, int),
	double (*ubound)(int, int, int, int),
	double (*bbound)(int, int, int, int));
void print_grid(double **u, int nx, int ny, int xs, int xe, int ys, int ye, int* coords, int *dims, MPI_Comm comm);
void print_row(double **u, int nx, int ny, int row, int xs, int xe, int ys, int ye, int* coords, int *dims, MPI_Comm comm);

#endif
