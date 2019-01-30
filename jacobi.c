#include "mpi.h"
#include "jacobi.h"

void exchange(double **x, int ys, int ye, int xs, int xe, MPI_Comm comm, int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Type coltype){
	MPI_Sendrecv(&x[xe][ys], (ye-ys+1), MPI_DOUBLE, nbrright, 0, 
		&x[xs-1][ys], (ye-ys+1), MPI_DOUBLE, nbrleft, 0, comm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(&x[xs][ys], (ye-ys+1), MPI_DOUBLE, nbrleft, 0, 
		&x[xe+1][ys], (ye-ys+1), MPI_DOUBLE, nbrright, 0, comm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(&x[xs][ye], 1, coltype, nbrup, 0, 
		&x[xs][ys-1], 1, coltype, nbrdown, comm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(&x[xs][ys], 1, coltype, nbrdown, 0, 
		&x[xs][ye+1], 1, coltype, nbrup, comm, MPI_STATUS_IGNORE);
}

void solve(double **uold, double **f, double **unew, int ys, int ye, int xs, int xe, int nx){
	double h;
	int i,j;

	h  = 1.0/(double)(nx+1);

	for(i=xs;i<=xe;i++)
		for(j=ys;j<ye;j++)
			unew[i][j] = 0.25*(uold[i-1][j]+uold[i+1][j]+uold[i][j+1]+uold[i][j-1]-h*h*f[i][j]);
}

double griddiff(double **x, double **y, int xs, int xe, int ys, int ye){
	double sum;
	double tmp;
	int i,j;

	for(i=xs;i<=xe;i++){
		for(j=ys;j<=ye;j++){
			tmp = (a[i][j] - b[i][j]);
			sum = sum + tmp*tmp;
		}
	}
	return sum;
}