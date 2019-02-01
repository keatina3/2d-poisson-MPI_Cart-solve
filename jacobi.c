#include <stdio.h>
#include <unistd.h>
#include "mpi.h"
#include "jacobi.h"

void exchange(double **x, int xs, int xe, int ys, int ye, MPI_Comm comm, int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Datatype coltype){
	
	MPI_Sendrecv(&x[xe][ys], (ye-ys+1), MPI_DOUBLE, nbrright, 0, 
			&x[xs-1][ys], (ye-ys+1), MPI_DOUBLE, nbrleft, 0, comm, MPI_STATUS_IGNORE);
	
	MPI_Sendrecv(&x[xs][ys], (ye-ys+1), MPI_DOUBLE, nbrleft, 0, 
			&x[xe+1][ys], (ye-ys+1), MPI_DOUBLE, nbrright, 0, comm, MPI_STATUS_IGNORE);
	
	MPI_Sendrecv(&x[xs][ye], 1, coltype, nbrup, 0, 
			&x[xs][ys-1], 1, coltype, nbrdown, 0, comm, MPI_STATUS_IGNORE);
	
	MPI_Sendrecv(&x[xs][ys], 1, coltype, nbrdown, 0, 
			&x[xs][ye+1], 1, coltype, nbrup, 0, comm, MPI_STATUS_IGNORE);
}

void solve(double **uold, double **f, double **unew, int xs, int xe, int ys, int ye, int nx){
	double h;
	int i,j;

	h  = 1.0/(double)(nx+1);

	for(i=xs;i<=xe;i++)
		for(j=ys;j<=ye;j++)
			unew[i][j] = 0.25*(uold[i-1][j]+uold[i+1][j]+uold[i][j+1]+uold[i][j-1]-h*h*f[i][j]);
}

void itersolve(double **unew, double **uold, double **f, int xs, int xe, int ys, int ye, int nx, int nbrleft, int nbrright, int nbrup, int nbrdown, int myid, MPI_Datatype coltype, MPI_Comm cart_comm){
    double ldiff = 1000.0, glob_diff = 0.0;
    int count = 0;

    while(count < MAX_ITERS){
        exchange(uold, xs, xe, ys, ye, cart_comm, nbrleft, nbrright, nbrup, nbrdown, coltype);
        solve(uold, f, unew, xs, xe, ys, ye, nx);

        exchange(unew, xs, xe, ys, ye, cart_comm, nbrleft, nbrright, nbrup, nbrdown, coltype);
        solve(unew, f, uold, xs, xe, ys, ye, nx);

        ldiff = griddiff(uold,unew,xs,xe,ys,ye);
        MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, cart_comm);

        if(glob_diff < ERR){
            if(myid==0){
                printf("\n\niterative solve converged (tol: %E, count: %d)\n\n",ERR, count);
            }
            break;
        }
        
        //fflush(stdout);
        usleep(500);
        count++;
        MPI_Barrier(cart_comm);
    }
}

double griddiff(double **x, double **y, int xs, int xe, int ys, int ye){
	double sum;
	double tmp;
	int i,j;
	
	sum = 0.0;

	for(i=xs;i<=xe;i++){
		for(j=ys;j<=ye;j++){
			tmp = (x[i][j] - y[i][j]);
			sum += tmp*tmp;
		}
	}
	return sum;
}
