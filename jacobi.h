#ifndef _JACOBI_H_
#define _JACOBI_H_

#define MAX_ITERS 2000
#define ERR 1.0E-07

void exchange(double **x, int xs, int xe, int ys, int ye, MPI_Comm comm, int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Datatype coltype);
void solve(double **uold, double **f, double **unew, int xs, int xe, int ys, int ye, int nx); 
void itersolve(double **unew, double **uold, double **f, int xs, int xe, int ys, int ye, int nx, int nbrleft, int nbrright, int nbrup, int nbrdown, int myid, MPI_Datatype coltype, MPI_Comm comm);
double griddiff(double **x, double **y, int xs, int xe, int ys, int ye);

#endif
