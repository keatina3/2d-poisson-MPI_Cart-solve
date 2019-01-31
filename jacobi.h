#ifndef _JACOBI_H_
#define _JACOBI_H_

void exchange(double **x, int ys, int ye, int xs, int xe, MPI_Comm comm, int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Datatype coltype);
void solve(double **uold, double **f, double **unew, int ys, int ye, int xs, int xe, int nx); 
double griddiff(double **x, double **y, int xs, int xe, int ys, int ye);

#endif
