#ifndef _POIS_H_
#define _POIS_H_

#define FD_LX 0.0
#define FD_UX 1.0
#define FD_LY 0.0
#define FD_UX 1.0

// functionality to evaluate Poisson PDE //

double fzero(int xind, int yind, int nx, int ny);
double fone(int xind, int yind, int nx, int ny);
double fiserles1(int xind, int yind, int nx, int ny);
double fiserles2(int xind, int yind, int nx, int ny);
double fiserles3(int xind, int yind, int nx, int ny);
double fiserlessoln(int xind, int yind, int nx, int ny);
int getxy(int nx, int ny, int xind, int yind, double* x, double* y);
int gethxy(int nx, int ny, double* hx, double* hy);

#endif
