#include "poisson.h"

double fzero(int xind, int yind, int nx, int ny){
    return 0.0;
}

double fone(int xind, int yind, int nx, int ny){
  return 1.0;
}

/*
  $\frac{1}{(1+x)^2 + 1}$
*/
double fiserles1(int xind, int yind, int nx, int ny){
  double x, y;
  double tmpv;

  getxy(nx, ny, xind, yind, &x, &y);
  tmpv = (1+x)*(1+x) + 1.0;
  
  return 1.0/tmpv;
}

/*
  $\frac{y}{y*2 + 1}$
*/
double fiserles2(int xind, int yind, int nx, int ny){
  double x, y;
  double tmpv;

  getxy(nx, ny, xind, yind, &x, &y);
  tmpv = 1 + y*y;
  
  return y/tmpv;
}

/*
  $\frac{y}{y*2 + 4}$
*/
double fiserles3(int xind, int yind, int nx, int ny){
  double x, y;
  double tmpv;

  getxy(nx, ny, xind, yind, &x, &y);
  tmpv = 4.0 + y*y;
  
  return y/tmpv;
}

/*
  $u(x,y)=\frac{y}{(1+x)^2 + y^2}$
*/
double fiserlessoln(int xind, int yind, int nx, int ny){
  double x, y;
  double fval;
  double tmpv;

  getxy(nx, ny, xind, yind, &x, &y);
  tmpv = (1+x)*(1+x) + y*y;
  fval = y/tmpv;
 
  return fval;

}

int getxy(int nx, int ny, int xind, int yind, double *x, double *y){
  double hx, hy;
  double lx, ly;

  gethxy(nx, ny, &hx, &hy);

  lx = FD_LX + xind*hx;
  ly = FD_LY + yind*hy;
  *x = lx;
  *y = ly;

  return 0;
}

int gethxy(int nx, int ny, double *hx, double *hy){

  *hx = (FD_UX - FD_LX)/( (double)(nx+1) );
  *hy = (FD_UY - FD_LY)/( (double)(ny+1) );

  return 0;
}
