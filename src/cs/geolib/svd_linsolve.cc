/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "geolib_math.h"
#include "geolib_methods.h"

/**
* Solve linear equation system by SVD
*/
namespace cseis_geolib {
void svd_linsolve(float** u, float* w, float** v, int nrows, int ncols, float* b, float* x, float vec_work[])
{
  int jj,i;
  float s;

  for( int j = 0; j <= ncols-1; j++ ) {
    s = 0.0;
    if (w[j]) {
      for( i = 0; i <= nrows-1; i++ ) s += u[i][j]*b[i];
      s /= w[j];
    }
    vec_work[j] = s;
  }
  for( int j = 0; j <= ncols-1; j++ ) {
    s = 0.0;
    for ( jj = 0; jj <= ncols-1; jj++ ) s += v[j][jj]*vec_work[jj];
    x[j] = s;
  }
}

void svd_linsolve(double** u, double* w, double** v, int nrows, int ncols, double* b, double* x, double vec_work[])
{
  int jj,i;
  double s;

  for( int j = 0; j <= ncols-1; j++ ) {
    s = 0.0;
    if (w[j]) {
      for( i = 0; i <= nrows-1; i++ ) s += u[i][j]*b[i];
      s /= w[j];
    }
    vec_work[j] = s;
  }
  for( int j = 0; j <= ncols-1; j++ ) {
    s = 0.0;
    for ( jj = 0; jj <= ncols-1; jj++ ) s += v[j][jj]*vec_work[jj];
    x[j] = s;
  }
}

} // namespace

