/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSplineInterpolation.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "csException.h"

using namespace cseis_geolib;

csSplineInterpolation::csSplineInterpolation() {
  myData_x1 = 0;
  myData_dx = 0;
  myData_x1 = 0;
  myData_dx = 0;
  myData_y  = NULL;
  myData_y2 = NULL;
}

csSplineInterpolation::csSplineInterpolation( int numValues, float x1, float dx ) {
  set( numValues, x1, dx );
}

csSplineInterpolation::~csSplineInterpolation() {
  if( myData_y != NULL ) {
    delete [] myData_y;
    myData_y = NULL;
  }
  if( myData_y2 != NULL ) {
    delete [] myData_y2;
    myData_y2 = NULL;
  }
}
void csSplineInterpolation::set( int numValues, float x1, float dx ) {
  if( numValues < 3 ) throw( csException("csSplineInterpolation(): Too few samples in input function: %d", numValues) );
  myNumValues = numValues;
  myData_x1 = x1;
  myData_dx = dx;
  myData_y  = new float[myNumValues];
  myData_y2 = new float[myNumValues];
}

void csSplineInterpolation::prepare( float const* y ) {
  memcpy( myData_y, y, myNumValues*sizeof(float) );
  int n = myNumValues;

  float* uArray = new float[n-1];
  uArray[0] = 0.0f;
  myData_y2[0] = 0.0f;

  for( int i = 1; i <= n - 2; i++ ) {
    float x_imin1  = myData_x1 + (float)(i-1)*myData_dx;
    float x_i      = myData_x1 + (float)(i)*myData_dx;
    float x_iplus1 = myData_x1 + (float)(i+1)*myData_dx;
    float sig = ( x_i - x_imin1 ) / ( x_iplus1 - x_imin1 );
    float p = sig * myData_y2[i-1] + 2.0f;
    myData_y2[i] = ( sig - 1.0f ) / p;
    uArray[i] = ( y[i+1] - y[i] ) / ( x_iplus1 - x_i ) - ( y[i] - y[i-1] ) / ( x_i - x_imin1 );
    uArray[i] = ( 6.0f * uArray[i] / ( x_iplus1 - x_imin1 ) - sig * uArray[i-1] ) / p;
  }
  
  float qn = 0.0f;
  float un = 0.0f;
  myData_y2[n-1] = ( un - qn * uArray[n-2] ) / ( qn * myData_y2[n-2] + 1.0f );
  
  for( int k = n - 2; k >= 0; k-- ) {
    myData_y2[k] = myData_y2[k] * myData_y2[k + 1] + uArray[k];
  }
  delete [] uArray;
}

//-----------------------------------------------------------

float csSplineInterpolation::compute( float xout ) {
  float result;
  int klo,khi,k;
  float h,b,a;

  klo = 0;
  khi = myNumValues-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    float x = myData_x1 + myData_dx*(float)k;
    if (x > xout) khi=k;
    else klo=k;
  }
  float xhi = myData_x1 + myData_dx * (float)khi;
  float xlo = myData_x1 + myData_dx * (float)klo;
  h = myData_dx * (float)( xhi - xlo );
  a = ( xhi - xout )/h;
  b = ( xout - xlo ) / h;
  result = a*myData_y[klo] + b*myData_y[khi] + ((a*a*a-a)*myData_y2[klo]+(b*b*b-b)*myData_y2[khi])*(h*h)/6.0f;
  return result;
}

