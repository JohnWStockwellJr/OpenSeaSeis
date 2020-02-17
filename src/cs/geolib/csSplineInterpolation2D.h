/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SPLINE_INTERPOLATION_2D_H
#define CS_SPLINE_INTERPOLATION_2D_H

namespace cseis_geolib {

class csSplineInterpolation;

/**
 * 2D spline interpolation
 */
class csSplineInterpolation2D {
 public:
  csSplineInterpolation2D( int numDim1, float x1, float dx1, int numDim2, float x2, float dx2 );
  ~csSplineInterpolation2D();

  void set( float x1, float dx1, float x2, float dx2 );
  void prepareAll( float const* yy );
  void prepareOne( float const* y, float valDim1 );
  float compute( float x1out, float x2out );
  void compute( float x1out, int numValuesOut, float const* x2out, float* yout );

 private:
  float myDim1_x0;
  float myDim1_dx;
  float myDim2_x0;
  float myDim2_dx;
  int    myNumDim1;
  int    myNumDim2;
  csSplineInterpolation* mySplinesDim1;
  csSplineInterpolation* mySplineDim2;
};


}

#endif

