/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SPLINE_INTERPOLATION_H
#define CS_SPLINE_INTERPOLATION_H

namespace cseis_geolib {

/**
 * 1D spline interpolation
 */
class csSplineInterpolation {
 public:
  csSplineInterpolation();
  csSplineInterpolation( int numSamples, float x1, float dx );
  ~csSplineInterpolation();

  void set( int numSamples, float x1, float dx );
  void prepare( float const* y );
  float compute( float x );

 private:
  float myData_x1;
  float myData_dx;
  float* myData_y;
  float* myData_y2;
  int    myNumValues;
};


}

#endif

