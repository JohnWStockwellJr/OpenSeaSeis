/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSplineInterpolation.h"
#include "csSplineInterpolation2D.h"
#include "csException.h"

using namespace cseis_geolib;

csSplineInterpolation2D::csSplineInterpolation2D( int numDim1, float x1, float dx1, int numDim2, float x2, float dx2 ) {
  myNumDim1 = numDim1;
  myNumDim2 = numDim2;
  myDim1_x0  = x1;
  myDim1_dx  = dx1;
  myDim2_x0  = x2;
  myDim2_dx  = dx2;
  if( myNumDim1 < 3 || myNumDim2 < 3 ) throw( csException("csSplineInterpolation2D(): Too few values in input data: %d x %d. Need to be at least 3 x 3 matrix", myNumDim1, myNumDim2) );
  mySplineDim2  = new csSplineInterpolation( myNumDim2, myDim2_x0, myDim2_dx );
  mySplinesDim1 = new csSplineInterpolation[myNumDim1];
  for( int i = 0; i < myNumDim1; i++ ) {
    mySplinesDim1[i].set( myNumDim2, x1, dx1 );
  }
}

csSplineInterpolation2D::~csSplineInterpolation2D() {
  if( mySplineDim2 != NULL ) {
    delete mySplineDim2;
    mySplineDim2 = NULL;
  }
  if( mySplinesDim1 != NULL ) {
    delete [] mySplinesDim1;
    mySplinesDim1 = NULL;
  }
}
void csSplineInterpolation2D::set( float x1, float dx1, float x2, float dx2 ) {
  myDim1_x0  = x1;
  myDim1_dx = dx1;
  myDim2_x0  = x2;
  myDim2_dx = dx2;
}

void csSplineInterpolation2D::prepareOne( float const* y, float valDim1 ) {
  int index = (int)( (valDim1 - myDim1_x0)/myDim1_dx );
  mySplinesDim1[index].prepare( y );
}

void csSplineInterpolation2D::prepareAll( float const* yy ) {
  for( int idim1 = 0; idim1 < myNumDim1; idim1++ ) {
    mySplinesDim1[idim1].prepare( &yy[idim1*myNumDim2] );
  }
}

float csSplineInterpolation2D::compute( float x1out, float x2out ) {
  float* values = new float[myNumDim1];
  for( int idim1 = 0; idim1 < myNumDim1; idim1++ ) {
    //    float dim1 = idim1 * myDim1_dx + myDim1_x0;
    values[idim1] = mySplinesDim1[idim1].compute( x1out );
  }
  mySplineDim2->prepare( values );
  float yout = mySplineDim2->compute( x2out );
  delete [] values;
  return yout;
}

void csSplineInterpolation2D::compute( float x1out, int numValuesOut, float const* x2out, float* yout ) {
  float* values = new float[myNumDim1];
  for( int idim1 = 0; idim1 < myNumDim1; idim1++ ) {
    //    float dim1 = idim1 * myDim1_dx + myDim1_x0;
    values[idim1] = mySplinesDim1[idim1].compute( x1out );
  }
  mySplineDim2->prepare( values );
  for( int i = 0; i < numValuesOut; i++ ) {
    yout[i] = mySplineDim2->compute( x2out[i] );
  }
  delete [] values;
}
