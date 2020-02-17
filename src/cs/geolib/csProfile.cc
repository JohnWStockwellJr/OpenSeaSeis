/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csProfile.h"
#include "csPoint3D.h"
#include "csLine3D.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>

using namespace std;
using namespace cseis_geolib;

//--------------------------------------------------
csProfile::csProfile() : csLine3D() {
  init(csProfile::EXTRAPOLATION_LINEAR);
}
csProfile::csProfile( csPoint3D const& p1, csPoint3D const& p2 ) : csLine3D() {
  setFromPoints( p1, p2 );
  init(csProfile::EXTRAPOLATION_LINEAR);
}
csProfile::csProfile( csPoint3D const& p1, csPoint3D const& p2, int extrapolationMethod ) : csLine3D() {
  setFromPoints( p1, p2 );
  init(extrapolationMethod);
}
void csProfile::init( int extrapolationMethod ) {
  myExtrapolationMethod = extrapolationMethod;
  myDepths   = NULL;
  myNumSteps = 0;
  myStep     = myOffset;
}
csProfile::~csProfile() {
  if( myDepths != NULL ) {
    delete [] myDepths;
    myDepths = NULL;
  }
}
csPoint3D csProfile::pointAt( double lamda ) const {
  csPoint3D p = csLine3D::pointAt( lamda );
  p.z = depthAt( lamda );
  return( p );
}
int csProfile::initDepths( double step ) {
  initDepths( (int)( myOffset/step ) + 1 );
  return myNumSteps;
}
//--------------------------------------------------------------------------------
void csProfile::initDepths( int numSteps, double const* depths ) {
  myNumSteps = numSteps;
  if( myNumSteps <= 0 ) myNumSteps = 1;
  myStep = myOffset / (double)(myNumSteps-1);
  if( myDepths != NULL ) {
    delete [] myDepths;
    myDepths = NULL;
  }
  myDepths = new double[myNumSteps];
  if( depths != NULL && numSteps >= 1 ) {
    memcpy( myDepths, depths, myNumSteps*sizeof(double) );
  }
  else {
    for( int i = 0; i < myNumSteps; i++ ) {
      myDepths[i] = 0.0;
    }
  }
}
//--------------------------------------------------------------------------------
void csProfile::setDepth( int index, double depth ) {
  if( index >= myNumSteps ) return;
  myDepths[index] = depth;
}
//--------------------------------------------------------------------------------
//
double csProfile::offsetAt( double lamda ) const {
  return( myOffset * lamda );
}
//--------------------------------------------------------------------------------
//
double csProfile::depthAt( double lamda ) const {
  if( myDepths == NULL ) return 0;

  double indexDouble = (double)(myNumSteps-1) * lamda;
  int index1 = (int)indexDouble;
  if( indexDouble < 0 ) {
    if( myExtrapolationMethod == csProfile::EXTRAPOLATION_CONSTANT ) return myDepths[0];
    index1 = 0;
  }
  else if( indexDouble >= myNumSteps-1 ) {
    if( myExtrapolationMethod == csProfile::EXTRAPOLATION_CONSTANT ) return myDepths[myNumSteps-1];
    index1 = myNumSteps-2;
  }
  int index2 = index1 + 1;
  double depth = myDepths[index1] + (indexDouble-(double)index1) * ( myDepths[index2] - myDepths[index1] );

  return depth;
}

