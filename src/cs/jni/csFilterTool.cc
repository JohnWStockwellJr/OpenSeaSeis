/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csFFTTools.h"
#include "csFilterTool.h"
#include <cmath>
#include <cstring>
#include <cstdio>

using namespace cseis_jni;
using namespace cseis_geolib;

csFilterTool::csFilterTool( int numSamples, float sampleInt )  : cseis_geolib::csFFTTools( numSamples, sampleInt ) {
  mySamples = new float[numSamples];
  myIsFilter = false;
}
csFilterTool::csFilterTool( int numSamplesIn, int numSamplesOut, float sampleIntIn, float sampleIntOut ) : cseis_geolib::csFFTTools(numSamplesIn, numSamplesOut, sampleIntIn, sampleIntOut) {
  mySamples = new float[numSamplesIn];
  myIsFilter = false;
}

csFilterTool::~csFilterTool() {
  if( mySamples != NULL ) {
    delete [] mySamples;
    mySamples = NULL;
  }
}
bool csFilterTool::isFilter() const {
  return myIsFilter;
}

void csFilterTool::setParam( float freqLowPass, float slopeLowPass, float freqHighPass, float slopeHighPass ) {
  int numLowPassFilters  = freqLowPass  > 0.0 ? 1 : 0;
  int numHighPassFilters = freqHighPass > 0.0 ? 1 : 0;
  myIsFilter = false;
  if( numLowPassFilters > 0 || numHighPassFilters != 0 ) {
    myIsFilter = true;
    prepareBandpassFilter( numLowPassFilters, &freqLowPass, &slopeLowPass, numHighPassFilters, &freqHighPass, &slopeHighPass );
  }
}

float* csFilterTool::retrieveSamplesPointer() {
  return mySamples;
}

void csFilterTool::applyFilter() {
  if( myIsFilter ) {
    applyBandpassFilter( mySamples, false );
  }
}
