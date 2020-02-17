/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <cstdio>
#include "csTaper.h"

using namespace std;
using namespace cseis_geolib;

csTaper::csTaper() {
  myTaperType = 0;
  myNormFlag  = 0;
  myNumSamplesIn = 0;
  myNumSamplesTaper = 0;
}
csTaper::csTaper( int numSamplesIn, int numSamplesTaper, int taperType, int normFlag ) {
  set( numSamplesIn, numSamplesTaper, taperType, normFlag );
}
void csTaper::set( int numSamplesIn, int numSamplesTaper, int taperType, int normFlag ) {
  myTaperType = taperType;
  myNormFlag  = normFlag;
  myNumSamplesIn = numSamplesIn;
  myNumSamplesTaper = numSamplesTaper;
}

csTaper::~csTaper() {
}

void csTaper::applyTaper( float* dataSamples, int numSamplesIn, int numSamplesTaper, int taperType, int normFlag, float& normScalar ) {
  set( numSamplesIn, numSamplesTaper, taperType, normFlag );
  apply( dataSamples, normScalar );
}

void csTaper::apply( float* dataSamples, float& normScalar ) {
  if( myTaperType == cseis_geolib::csTaper::COSINE || myTaperType == cseis_geolib::csTaper::HANNING ) {
    for( int i = 0; i < myNumSamplesTaper; i++ ) {
      float scalar = 0.5 * ( 1.0 + cos( M_PI * (float)(myNumSamplesTaper-i)/(float)myNumSamplesTaper ) );
      dataSamples[i] *= scalar;
    }
    int temp = myNumSamplesTaper - myNumSamplesIn + 1;
    for( int i = myNumSamplesIn-myNumSamplesTaper; i < myNumSamplesIn; i++ ) {
      float scalar = 0.5 * (1.0 + cos( M_PI * (float)(temp+i)/(float)myNumSamplesTaper ) );
      dataSamples[i] *= scalar;
    }
    if( myNormFlag != cseis_geolib::csTaper::NORM_NO ) {
      float weightSum = 0.0;
      int weightCounter = 0;
      for( int i = 0; i < myNumSamplesTaper; i++ ) {
        if( myNormFlag == cseis_geolib::csTaper::NORM_ZERO && dataSamples[i] == 0 ) continue;
        float scalar = cos( M_PI_2 * (float)(myNumSamplesTaper-i)/(float)myNumSamplesTaper );
        weightSum += scalar;
        weightCounter += 1;
      }
      for( int i = myNumSamplesIn-myNumSamplesTaper; i < myNumSamplesIn; i++ ) {
        if( myNormFlag == cseis_geolib::csTaper::NORM_ZERO && dataSamples[i] == 0 ) continue;
        float scalar = cos( M_PI_2 * (float)(myNumSamplesTaper-myNumSamplesIn+i+1)/(float)myNumSamplesTaper );
        weightSum += scalar;
        weightCounter += 1;
      }
      normScalar = 1.0 / (weightSum+(1.0/normScalar)-weightCounter);
    }
  }
  else if( myTaperType == cseis_geolib::csTaper::BLACKMAN ) {
    float alpha = 0.16;
    float a0 = 0.5 * (1.0 - alpha);
    float a1 = 0.5;
    float a2 = 0.5 * alpha;
    for( int i = 0; i < myNumSamplesIn; i++ ) {
      float piFactor = (2.0 * M_PI) * (float)i / (float)(myNumSamplesIn - 1);
      float weight = a0 - a1*cos( piFactor ) + a2*cos( 2 * piFactor );
      dataSamples[i] *= weight;
    }
    if( myNormFlag != cseis_geolib::csTaper::NORM_NO ) {
      float weightSum = 0.0;
      int weightCounter = 0;
      for( int i = 0; i < myNumSamplesIn; i++ ) {
        if( myNormFlag == cseis_geolib::csTaper::NORM_ZERO && dataSamples[i] == 0 ) continue;
        float piFactor = (2.0 * M_PI) * (float)i / (float)(myNumSamplesIn - 1);
        float weight = a0 - a1*cos( piFactor ) + a2*cos( 2 * piFactor );
        weightSum += weight;
        weightCounter += 1;
      }
      normScalar = 1.0 / (weightSum+(1.0/normScalar)-weightCounter);
    }
  }
}
