/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <cstring>
#include "csTimeStretch.h"
#include "csInterpolation.h"
#include "geolib_methods.h"
#include "csException.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace cseis_geolib;

csTimeStretch::csTimeStretch( double sampleInt_ms, int numSamples ) {
  init( sampleInt_ms, numSamples, csTimeStretch::SAMPLE_INTERPOL_SINC );
}
csTimeStretch::csTimeStretch( double sampleInt_ms, int numSamples, int methodSampleInterpolation ) {
  init( sampleInt_ms, numSamples, methodSampleInterpolation );
}
void csTimeStretch::init( double sampleInt_ms, int numSamples, int methodSampleInterpolation ) {
  myIsBottomStretch = false;
  myLayerInterpolationMethod = csTimeStretch::LAYER_INTERPOL_LIN;
  mySampleInt           = sampleInt_ms;
  myNumSamples          = numSamples;
  myTimeOfLastSample    = (myNumSamples-1)*mySampleInt;
  mySampleInterpolationMethod = methodSampleInterpolation;

  myIndexBufferOut = NULL;
  myTimeDepthFunc  = NULL;
  myInterpol_timeDepthConversion = NULL;

  if( mySampleInterpolationMethod == csTimeStretch::SAMPLE_INTERPOL_SINC ) {
    mySampleInterpolation = new csInterpolation( numSamples, (float)sampleInt_ms );
  }
  else {
    mySampleInterpolation = NULL;
  }
}

csTimeStretch::~csTimeStretch() {
  if( mySampleInterpolation != NULL ) {
    delete mySampleInterpolation;
    mySampleInterpolation = NULL;
  }
  if( myIndexBufferOut != NULL ) {
    delete [] myIndexBufferOut; myIndexBufferOut = NULL;
  }
  if( myTimeDepthFunc != NULL ) {
    delete [] myTimeDepthFunc; myTimeDepthFunc = NULL;
  }
  if( myInterpol_timeDepthConversion != NULL ) {
    delete myInterpol_timeDepthConversion; myInterpol_timeDepthConversion = NULL;
  }
}

void csTimeStretch::setLayerInterpolationMethod( int method ) {
  myLayerInterpolationMethod = method;
}
//--------------------------------------------------------------------
//
//
void csTimeStretch::applyStretchFunction( float const* samplesIn, float const* tIn_ms, float const* stretch_ms, int numTimes, float* samplesOut )
{
  int numLayers = numTimes - 1;
  float* tOut_ms = new float[numTimes];
  tOut_ms[0] = tIn_ms[0];
  double stretchSumDbl = 0.0;
  for( int ilay = 0; ilay < numLayers; ilay++ ) {
    stretchSumDbl += (double)stretch_ms[ilay];
    tOut_ms[ilay+1] = tIn_ms[ilay+1] + (float)stretchSumDbl;
  }
  //  for( int ilay = 0; ilay < numLayers+1; ilay++ ) {
  //   fprintf(stderr,"STRETCH   %.4f %.4f\n", tIn_ms[ilay], tOut_ms[ilay]);
  //  }
  applyTimeInterval( samplesIn, tIn_ms, tOut_ms, numTimes, samplesOut );
  delete [] tOut_ms;
}

void csTimeStretch::applyTimeInterval( float const* samplesIn, float const* tIn, float const* tOut, int numTimes, float* samplesOut )
{
  int numLayers = numTimes - 1;
  int ilay = 0;
  // Traces always start at time = 0.0. First layer to stretch may start further down. If that is the case, don't stretch the top, just copy.
  float tTopIn  = 0.0;
  float tTopOut = 0.0;
  float tBotIn  = tIn[ilay];
  float tBotOut = tOut[ilay];
  float dtIn  = tBotIn  -  tTopIn;  // Input layer thickness in [ms]
  float dtOut = tBotOut - tTopOut; // Output layer thickness in [ms]
  float timeLast = (float)((myNumSamples-1) * mySampleInt);

  //  fprintf(stdout,"TIMEIO  START %f %f   %f %f  %d\n", tBotIn, tBotOut,dtIn,dtOut,ilay);
  for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
    float timeOut = (float)(isamp * mySampleInt);
    float timeIn  = timeOut;
    while( ilay <= numLayers && timeOut > tBotOut ) {
      ilay += 1;
      tTopIn = tBotIn;
      tTopOut = tBotOut;
      if( ilay <= numLayers ) {
        tBotIn  = tIn[ilay];
        tBotOut = tOut[ilay];
        dtIn  = tBotIn  -  tTopIn;  // Input layer thickness in [ms]
        dtOut = tBotOut - tTopOut; // Output layer thickness in [ms]
      }
      else {
        // For bottom of data beyond bottom of last layer: Make stretch ratio = dtIN/dtOut = 1.0 = no stretching, just static shift
        if( !myIsBottomStretch ) {
          dtOut = dtIn;
        }
        // ..otherwise, keep stretch factor from last specified layer
      }
    }
    if( dtOut != 0.0 ) {
      timeIn = tTopIn + (dtIn/dtOut) * ( timeOut - tTopOut );
      if( timeIn < 0.0 ) timeIn = 0.0;
      else if( timeIn > timeLast ) timeIn = timeLast;
    }
    //    fprintf(stdout,"TIMEIO  %f %f   %f %f  %d   %f %f\n", timeIn, timeOut,dtIn,dtOut,ilay,tBotIn,tBotOut);
    if( mySampleInterpolationMethod == csTimeStretch::SAMPLE_INTERPOL_SINC ) {
      samplesOut[isamp] = mySampleInterpolation->valueAt( timeIn, samplesIn );
    }
    else if( mySampleInterpolationMethod == csTimeStretch::SAMPLE_INTERPOL_QUAD ) {
      samplesOut[isamp] = getQuadAmplitudeAtSample( samplesIn, timeIn/mySampleInt, myNumSamples );
    }
    else {
      samplesOut[isamp] = getLinAmplitudeAtSample( samplesIn, timeIn/mySampleInt, myNumSamples );
    }
  }
}

//void csTimeStretch::applyDepthTimeConversion( float* samplesInOut, float const* velFunc, float sampleInt_velFunc_m, int numSamp_velFunc, float sampleInt_out ) {

void csTimeStretch::initialize_timeDepthConversion( int numSamplesOut, float sampleIntOut ) {
  if( myIndexBufferOut != NULL ) {
    delete [] myIndexBufferOut;
  }
  if( myTimeDepthFunc != NULL ) {
    delete [] myTimeDepthFunc;
  }
  myNumSamplesOut  = numSamplesOut;
  mySampleIntOut   = sampleIntOut;
  myIndexBufferOut = new float[myNumSamplesOut];
  myTimeDepthFunc  = new float[myNumSamplesOut];
  
  if( myInterpol_timeDepthConversion != NULL ) {
    delete myInterpol_timeDepthConversion;
  }
  myInterpol_timeDepthConversion = new csInterpolation( myNumSamples, csInterpolation::METHOD_SINC, mySampleInt );
}

void csTimeStretch::apply_timeDepthConversion( float const* samplesIn, float const* velFunc, float sampleIntVel_m, int numSamplesVel, float* samplesOut, bool isTime2Depth ) {
  if( myIndexBufferOut == NULL ) {
    throw( csException("csTimeStretch::apply_timeDepthConversion(%d):: Function not initialized. This is a bug in the calling function", isTime2Depth) );
  }
  float* depthTimeFunc = new float[numSamplesVel];
  // 1) Compute time-depth function from velocity-depth function
  depthTimeFunc[0] = 0;
  for( int isamp = 1; isamp < numSamplesVel; isamp++ ) {
    depthTimeFunc[isamp] = depthTimeFunc[isamp-1] + 2.0 * 1000.0 * sampleIntVel_m / velFunc[isamp-1];   // 2x for TWT
  }

  csInterpolation interpol( numSamplesVel, csInterpolation::METHOD_LIN, sampleIntVel_m );
  if( isTime2Depth ) {
    for( int isampOut = 0; isampOut < myNumSamplesOut; isampOut++ ) {
      float time = interpol.valueAt( isampOut*mySampleIntOut, depthTimeFunc );
      myIndexBufferOut[isampOut] = time/mySampleInt;
    }
  }
  else {
    interpol.timeAt_lin( depthTimeFunc, mySampleIntOut, myNumSamplesOut, myTimeDepthFunc );
    for( int isampOut = 0; isampOut < myNumSamplesOut; isampOut++ ) {
      myIndexBufferOut[isampOut] = myTimeDepthFunc[isampOut]/mySampleInt;
    }
  }

  myInterpol_timeDepthConversion->process( 1.0, 0.0, samplesIn, myIndexBufferOut, samplesOut, myNumSamplesOut );

  delete [] depthTimeFunc;
}
/*
void csTimeStretch::convert_depth2time( float const* samplesIn, float const* velFunc, float sampleIntVel_m, int numSamplesVel, float* samplesOut ) {
  if( myIndexBufferOut == NULL ) {
    throw( csException("csTimeStretch::convert_depth2time:: Function not initialized. This is a bug in the calling function") );
  }
  float* depthTimeFunc = new float[numSamplesVel];
  // 1) Compute time-depth function from velocity-depth function
  depthTimeFunc[0] = 0;
  for( int isamp = 1; isamp < numSamplesVel; isamp++ ) {
    depthTimeFunc[isamp] = depthTimeFunc[isamp-1] + 2.0 * 1000.0 * sampleIntVel_m / velFunc[isamp-1];   // 2x for TWT
  }

  csInterpolation interpol_timeDepth( numSamplesVel, csInterpolation::METHOD_LIN, sampleIntVel_m );
  interpol_timeDepth.timeAt_lin( depthTimeFunc, mySampleIntOut, myNumSamplesOut, myTimeDepthFunc );

  for( int isampOut = 0; isampOut < myNumSamplesOut; isampOut++ ) {
    myIndexBufferOut[isampOut] = myTimeDepthFunc[isampOut]/mySampleInt;
  }
  myInterpol_timeDepthConversion->process( 1.0, 0.0, samplesIn, myIndexBufferOut, samplesOut, myNumSamplesOut );

  delete [] depthTimeFunc;
}


void csTimeStretch::convert_time2depth( float const* samplesIn, float const* velFunc, float sampleIntVel_m, int numSamplesVel, float* samplesOut ) {
  if( myIndexBufferOut == NULL ) {
    throw( csException("csTimeStretch::convert_time2depth:: Function not initialized. This is a bug in the calling function") );
  }
  float* depthTimeFunc = new float[numSamplesVel];

  // 1) Compute time-depth function from velocity-depth function
  depthTimeFunc[0] = 0;
  for( int isamp = 1; isamp < numSamplesVel; isamp++ ) {
    depthTimeFunc[isamp] = depthTimeFunc[isamp-1] + 2.0 * 1000.0 * sampleIntVel_m / velFunc[isamp-1];   // 2x for TWT
  }

  csInterpolation interpol_depthTime( numSamplesVel, csInterpolation::METHOD_LIN, sampleIntVel_m );
  for( int isampOut = 0; isampOut < myNumSamplesOut; isampOut++ ) {
    float time = interpol_depthTime.valueAt( isampOut*mySampleIntOut, depthTimeFunc );
    myIndexBufferOut[isampOut] = time/mySampleInt;
  }

  myInterpol_timeDepthConversion->process( 1.0, 0.0, samplesIn, myIndexBufferOut, samplesOut, myNumSamplesOut );

  delete [] depthTimeFunc;
}
*/
