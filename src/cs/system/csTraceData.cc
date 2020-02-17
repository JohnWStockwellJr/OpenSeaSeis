/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csTraceData.h"
#include "csException.h"
#include <string>
#include <cstdio>
#include <cstring>

using namespace cseis_system;

csTraceData::csTraceData() {
  myNumSamples = 0;
  myNumAllocatedSamples = 0;
  myDataSamples = NULL;
  myDoTrimOnNextCall = false;
}
csTraceData::csTraceData( int numSamples ) {
  myNumSamples = numSamples;
  myNumAllocatedSamples = numSamples;
  myDataSamples = new float[myNumAllocatedSamples];
  myDoTrimOnNextCall = false;
}
csTraceData::~csTraceData() {
  if( myDataSamples ) {
    delete [] myDataSamples;
    myDataSamples = NULL;
  }
}
//---------------------------------------------------------------------------
//
void csTraceData::setData( csTraceData const* data ) {
  if( myNumSamples != data->myNumSamples ) {
    // BUGFIX 080630: Previously, no check was made whether this data object had the same number of samples. This lead to data objects with 0 numSamples etc.
    set( data->myNumSamples );
  }
  setData( data->myDataSamples, data->myNumSamples );
}
void csTraceData::setData( float const* samples, int nSamples ) {
  memcpy( myDataSamples, samples, std::min(nSamples,myNumSamples)*sizeof(float) );
}
void csTraceData::setData( char const* samples, int nSamples ) {
  memcpy( myDataSamples, samples, std::min(nSamples,myNumSamples)*sizeof(float) );
}
void csTraceData::trim() {
  //  printf("---Trimmed from %d to %d samples\n", myNumAllocatedSamples, myNumSamples );
  myDoTrimOnNextCall = true;
}
void csTraceData::set( int numSamplesNew, int firstLiveSample ) {
  //  if( myDoTrimOnNextCall ) printf("Trimmed from %d to %d to %d samples\n", myNumAllocatedSamples, myNumSamples, numSamplesNew );
  if( numSamplesNew > myNumAllocatedSamples || myDoTrimOnNextCall ) {
    float* dataNew = NULL;
    try {
      dataNew = new float[numSamplesNew];
    }
    catch(...) {
      throw( cseis_geolib::csException("csTraceData::set: Unable to allocate new trace data buffer. Out of memory.") );
    }
    if( myDataSamples ) {
      memcpy( dataNew, myDataSamples, std::min( myNumSamples, numSamplesNew )*sizeof(float) );
      delete [] myDataSamples;
    }
    myDataSamples = dataNew;
    myNumAllocatedSamples = numSamplesNew;
    myNumSamples = numSamplesNew;
    myDoTrimOnNextCall = false;
  }
  if( numSamplesNew > myNumSamples ) {
    memset( &myDataSamples[myNumSamples], 0, (numSamplesNew-myNumSamples)*sizeof(float) );
    myNumSamples = numSamplesNew;
  }
  if( firstLiveSample > 0 ) {
    if( firstLiveSample > myNumSamples ) {
      firstLiveSample = myNumSamples;
    }
    memset( &myDataSamples[0], 0, firstLiveSample*sizeof(float) );
  }
}

