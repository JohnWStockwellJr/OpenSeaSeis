/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TRACE_DATA_H
#define CS_TRACE_DATA_H

#include <cstdlib>
#include "csException.h"

namespace cseis_system {

/**
* Trace samples/trace data
*
* Manages seismic trace samples for one trace
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTraceData {
public:
  csTraceData();
  csTraceData( int numSamples );
  ~csTraceData();
  /// Return data samples
  inline float* getSamples() { return myDataSamples; }
  /// Return number of samples
  inline int numSamples() const { return myNumSamples; }
  /// Return number of samples
  inline int getNumAllocatedSamples() const { return myNumAllocatedSamples; }
  /// Set data samples
  void setData( csTraceData const* data );
  inline void setData( float const* samples ) {
    setData( samples, myNumSamples );
  }
  void setData( float const* samples, int nSamples );
  void setData( char const* samples, int nSamples );
  inline float& operator [] ( int index ) {  // May throw exception
    if( index >= 0 && index < myNumSamples ) {
      return myDataSamples[index];
    }
    throw cseis_geolib::csException("Wrong sample index passed to trace");
  }
  void trim();
  friend class csMemoryPoolManager;
  friend class csModule;
  friend class csTrace;
private:
  float* myDataSamples;
  int myNumSamples;
  int myNumAllocatedSamples;
  bool myDoTrimOnNextCall;

  /// Set number of samples to maximum between numSamples passed as argument and numSamples as currently set
  inline void setMax( int numSamplesNew ) {
    set( std::max(numSamplesNew,myNumSamples), 0 );
  }
  inline void set( int numSamplesNew ) {
    set( numSamplesNew, 0 );
  }
  void set( int numSamplesNew, int firstLiveSample );
};

} // namespace
#endif

