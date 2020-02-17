/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGY_TRACE_H
#define CS_SEGY_TRACE_H

#include <string>
#include "csSegyDefines.h"
#include "csByteHeader.h"

namespace cseis_geolib {
  template<typename T> class csVector;

/**
 * Segy trace
 *
 * @author Bjorn Olofsson
 * @date 2008
 */
class csSegyTrace {
public:
  csSegyTrace( int numSamples ) {
    myNumSamples = numSamples;
    myDataBuffer = new float[myNumSamples];
    myHdr = new csByteHeader( cseis_segy::SIZE_TRCHDR );
  }
  float const* getSamples() {
    return myDataBuffer;
  }
  csByteHeader const* getHeader() const {
    return myHdr;
  }
  int getInt( int bytePos, int byteSize );
  float getFloat( int bytePos, int byteSize );
  double getDouble( int bytePos, int byteSize );
  std::string getString( int bytePos, int byteSize );

private:
  int myNumSamples;
  float* myDataBuffer;
  csByteHeader* myHdr;
};

} // end namespace
#endif

