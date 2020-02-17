/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdlib>
#include <cstring>
#include "cseis_defines.h"
#include "csTraceHeaderData.h"
#include "csTraceHeaderDef.h"
#include "csTraceHeaderInfo.h"
#include "csGeolibUtils.h"
#include "csVector.h"
#include "csException.h"

using namespace cseis_system;

int csTraceHeaderData::counter = 0;

int csTraceHeaderData::NUM_ADD_HEADERS = 5;
int csTraceHeaderData::NUM_ADD_BYTES = 20;

csTraceHeaderData::csTraceHeaderData() {
  csTraceHeaderData::counter += 1;
  myIndex = csTraceHeaderData::counter;
  
  myNumHeaders = 0;
  myNumBytes   = 0;
  myNumAllocatedBytes   = 0;

  myValueBlock   = NULL;
  myByteLocationPtr = NULL;

  reallocateBytes( NUM_ADD_BYTES );
}
//----------------------------------------------------------------------
//
csTraceHeaderData::~csTraceHeaderData() {
  if( myValueBlock ) {
    delete [] myValueBlock;
    myValueBlock = NULL;
  }
}
//------------------------------------------------------------
//
void csTraceHeaderData::setHeaders( csTraceHeaderDef const* hdef, int inPort ) {
  int totalNumBytes = hdef->getTotalNumBytes();
  if( myNumBytes < totalNumBytes ) {
    reallocateBytes( totalNumBytes );    
  }
  myNumBytes = totalNumBytes;

  myByteLocationPtr = hdef->getHandleByteLocation();
  myNumHeaders      = hdef->numHeaders();
}
//----------------------------------------------------------------------
//
void csTraceHeaderData::deleteHeaders( csTraceHeaderDef const* hdef ) {
  cseis_geolib::csVector<int> const* indexOfHeadersToDel = hdef->getIndexOfHeadersToDel();
  int numHdrToDel = indexOfHeadersToDel->size();
  int const* byteLocation = hdef->getHandleByteLocation();

  int sumNumBytesToDel = 0;
  for( int iHdrDel = 0; iHdrDel < numHdrToDel; iHdrDel++ ) {
    int indexHdrDel = indexOfHeadersToDel->at(iHdrDel);
    csTraceHeaderInfo const* info = hdef->headerInfo(indexHdrDel);
    int numBytesHdrDel = cseis_geolib::csGeolibUtils::numBytes(info->type);  // or info->numBytes
    sumNumBytesToDel += numBytesHdrDel;
  }
  int numOldBytes = myNumBytes;
  int numNewBytes = myNumBytes-sumNumBytesToDel;

  char* temp_buffer = new char[numNewBytes];
  int bytePos = 0;
  int bytePosNew = 0;
  int bytePosCurrent = 0;

  for( int iHdrDel = 0; iHdrDel < numHdrToDel; iHdrDel++ ) {
    int indexHdrDel = indexOfHeadersToDel->at(iHdrDel);
    csTraceHeaderInfo const* info = hdef->headerInfo(indexHdrDel);
    int numBytesHdrDel = cseis_geolib::csGeolibUtils::numBytes(info->type);  // or info->numBytes
    bytePosNew = byteLocation[indexHdrDel];
    memcpy( &temp_buffer[bytePosCurrent], &myValueBlock[bytePos], bytePosNew-bytePos );
    bytePosCurrent += bytePosNew-bytePos;
    bytePos = bytePosNew + numBytesHdrDel;
  }
  if( numOldBytes-bytePos+bytePosCurrent != numNewBytes ) {
    throw cseis_geolib::csException("Error: Unequal number of bytes, %d != %d  (%d %d %d)\n",
                                    numOldBytes-bytePos+bytePosCurrent, numNewBytes,
                                    numOldBytes, bytePos, bytePosCurrent );
  }

  memcpy( &temp_buffer[bytePosCurrent], &myValueBlock[bytePos], numOldBytes-bytePos );
  memcpy( myValueBlock, temp_buffer, numNewBytes );
  delete [] temp_buffer;
}
//----------------------------------------------------------------------
//
void csTraceHeaderData::copyFrom( csTraceHeaderData const* hdr ) {
  myByteLocationPtr = hdr->myByteLocationPtr;
  myNumHeaders = hdr->myNumHeaders;
  int newNumBytes = hdr->myNumBytes;

  if( newNumBytes > myNumAllocatedBytes ) {
    reallocateBytes( newNumBytes );
  }
  myNumBytes = newNumBytes;
  memcpy( myValueBlock, hdr->myValueBlock, myNumBytes );
}
//----------------------------------------------------------------------
//
void csTraceHeaderData::reallocateBytes( int numBytesToAllocate ) {
  myNumAllocatedBytes  = numBytesToAllocate;
  char* newValueBlock = NULL;
  newValueBlock     = new char[myNumAllocatedBytes];
  memset( newValueBlock, 0, myNumAllocatedBytes );
  if( myNumBytes != 0 ) {
    memcpy( newValueBlock, myValueBlock, myNumBytes );
  }

  if( myValueBlock ) delete [] myValueBlock;
  myValueBlock = newValueBlock;
}
//---------------------------------------------------------------------
//
void csTraceHeaderData::clear() {
  myNumBytes = 0;
}
void csTraceHeaderData::clearMemory() {
  myNumBytes = 0;
}
//---------------------------------------------------------------------
//
void csTraceHeaderData::setIntValue( int index, int value ) {
#ifdef CS_DEBUG
  if( myByteLocationPtr[index] >= myNumBytes ) throw cseis_geolib::csException("csTraceHeaderData::setIntValue: myByteLocation is wrong...");
#endif
#ifdef ARCHITECTURE_ITANIUM
  memcpy( &myValueBlock[myByteLocationPtr[index]], &value, 4 );
#else
  *((int*)&myValueBlock[myByteLocationPtr[index]]) = value;
#endif

}
//---------------------------------------------------------------------
//
void csTraceHeaderData::setInt64Value( int index, csInt64_t value ) {
#ifdef CS_DEBUG
  if( myByteLocationPtr[index] >= myNumBytes ) throw cseis_geolib::csException("csTraceHeaderData::setInt64Value: myByteLocation is wrong...");
#endif
#ifdef ARCHITECTURE_ITANIUM
  memcpy( &myValueBlock[myByteLocationPtr[index]], &value, 8 );
#else
  *((csInt64_t*)&myValueBlock[myByteLocationPtr[index]]) = value;
#endif
}
//-----------------------------------------------------------
//
void csTraceHeaderData::setFloatValue( int index, float value ) {
#ifdef CS_DEBUG
  if( myByteLocationPtr[index] >= myNumBytes ) throw cseis_geolib::csException("csTraceHeaderData::setFloatValue: myByteLocation is wrong...");
#endif
#ifdef ARCHITECTURE_ITANIUM
  memcpy( &myValueBlock[myByteLocationPtr[index]], &value, 4 );
#else
  *((float*)&myValueBlock[myByteLocationPtr[index]]) = value;
#endif
}
//-----------------------------------------------------------
//
void csTraceHeaderData::setDoubleValue( int index, double value ) {
#ifdef CS_DEBUG
  if( myByteLocationPtr[index]/4 >= myNumBytes ) throw cseis_geolib::csException("csTraceHeaderData::setDoubleValue: myByteLocation is wrong...");
#endif
#ifdef ARCHITECTURE_ITANIUM
  memcpy( &myValueBlock[myByteLocationPtr[index]], &value, 8 );
#else
  *((double*)&myValueBlock[myByteLocationPtr[index]]) = value;  // ATTENTION: THIS WON'T WORK ON SGI!
#endif
}
//-----------------------------------------------------------
//
void csTraceHeaderData::setStringValue( int index, std::string value ) {
  int oldLength;
  int newLength = (int)value.length();
  if( index == myNumHeaders-1 ) {
    oldLength = myNumBytes - myByteLocationPtr[index];
  }
  else {
    oldLength = myByteLocationPtr[index+1] - myByteLocationPtr[index];
  }
  if( oldLength < newLength ) {
    newLength = oldLength;
  }
  memcpy( &myValueBlock[myByteLocationPtr[index]], value.c_str(), newLength );
}
std::string csTraceHeaderData::stringValue( int index ) const {
  int size;
  if( index == myNumHeaders-1 ) {
    size = myNumBytes - myByteLocationPtr[index];
  }
  else {
    size = myByteLocationPtr[index+1] - myByteLocationPtr[index];
  }
  std::string text( &myValueBlock[myByteLocationPtr[index]], size );
  return text;
}

void csTraceHeaderData::setVectorValueX( int index, double value ) {
  *((double*)&myValueBlock[myByteLocationPtr[index]])    = value;
}
void csTraceHeaderData::setVectorValueY( int index, double value ) {
  *((double*)&myValueBlock[myByteLocationPtr[index+8]])    = value;
}
void csTraceHeaderData::setVectorValueZ( int index, double value ) {
  *((double*)&myValueBlock[myByteLocationPtr[index+16]])   = value;
}
void csTraceHeaderData::setVectorValue( int index, double value, cseis_geolib::type_t type ) {
  int add = 0;
  if( type == cseis_geolib::TYPE_VECTOR_Y ) add = 8;
  else if( type == cseis_geolib::TYPE_VECTOR_Z ) add = 16;
  *((double*)&myValueBlock[myByteLocationPtr[index]+add])   = value;
}

void csTraceHeaderData::setVectorValue( int index, cseis_geolib::csPoint3D value ) {
  *((double*)&myValueBlock[myByteLocationPtr[index]])    = value.x;  // ATTENTION: THIS WON'T WORK ON SGI!
  *((double*)&myValueBlock[myByteLocationPtr[index]+8])  = value.y;
  *((double*)&myValueBlock[myByteLocationPtr[index]+16]) = value.z;
}
cseis_geolib::csPoint3D csTraceHeaderData::vectorValue( int index ) const {
  cseis_geolib::csPoint3D value;
  memcpy( &value.x, &myValueBlock[myByteLocationPtr[index]], 8 );
  memcpy( &value.y, &myValueBlock[myByteLocationPtr[index]+8], 8 );
  memcpy( &value.z, &myValueBlock[myByteLocationPtr[index]+16], 8 );
  return value;
}
double csTraceHeaderData::vectorValue( int index, cseis_geolib::type_t type ) const {
  double value;
  int add = 0;
  if( type == cseis_geolib::TYPE_VECTOR_Y ) add = 8;
  else if( type == cseis_geolib::TYPE_VECTOR_Z ) add = 16;
  memcpy( &value, &myValueBlock[myByteLocationPtr[index]+add], 8 );
  return value;
}
