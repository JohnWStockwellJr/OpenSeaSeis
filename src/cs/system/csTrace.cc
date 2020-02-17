/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csTrace.h"
#include "csTraceHeader.h"
#include "csTraceData.h"
#include "csTracePool.h"
#include "csException.h"
#include <cstring>

using namespace cseis_system;

int csTrace::myIdentCounter = 0;

//---------------------------------------------------------
csTrace::csTrace( csTracePool* const tracePoolPtr ) :
  myTracePoolPtr( tracePoolPtr ),
  myIdentNumber( myIdentCounter++ )
{
  myTraceHeader = new csTraceHeader();
  myData        = new csTraceData();
}
//---------------------------------------------------------
csTrace::csTrace() :
  myTracePoolPtr( NULL ),
  myIdentNumber( myIdentCounter++ )
{
  myTraceHeader = new csTraceHeader();
  myData        = new csTraceData();
}
//---------------------------------------------------------
csTrace::~csTrace() {
  if( myTraceHeader ) {
    delete myTraceHeader;
    myTraceHeader = NULL;
  }
  if( myData ) {
    delete myData;
    myData = NULL;
  }
}
//--------------------------------------------------
csTraceHeader* csTrace::getTraceHeader() {
  return myTraceHeader;
}
csTraceHeader const* csTrace::getTraceHeader() const {
  return myTraceHeader;
}
csTraceData const* csTrace::getTraceDataObject() const {
  return myData;
}
csTraceData* csTrace::getTraceDataObject() {
  return myData;
}
float* csTrace::getTraceSamples() {
  return myData->getSamples();
}
float const* csTrace::getTraceSamples() const {
  return myData->getSamples();
}
int csTrace::numSamples() const {
  return myData->numSamples();
}
int csTrace::numHeaders() const {
  return myTraceHeader->numHeaders();
}
//--------------------------------------------------
void csTrace::free() {
  if( myTracePoolPtr != NULL ) {
    myTracePoolPtr->freeTrace( this );
    myTraceHeader->clear();
  }
  else {  // TEMP
    throw cseis_geolib::csException("csTrace: ERROR, pool pointer is NULL");
  }
}
void csTrace::trim() {
  return myData->trim();
}


int csTrace::mpi_getByteSize() const {
  char* dataBuffer = new char[10000];
  int byteSize = mpi_compress( dataBuffer );
  delete [] dataBuffer;
  return byteSize;
}
int csTrace::mpi_compress( char* data ) const {
  /*
 csTraceData:
    *   float* myDataSamples;
    *   int myNumSamples;
    --  int myNumAllocatedSamples;
    --  bool myDoTrimOnNextCall;
 csTraceHeader:
    /// The actual trace header values:
    csTraceHeaderData* myTraceHeaderData;
    /// Constant pointer to trace header definition (managed somewhere else)
    -- csTraceHeaderDef const* myHeaderDefPtr;
csTraceHeaderData:
  /// Buffer that holds all trace header values in one chunk of memory
  char* myValueBlock;
  /// Maps the sequential header index 0,1,2... to the byte location index in the char* 'value block'. Pointer to object, kept in csTraceHeaderDef
  --   int const*  myByteLocationPtr;
  int   myNumBytes;
  int   myNumAllocatedBytes;
  int   myNumHeaders;
  // Helper method that sets header data according header definition object (Adds new headers = allocates enough memory)
  int myIndex;
  --  static int counter;

  */

  int numSamples = myData->numSamples();
  int numBytesHdrValueBlock = myTraceHeader->getHeaderDefPtr()->getTotalNumBytes();
    
  int byteLoc = 0;
  memcpy( &data[byteLoc], &numSamples, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], &numBytesHdrValueBlock, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], myData->myDataSamples, numSamples*sizeof(float) );
  byteLoc += numSamples*sizeof(float);
  memcpy( &data[byteLoc], myTraceHeader->getTraceHeaderValueBlock(), numBytesHdrValueBlock*sizeof(char) );
  byteLoc += numBytesHdrValueBlock*sizeof(char);

  return byteLoc;
}

void csTrace::mpi_decompress( char const* data ) {
  int numSamples;
  int numBytesHdrValueBlock;

  int byteLoc = 0;
  memcpy( &numSamples, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &numBytesHdrValueBlock, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);

  myData->setData( &data[byteLoc], numSamples );
  byteLoc += numSamples*sizeof(float);
  myTraceHeader->setTraceHeaderValueBlock( &data[byteLoc], numBytesHdrValueBlock );
  byteLoc += numBytesHdrValueBlock*sizeof(char);
}
