/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSUTraceManager.h"
#include "csSegyTraceHeader.h"
#include "csSegyHdrMap.h"
#include "geolib_endian.h"
#include <cstdio>
#include <cstring>
#include <stdarg.h>

using namespace cseis_su;
using namespace std;
 
csSUTraceManager::csSUTraceManager() {
  myStatus     = STATUS_EMPTY;
  myIsEOF      = false;
  myIsError    = false;
  myTrace      = -1;
  myLogFilePtr = NULL;
  mySleepMilliSec = 2.0f;
  mySUDoc      = "";
  myBuffer     = NULL;
  myIsDocRequestOnly = false;
  reallocateBuffer( SU_NFLTS );
}

csSUTraceManager::~csSUTraceManager() {
  if( myBuffer != NULL ) {
    delete [] myBuffer;
    myBuffer = NULL;
  }
}

void csSUTraceManager::reallocateBuffer( int numSamples ) {
  if( myBuffer != NULL ) delete [] myBuffer;
  myNumSamples = numSamples;
  myBuffer = new unsigned char[cseis_su::HDRBYTES+myNumSamples*4];
  for( int i = 0; i < (cseis_su::HDRBYTES+myNumSamples*4); i++ ) {
    myBuffer[i] = 0;
  }
}
void csSUTraceManager::setLogFile( FILE* logFile ) {
  myLogFilePtr = logFile;
}
void csSUTraceManager::setDocRequestOnly() {
  myIsDocRequestOnly = true;
}
bool csSUTraceManager::isDocRequestOnly() const {
  return myIsDocRequestOnly;
}

int csSUTraceManager::putTrace( cseis_geolib::csSegyTraceHeader const* suTrcHdr, float const* samplesPtr, int numSamples ) {
  while( myStatus == STATUS_TRACE_WAITING && !myIsEOF ) {
    csSUTraceManager::sleep_ms( mySleepMilliSec );
  }
  if( myIsEOF ) return SU_FALSE;

  suTrcHdr->writeHeaderValues(myBuffer,false,true);   // swapEndian=false, isAutoScaleHeaders=true
  memcpy( &myBuffer[240], samplesPtr, numSamples*sizeof(float) );
  // End: copy trace

  myStatus = STATUS_TRACE_WAITING;
  // fprintf(stdout,"Manager: Put trace ...\n");
  return SU_TRUE;
}


int csSUTraceManager::getTracePtr( unsigned char const** bufferPtr ) {
  while( myStatus == STATUS_EMPTY && !myIsEOF ) {
    csSUTraceManager::sleep_ms( mySleepMilliSec );
  }
  if( myStatus == STATUS_EMPTY && myIsEOF ) return SU_FALSE;
  //  fprintf(stdout,"Pulling the buffer... %x\n", myBuffer);
  *bufferPtr = myBuffer;
  return SU_TRUE;
}
void csSUTraceManager::freeTrace() {
  myStatus = STATUS_EMPTY;
}

int csSUTraceManager::putTrace( segy const* trace ) {
  while( myStatus == STATUS_TRACE_WAITING && !myIsEOF ) {
    csSUTraceManager::sleep_ms( mySleepMilliSec );
  }
  if( myIsEOF ) return SU_FALSE;

  memcpy( myBuffer, trace, 240+myNumSamples*sizeof(float) );

  myStatus = STATUS_TRACE_WAITING;
  //  fprintf(stdout,"Manager: Put trace...\n");
  return SU_TRUE;
}

int csSUTraceManager::getTrace( segy* trace ) {
  while( myStatus == STATUS_EMPTY && !myIsEOF ) {
    csSUTraceManager::sleep_ms( mySleepMilliSec );
  }
  if( myStatus == STATUS_EMPTY && myIsEOF ) return SU_FALSE;
  memcpy( trace, myBuffer, 240+myNumSamples*sizeof(float) );

  myStatus = STATUS_EMPTY;
  return SU_TRUE;
}

int csSUTraceManager::getTraceMaybe( int* trace ) {
  if( myStatus == STATUS_EMPTY ) return SU_FALSE;
  myStatus = STATUS_EMPTY;
  *trace = myTrace;
  //  fprintf(stdout,"Manager: Get trace %d\n", *trace);
  return SU_TRUE;
}

void csSUTraceManager::setEOF() {
  myIsEOF = true;
}

void csSUTraceManager::setError( char const* text, ... ) {
  if( myLogFilePtr ) {
    va_list argList;
    va_start( argList, text );
    vfprintf( myLogFilePtr, text, argList );
    fprintf( myLogFilePtr, "\n" );
  }
  myIsError = true;
  setEOF();
}

bool csSUTraceManager::isError() const {
  return myIsError;
}

bool csSUTraceManager::isEmpty() const {
  return( myStatus == STATUS_EMPTY );
}

void csSUTraceManager::setSUDoc( std::string& sdoc ) {
  mySUDoc = sdoc;
}
char const* csSUTraceManager::getSUDoc() const {
  return mySUDoc.c_str();
}

int csSUTraceManager::numSamples() const {
  if( isEmpty() ) return -1;
  segy const* trace = (segy*)myBuffer;
  return trace->ns;
}

float csSUTraceManager::sampleInt() const {
  if( isEmpty() ) return -1;
  segy const* trace = (segy*)myBuffer;
  return ( (float)trace->dt / 1000.0 );
}
