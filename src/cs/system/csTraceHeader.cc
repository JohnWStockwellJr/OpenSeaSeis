/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include "csTraceHeader.h"
#include "csException.h"
#include "csGeolibUtils.h"
#include "geolib_defines.h"
#include <cstring>

using namespace cseis_system;

csTraceHeader::csTraceHeader() :
  myTraceHeaderData( NULL ),
  myHeaderDefPtr( NULL )
{
  myTraceHeaderData = new csTraceHeaderData();
}
//----------------------------------------------------------------------
//
void csTraceHeader::setHeaders( csTraceHeaderDef const* traceHeaderDef, int inPort ) {
  myHeaderDefPtr = traceHeaderDef;
  myTraceHeaderData->setHeaders( traceHeaderDef, inPort );
}
void csTraceHeader::copyFrom( csTraceHeader const* hdr ) {
  myTraceHeaderData->copyFrom( hdr->myTraceHeaderData );
  myHeaderDefPtr = hdr->myHeaderDefPtr;  // BUGFIX 080630: Added this line. Previously, hdef pointer was not set correctly (how did this ever work..?)
}
//----------------------------------------------------------------------
//
csTraceHeader::~csTraceHeader() {
  if( myTraceHeaderData != NULL ) {
    delete myTraceHeaderData;
    myTraceHeaderData = NULL;
  }
}
//---------------------------------------------------------------------
//
void csTraceHeader::clear() {
  if( myTraceHeaderData != NULL ) {
    myTraceHeaderData->clear();
  }
}
void csTraceHeader::clearMemory() {
  if( myTraceHeaderData != NULL ) {
    myTraceHeaderData->clearMemory();
  }
}
//------------------------------------------------------
char const* csTraceHeader::getTraceHeaderValueBlock() const {
  return myTraceHeaderData->myValueBlock;
}
//------------------------------------------------------
char* csTraceHeader::getTraceHeaderValueBlockHandle() {
  return myTraceHeaderData->myValueBlock;
}
void csTraceHeader::setTraceHeaderValueBlock( char const* hdrValueBlock, int byteSize ) {
  if( byteSize > myTraceHeaderData->myNumBytes ) {
    throw( cseis_geolib::csException("csTraceHeader::setTraceHeaderValueBlock: Program bug. Incorrect byte size (%d != %d).", byteSize, myTraceHeaderData->myNumBytes) );
  }
  memcpy( myTraceHeaderData->myValueBlock, hdrValueBlock, byteSize );

  // TEMP:
  /*
  int value1;
  int value2;
  int index = 2;
  memcpy( &value1, &hdrValueBlock[myTraceHeaderData->myByteLocationPtr[index]], 4 );
  memcpy( &value2, &myTraceHeaderData->myValueBlock[myTraceHeaderData->myByteLocationPtr[index]], 4 );
  fprintf(stderr,"csTraceHeader::setTraceHeaderValueBlock(): numHdr: %d, numBytes: %d; trcno values: %d ?= %d\n",
          myTraceHeaderData->numHeaders(), myTraceHeaderData->numBytes(), value1, value2 );
  for( int i = 0; i < myTraceHeaderData->numHeaders(); i++ ) {
    memcpy( &value1, &hdrValueBlock[myTraceHeaderData->myByteLocationPtr[i]], 4 );
    fprintf(stderr,"  byteLocPtr #%d: %d - INTvalue: %d\n", i, myTraceHeaderData->myByteLocationPtr[i], value1 );
  }
  */
}
void csTraceHeader::writeTraceHeaderValueBlock( char const* hdrValueBlock_in, int const* byteMap_in, int const* hdrMap_out, int numHeaders_in ) {
  char* hdrValueBlock = myTraceHeaderData->myValueBlock;
  for( int ihdr = 0; ihdr < numHeaders_in; ihdr++ ) {
    int byteLoc_in   = byteMap_in[ihdr];
    int hdrIndex_out = hdrMap_out[ihdr];
    memcpy( &hdrValueBlock[myHeaderDefPtr->getByteLocation(hdrIndex_out)], &hdrValueBlock_in[byteLoc_in],
            cseis_geolib::csGeolibUtils::numBytes( myHeaderDefPtr->headerType(hdrIndex_out) ) );
  }
}
void csTraceHeader::readTraceHeaderValueBlock( char* hdrValueBlock_out, int const* byteMap_out, int const* hdrMap_in, int numHeaders_out ) const {
  char* hdrValueBlock = myTraceHeaderData->myValueBlock;
  for( int ihdr = 0; ihdr < numHeaders_out; ihdr++ ) {
    int byteLoc_out  = byteMap_out[ihdr];
    int hdrIndex_in  = hdrMap_in[ihdr];
    memcpy( &hdrValueBlock_out[byteLoc_out], &hdrValueBlock[myHeaderDefPtr->getByteLocation(hdrIndex_in)],
            cseis_geolib::csGeolibUtils::numBytes( myHeaderDefPtr->headerType(hdrIndex_in) ) );
  }
}
//--------------------------------------------------
//
void csTraceHeader::dump( std::FILE* stream ) const {
  int nHeaders = numHeaders();
  if( stream == NULL ) stream = stdout;
  for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
    cseis_geolib::type_t theType = type(ihdr);
    fprintf( stream, "#%03d %20s %d ", ihdr+1, name(ihdr).c_str(), theType );
    switch( theType ) {
    case cseis_geolib::TYPE_INT:
      fprintf( stream, "%d\n", intValue(ihdr) );
      break;
    case cseis_geolib::TYPE_INT64:
      fprintf( stream, "%lld\n", int64Value(ihdr) );
      break;
    case cseis_geolib::TYPE_FLOAT:
      fprintf( stream, "%f\n", floatValue(ihdr) );
      break;
    case cseis_geolib::TYPE_DOUBLE:
      fprintf( stream, "%f\n", doubleValue(ihdr) );
      break;
    case cseis_geolib::TYPE_CHAR:
    case cseis_geolib::TYPE_STRING:
      fprintf( stream, "%s\n", stringValue(ihdr).c_str() );
      break;
    case cseis_geolib::TYPE_VECTOR:
      fprintf( stream, "%s\n", vectorValue(ihdr).toString().c_str() );
      break;
    }
  }
}

