/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csGeneralSeismicReader.h"
#include "csSeismicReader_ver.h"
#include "csSeismicIOConfig.h"
#include "csHeaderInfo.h"
#include "csException.h"
#include "csGeolibUtils.h"

#include "csFlexHeader.h"
#include "csFlexNumber.h"
#include "csIOSelection.h"

#include <cstring>

using namespace cseis_io;

csGeneralSeismicReader::csGeneralSeismicReader( std::string filename, bool enableRandomAccess, int numTracesToBuffer ) {
  try {
    myReader = csSeismicReader_ver::createReaderObject( filename, enableRandomAccess, numTracesToBuffer );
  }
  catch( cseis_geolib::csException& e ) {
    throw( cseis_geolib::csException( e.getMessage()) );
  }
  myConfig = new csSeismicIOConfig();
  myHdrValueBlock = NULL;
  myTraceBuffer   = NULL;

  myIntPtr   = NULL;
  myFloatPtr = NULL;

  myByteLocation     = NULL;
  myIsFileHeaderRead = false;
  
  myHdrCheckByteOffset = 0;
  myHdrCheckType       = cseis_geolib::TYPE_UNKNOWN;
  myHdrCheckByteSize   = 0;
  myHdrCheckBuffer     = NULL;
  myIOSelection        = NULL;
}

csGeneralSeismicReader::~csGeneralSeismicReader() {
  if( myIOSelection != NULL ) {
    delete myIOSelection;
    myIOSelection = NULL;
  }
  if( myByteLocation != NULL ) {
    delete [] myByteLocation;
    myByteLocation = NULL;
  }
  if( myHdrValueBlock != NULL ) {
    delete [] myHdrValueBlock;
    myHdrValueBlock = NULL;
  }
  if( myReader != NULL ) {
    delete myReader;
    myReader = NULL;
  }
  if( myConfig != NULL ) {
    delete myConfig;
    myConfig = NULL;
  }
  if( myTraceBuffer != NULL ) {
    delete [] myTraceBuffer;
    myTraceBuffer = NULL;
  }
  if( myHdrCheckBuffer != NULL ) {
    delete [] myHdrCheckBuffer;
    myHdrCheckBuffer = NULL;
  }    
}

//--------------------------------------------------------------------
void csGeneralSeismicReader::closeFile() {
  myReader->closeFile();
}

//--------------------------------------------------------------------
bool csGeneralSeismicReader::readFileHeader() {

  if( myIsFileHeaderRead ) return false;
  
  bool success = myReader->readFileHeader( myConfig );
  
  if( !success ) return false;

  myHdrValueBlock = new char[myConfig->byteSizeHdrValueBlock];
  myIntPtr   = reinterpret_cast<int*>( myHdrValueBlock );
  myFloatPtr = reinterpret_cast<float*>( myHdrValueBlock );

  setByteLocation();

  myIsFileHeaderRead = true; 
  
  return true;
}

void csGeneralSeismicReader::setByteLocation() {
  if( myByteLocation != NULL ) {
    delete [] myByteLocation;
    myByteLocation = NULL;
  }
  int nHeaders = numTraceHeaders();
  myByteLocation = new int[nHeaders];
  myByteLocation[0] = 0;
  int numBytes = 0;
  for( int ihdr = 1; ihdr < nHeaders; ihdr++ ) {
    cseis_geolib::csHeaderInfo const* info = myConfig->headerInfo( ihdr-1 );
    if( info->type != cseis_geolib::TYPE_STRING ) {
      numBytes += cseis_geolib::csGeolibUtils::numBytes( info->type );
    }
    else { // if( type == TYPE_STRING ) {
      numBytes += info->nElements;
    }
    myByteLocation[ihdr] = numBytes;
  }
}

int csGeneralSeismicReader::numTraceHeaders() const {
  return myConfig->numTrcHeaders();
}

int csGeneralSeismicReader::headerIndex( std::string const& headerName ) const {
  return myConfig->headerIndex( headerName );
}

std::string csGeneralSeismicReader::headerName( int hdrIndex ) const {
  if( hdrIndex < 0 || hdrIndex >= myConfig->numTrcHeaders() ) return NULL;
  //  cseis_geolib::csHeaderInfo const* info = myConfig->headerInfo( hdrIndex )->name;
  return myConfig->headerInfo( hdrIndex )->name;
}

std::string csGeneralSeismicReader::headerDesc( int hdrIndex ) const {
  if( hdrIndex < 0 || hdrIndex >= myConfig->numTrcHeaders() ) return NULL;
  //  cseis_geolib::csHeaderInfo const* info = myConfig->headerInfo( hdrIndex );
  //  return info->description;
  return myConfig->headerInfo( hdrIndex )->description;
}

cseis_geolib::type_t csGeneralSeismicReader::headerType( int hdrIndex ) const {
  if( hdrIndex < 0 || hdrIndex >= myConfig->numTrcHeaders() ) return cseis_geolib::TYPE_UNKNOWN;
  cseis_geolib::csHeaderInfo const* info = myConfig->headerInfo( hdrIndex );
  return info->type;
}

int csGeneralSeismicReader::headerElements( int hdrIndex ) const {
  if( hdrIndex < 0 || hdrIndex >= myConfig->numTrcHeaders() ) {
    return 0;
  }
  else {  
    return myConfig->headerInfo( hdrIndex )->nElements;
  }  
}

int csGeneralSeismicReader::hdrIntValue( int hdrIndex ) const {
  int value;
  memcpy( &value, &myHdrValueBlock[myByteLocation[hdrIndex]], sizeof(int) );
  return value;
  //  return *(reinterpret_cast<int*>( &myHdrValueBlock[myByteLocation[hdrIndex]] ));
}

float csGeneralSeismicReader::hdrFloatValue( int hdrIndex ) const {
  float value;
  memcpy( &value, &myHdrValueBlock[myByteLocation[hdrIndex]], sizeof(float) );
  return value;
  //  return *(reinterpret_cast<float*>( &myHdrValueBlock[myByteLocation[hdrIndex]] ));
}

double csGeneralSeismicReader::hdrDoubleValue( int hdrIndex ) const {
   double value;
   memcpy( &value, &myHdrValueBlock[myByteLocation[hdrIndex]], sizeof(double) );
   return value;
  //  return *(reinterpret_cast<double*>( &myHdrValueBlock[myByteLocation[hdrIndex]] ));
}

csInt64_t csGeneralSeismicReader::hdrInt64Value( int hdrIndex ) const {
  //return 0;
  csInt64_t value;
  memcpy( &value, &myHdrValueBlock[myByteLocation[hdrIndex]], sizeof(csInt64_t) );
  return value;
  //  return *(reinterpret_cast<csInt64_t*>( &myHdrValueBlock[myByteLocation[hdrIndex]] ));
}

std::string csGeneralSeismicReader::hdrStringValue( int hdrIndex ) const {
  int size;
  if( hdrIndex == numTraceHeaders() - 1 ) {
    size = myConfig->byteSizeHdrValueBlock - myByteLocation[hdrIndex];
  }
  else {
    size = myByteLocation[hdrIndex + 1] - myByteLocation[hdrIndex];
  }
  std::string text( &myHdrValueBlock[myByteLocation[hdrIndex]], size );
  return text;
  //  return "";
}

bool csGeneralSeismicReader::moveToTrace( int traceIndex ) {
  return myReader->moveToTrace( traceIndex );
}

bool csGeneralSeismicReader::moveToTrace( int traceIndex, int numTracesToRead ) {
  return myReader->moveToTrace( traceIndex, numTracesToRead );
}

float const* csGeneralSeismicReader::readTraceReturnPointer() {
  if( myTraceBuffer == NULL ) {
    myTraceBuffer = new float[myConfig->numSamples];
  }
  if( myIOSelection ) {
    if( !performIOSelection() ) return NULL;
  }
  bool success = readTrace( myTraceBuffer );
  if( success ) {
    return myTraceBuffer;
  }
  else {
    return NULL;
  }
}

bool csGeneralSeismicReader::readTrace( float* samples ) {
  if( myIOSelection ) {
    if( !performIOSelection() ) return false;
  }
  return myReader->readTrace( samples, myHdrValueBlock );
}

bool csGeneralSeismicReader::performIOSelection() {
  bool success = false;
  int traceIndex = myIOSelection->getNextTraceIndex();
  if( traceIndex >= 0 ) {
    // If traceIndex progressed more than 1 --> Move to new trace
    if( traceIndex-1 != myIOSelection->getCurrentTraceIndex() ) {
      success = myReader->moveToTrace( traceIndex );
    }
  }
  return( success );
}

//--------------------------------------------------------------------
int csGeneralSeismicReader::numTraces() const {
 return myReader->numTraces();
}   

bool csGeneralSeismicReader::setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType ) {
  bool success = setHeaderToPeek( headerName );
  headerType = myHdrCheckType;
  return success;
}
bool csGeneralSeismicReader::setHeaderToPeek(std::string const& headerName) {
  int idx = headerIndex(headerName);    // returns -1 on failure

  if (idx < 0) {
    return false;
  }
  else {
    myReader->revertFromPeekPosition();
    myHdrCheckType       = headerType(idx);
    myHdrCheckByteOffset = myByteLocation[idx]; 
    int hdrCheckByteSize   = headerElements(idx) * cseis_geolib::csGeolibUtils::numBytes(myHdrCheckType);      
    
    if( hdrCheckByteSize != myHdrCheckByteSize && myHdrCheckBuffer != NULL ) {
      delete [] myHdrCheckBuffer;
      myHdrCheckBuffer = NULL;
    }
    if( myHdrCheckBuffer == NULL ) {
      myHdrCheckByteSize = hdrCheckByteSize;
      myHdrCheckBuffer = new char[myHdrCheckByteSize];
    }
    return true;
  }      
}

bool csGeneralSeismicReader::peekHeaderValue(cseis_geolib::csFlexHeader* hdrValue, int traceIndex) {

  if (myHdrCheckBuffer == NULL) {
    throw(cseis_geolib::csException("csGeneralSeismicReader::peekHeaderValue: No header has been set for checking. This is a program bug in the calling function"));
  }
  
  bool success = true;
  if (traceIndex < 0) {
    success = myReader->peek(myHdrCheckByteOffset, myHdrCheckByteSize, myHdrCheckBuffer);
  }
  else {
    success = myReader->peek(myHdrCheckByteOffset, myHdrCheckByteSize, myHdrCheckBuffer, traceIndex);
  }
  
  if (success) {
    if( myHdrCheckType == cseis_geolib::TYPE_FLOAT ) {
      float value = *(reinterpret_cast<float*>( myHdrCheckBuffer ));    
      hdrValue->setFloatValue( value );
    }
    else if( myHdrCheckType == cseis_geolib::TYPE_INT ) {
      int value = *(reinterpret_cast<int*>( myHdrCheckBuffer ));      
      hdrValue->setIntValue( value );
    }
    else if( myHdrCheckType == cseis_geolib::TYPE_DOUBLE ) {
      double value = *(reinterpret_cast<double*>( myHdrCheckBuffer ));
      hdrValue->setDoubleValue( value );
    }
  }
  
  return success;
}

bool csGeneralSeismicReader::setSelection( std::string const& hdrValueSelectionText,
                                           std::string const& headerName,
                                           int sortOrder,
                                           int sortMethod )
{
  myIOSelection = new cseis_geolib::csIOSelection( headerName, sortOrder, sortMethod );
  return myIOSelection->initialize( this, hdrValueSelectionText );
}
void csGeneralSeismicReader::setSelectionStep1( std::string const& hdrValueSelectionText,
                                                std::string const& headerName,
                                                int sortOrder,
                                                int sortMethod )
{
  myIOSelection = new cseis_geolib::csIOSelection( headerName, sortOrder, sortMethod );
  myIOSelection->step1( this, hdrValueSelectionText );
}
bool csGeneralSeismicReader::setSelectionStep2( int& numTracesToRead )
{
  return myIOSelection->step2( this, numTracesToRead );
}
bool csGeneralSeismicReader::setSelectionStep3()
{
  return myIOSelection->step3( this );
}


int csGeneralSeismicReader::getNumSelectedTraces() const {
  return myIOSelection->getNumSelectedTraces();
}
cseis_geolib::csFlexNumber const* csGeneralSeismicReader::getSelectedValue( int traceIndex ) const {
  return myIOSelection->getSelectedValue( traceIndex );
}
int csGeneralSeismicReader::getSelectedIndex( int traceIndex ) const {
  return myIOSelection->getSelectedIndex( traceIndex );
}
int csGeneralSeismicReader::getCurrentTraceIndex() const {
  return myReader->currentTraceIndex();
}
