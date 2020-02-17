/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csExternalHeader.h"
#include "csNavHeader.h"
#include "csGCS90Header.h"
#include "csNavInterface.h"

using namespace cseis_segd;
using std::endl;

csExternalHeader::csExternalHeader( int numBlocks, int navSystemID, int navInterfaceID ) :
  csBaseExternalHeader(numBlocks)
{
  myNavHeader    = new csNavHeader( navSystemID );
  myGCS90Header  = new csGCS90Header();
  myNavInterface = new csNavInterface( navInterfaceID );
}

csExternalHeader::~csExternalHeader() {
  if( myNavInterface ) {
    delete myNavInterface; myNavInterface = NULL;
  }
  if( myNavHeader ) {
    delete myNavHeader; myNavHeader = NULL;
  }
  if( myGCS90Header ) {
    delete myGCS90Header; myGCS90Header = NULL;
  }
}
void csExternalHeader::extractHeaders( byte const* buffer ) {
  
}
void csExternalHeader::extractHeaders( byte const* buffer, int totalNumBytes ) {
  // int totalNumBytes = 0;
  int bytePos = 0;

  myNavHeader->extractHeaders( &buffer[bytePos], totalNumBytes-bytePos );
  bytePos += myNavHeader->size();
  if( bytePos >= totalNumBytes ) return;

  myGCS90Header->extractHeaders( &buffer[bytePos], totalNumBytes-bytePos );
  bytePos += myGCS90Header->totalSize();
  if( bytePos >= totalNumBytes ) return;

  myNavInterface->extractHeaders( &buffer[bytePos], totalNumBytes-bytePos );
}
void csExternalHeader::setCommonRecordHeaders( commonRecordHeaderStruct& comRecHdr ) {
  if( myNavHeader->navHeaderID() != NAV_HEADER_NONE ) {
    myNavHeader->setCommonRecordHeaders( comRecHdr );
  }
  else {
    if( myGCS90Header->isOK() ) {
      comRecHdr.shotNum    = (int)myGCS90Header->shotPoint;
      //      comRecHdr.seq        = (int)myGCS90Header.currentSeqNumber;
    }
  }
}
void csExternalHeader::dump( std::ostream& os )
{
  os <<
    header("ExternalHeader start") << std::endl;;
  myNavHeader->dump( os );
  myGCS90Header->dump( os );
  myNavInterface->dump( os );
  //"...external header format not supported..." << std::endl <<
  os << header("ExternalHeader end  ") << std::endl;
}


