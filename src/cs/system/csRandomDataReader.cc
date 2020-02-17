/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csRandomDataReader.h"
#include "csSeismicReader.h"
#include "csSegyReader.h"
#include "csSegyHeaderInfo.h"
#include "csSegyTraceHeader.h"
#include "csRSFReader.h"
#include "csRSFHeader.h"
#include "csIReader.h"
#include "csException.h"
#include "csSuperHeader.h"
#include "csTraceHeaderDef.h"
#include "csTraceHeader.h"
#include "csStandardHeaders.h"

using namespace cseis_geolib;
using namespace cseis_system;

csRandomDataReader::csRandomDataReader( std::string filename, bool enableRandomAccess, int numTracesToBuffer ) {
  mySeismicReader = NULL;
  mySegyReader    = NULL;
  myRSFReader     = NULL;

  myHdrValueBlock = NULL;
  myHdrValueBlockSize = 0;
  myHdrIndex = NULL;
  myHdrType = NULL;
  
  // Determine file format:
  int length = (int)filename.length();
  //  fprintf(stderr,"FILE EXTENSION: %s\n",filename.substr(length-6).c_str());
  if( length >= 6 && filename.substr(length-6).compare(".cseis") == 0 ) {
    myFileFormat = csRandomDataReader::FORMAT_CSEIS;
    mySeismicReader = new cseis_system::csSeismicReader( filename, enableRandomAccess, numTracesToBuffer );
    myIReader = mySeismicReader;
  }
  else if( length >= 4 && filename.substr(length-4).compare(".rsf") == 0 ) {
    myFileFormat = csRandomDataReader::FORMAT_RSF;
    bool reverseByteOrder = false;
    myRSFReader = new cseis_io::csRSFReader( filename, numTracesToBuffer, reverseByteOrder );
    myIReader = myRSFReader;
  }
  else {
    myFileFormat = csRandomDataReader::FORMAT_SEGY;
    cseis_geolib::csSegyReader::SegyReaderConfig config;
    config.numTracesBuffer   = numTracesToBuffer;
    config.segyHeaderMapping = csSegyHdrMap::SEGY_STANDARD;
    config.reverseByteOrderData = false;
    config.reverseByteOrderHdr  = false;
    config.enableRandomAccess = enableRandomAccess;
    config.isSUFormat = false;
    if( length >= 3 && filename.substr(length-6).compare(".su") == 0 ) {
      config.isSUFormat = true;
    }
    mySegyReader = new cseis_geolib::csSegyReader( filename, config );
    myIReader = mySegyReader;
  }
}
csRandomDataReader::~csRandomDataReader() { 
  if( mySegyReader != NULL ) {
    delete mySegyReader;
    mySegyReader = NULL;
  }
  if( myRSFReader != NULL ) {
    delete myRSFReader;
    myRSFReader = NULL;
  }
  if( mySeismicReader != NULL ) {
    delete mySeismicReader;
    mySeismicReader = NULL;
  }
  if( myHdrValueBlock != NULL ) {
    delete [] myHdrValueBlock;
    myHdrValueBlock = NULL;
  }
  if( myHdrIndex != NULL ) {
    delete [] myHdrIndex;
    myHdrIndex = NULL;
  }
  if( myHdrType != NULL ) {
    delete [] myHdrType;
    myHdrType = NULL;
  }
}

bool csRandomDataReader::readFileHeader( cseis_system::csSuperHeader* shdr, cseis_system::csTraceHeaderDef* hdef, std::FILE* errorStream ) {
  bool retValue = true;
  if( myFileFormat == csRandomDataReader::FORMAT_CSEIS ) {
    retValue = mySeismicReader->readFileHeader( shdr, hdef, &myHdrValueBlockSize, errorStream );
    myHdrValueBlock = new char[myHdrValueBlockSize];
    //    fprintf(stderr,"value block size: %d\n", myHdrValueBlockSize);
  }
  else if( myFileFormat == csRandomDataReader::FORMAT_SEGY ) {
    try {
      mySegyReader->initialize();
    }
    catch( cseis_geolib::csException& e ) {
      fprintf(errorStream, "Error when initializing SEGY reader object.\nSystem message: %s\n", e.getMessage() );
      return false;
    }
    shdr->numSamples = mySegyReader->numSamples();
    shdr->sampleInt  = mySegyReader->sampleIntMS();
    setSEGYTraceHeaders( hdef );
  }
  else if( myFileFormat == csRandomDataReader::FORMAT_RSF ) {
    cseis_io::csRSFHeader rsfHdr;
    try {
      myRSFReader->initialize( &rsfHdr );
    }
    catch( cseis_geolib::csException& e ) {
      fprintf(errorStream,"Error when initializing RSF reader object.\nSystem message: %s\n", e.getMessage() );
      return false;
    }
    shdr->numSamples = myRSFReader->numSamples();
    shdr->sampleInt  = myRSFReader->sampleInt();
    setRSFTraceHeaders( hdef );
  }
  return retValue;
}

int csRandomDataReader::numTraces() const {
  return myIReader->numTraces();
}
int csRandomDataReader::getCurrentTraceIndex() const {
  return myIReader->getCurrentTraceIndex();
}
bool csRandomDataReader::moveToTrace( int traceIndex ) {
  return myIReader->moveToTrace( traceIndex );
}
bool csRandomDataReader::setHeaderToPeek( std::string const& headerName ) {
  return myIReader->setHeaderToPeek( headerName );
}
bool csRandomDataReader::setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType ) {
  return myIReader->setHeaderToPeek( headerName, headerType );
}
bool csRandomDataReader::peekHeaderValue( cseis_geolib::csFlexHeader* value, int traceIndex ) {
  return myIReader->peekHeaderValue( value, traceIndex );
}
bool csRandomDataReader::readTrace( float* samples, int numSamples, cseis_system::csTraceHeader* trcHdr ) {
  bool retValue = false;
  if( myFileFormat == csRandomDataReader::FORMAT_CSEIS ) {
    //    fprintf(stderr,"Read %d samples...\n", numSamples);
    retValue = mySeismicReader->readTrace( samples, myHdrValueBlock, numSamples );
    //    fprintf(stderr,"Read %d samples... DONE\n", numSamples);
    if( trcHdr != NULL ) trcHdr->setTraceHeaderValueBlock( myHdrValueBlock, myHdrValueBlockSize );
  }
  else if( myFileFormat == csRandomDataReader::FORMAT_SEGY ) {
    retValue = mySegyReader->getNextTrace( (cseis_geolib::byte_t*)samples, numSamples );
    if( trcHdr != NULL ) readSegyTraceHeaders( trcHdr );
  }
  else if( myFileFormat == csRandomDataReader::FORMAT_RSF ) {
    retValue = myRSFReader->getNextTrace( (cseis_geolib::byte_t*)samples, numSamples );
    if( trcHdr != NULL ) readRSFTraceHeaders( trcHdr );
  }
  return retValue;
}
void csRandomDataReader::setTraceHeaderValues( cseis_system::csTraceHeader* trcHdr ) {
  if( myFileFormat == csRandomDataReader::FORMAT_CSEIS ) {
    trcHdr->setTraceHeaderValueBlock( myHdrValueBlock, myHdrValueBlockSize );
  }
  else if( myFileFormat == csRandomDataReader::FORMAT_SEGY ) {
    readSegyTraceHeaders( trcHdr );
  }
  else if( myFileFormat == csRandomDataReader::FORMAT_RSF ) {
    void readRSFTraceHeaders( cseis_system::csTraceHeader* trcHdr );
  }

}

int csRandomDataReader::fileFormat() const {
  return myFileFormat;
}

void csRandomDataReader::setSEGYTraceHeaders( cseis_system::csTraceHeaderDef* hdef ) {
  cseis_geolib::csSegyHdrMap hdrMap( cseis_geolib::csSegyHdrMap::SEGY_STANDARD );
  int nHeaders = mySegyReader->numTraceHeaders();
  myHdrIndex = new int[nHeaders];
  myHdrType  = new cseis_geolib::type_t[nHeaders];
  for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
    cseis_geolib::csSegyHeaderInfo const* info = mySegyReader->header(ihdr);
    myHdrType[ihdr]  = info->outType;
    if( myHdrType[ihdr] != TYPE_STRING ) {
      myHdrIndex[ihdr] = hdef->addHeader( myHdrType[ihdr], info->name, info->description );
    }
    else { // TYPE_STRING
      myHdrIndex[ihdr] = hdef->addHeader( myHdrType[ihdr], info->name, info->description, info->byteSize );
    }
  }
}

void csRandomDataReader::setRSFTraceHeaders( cseis_system::csTraceHeaderDef* hdef ) {
  int nHeaders = 3;
  myHdrIndex = new int[nHeaders];
  myHdrType  = new cseis_geolib::type_t[nHeaders];
  myHdrType[0] = cseis_geolib::TYPE_DOUBLE;
  myHdrType[1] = cseis_geolib::TYPE_DOUBLE;
  myHdrType[2] = cseis_geolib::TYPE_DOUBLE;
  
  myHdrIndex[0] = hdef->addHeader( myHdrType[0], "dim2", "RSF data dimension 2" );
  myHdrIndex[1] = hdef->addHeader( myHdrType[1], "dim3", "RSF data dimension 3" );
  myHdrIndex[2] = hdef->addStandardHeader( cseis_geolib::HDR_TRCNO.name );
}

void csRandomDataReader::readRSFTraceHeaders( cseis_system::csTraceHeader* trcHdr ) {
  int traceIndex = myRSFReader->getCurrentTraceIndex();
  trcHdr->setDoubleValue( myHdrIndex[0], myRSFReader->computeDim2(traceIndex) );
  trcHdr->setDoubleValue( myHdrIndex[1], myRSFReader->computeDim3(traceIndex) );
  trcHdr->setIntValue( myHdrIndex[2], traceIndex+1 );
}

void csRandomDataReader::readSegyTraceHeaders( cseis_system::csTraceHeader* trcHdr ) {
  cseis_geolib::csSegyTraceHeader const* segyTrcHdr = mySegyReader->getTraceHeader();
 
  int nHeaders = segyTrcHdr->numHeaders();
  for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
    int hdrIdOut = myHdrIndex[ihdr];
    switch( myHdrType[ihdr] ) {
    case TYPE_FLOAT:
      trcHdr->setFloatValue( hdrIdOut, segyTrcHdr->floatValue(ihdr) );
      break;
    case TYPE_DOUBLE:
      trcHdr->setDoubleValue( hdrIdOut, segyTrcHdr->doubleValue(ihdr) );
      break;
    case TYPE_INT:
      trcHdr->setIntValue( hdrIdOut, segyTrcHdr->intValue(ihdr) );
      break;
    case TYPE_INT64:
      trcHdr->setInt64Value( hdrIdOut, (csInt64_t)segyTrcHdr->intValue(ihdr) );
      break;
    case TYPE_STRING:
      trcHdr->setStringValue( hdrIdOut, segyTrcHdr->stringValue(ihdr) );
      break;
    }
  }
}
