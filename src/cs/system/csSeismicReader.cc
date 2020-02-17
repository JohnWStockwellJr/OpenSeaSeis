/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSeismicWriter.h"
#include "csSeismicReader.h"
#include "csSeismicReader_ver.h"
#include "csSeismicIOConfig.h"
#include "csSuperHeader.h"
#include "csTraceHeaderDef.h"
#include "csHeaderInfo.h"
#include "csStandardHeaders.h"
#include "csGeolibUtils.h"
#include "csException.h"
#include "csTraceHeaderInfo.h"
#include "csFlexHeader.h"
#include "csIOSelection.h"
#include <string>

using namespace cseis_system;

csSeismicReader::csSeismicReader( std::string filename, int numTraces ) {
  bool enableRandomAccess = false;
  myReader = cseis_io::csSeismicReader_ver::createReaderObject( filename, enableRandomAccess, numTraces );
  init();
}

csSeismicReader::csSeismicReader( std::string filename, bool enableRandomAccess, int numTraces ) {
  myReader = cseis_io::csSeismicReader_ver::createReaderObject( filename, enableRandomAccess, numTraces );
  init();
}

void csSeismicReader::init() {
  myNumTraces = 0;
  myHdrCheckByteOffset = 0;
  myHdrCheckType       = cseis_geolib::TYPE_UNKNOWN;
  myHdrCheckByteSize   = 0;
  myHdrCheckBuffer     = NULL;
  myTrcHdrDef          = NULL;
  myIOSelection        = NULL;
}

csSeismicReader::~csSeismicReader() {
  if( myIOSelection != NULL ) {
    delete myIOSelection;
    myIOSelection = NULL;
  }
  if( myReader != NULL ) {
    delete myReader;
    myReader = NULL;
  }
  if( myHdrCheckBuffer != NULL ) {
    delete [] myHdrCheckBuffer;
    myHdrCheckBuffer = NULL;
  }
}
//--------------------------------------------------------------------
bool csSeismicReader::readFileHeader( csSuperHeader* shdr, csTraceHeaderDef* hdef, int* hdrValueBlockSize, std::FILE* stream ) {
  cseis_io::csSeismicIOConfig config;
  bool success = myReader->readFileHeader( &config );
  if( !success ) return false;
  
  shdr->domain      = config.domain;
  shdr->fftDataType  = config.fftDataType;
  shdr->numSamplesXT  = config.numSamplesXT;
  shdr->sampleIntXT  = config.sampleIntXT;
  shdr->numSamples = config.numSamples;
  shdr->sampleInt  = config.sampleInt;
  shdr->grid_orig_x      = config.grid_orig_x;
  shdr->grid_orig_y      = config.grid_orig_y;
  shdr->grid_orig_il     = config.grid_orig_il;
  shdr->grid_orig_xl     = config.grid_orig_xl;
  shdr->grid_binsize_il  = config.grid_binsize_il;
  shdr->grid_binsize_xl  = config.grid_binsize_xl;
  shdr->grid_azim_il     = config.grid_azim_il;
  shdr->grid_azim_xl     = config.grid_azim_xl;
  int numEnsKeys   = config.ensKeyNames.size();
  for( int ikey = 0; ikey < numEnsKeys; ikey++ ) {
    shdr->setEnsembleKey( config.ensKeyNames.at(ikey), ikey );
  }
//  config.byteSizeSamples = config.numSamples * 4;  // Get sample size dynamically
//  config.byteSizeHdrValueBlock = hdef->getTotalNumBytes();

  if( hdef != NULL ) {
    int numTrcHdr = config.numTrcHeaders();
    for( int ihdr = 0; ihdr < numTrcHdr; ihdr++ ) {
      cseis_geolib::csHeaderInfo const* info = config.headerInfo( ihdr );
      cseis_geolib::csHeaderInfo const* info_std = cseis_geolib::csStandardHeaders::get( info->name );
      if( info_std == NULL ) {
        hdef->addHeader( info->type, info->name, info->description, info->nElements );
      }
      else if( info_std->type == info->type ) {
        hdef->addStandardHeader( info->name );
      }
      else {
        std::string changedName = "_";
        changedName.append(info->name);
        hdef->addHeader( info->type, changedName, info->description, info->nElements );
        if( stream != NULL ) {
          fprintf(stream,"Input file trace header '%s' is a defined standard header, however of different type. Actual: %s, standard: %s.\n",
            info->name.c_str(), cseis_geolib::csGeolibUtils::typeText(info->type), cseis_geolib::csGeolibUtils::typeText(info_std->type) );
          fprintf(stream,"This may be due to a update of the standard header definitions since the data file was written.\n");
          fprintf(stream,"To avoid conflicts with the standard definition, the trace header will be renamed to '%s'.\n", changedName.c_str() );
        }
      }
    }
  }

  *hdrValueBlockSize = config.byteSizeHdrValueBlock;
  myNumTraces = config.numTraces;

  myTrcHdrDef = hdef;

  return true;
}

bool csSeismicReader::readTrace( float* samples, char* hdrValueBlock, int numSamples ) {
  if( myIOSelection ) {
    if( !performIOSelection() ) return false;
  }
  return myReader->readTrace( samples, hdrValueBlock, numSamples );
}
bool csSeismicReader::readTrace( float* samples, char* hdrValueBlock ) {
  if( myIOSelection ) {
    if( !performIOSelection() ) return false;
  }
  return myReader->readTrace( samples, hdrValueBlock );
}
bool csSeismicReader::performIOSelection() {
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
bool csSeismicReader::setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType ) {
  bool success = setHeaderToPeek( headerName );
  headerType = myHdrCheckType;
  return success;
}
bool csSeismicReader::setHeaderToPeek( std::string const& headerName ) {
  if( !myTrcHdrDef->headerExists( headerName ) ) return false;
  int index = myTrcHdrDef->headerIndex( headerName );
  csTraceHeaderInfo const* info = myTrcHdrDef->headerInfo( index );

  myHdrCheckByteOffset   = myTrcHdrDef->getByteLocation( index );
  myHdrCheckType         = info->type;
  int hdrCheckByteSize   = info->nElements * cseis_geolib::csGeolibUtils::numBytes(myHdrCheckType);

  if( hdrCheckByteSize != myHdrCheckByteSize && myHdrCheckBuffer != NULL ) {
    delete [] myHdrCheckBuffer;
    myHdrCheckBuffer = NULL;
  }
  if( myHdrCheckBuffer == NULL ) {
    myHdrCheckByteSize =  hdrCheckByteSize;
    myHdrCheckBuffer = new char[myHdrCheckByteSize];
  }
  return true;
}
int csSeismicReader::getPeekByteOffset() const {
  return myHdrCheckByteOffset;
}
bool csSeismicReader::peekHeaderValue( cseis_geolib::csFlexHeader* hdrValue, int traceIndex ) {
  if( myHdrCheckBuffer == NULL ) {
    throw( cseis_geolib::csException("csSeismicReader::checkHeaderValue: No header has been set for checking. This is a program bug in the calling function") );
  }
  bool success = true;
  if( traceIndex < 0 ) {
    success = myReader->peek( myHdrCheckByteOffset, myHdrCheckByteSize, myHdrCheckBuffer );
  }
  else {
    //    fprintf(stdout,"Peek: %d %d %d\n", myHdrCheckByteOffset, myHdrCheckByteSize, traceIndex );
    success = myReader->peek( myHdrCheckByteOffset, myHdrCheckByteSize, myHdrCheckBuffer, traceIndex );
  }
  if( !success ) return false;
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
  else if( myHdrCheckType == cseis_geolib::TYPE_STRING ) {
    std::string text = std::string( myHdrCheckBuffer );
    hdrValue->setStringValue( text );
  }
  return true;
}
int csSeismicReader::numTracesCapacity() const {
  if( myReader == NULL ) return 0;
  return myReader->numTracesCapacity();
}
bool csSeismicReader::moveToTrace( int traceIndex ) {
  if( myReader == NULL ) return false;
  return myReader->moveToTrace( traceIndex );
}

bool csSeismicReader::moveToTrace( int traceIndex, int numTracesToRead ) {
  if( myReader == NULL ) return false;
  return myReader->moveToTrace( traceIndex, numTracesToRead );
}
int csSeismicReader::numSamples() const {
  if( myReader == NULL ) return 0;
  return myReader->numSamples();
}

bool csSeismicReader::setSelection( std::string const& hdrValueSelectionText, std::string const& headerName,
                                    int sortOrder, int sortMethod )
{
  myIOSelection = new cseis_geolib::csIOSelection( headerName, sortOrder, sortMethod );
  return myIOSelection->initialize( this, hdrValueSelectionText );
}
int csSeismicReader::getCurrentTraceIndex() const {
  return myReader->currentTraceIndex();
}
