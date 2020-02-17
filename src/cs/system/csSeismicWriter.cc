/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSeismicWriter.h"
#include "csSeismicWriter_ver.h"
#include "csSeismicIOConfig.h"
#include "csSuperHeader.h"
#include "csTraceHeaderDef.h"
#include "csHeaderInfo.h"
#include "csTraceHeaderInfo.h"
#include "csGeolibUtils.h"
#include "csException.h"
#include <cstring>

using namespace cseis_system;

csSeismicWriter::csSeismicWriter( std::string filename, int numTracesBuffer, int sampleByteSize, bool overwrite ) {
  myWriter = new cseis_io::csSeismicWriter_ver( filename, numTracesBuffer, sampleByteSize, overwrite );
  myHdrTempBuffer    = NULL;
  myHdef = NULL;
}
csSeismicWriter::~csSeismicWriter() {
  if( myWriter != NULL ) {
    delete myWriter;
    myWriter = NULL;
  }
  if( myHdrTempBuffer != NULL ) {
    delete [] myHdrTempBuffer;
    myHdrTempBuffer = NULL;
  }
}
//--------------------------------------------------------------------
bool csSeismicWriter::writeFileHeader( csSuperHeader const* shdr, csTraceHeaderDef const* hdef ) {
  myHdef = hdef;

  cseis_io::csSeismicIOConfig config;
  config.domain     = shdr->domain;
  config.fftDataType  = shdr->fftDataType;
  config.numSamplesXT = shdr->numSamplesXT;
  config.sampleIntXT  = shdr->sampleIntXT;
  config.numSamples = shdr->numSamples;
  config.sampleInt  = shdr->sampleInt;
  config.grid_orig_x     = shdr->grid_orig_x;
  config.grid_orig_y     = shdr->grid_orig_y;
  config.grid_orig_il    = shdr->grid_orig_il;
  config.grid_orig_xl    = shdr->grid_orig_xl;
  config.grid_binsize_il = shdr->grid_binsize_il;
  config.grid_binsize_xl = shdr->grid_binsize_xl;
  config.grid_azim_il    = shdr->grid_azim_il;
  config.grid_azim_xl    = shdr->grid_azim_xl;
  int numEnsKeys = shdr->numEnsembleKeys();
  for( int ikey = 0; ikey < numEnsKeys; ikey++ ) {
    std::string name = *shdr->ensembleKey(ikey);
    config.ensKeyNames.insertEnd( name );
  }

  config.byteSizeSamples       = config.numSamples * 4;  // Get sample size dynamically
  config.byteSizeHdrValueBlock = hdef->getTotalNumBytes();

  int numTrcHdr = hdef->numHeaders();
  for( int ihdr = 0; ihdr < numTrcHdr; ihdr++ ) {
    cseis_geolib::csHeaderInfo info;
    info.name        = hdef->headerName(ihdr);
    info.description = hdef->headerDesc(ihdr);
    info.type        = hdef->headerType(ihdr);
    info.nElements   = hdef->numElements(ihdr);
    config.addHeader( &info );
  }

  // Check whether any trace headers should be deleted
  int numHdrBytes = 0;
  for( int ihdr = 0; ihdr < numTrcHdr; ihdr++ ) {
    cseis_geolib::csHeaderInfo const* info = config.headerInfo( ihdr );
    if( info->type != cseis_geolib::TYPE_STRING ) {
      numHdrBytes += cseis_geolib::csGeolibUtils::numBytes( info->type );
    }
    else {
      numHdrBytes += info->nElements;
    }
  }
  if( numHdrBytes != config.byteSizeHdrValueBlock ) {
    config.byteSizeHdrValueBlock = numHdrBytes;
    myHdrTempBuffer = new char[config.byteSizeHdrValueBlock];
  }

  return myWriter->writeFileHeader( &config );
}
bool csSeismicWriter::writeTrace( float* samples, char const* hdrValueBlock ) {
  if( myHdrTempBuffer == NULL ) {
    return myWriter->writeTrace( samples, hdrValueBlock );
  }
  else {
    int numTrcHdr = myHdef->numHeaders();
    int counterBytes = 0;
    int const* byteLocation = myHdef->getHandleByteLocation();
    for( int ihdr = 0; ihdr < numTrcHdr; ihdr++ ) {
      cseis_system::csTraceHeaderInfo const* info = myHdef->headerInfo( ihdr );
      int numHdrBytes = 0;
      if( info->type != cseis_geolib::TYPE_STRING ) {
        numHdrBytes = cseis_geolib::csGeolibUtils::numBytes( info->type );
      }
      else {
        numHdrBytes = info->nElements;
      }
      memcpy( &myHdrTempBuffer[counterBytes], &hdrValueBlock[byteLocation[ihdr]], numHdrBytes );
      counterBytes += numHdrBytes;
    }
    return myWriter->writeTrace( samples, myHdrTempBuffer );
  }
}
