/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSeismicReader_ver02.h"
#include "csSeismicIOConfig.h"
#include "csException.h"
#include "csHeaderInfo.h"
#include "csIODefines.h"
#include "csFileUtils.h"
#include "csFlexHeader.h"
#include <cstring>
#include <limits>

extern "C" {
  #include <sys/stat.h>
  #include <unistd.h>
  #include <stdio.h>
}

using namespace cseis_io;

csSeismicReader_ver02::csSeismicReader_ver02( std::string filename, bool enableRandomAccess, int numTracesBuffer ) :
  csSeismicReader_ver( filename, enableRandomAccess, numTracesBuffer ) {
  cseis_io::csIODefines::createVersionString( VERSION_SEISMIC_READER, myVersionText );

  int versionNumber = VERSION_SEISMIC_READER % 100;
  myVersionMajor = (short int)(versionNumber/10);
  myVersionMinor = (short int)(versionNumber - myVersionMajor*10);

  initialize();
}
//----------------------------------------------------------------
csSeismicReader_ver02::~csSeismicReader_ver02() {
}
//--------------------------------------------------------------------
bool csSeismicReader_ver02::readFileHeader( csSeismicIOConfig* config ) {
  if( myFile == NULL ) return false;
  if( myIsReadFileHeader ) throw( cseis_geolib::csException("csSeismicReader_ver00::readFileHeader: Attempt to re-read file header. This is probably a program bug in the calling function") );
  myIsReadFileHeader = true;
  
  char byteSizeChar[4];
  myFile->read( byteSizeChar, 4 );
  int byteSize = 0;
  memcpy( &byteSize, byteSizeChar, 4 );
  if( myFile->fail() ) {
    throw cseis_geolib::csException("Unexpected error occurred when reading SeaSeis header");
  }
  myHeaderByteSize += 4;

  myTempBuffer = new char[byteSize];
  myFile->read( myTempBuffer, byteSize );
  if( myFile->fail() ) {
    throw cseis_geolib::csException("Unexpected error occurred when reading SeaSeis header");
  }
  myHeaderByteSize += byteSize;

  //  fprintf(stdout,"Byte size read: %d\n", byteSize);
  // Set super header
  int byteLoc = 0;
  memcpy( &config->numSamples, &myTempBuffer[byteLoc], 4 );
  memcpy( &config->sampleInt, &myTempBuffer[byteLoc+4], 4 );
  memcpy( &config->domain, &myTempBuffer[byteLoc+8], 4 );
  memcpy( &config->grid_orig_x, &myTempBuffer[byteLoc+12], 8 );
  memcpy( &config->grid_orig_y, &myTempBuffer[byteLoc+20], 8 );
  memcpy( &config->grid_orig_il, &myTempBuffer[byteLoc+28], 4 );
  memcpy( &config->grid_orig_xl, &myTempBuffer[byteLoc+32], 4 );
  memcpy( &config->grid_binsize_il, &myTempBuffer[byteLoc+36], 8 );
  memcpy( &config->grid_binsize_xl, &myTempBuffer[byteLoc+44], 8 );
  memcpy( &config->grid_azim_il, &myTempBuffer[byteLoc+52], 8 );
  memcpy( &config->grid_azim_xl, &myTempBuffer[byteLoc+60], 8 );
  myNumSamples = config->numSamples;
  int numEnsKeys = 0;
  memcpy( &numEnsKeys, &myTempBuffer[byteLoc+68], 4 );
  byteLoc += 72;

  //  fprintf(stdout,"Num ens keys: %d\n", numEnsKeys);

  for( int ikey = 0; ikey < numEnsKeys; ikey++ ) {
    int sizeText;
    memcpy( &sizeText, &myTempBuffer[byteLoc], 4 );
    char* text = new char[sizeText+1];
    text[sizeText] = '\0';
    memcpy( text, &myTempBuffer[byteLoc+4], sizeText );
    config->ensKeyNames.insertEnd( text );
    delete [] text;
    byteLoc += 4 + sizeText;
  }


  memcpy( &config->byteSizeSamples, &myTempBuffer[byteLoc], 4 );
  memcpy( &config->byteSizeHdrValueBlock, &myTempBuffer[byteLoc+4], 4 );
  //  fprintf(stdout,"Num bytes: %d %d\n", config->byteSizeSamples, config->byteSizeHdrValueBlock );

  int numTrcHdrs = 0;
  memcpy( &numTrcHdrs, &myTempBuffer[byteLoc+8], 4 );
  byteLoc += 12;
  //  fprintf(stdout,"Num trc hdrs: %d\n", numTrcHdrs);

  for( int ihdr = 0; ihdr < numTrcHdrs; ihdr++ ) {
    cseis_geolib::type_t type = (cseis_geolib::type_t)myTempBuffer[byteLoc];
    //fprintf(stdout,"Header %d: Type: %d, ", ihdr, type );
    byteLoc += 1;

    int nElements;
    memcpy( &nElements, &myTempBuffer[byteLoc], 4 );
  //  fprintf(stdout,"nElements: %d, ", nElements );
    byteLoc += 4;

    int sizeName;
    memcpy( &sizeName, &myTempBuffer[byteLoc], 4 );
    //    fprintf(stdout,"SizeName: %d, ", sizeName );
    byteLoc += 4;
    char* name = new char[sizeName+1];
    name[sizeName] = '\0';
    memcpy( name, &myTempBuffer[byteLoc], sizeName );
    byteLoc += sizeName;
    //fprintf(stdout,"Name: '%s', ", name);
    
    int sizeDesc;
    memcpy( &sizeDesc, &myTempBuffer[byteLoc], 4 );
    //fprintf(stdout,"SizeDesc: %d, ", sizeDesc );
    byteLoc += 4;
    char* desc = new char[sizeDesc+1];
    desc[sizeDesc] = '\0';
    memcpy( desc, &myTempBuffer[byteLoc], sizeDesc );
    byteLoc += sizeDesc;
    //fprintf(stdout,"Desc: '%s', ", desc);
    //fprintf(stdout,"'\n");

    config->addHeader( type, name, desc, nElements );

    delete [] name;
    delete [] desc;
  }

  myByteSizeSamples       = config->byteSizeSamples;
  myByteSizeHdrValueBlock = config->byteSizeHdrValueBlock;
  myTraceByteSize         = myByteSizeSamples + myByteSizeHdrValueBlock;

  if( myBufferCapacityNumTraces <= 0 ) {  // Set number of buffered traces if it wasn't set explicitely in constructor
    myBufferCapacityNumTraces = DEFAULT_BUFFERED_SAMPLES / config->numSamples;
    if( myBufferCapacityNumTraces <= 0 ) myBufferCapacityNumTraces = 1;
    else if( myBufferCapacityNumTraces > 20 ) myBufferCapacityNumTraces = 20;
  }

  // Determine number of traces in file
  if( myFileSize != cseis_geolib::csFileUtils::FILESIZE_UNKNOWN ) {
    csInt64_t sizeTraces  = myFileSize - (csInt64_t)myHeaderByteSize;
    myNumTraces       = (int)(sizeTraces/(csInt64_t)(myByteSizeHdrValueBlock+myByteSizeSamples));  // Total number of full traces in input file
    config->numTraces = myNumTraces;
    myLastTraceIndex  = myNumTraces-1;
    if( myBufferCapacityNumTraces > myNumTraces ) myBufferCapacityNumTraces = myNumTraces;
  }
  else {
    config->numTraces = 0;
    myLastTraceIndex = 0;
  }

  if( myBufferCapacityNumTraces > 1 ) {
    resizeDataBuffer( myBufferCapacityNumTraces );
  }

  //  fprintf(stderr,"IN:  Byte size: %d %d %d   %d, num trchdrs: %d\n", myByteSizeSamples, myByteSizeHdrValueBlock, myHeaderByteSize, byteSize, numTrcHdrs );

  return true;
}

