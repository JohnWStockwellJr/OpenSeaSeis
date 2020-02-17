/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSeismicWriter_ver.h"
#include "csSeismicIOConfig.h"
#include "csException.h"
#include "csGeolibUtils.h"
#include "csHeaderInfo.h"
#include "csIODefines.h"
#include <cstring>
#include <limits>
#include <cmath>

using namespace cseis_io;
using namespace std;

csSeismicWriter_ver::csSeismicWriter_ver( std::string filename, int numTracesBuffer, int sampleByteSize, bool overwrite ) {
  cseis_io::csIODefines::createVersionString( VERSION_SEISMIC_WRITER, myVersionText );

  myFileName = filename;
  myByteSizeOneSample = sampleByteSize;
  if( myByteSizeOneSample != 4 ) {
    myByteSizeCompression = 8;  // 8 additional bytes needed for data compression, in order to store minValue & rangeValue
  }
  else {
    myByteSizeCompression = 0;
  }
  myCompressedSampleBuffer = NULL;

  myTempByteSize = 0;
  myTempBuffer   = NULL;
  myByteLoc      = 0;

  myDataBufferSize   = 0;
  myDataBuffer       = NULL;
  myCurrentDataBufferSize = 0;
  myNumBufferTraces  = numTracesBuffer;

  open( overwrite );
}
//----------------------------------------------------------------
csSeismicWriter_ver::~csSeismicWriter_ver() {
  close();
  if( myTempBuffer != NULL ) {
    delete [] myTempBuffer;
    myTempBuffer = NULL;
  }
  if( myDataBuffer ) {
    delete [] myDataBuffer;
    myDataBuffer = NULL;
  }
  if( myCompressedSampleBuffer != NULL ) {
    delete [] myCompressedSampleBuffer;
    myCompressedSampleBuffer = NULL;
  }
}
//----------------------------------------------------------------
void csSeismicWriter_ver::open( bool overwrite ) {
  if( !overwrite ) {
    try {
      myFile = fopen( myFileName.c_str(), "rb" );
    }
    catch( ... ) {
      myFile = NULL;
      throw( cseis_geolib::csException("Cannot open file '%s'", myFileName.c_str() ) );
    }
    if( myFile != NULL ) {
      throw( cseis_geolib::csException("File exists and shall NOT be overwritten: '%s'", myFileName.c_str() ) );
    }
  }
  try {
    myFile = fopen( myFileName.c_str(), "wb" );
  }
  catch( ... ) {
    myFile = NULL;
    throw( cseis_geolib::csException("Cannot open file '%s'", myFileName.c_str() ) );
  }
  if( myFile == NULL ) {
    throw( cseis_geolib::csException("Error occurred when opening file '%s'", myFileName.c_str() ) );
  }
}
//----------------------------------------------------------------
void csSeismicWriter_ver::close() {
  if( myFile != NULL ) {
    if( myCurrentDataBufferSize != 0 ) {
      // Some traces are still buffered and haven't been flushed yet --> Write them out now
      writeCurrentDataBuffer();
    }
    fclose( myFile );
    myFile = NULL;
  }
}
bool csSeismicWriter_ver::writeCurrentDataBuffer() {
  int sizeWrite = (int)fwrite( myDataBuffer, myCurrentDataBufferSize, 1, myFile );
  bool retValue = (sizeWrite == 1);
  myCurrentDataBufferSize = 0;
  return retValue;
}
//----------------------------------------------------------------
void csSeismicWriter_ver::initialize() {
  std::string text;
  text.append( ID_TEXT_CSEIS );
  text.append( myVersionText );
  int size = (int)text.length();
  int sizeWrite = 0;
  if( (sizeWrite = (int)fwrite( text.c_str(), size, 1, myFile ) ) != 1 ) return;
}
//--------------------------------------------------------------------
bool csSeismicWriter_ver::writeFileHeader( csSeismicIOConfig const* config ) {
  if( myFile == NULL ) return false;
  initialize();

  myByteLoc = 0;
  myTempByteSize = 200;
  myTempBuffer = new char[myTempByteSize];

  myNumSamples            = config->numSamples;
  myByteSizeSamples       = config->numSamples * myByteSizeOneSample;
  myByteSizeHdrValueBlock = config->byteSizeHdrValueBlock;

  if( myByteSizeOneSample != 4 ) {
    myCompressedSampleBuffer = new char[myByteSizeSamples];
  }

  int numTrcHdrs = config->numTrcHeaders();

  if( myNumBufferTraces <= 0 ) {  // Set number of buffered traces if it wasn't set explicitely in constructor
    myNumBufferTraces = DEFAULT_BUFFERED_BYTES / myByteSizeSamples;
    if( myNumBufferTraces <= 0 ) myNumBufferTraces = 1;
    else if( myNumBufferTraces > 20 ) myNumBufferTraces = 20;
  }

  resizeDataBuffer( myNumBufferTraces * (myByteSizeSamples + myByteSizeHdrValueBlock + myByteSizeCompression) );

// Set super header
  appendInt(    config->numSamples );
  appendFloat(  config->sampleInt );
  appendInt(    config->domain );
  appendInt(    config->fftDataType );
  appendInt(    config->numSamplesXT );
  appendFloat(    config->sampleIntXT );
  appendDouble( config->grid_orig_x );
  appendDouble( config->grid_orig_y );
  appendInt( config->grid_orig_il );
  appendInt( config->grid_orig_xl );
  appendDouble( config->grid_binsize_il );
  appendDouble( config->grid_binsize_xl );
  appendDouble( config->grid_azim_il );
  appendDouble( config->grid_azim_xl );
  int numEnsKeys = config->ensKeyNames.size();
  appendInt(   numEnsKeys );
  for( int ikey = 0; ikey < numEnsKeys; ikey++ ) {
    std::string name = config->ensKeyNames.at(ikey);
    int sizeName = (int)name.length();
    appendInt( sizeName );
    appendString( name.c_str(), sizeName );
  }

  appendInt( myByteSizeOneSample );
// Set trace header
  appendInt( config->byteSizeHdrValueBlock );

  appendInt( numTrcHdrs );
  char randomName[4];
  int randomCounter = 0;
  for( int ihdr = 0; ihdr < numTrcHdrs; ihdr++ ) {
    cseis_geolib::csHeaderInfo const* info;
    info = config->headerInfo( ihdr );
    appendChar( info->type );
    appendInt( info->nElements );
    int sizeName = (int)info->name.length();
    int sizeDesc = (int)info->description.length();
    if( sizeName > 0 ) {
      appendInt( sizeName );
      appendString( info->name.c_str(), sizeName );
    }
    else {
      // ??? Zero length header name. This must be a bug in CSEIS. Give this header a random name...
      if( randomCounter >= 99 ) randomCounter = 0;
      sprintf(randomName,"A%-2d", randomCounter++);
      appendInt( 3 );
      appendString( randomName, 3 );
    }
    appendInt( sizeDesc );
    if( sizeDesc > 0 ) {
      //      fprintf(stdout,"Description length: %d, '%s'\n", sizeDesc, info->description.c_str());
      appendString( info->description.c_str(), sizeDesc );
    }
  }

  int sizeWrite = 0;
  if( (sizeWrite = (int)fwrite( &myByteLoc, 4, 1, myFile ) ) != 1 ) {
  }
  if( (sizeWrite = (int)fwrite( myTempBuffer, myByteLoc, 1, myFile ) ) != 1 ) {
  }

  // fprintf(stderr,"OUT: Byte size: %d %d %d, numTrcHdrs: %d, numbytes orig:\n", myByteSizeSamples, myByteSizeHdrValueBlock, myByteLoc, numTrcHdrs );

  if( myTempBuffer != NULL ) {
    delete [] myTempBuffer;
    myTempBuffer = NULL;
  }
  return true;
}

//----------------------------------------------------------------
void csSeismicWriter_ver::resizeDataBuffer( int newSize ) {
  if( newSize > myDataBufferSize ) {
    char* newBuffer = NULL;
    try {
      newBuffer = new char[newSize];
    }
    catch(...) {
      throw( cseis_geolib::csException("csSeismicWriter_ver::resizeDataBuffer: Cannot allocate data buffer. Try to reduce number of buffered output traces.") );
    }
    if( myDataBuffer ) {
      delete [] myDataBuffer;
      myDataBuffer = NULL;
    }
    myDataBuffer     = newBuffer;
    myDataBufferSize = newSize;
  }
}
void csSeismicWriter_ver::checkBufferSize( int sizeAdd ) {
  if( myByteLoc+sizeAdd >= myTempByteSize ) {
    int newSize = std::max( myByteLoc + sizeAdd, 2*myTempByteSize );
    char* newBuffer = new char[newSize];
    memcpy( newBuffer, myTempBuffer, myByteLoc );
    delete [] myTempBuffer;
    myTempBuffer = newBuffer;
    myTempByteSize = newSize;
  }
}
void csSeismicWriter_ver::appendChar( char value ) {
  checkBufferSize( 1 );
  myTempBuffer[myByteLoc] = value;
  myByteLoc += 1;
}
void csSeismicWriter_ver::appendInt( int value ) {
  checkBufferSize( 4 );
  memcpy( &myTempBuffer[myByteLoc], &value, 4 );
  myByteLoc += 4;
}
void csSeismicWriter_ver::appendFloat( float value ) {
  checkBufferSize( 4 );
  memcpy( &myTempBuffer[myByteLoc], &value, 4 );
  myByteLoc += 4;
}
void csSeismicWriter_ver::appendDouble( double value ) {
  checkBufferSize( 8 );
  memcpy( &myTempBuffer[myByteLoc], &value, 8 );
  myByteLoc += 8;
}
void csSeismicWriter_ver::appendString( char const* value, int size ) {
  checkBufferSize( size );
  memcpy( &myTempBuffer[myByteLoc], value, size );
  myByteLoc += size;
}
//----------------------------------------------------------------
bool csSeismicWriter_ver::writeTrace( float* samples, char const* hdrValueBlock ) {
  if( myFile != NULL ) {
    memcpy( &(myDataBuffer[myCurrentDataBufferSize]), hdrValueBlock, myByteSizeHdrValueBlock );
    myCurrentDataBufferSize += myByteSizeHdrValueBlock;

    if( myByteSizeOneSample == 4 ) {
      memcpy( &(myDataBuffer[myCurrentDataBufferSize]), samples, myByteSizeSamples );
    }
    else {
      float minValue;
      float rangeValue;
      compressData( samples, myCompressedSampleBuffer, minValue, rangeValue );
      memcpy( &(myDataBuffer[myCurrentDataBufferSize]), &minValue, sizeof(float) );
      memcpy( &(myDataBuffer[myCurrentDataBufferSize+sizeof(float)]), &rangeValue, sizeof(float) );
      memcpy( &(myDataBuffer[myCurrentDataBufferSize+myByteSizeCompression]), myCompressedSampleBuffer, myByteSizeSamples );
    }
    myCurrentDataBufferSize += myByteSizeSamples + myByteSizeCompression;

    if( myCurrentDataBufferSize == myDataBufferSize ) {
      return writeCurrentDataBuffer();
    }
    return true;
  }
  else {
    return false;
  }
}

void csSeismicWriter_ver::computeCompressionValues( float const* samples, float& minValue, float& rangeValue ) {
  minValue = samples[0];
  float maxVal = samples[0];
  for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
    float value = samples[isamp];
    if( value < minValue ) {
      minValue = value;
    }
    if( value > maxVal ) {
      maxVal = value;
    }
  }
  rangeValue = maxVal - minValue;
}

void csSeismicWriter_ver::compressData( float const* samplesIn, char* samplesOut, float& minValue, float& rangeValue ) {
  computeCompressionValues( samplesIn, minValue, rangeValue );
  //  fprintf(stderr,"compressData, minValue = %.6e, rangeValue = %.6e\n", minValue, rangeValue);

  if( myByteSizeOneSample == 2 ) { // 16 bit
    int max16bitValue = (int)numeric_limits<unsigned short>::max();
    float stepValue   = rangeValue / (float)max16bitValue;
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      unsigned short newValue = (unsigned short)round( ( samplesIn[isamp] - minValue ) / stepValue );
      memcpy( &samplesOut[isamp*myByteSizeOneSample], &newValue, myByteSizeOneSample );
    }
  }
  else { // 8 bit
    int max8bitValue = (int)numeric_limits<unsigned char>::max();
    float stepValue = rangeValue / (float)max8bitValue;
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      unsigned char newValue = (unsigned char)round( ( samplesIn[isamp] - minValue ) / stepValue );
      memcpy( &samplesOut[isamp*myByteSizeOneSample], &newValue, myByteSizeOneSample );
    }
  }

}

