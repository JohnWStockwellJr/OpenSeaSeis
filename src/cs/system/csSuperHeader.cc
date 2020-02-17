/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_defines.h"
#include "csSuperHeader.h"
#include "csVector.h"
#include "csException.h"
#include "csGeolibUtils.h"
#include <cstring>

using namespace cseis_system;

csSuperHeader::csSuperHeader() {
  numSamples   = 0;
  numSamplesXT = 0;
  sampleInt    = 1.0;
  sampleIntXT  = 1.0;
  domain       = cseis_geolib::DOMAIN_XT;
  fftDataType  = cseis_geolib::FX_NONE;

  grid_orig_x     = 0.0;
  grid_orig_y     = 0.0;
  grid_orig_il    = 1;
  grid_orig_xl    = 1;
  grid_binsize_il = 0.0;
  grid_binsize_xl = 0.0;
  grid_azim_il    = 0.0;
  grid_azim_xl    = 90.0;

  myEnsKeyInfo = NULL;
  myNumEnsKeys = 0;
  myNumAllocatedEnsKeys = 0;
}

csSuperHeader::~csSuperHeader() {
  if( myEnsKeyInfo != NULL ) {
    delete [] myEnsKeyInfo;
    myEnsKeyInfo = NULL;
  }
}
int csSuperHeader::mpi_getByteSize() const {
  char* dataBuffer = new char[10000];
  int byteSize = mpi_compress( dataBuffer );
  delete [] dataBuffer;
  return byteSize;
}
int csSuperHeader::mpi_compress( char* data ) const {
  int byteLoc = 0;
  memcpy( &data[byteLoc], &numSamples, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], &sampleInt, sizeof(float) );
  byteLoc += sizeof(float);
  memcpy( &data[byteLoc], &domain, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], &numSamplesXT, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], &sampleIntXT, sizeof(float) );
  byteLoc += sizeof(float);
  memcpy( &data[byteLoc], &fftDataType, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], &grid_orig_x, sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &data[byteLoc], &grid_orig_y, sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &data[byteLoc], &grid_orig_il, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], &grid_orig_xl, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], &grid_binsize_il, sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &data[byteLoc], &grid_binsize_xl, sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &data[byteLoc], &grid_azim_il, sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &data[byteLoc], &grid_azim_xl, sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &data[byteLoc], &myNumEnsKeys, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], &myNumAllocatedEnsKeys, sizeof(int) );
  byteLoc += sizeof(int);

  for( int i = 0; i < myNumEnsKeys; i++ ) {
    keyInfoStruct info = myEnsKeyInfo[i];
    int strLength = info.name.length();
    memcpy( &data[byteLoc], &strLength, sizeof(int) );
    byteLoc += sizeof(int);
    memcpy( &data[byteLoc], info.name.c_str(), strLength*sizeof(char) );
    byteLoc += strLength*sizeof(char);
    memcpy( &data[byteLoc], &info.hdrType, sizeof(int) );
    byteLoc += sizeof(int);
    memcpy( &data[byteLoc], &info.hdrIndex, sizeof(int) );
    byteLoc += sizeof(int);
  }
  return byteLoc;
}

void csSuperHeader::mpi_decompress( char const* data ) {
  int byteLoc = 0;
  memcpy( &numSamples, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &sampleInt, &data[byteLoc], sizeof(float) );
  byteLoc += sizeof(float);
  memcpy( &domain, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &numSamplesXT, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &sampleIntXT, &data[byteLoc], sizeof(float) );
  byteLoc += sizeof(float);
  memcpy( &fftDataType, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &grid_orig_x, &data[byteLoc], sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &grid_orig_y, &data[byteLoc], sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &grid_orig_il, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &grid_orig_xl, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &grid_binsize_il, &data[byteLoc], sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &grid_binsize_xl, &data[byteLoc], sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &grid_azim_il, &data[byteLoc], sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &grid_azim_xl, &data[byteLoc], sizeof(double) );
  byteLoc += sizeof(double);
  memcpy( &myNumEnsKeys, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &myNumAllocatedEnsKeys, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);

  if( myEnsKeyInfo != NULL ) delete [] myEnsKeyInfo;
  myEnsKeyInfo = NULL;
  if( myNumAllocatedEnsKeys > 0 ) myEnsKeyInfo = new keyInfoStruct[myNumAllocatedEnsKeys];
  
  for( int i = 0; i < myNumEnsKeys; i++ ) {
    keyInfoStruct info;
    int strLength;
    memcpy( &strLength, &data[byteLoc], sizeof(int) );
    byteLoc += sizeof(int);
    char* nameStr = new char[strLength];
    memcpy( nameStr, &data[byteLoc], strLength*sizeof(char) );
    byteLoc += strLength*sizeof(char);
    myEnsKeyInfo[i].name = nameStr;
    delete [] nameStr;
    memcpy( &myEnsKeyInfo[i].hdrType, &data[byteLoc], sizeof(int) );
    byteLoc += sizeof(int);
    memcpy( &myEnsKeyInfo[i].hdrIndex, &data[byteLoc], sizeof(int) );
    byteLoc += sizeof(int);
  }
}

void csSuperHeader::set( csSuperHeader const* superHeader ) {
  numSamples   = superHeader->numSamples;
  numSamplesXT = superHeader->numSamplesXT;
  sampleInt    = superHeader->sampleInt;
  sampleIntXT  = superHeader->sampleIntXT;
  domain       = superHeader->domain;
  fftDataType  = superHeader->fftDataType;

  grid_orig_x     = superHeader->grid_orig_x;
  grid_orig_y     = superHeader->grid_orig_y;
  grid_orig_il    = superHeader->grid_orig_il;
  grid_orig_xl    = superHeader->grid_orig_xl;
  grid_binsize_il = superHeader->grid_binsize_il;
  grid_binsize_xl = superHeader->grid_binsize_xl;
  grid_azim_il    = superHeader->grid_azim_il;
  grid_azim_xl    = superHeader->grid_azim_xl;

  clearEnsembleKeys();
  for( int i = 0; i < superHeader->numEnsembleKeys(); i++ ) {
    setEnsembleKey( superHeader->myEnsKeyInfo[i].name, superHeader->myEnsKeyInfo[i].hdrIndex,
      superHeader->myEnsKeyInfo[i].hdrType, i );
  }
}

int csSuperHeader::numEnsembleKeys() const {
  return myNumEnsKeys;
}

void csSuperHeader::clearEnsembleKeys() {
  myNumEnsKeys = 0;
}
void csSuperHeader::reallocate( int newNumAllocatedKeys ) {
  if( myNumAllocatedEnsKeys < newNumAllocatedKeys ) {
    keyInfoStruct* keyInfo = new keyInfoStruct[newNumAllocatedKeys];
    for( int i = 0; i < myNumAllocatedEnsKeys; i++ ) {
      keyInfo[i].name = myEnsKeyInfo[i].name;
      keyInfo[i].hdrIndex = myEnsKeyInfo[i].hdrIndex;
      keyInfo[i].hdrType = myEnsKeyInfo[i].hdrType;
    }
    delete [] myEnsKeyInfo;
    myEnsKeyInfo = keyInfo;
  }
  myNumAllocatedEnsKeys = newNumAllocatedKeys;
}
void csSuperHeader::setEnsembleKey( std::string const& keyName, int keyIndex ) {
  setEnsembleKey( keyName, 0, 0, keyIndex );  // TEMP Workaround
}
void csSuperHeader::setEnsembleKey( std::string const& keyName, int hdrIndex, cseis_geolib::type_t type, int keyIndex ) {
  if( keyIndex == myNumAllocatedEnsKeys ) {
    reallocate(myNumAllocatedEnsKeys+1);
  }
  if( keyIndex == myNumEnsKeys ) {
    myNumEnsKeys += 1;
  }
  if( keyIndex < myNumAllocatedEnsKeys ) {
    myEnsKeyInfo[keyIndex].hdrIndex = hdrIndex;
    myEnsKeyInfo[keyIndex].name     = keyName;
    myEnsKeyInfo[keyIndex].hdrType  = type;
  }
  else {
    throw cseis_geolib::csException("Cannot set key. Key index higher than number of keys.");
  }
//  printf("Super header set: %s  =  %s\n", keyName.c_str(), myEnsKeyInfo[keyIndex].name.c_str() );
}
std::string const* csSuperHeader::ensembleKey( int keyIndex ) const {
  return &myEnsKeyInfo[keyIndex].name;
}
void csSuperHeader::removeEnsembleKey( int keyIndex ) {
  for( int i = keyIndex+1; i < myNumEnsKeys; i++ ) {
    myEnsKeyInfo[i-1].hdrIndex = myEnsKeyInfo[i].hdrIndex;
    myEnsKeyInfo[i-1].name     = myEnsKeyInfo[i].name;
    myEnsKeyInfo[i-1].hdrType  = myEnsKeyInfo[i].hdrType;
  }
  if( myNumEnsKeys > 0 && keyIndex < myNumEnsKeys && keyIndex >= 0 ) {
    myNumEnsKeys -= 1;
  }
}
void csSuperHeader::dump( FILE* fout ) const {
  fprintf( fout, "*** START Super header ***\n" );
  fprintf( fout, "Sample interval:        %f %s\n", sampleInt, cseis_geolib::csGeolibUtils::domain2UnitText( domain ) );
  fprintf( fout, "Number of samples:      %d\n", numSamples );
  fprintf( fout, "Domain code:            %d (= %s domain)\n", domain, cseis_geolib::csGeolibUtils::domain2Text( domain ) );
  fprintf( fout, "FFT data type:          %d (= %s)\n", fftDataType, cseis_geolib::csGeolibUtils::fftDataType2Text( fftDataType ) );
  fprintf( fout, "   FFT parameters (if applicable): XT Number of samples/sampleInt: %d / %f\n", numSamplesXT, sampleIntXT );
  fprintf( fout, "# ensemble keys:        %d\n", myNumEnsKeys );
  for( int ikey = 0; ikey < myNumEnsKeys; ikey++ ) {
    fprintf( fout, "  Key #%d:  '%s'\n", ikey+1, myEnsKeyInfo[ikey].name.c_str() );
  }
  fprintf( fout, "Grid origin row/col:    %16d     %16d\n", grid_orig_il, grid_orig_xl );
  fprintf( fout, "Grid origin XY coord.:  %16.5f m   %16.5f m\n", grid_orig_x, grid_orig_y );
  fprintf( fout, "Grid bin size (il/xl):  %16.5f m   %16.5f m  (corresponding to row/col interval = 1)\n", grid_binsize_il, grid_binsize_xl );
  fprintf( fout, "Inline/xl direction:    %16.5f deg %16.5f deg  (clock-wise from North)\n", grid_azim_il, grid_azim_xl );
  fprintf( fout, "*** END Super header ***\n" );
}


