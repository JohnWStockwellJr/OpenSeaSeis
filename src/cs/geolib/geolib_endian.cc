/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "geolib_endian.h"
#include <cstring>

bool cseis_geolib::isPlatformLittleEndian() {
  union {
    unsigned char  cc[2];
    unsigned short int si;
  } aa;
  aa.si = 258;

  // Intel: aa.cc[0]=2  // LittleEndian
  // SPARC: aa.cc[0]=1  // BigEndian

  return( aa.cc[0] == 2 );
}

void cseis_geolib::swapEndian4( char* array, int size ) {
  char tmp;
  if( size == 4 ) {
    tmp = array[0];
    array[0] = array[3];
    array[3] = tmp;
    tmp = array[1];
    array[1] = array[2];
    array[2] = tmp;
  }
  else {
    for( int i = 0; i < size; i+=4 ) {
      tmp        = array[i+3];
      array[i+3] = array[i];
      array[i]   = tmp;
      tmp        = array[i+2];
      array[i+2] = array[i+1];
      array[i+1] = tmp;
    }
  }
}

void cseis_geolib::swapEndian2( char* array, int size ) {
  if( size == 2 ) {
    char c2[2];
    memcpy( c2, array, 2  );
    array[0] = c2[1];
    array[1] = c2[0];
  }
  else {
    char* arrayNew = new char[size];
    memcpy( arrayNew, array, size );
    for( int i = 0; i < size; i+=2 ) {
      array[i]   = arrayNew[i+1];
      array[i+1] = arrayNew[i];
    }
    delete [] arrayNew;
  }
}

void cseis_geolib::swapEndian( char* array, int numBytesSwap, int numTotalBytes ) {
  char* arrayNew = new char[numTotalBytes];
  memcpy( arrayNew, array, numTotalBytes );
  for( int startPos = 0; startPos < numTotalBytes; startPos += numBytesSwap ) {
    int endPos = startPos+numBytesSwap-1;
    for( int ib = 0; ib < numBytesSwap; ib++ ) {
      array[startPos+ib] = arrayNew[endPos-ib];
    }
  }
  delete [] arrayNew;
}
void cseis_geolib::swapEndian8( char* array, int size ) {
  char* arrayNew = new char[size];
  memcpy( arrayNew, array, size );
  for( int i = 0; i < size; i+=8 ) {
    array[i]   = arrayNew[i+7];
    array[i+1] = arrayNew[i+6];
    array[i+2] = arrayNew[i+5];
    array[i+3] = arrayNew[i+4];
    array[i+4] = arrayNew[i+3];
    array[i+5] = arrayNew[i+2];
    array[i+6] = arrayNew[i+1];
    array[i+7] = arrayNew[i];
  }
  delete [] arrayNew;
}

