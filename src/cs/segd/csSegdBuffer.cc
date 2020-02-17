/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegdBuffer.h"
#include <cstdlib>

using namespace cseis_segd;

csSegdBuffer::csSegdBuffer() {
  myBuffer = NULL;
  mySize = 0;
}
csSegdBuffer::~csSegdBuffer() {
  clear();
}
byte* csSegdBuffer::setNumBytes( int numBytes ) {
  if( numBytes <= mySize ) return myBuffer;
  clear();
  mySize   = numBytes;
  myBuffer = new byte[mySize];
  return myBuffer;
}
void csSegdBuffer::clear() {
  if( myBuffer ) {
    delete [] myBuffer;
    myBuffer = NULL;
  }
}

