/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGD_BUFFER_H
#define CS_SEGD_BUFFER_H

#include "csSegdDefines.h"

namespace cseis_segd {

class csSegdBuffer {
public:
  csSegdBuffer();
  ~csSegdBuffer();
  byte* setNumBytes( int numBytes );
  byte* buffer() { return myBuffer; }
  int size() { return mySize; }

private:  
  void clear();
  byte* myBuffer;
  int   mySize;
};

} // end namespace
#endif

