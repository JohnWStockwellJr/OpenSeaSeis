/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TRACE_HEADER_DATA_H
#define CS_TRACE_HEADER_DATA_H

#include <string>
#include "geolib_defines.h"
#include "csPoint3D.h"

namespace cseis_geolib {
  class csAbsoluteTime;
}

namespace cseis_system {

  class csTraceHeaderDef;
static int const DEFAULT_STRING_LENGTH = 20;

/**
* Trace header data
*
* Manages values for all trace headers of one seismic trace
* Header values are stored as a continuous char buffer
*
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTraceHeaderData {
public:
  csTraceHeaderData();
  ~csTraceHeaderData();
  /// Set/create trace header data from trace header definition
  void setHeaders( csTraceHeaderDef const* hdef, int inPort = -1 );
  /// Make exact copy of input object
  void copyFrom( csTraceHeaderData const* hdr );
  //------------------------------------------------------
  void setIntValue( int index, int value );
  void setInt64Value( int index, csInt64_t value );
  void setFloatValue( int index, float value );
  void setDoubleValue( int index, double value );
  void setStringValue( int index, std::string value );
  void setVectorValue( int index, cseis_geolib::csPoint3D point );
  void setVectorValue( int index, double value, cseis_geolib::type_t type );
  void setVectorValueX( int index, double value );
  void setVectorValueY( int index, double value );
  void setVectorValueZ( int index, double value );
  //------------------------------------------------------
  inline int intValue( int index ) const {
#ifndef ARCHITECTURE_ITANIUM
    return *(reinterpret_cast<int*>( &myValueBlock[myByteLocationPtr[index]] ));
#else
    int value;
    memcpy( &value, &myValueBlock[myByteLocationPtr[index]], 4 );
    return value;
#endif
  }
  inline csInt64_t int64Value( int index ) const {
#ifndef ARCHITECTURE_ITANIUM
    return *(reinterpret_cast<csInt64_t*>( &myValueBlock[myByteLocationPtr[index]] ));
#else
    csInt64_t value;
    memcpy( &value, &myValueBlock[myByteLocationPtr[index]], 8 );
    return value;
#endif
  }
  inline float floatValue( int index ) const {
#ifndef ARCHITECTURE_ITANIUM
    return *(reinterpret_cast<float*>( &myValueBlock[myByteLocationPtr[index]] ));
#else
    float value;
    memcpy( &value, &myValueBlock[myByteLocationPtr[index]], 4 );
    return value;
#endif
  }
  inline double doubleValue( int index ) const {
#ifndef ARCHITECTURE_ITANIUM
    return *(reinterpret_cast<double*>( &myValueBlock[myByteLocationPtr[index]] ));
#else
    double value;
    memcpy( &value, &myValueBlock[myByteLocationPtr[index]], 8 );
    return value;
#endif
  }
  cseis_geolib::csPoint3D vectorValue( int index ) const;
  double vectorValue( int index, cseis_geolib::type_t type ) const;
  std::string stringValue( int index ) const;
  //------------------------------------------------------
  inline int numBytes() const { return myNumBytes; }
  //------------------------------------------------------
  inline int numHeaders() const { return myNumHeaders; }
  //------------------------------------------------------
  void clear();
  void clearMemory();
  friend class csTraceHeader;

private:
  void deleteHeaders( csTraceHeaderDef const* hdef );
  /// Buffer that holds all trace header values in one chunk of memory
  char* myValueBlock;
  /// Maps the sequential header index 0,1,2... to the byte location index in the char* 'value block'
  /// Pointer to object, kept in csTraceHeaderDef
  int const*  myByteLocationPtr;

  int   myNumBytes;
  int   myNumAllocatedBytes;
  int   myNumHeaders;

  void reallocateBytes( int numAllocatedBytes );
  static int NUM_ADD_HEADERS;
  static int NUM_ADD_BYTES;

  // Helper method that sets header data according header definition object (Adds new headers = allocates enough memory)
  int myIndex;
  static int counter;
};

} // namespace

#endif
