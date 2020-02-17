/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGY_TRACE_HEADER_H
#define CS_SEGY_TRACE_HEADER_H

#include "geolib_defines.h"
#include "csHeaderInfo.h"

namespace cseis_geolib {

template<typename T> class csVector;
class csFlexHeader;
class csSegyHdrMap;

/**
 * Segy trace header
 * @author Bjorn Olofsson
 * @date 2008
 */
class csSegyTraceHeader {

public:
  csSegyTraceHeader( csSegyHdrMap const* hdrMap );
  ~csSegyTraceHeader();
  int numHeaders() const;
  cseis_geolib::csHeaderInfo const* info( int hdrIndex ) const;

  int byteLoc( int hdrIndex ) const;
  int bytePos( int hdrIndex ) const;
  char const* headerDesc( int hdrIndex ) const;
  char const* headerName( int hdrIndex ) const;
  int headerIndex( std::string name ) const;
  type_t headerType( int hdrIndex ) const;
    
  int intValue( int hdrIndex ) const;
  float floatValue( int hdrIndex ) const;
  double doubleValue( int hdrIndex ) const;
  std::string stringValue( int hdrIndex ) const;
  
  void setIntValue( int hdrIndex, int value );
  void setFloatValue( int hdrIndex, float value );
  void setDoubleValue( int hdrIndex, double value );
  void setStringValue( int hdrIndex, std::string const& value );

  /**
   * Dump contents of 240 block as 2-byte and 4-byte fields
   * @param buffer       (i) 240 byte block
   * @param doSwapEndian (i) True if endian byte order needs to be swapped
   * @param fout         (i/o) Output stream
   */
  static void dump( byte_t const* buffer, bool doSwapEndian, FILE* fout );

  /*
    int getInt( int bytePos, int byteSize ) const;
    float getFloat( int bytePos, int byteSize ) const;
    double getDouble( int bytePos, int byteSize ) const;
    std::string getString( int bytePos, int byteSize ) const;
  */

  friend class csSegyReader;
  friend class csSegyWriter;  
  
  void readHeaderValues( byte_t const* buffer, bool doSwapEndian, bool isAutoScaleHeaders );
  void writeHeaderValues( byte_t* buffer, bool doSwapEndian, bool isAutoScaleHeaders ) const;

 protected:
  cseis_geolib::csFlexHeader* getHandleHdrValues();

private:
  void init();
  csSegyTraceHeader();
  cseis_geolib::csFlexHeader* myHdrValues;
  csSegyHdrMap const* myHdrMapPtr;
};

} // end namespace
#endif

