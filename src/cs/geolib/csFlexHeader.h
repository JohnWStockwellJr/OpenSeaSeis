/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_FLEX_HEADER_H
#define CS_FLEX_HEADER_H

#include "geolib_defines.h"
#include <algorithm>
#include <string>
#include <cstring>

namespace cseis_geolib {

  /**
   * Flexible header class.
   *
   * This class stores header value in a memory efficient way.
   *
   * Supported formats are
   *  int    (32bit)
   *  float  (32bit)
   *  double (64bit)
   *  int64  (64bit)
   *  char   (8bit)
   *  string (variable size)
   *
   * @author Bjorn Olofsson
   * @date 2007
   */

class csFlexHeader {
 public:
  csFlexHeader();
  csFlexHeader( double d );
  csFlexHeader( float f );
  csFlexHeader( int i );
  csFlexHeader( char c );
  csFlexHeader( csInt64_t i );
  csFlexHeader( std::string& text );
  csFlexHeader( csFlexHeader const& obj );
  csFlexHeader( cseis_geolib::type_t type, std::string value );
  void init();
  void init( int byteSize );
  int byteSize() const { return myByteSize; }
  ~csFlexHeader();
  csFlexHeader& operator=( const csFlexHeader& obj );
  csFlexHeader& operator=( double d );
  csFlexHeader& operator=( int i );
  csFlexHeader& operator=( char c );
  csFlexHeader& operator=( csInt64_t i );
  csFlexHeader& operator=( std::string& text );
  /// @return  of number
  inline cseis_geolib::type_t type() const { return myType; }
  void setValueFromString( cseis_geolib::type_t type, std::string value );
  inline void setStringValue( std::string& value ) {
    int newByteSize = std::max( MIN_BYTE_SIZE, (int)value.length()+1 );
    if( myByteSize < newByteSize ) {
      freeMemory();
      myValue = new char[newByteSize];
    }
    myByteSize = newByteSize;
    memcpy( myValue, value.c_str(), value.length() );
    myValue[value.length()] = '\0';
    myType = TYPE_STRING;
  }
  inline void setFloatValue( float value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    *((float*)myValue) = value;
    myType = TYPE_FLOAT;
  }
  inline void setDoubleValue( double value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    *((double*)myValue) = value;
    myType = TYPE_DOUBLE;
  }
  inline void setIntValue( int value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    *((int*)myValue) = value;
    myType = TYPE_INT;
  }
  inline void setInt64Value( csInt64_t value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    *((csInt64_t*)myValue) = value;
    myType = TYPE_INT64;
  }
  inline void setCharValue( char value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    myValue[0] = value;
    myValue[1] = '\0';
    myType = TYPE_CHAR;
  }

  inline char charValue() const {
    switch( myType ) {
    case TYPE_CHAR:
    case TYPE_STRING:
      return myValue[0];
    default:
      return ' ';
    }
  }
  inline int intValue() const {
    switch( myType ) {
    case TYPE_INT:
      return *((int*)myValue);
    case TYPE_FLOAT:
      return (int)( *((float*)myValue) );
    case TYPE_DOUBLE:
      return (int)( *((double*)myValue) );
    case TYPE_INT64:
      return (int)( *((csInt64_t*)myValue) );
    default:
      return 0;
    }
  }
  inline csInt64_t int64Value() const {
    switch( myType ) {
    case TYPE_INT64:
      return *((csInt64_t*)myValue);
    case TYPE_INT:
      return (csInt64_t)( *((int*)myValue) );
    case TYPE_FLOAT:
      return (csInt64_t)( *((float*)myValue) );
    case TYPE_DOUBLE:
      return (csInt64_t)( *((double*)myValue) );
    default:
      return 0;
    }
  }
  inline float floatValue() const {
    switch( myType ) {
    case TYPE_FLOAT:
      return *((float*)myValue);
    case TYPE_DOUBLE:
      return (float)( *((double*)myValue) );
    case TYPE_INT:
      return (float)( *((int*)myValue) );
    case TYPE_INT64:
      return (float)( *((csInt64_t*)myValue) );
    default:
      return 0.0;
    }
  }
  inline double doubleValue() const {
    switch( myType ) {
    case TYPE_DOUBLE:
      return *((double*)myValue);
    case TYPE_FLOAT:
      return (double)( *((float*)myValue) );
    case TYPE_INT:
      return (double)( *((int*)myValue) );
    case TYPE_INT64:
      return (double)( *((csInt64_t*)myValue) );
    default:
      return 0.0;
    }
  }
  inline char const* stringValue() const {
    switch( myType ) {
    case TYPE_STRING:
    case TYPE_CHAR:
      return myValue;
    default:
      return "";
    }
  }
  int stringSize();
  
  inline bool isUndefined() const { return( myType == TYPE_UNKNOWN ); }

  bool operator==( const csFlexHeader& obj ) const;
  bool operator!=( const csFlexHeader& obj ) const;
  bool operator>=( const csFlexHeader& obj ) const;
  bool operator<=( const csFlexHeader& obj ) const;
  bool operator>( const csFlexHeader& obj ) const;
  bool operator<( const csFlexHeader& obj ) const;

  void dump() const;
  std::string toString() const;
private:
  static const int MIN_BYTE_SIZE;
  char*    myValue;
  cseis_geolib::type_t   myType;
  int      myByteSize;

// Internal methods, without any checks
  inline int intValue_internal() const {
    return *((int*)myValue);
  }
  inline csInt64_t int64Value_internal() const {
    return *((csInt64_t*)myValue);
  }
  inline float floatValue_internal() const {
    return *((float*)myValue);
  }
  inline double doubleValue_internal() const {
    return *((double*)myValue);
  }
  inline char const* stringValue_internal() const {
    return myValue;
  }
  inline void freeMemory() {
    if( myValue != NULL ) {
      delete [] myValue;
      myValue = NULL;
    }
  }

};

}

#endif

