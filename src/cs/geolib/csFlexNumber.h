/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_FLEX_NUMBER_H
#define CS_FLEX_NUMBER_H

#include "csException.h"
#include "geolib_defines.h"
#include <string>

namespace cseis_geolib {
  class csFlexHeader;

  /**
   * Flexible number class
   *
   * Memory-efficient number class for 4-byte and 8-byte floating point and integer numbers
   * Quite well optimised for speed as well..
   *
   * Differences to csNumber:
   * - No matter what number type the value is, this class always returns a valid number
   * - This means, csFlexNumber automatically converts into the number type that is requested
   *
   * @author Bjorn Olofsson
   * @date 2006
   */
  
class csFlexNumber {
public:
  csFlexNumber();
  csFlexNumber( double d );
  csFlexNumber( float f );
  csFlexNumber( int i );
  csFlexNumber( csInt64_t i );
  csFlexNumber( csFlexNumber const& obj );
  csFlexNumber( cseis_geolib::type_t type, std::string value );
  csFlexNumber( csFlexHeader const* flexHdr, bool invertPolarity = false );
  csFlexNumber( csFlexNumber const* flexNum, bool invertPolarity = false );
  void init();
  ~csFlexNumber();
  void setValueFromString( cseis_geolib::type_t type, std::string value );
  csFlexNumber& operator=( csFlexNumber const& obj );
  csFlexNumber& operator=( double d );
  csFlexNumber& operator=( int i );
  csFlexNumber& operator=( csInt64_t i );
  /**
  * Convert text string into number
  * @param text Input text string
  * @return true if text string contains a valid number, false otherwise
  */
  bool convertToNumber( std::string text );
  /// @return type of number
  inline cseis_geolib::type_t type() const { return myType; }

  inline void setFloatValue( float value ) {
    *((float*)myValue) = value;
    myType = TYPE_FLOAT;
  }
  inline void setDoubleValue( double value ) {
    *((double*)myValue) = value;
    myType = TYPE_DOUBLE;
  }
  inline void setIntValue( int value ) {
    *((int*)myValue) = value;
    myType = TYPE_INT;
  }
  inline void setInt64Value( csInt64_t value ) {
    *((csInt64_t*)myValue) = value;
    myType = TYPE_INT64;
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
      return 0;
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
      return 0;
    }
  }
  inline csInt64_t int64Value() const {
    switch( myType ) {
    case TYPE_INT64:
      return *((csInt64_t*)myValue);
    case TYPE_DOUBLE:
      return (csInt64_t)( *((double*)myValue) );
    case TYPE_FLOAT:
      return (csInt64_t)( *((float*)myValue) );
    case TYPE_INT:
      return (csInt64_t)( *((int*)myValue) );
    default:
      return 0;
    }
  }
  
  inline bool isUndefined() const { return( myType == TYPE_UNKNOWN ); }
  
  csFlexNumber& operator*=( const csFlexNumber& obj );
  csFlexNumber operator-( const csFlexNumber& obj ) const;
  csFlexNumber operator+( const csFlexNumber& obj ) const;
  csFlexNumber operator*( const csFlexNumber& obj ) const;
  csFlexNumber operator/( const csFlexNumber& obj ) const;
  csFlexNumber operator%( const csFlexNumber& obj ) const;

  bool operator==( const csFlexNumber& obj ) const;
  bool operator!=( const csFlexNumber& obj ) const;
  bool operator>=( const csFlexNumber& obj ) const;
  bool operator<=( const csFlexNumber& obj ) const;
  bool operator>( const csFlexNumber& obj ) const;
  bool operator<( const csFlexNumber& obj ) const;

  void dump() const;
private:
  static const int MAX_SIZE = 8;
  unsigned char*  myValue;
  cseis_geolib::type_t myType;

// Internal methods, without any checks
  inline csInt64_t int64Value_internal() const {
    return *((csInt64_t*)myValue);
  }
  inline int intValue_internal() const {
    return *((int*)myValue);
  }
  inline float floatValue_internal() const {
    return *((float*)myValue);
  }
  inline double doubleValue_internal() const {
    return *((double*)myValue);
  }

};

class csFlexNumberException : public csException {
public:
  csFlexNumberException( std::string message );
  ~csFlexNumberException();
private:
  csFlexNumberException();
};

}

#endif

