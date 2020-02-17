/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_GEOLIB_UTILS_H
#define CS_GEOLIB_UTILS_H

#include "geolib_defines.h"
#include "csException.h"
#include <string>
#include <cstdio>

namespace cseis_geolib {

  /**
   * Various utilities
   *
   * @author Bjorn Olofsson
   * @date 2007
   */
class csGeolibUtils {
public:
  /// Return number of bytes of fundamental type
  static inline int numBytes( type_t type ) {
    switch( type ) {
    case TYPE_INT:
    case TYPE_ARRAY_INT:
       return (int)sizeof(int);   
    case TYPE_FLOAT:
    case TYPE_FLOAT_IBM:
    case TYPE_ARRAY_FLOAT:
       return (int)sizeof(float);  
    case TYPE_DOUBLE:
    case TYPE_ARRAY_DOUBLE:
       return (int)sizeof(double);
    case TYPE_VECTOR:
       return 3*(int)sizeof(double);
    case  TYPE_INT64:
      return (int)sizeof(csInt64_t);  
    case TYPE_CHAR:
    case TYPE_STRING:
      return (int)sizeof(char); 
    }
    return 0;
  }

  static char const* domain2UnitText( int domain ) {
    switch( domain ) {
    case DOMAIN_UNKNOWN:
      return( "unknown" );
    case DOMAIN_XT:
    case DOMAIN_KT:
      return( "ms" );
    case DOMAIN_XD:
      return( "m" );
    case DOMAIN_FX:
    case DOMAIN_FK:
      return( "Hz" );
    default:
      throw( csException("csGeolibUtils::domain2UnitText: Undefined domain type: %d", domain) );
    }
  }

  static char const* domain2Text( int domain ) {
    switch( domain ) {
    case DOMAIN_UNKNOWN:
      return( "unknown" );
    case DOMAIN_XT:
      return( "x-time" );
    case DOMAIN_KT:
      return( "k-time" );
    case DOMAIN_XD:
      return( "x-depth" );
    case DOMAIN_FX:
      return( "x-frequency" );
    case DOMAIN_FK:
      return( "k-frequency" );
    default:
      throw( csException("csGeolibUtils::domain2Text: Undefined domain type: %d", domain) );
    }
  }

  static char const* fftDataType2Text( int fftDataType ) {
    switch( fftDataType ) {
    case 0:
    case FX_NONE:
      return( "None" );
    case FX_PSD:
      return( "Power spectrum" );
    case FX_AMP:
      return( "Amplitude spectrum" );
    case FX_PHASE:
      return( "Phase spectrum" );
    case FX_AMP_PHASE:
      return( "Amplitude/phase spectrum" );
    case FX_REAL_IMAG:
      return( "Real/imaginary spectrum" );
    case FK_PSD:
      return( "Power k-spectrum" );
    case FK_AMP:
      return( "Amplitude k-spectrum" );
    case FK_AMP_PHASE:
      return( "Amplitude/phase k-spectrum" );
    case FK_REAL_IMAG:
      return( "Real/imaginary k-spectrum" );
    default:
      throw( csException("csGeolibUtils::fftDataType2Text: Undefined data type: %d", fftDataType) );
    }
  }

  static inline int numBytes( type_t type, int nElements ) {
    return nElements*numBytes( type );
  }
  
  static char const* typeText( type_t type ) {
    switch( type ) {
      case TYPE_UNKNOWN:
        return( "unknown" );
      case TYPE_EMPTY:
        return( "empty" );
      case TYPE_INT:
        return( "int" );
      case TYPE_INT64:
        return( "int64" );
      case TYPE_FLOAT:
        return( "float" );
      case TYPE_FLOAT_IBM:
        return( "float_ibm" );
      case TYPE_DOUBLE:
        return( "double" );
      case TYPE_CHAR:
        return( "char" );
      case TYPE_STRING:
        return( "string" );
      case TYPE_OPTION:
        return( "option" );
      case TYPE_FLOAT_DOUBLE:
        return( "float/double" );
      case TYPE_SHORT:
        return( "short" );
      case TYPE_USHORT:
        return( "ushort" );
      default:
        throw( csException("csGeolibUtils::typeText: Undefined type id: %d", (int)type) );
    }
  }

  static type_t text2Type( std::string text );

  /*  static type_t const TYPE_UNKNOWN = 255;
  static type_t const TYPE_EMPTY   = 0;
  static type_t const TYPE_INT     = 1;  // 32bit
  static type_t const TYPE_FLOAT   = 2;  // 32bit
  static type_t const TYPE_DOUBLE  = 3;  // 64bit
  static type_t const TYPE_CHAR    = 4;  // 8bit
  static type_t const TYPE_STRING  = 5;  // = ARRAY_CHAR
  static type_t const TYPE_INT64   = 6;  // 64bit integer
  static type_t const TYPE_SHORT   = 11;  // 16bit signed short
  static type_t const TYPE_USHORT  = 12;  // 16bit unsigned short
  static type_t const TYPE_FLOAT_IBM = 13;  // 32bit IBM float
  */

  static csInt64_t date2UNIXmsec( int year, int day, int hour, int min, int sec, int msec );
  static csInt64_t date2UNIXmsec( int year, int month, int day_of_month, int hour, int min, int sec, int msec );
  static int date2UNIXsec( int year, int day, int hour, int min, int sec );
  static int date2UNIXsec( int year, int month, int day_of_month, int hour, int min, int sec );

  static int date2JulianDay( int year, int month, int day );

  static void UNIXsec2date( int& year, int& julianDay, int& hour, int& min, int& sec, int unixSec );
  static void UNIXsec2date( int& year, int& month, int& day, int& hour, int& min, int& sec, int unixSec );

  static void julianDay2date( int& month, int& day, int year, int julianDay );
  static void printDate( FILE* stream, int unixSec );
  static std::string UNIXsec2dateString( int unixSec );
};
 
} // namespace
#endif

