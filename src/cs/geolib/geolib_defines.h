/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef GEOLIB_DEFINES_H
#define GEOLIB_DEFINES_H

#include <cstdlib>
#include "geolib_platform_dependent.h"
#include <iostream>
#include <sys/types.h>

// Define internal 64bit integer type because of missing standard on supported platforms (Windows/Linux)
typedef signed long long csInt64_t;

namespace cseis_geolib {

 typedef unsigned char type_t;
 typedef unsigned char byte_t;

/// Trace and parameter types
 static type_t const TYPE_UNKNOWN = 255;
 static type_t const TYPE_EMPTY   = 0;
 static type_t const TYPE_INT     = 1;  // 32bit
 static type_t const TYPE_FLOAT   = 2;  // 32bit
 static type_t const TYPE_DOUBLE  = 3;  // 64bit
 static type_t const TYPE_CHAR    = 4;  // 8bit
 static type_t const TYPE_STRING  = 5;  // = ARRAY_CHAR
 static type_t const TYPE_INT64   = 6;  // 64bit integer
 static type_t const TYPE_SHORT   = 11;  // 16bit signed short
 static type_t const TYPE_USHORT  = 12;  // 16bit unsigned short
 static type_t const TYPE_FLOAT_IBM  = 13;  // 32bit IBM float

 static type_t const TYPE_ARRAY_INT    = 20;
 static type_t const TYPE_ARRAY_FLOAT  = 21;
 static type_t const TYPE_ARRAY_DOUBLE = 22;
 static type_t const TYPE_VECTOR   = 23;
 static type_t const TYPE_VECTOR_X = 24;
 static type_t const TYPE_VECTOR_Y = 25;
 static type_t const TYPE_VECTOR_Z = 26;
 static type_t const TYPE_OPTION       = 30;
 static type_t const TYPE_FLOAT_DOUBLE = 40;

 static int const FX_REAL_IMAG = 501;
 static int const FX_AMP_PHASE = 502;
 static int const FX_AMP       = 503;
 static int const FX_PSD       = 504;
 static int const FX_PHASE     = 505;
 static int const FX_NONE      = 500;

 static int const FK_NONE      = FX_NONE;
 static int const FK_REAL_IMAG = 601;
 static int const FK_AMP_PHASE = 602;
 static int const FK_AMP       = 603;
 static int const FK_PSD       = 604;
 static int const FK_COMPLEX   = 605;

 /// Data domain types
 // Never change these numbers, for backward compatilbility
 static int const DOMAIN_UNKNOWN = -1;
 /// Time domain
 static int const DOMAIN_XT = 110;
 /// Frequency-trace domain
 static int const DOMAIN_FX = 111;
 /// Frequency-K domain
 static int const DOMAIN_FK = 112;
 /// Time series-K domain
 static int const DOMAIN_KT = 113;
 /// Depth domain
 static int const DOMAIN_XD = 114;

 /// Vertical domain
 static int const VERTICAL_DOMAIN_TIME  = 1001;
 static int const VERTICAL_DOMAIN_DEPTH = 1002;
 static int const VERTICAL_DOMAIN_FREQ  = 1003;

 // Matlab variables
 static int  const    MATLAB_DAYS_01JAN1970 = 719529;         // Matlab 'datenum' for 01-Jan-1970 (Number of days since 01-Jan-0000)
 static csInt64_t const MATLAB_SEC_01JAN1970  = 62167305600LL;  // Matlab 'datenum' for 01-Jan-1970 (Number of seconds since 01-Jan-0000)

} // end namespace
#endif

