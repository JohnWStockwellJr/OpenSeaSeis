/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef METHODS_NUMBER_CONVERSIONS_H
#define METHODS_NUMBER_CONVERSIONS_H

namespace cseis_geolib {

  void convertInt2Float( int* array, int nValues );

  void convertShort2Float( short* arrayShort, float* arrayFloat, int nValues );

  /**
  * Convert from EBCDIC code to C char
  */
  char ebcdic2char( short c );

  /**
  * Convert from C char to EBCDIC code
  */
  short char2ebcdic( char c );

  void ibm2ieee( unsigned char* values, int numValues );
  void ieee2ibm( unsigned char* values, int numValues );

} // end namespace

#endif

