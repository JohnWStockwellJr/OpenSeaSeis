/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "methods_number_conversions.h"
#include "geolib_endian.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

namespace{
const char e2a[] = {
  '@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@',
  '@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@',
  '@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@','@',
  '@','@','@','@',' ','@','@','@','@','@','@','@','@','@','[','.','<','(','+','|',
  '&','@','@','@','@','@','@','@','@','@',']','$','*',')',';','^','-','/','@','@',
  '@','@','@','@','@','@','|',',','%','_','>','?','@','@','@','@','@','@','@','@',
  '@','.',':','#','@','.','=','"',' ','a','b','c','d','e','f','g','h','i',' ',' ',
  ' ',' ',' ',' ',' ','j','k','l','m','n','o','p','q','r',' ',' ',' ',' ',' ',' ',
  ' ','~','s','t','u','v','w','x','y','z',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
  '@','@','@','@','@','@','@','@','@','@','@','@','{','A','B','C','D','E','F','G',
  'H','I',' ',' ',' ',' ',' ',' ','}','J','K','L','M','N','O','P','Q','R',' ',' ',
  ' ',' ',' ',' ','/',' ','S','T','U','V','W','X','Y','Z',' ',' ',' ',' ',' ',' ',
  '0','1','2','3','4','5','6','7','8','9',' ',' ',' ',' ',' ',' '
};

const short a2e[] = {
  32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
  32,32,32,32,32,32,32,32,32,32,32,32,64,32,127,123,91,108,80,32,
  77,93,92,78,107,96,75,97,240,241,242,243,244,245,246,247,248,249,122,94,
  76,126,110,111,0,193,194,195,196,197,198,199,200,201,209,210,211,212,213,214,
  215,216,217,226,227,228,229,230,231,232,233,74,32,90,95,109,32,129,130,131,
  132,133,134,135,136,137,145,146,147,148,149,150,151,152,153,162,163,164,165,166,
  167,168,169,192,79,208,161,32,32,32,32,32,32,32,32,32,32,32,32,32,
  32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
  32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
  32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
  32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
  32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
  32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32
};
}

/**
* Convert from 32bit integer to 32bit float
*/

void cseis_geolib::convertInt2Float( int* array, int nValues ) {
  float valueFloat;
  float* arrayFloatPtr = reinterpret_cast<float*>(array);

  for( int i = 0; i < nValues; i++ ) {
    valueFloat = (float)array[i];
    arrayFloatPtr[i] = valueFloat;
  }
}

/**
* Convert from 16bit integer to 32bit float
*/

void cseis_geolib::convertShort2Float( short* arrayShort, float* arrayFloat, int nValues ) {
  for( int i = 0; i < nValues; i++ ) {
    arrayFloat[i] = (float)arrayShort[i];
  }
}

/**
 * Convert from EBCDIC code to C char
 */
char cseis_geolib::ebcdic2char( short c ) {
  if( c < 0 ) {
    c += 256;
  }
  if( c < 0 || c > 255 ) {
    return ' ';
  }
  else {
    return e2a[c];
  }
}

/**
 * Convert from C char to EBCDIC code
 */
short cseis_geolib::char2ebcdic( char c ) {
  return a2e[ (int)c ];
}

void cseis_geolib::ibm2ieee( unsigned char* values, int numValues ) {
  register unsigned fraction;
  register int exponent;
  register int signum;

  for( int i = 0; i < numValues; i++ ) {
    memcpy( &fraction, &values[i*4], 4 );

    signum = fraction >> 31;
    fraction <<= 1;
    exponent = fraction >> 25;
    fraction <<= 7;
    
    if( fraction == 0 ) {
      exponent = 0;
    }
    else {
      exponent = (exponent << 2) - 130;
    
      while (fraction < 0x80000000) {
        --exponent;
        fraction <<= 1;
      }
    
      if( exponent <= 0 ) {
        if( exponent < -24 ) {
          fraction = 0;
        }
        else {
          fraction >>= -exponent;
        }
        exponent = 0;
      }
      else if( exponent >= 255 ) {
        fraction = 0;
        exponent = 255;
      }
      else {
        fraction <<= 1;
      }
    }
    
    fraction = (fraction >> 9) | (exponent << 23) | (signum << 31);
    memcpy( &values[i*4], &fraction, 4 );
  }
}

void cseis_geolib::ieee2ibm( unsigned char* values, int numValues ) {
  register unsigned fraction;
  register int exponent;
  register int signum;

  for( int i = 0; i < numValues; i++ ) {
    memcpy( &fraction, &values[i*4], 4 );
    signum = fraction >> 31;
    fraction <<= 1;
    exponent = fraction >> 24;
    fraction <<= 8;

    if( exponent > 0 && exponent != 255 ) {
      fraction = (fraction >> 1) | 0x80000000;
      exponent += 130;
      fraction >>= -exponent & 3;
      exponent = (exponent + 3) >> 2;
    
      while (fraction < 0x10000000) {
        --exponent;
        fraction <<= 4;
      }
    }
    else { // fraction == 0 || fraction == 255
      if( exponent == 255 ) {
        fraction = 0xffffff00;
        exponent = 0x7f;
      }
    }
    
    fraction = (fraction >> 8) | (exponent << 24) | (signum << 31);
    memcpy( &values[i*4], &fraction, 4 );
  }

}
