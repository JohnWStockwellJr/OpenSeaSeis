/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_BYTE_CONVERSIONS_H
#define CS_BYTE_CONVERSIONS_H

#include "geolib_defines.h"
#include "methods_number_conversions.h"
#include <cmath>
#include <cstring>

namespace cseis_geolib {

inline int byte2Short_SWAP( byte_t const* ptr ) {
  return( (short)( (ptr[0] << 8) + ptr[1] ) );
}
inline int byte2Short( byte_t const* ptr ) {
  return( (short)( (ptr[1] << 8) + ptr[0] ) );
}

inline unsigned short byte2UShort_SWAP( byte_t const* ptr ) {
  return( (unsigned short)( (ptr[0] << 8) + ptr[1] ) );
}
inline unsigned short byte2UShort( byte_t const* ptr ) {
  return( (unsigned short)( (ptr[1] << 8) + ptr[0] ) );
}

inline int byte2Int( byte_t const* ptr ) {
  return( (int)( (ptr[3] << 24) + (ptr[2] << 16) + (ptr[1] << 8) + ptr[0] ) );
}
inline csInt64_t byte2Int64( byte_t const* ptr ) {
  csInt64_t value;
  memcpy( &value, ptr, 8 );
  return value;
}
inline int byte2Int_SWAP( byte_t const* ptr ) {
  return( (int)( (ptr[0] << 24) + (ptr[1] << 16) + (ptr[2] << 8) + ptr[3] ) );
}

inline float byte2Float( byte_t const* ptr ) {
  return( (float)( (ptr[3] << 24) + (ptr[2] << 16) + (ptr[1] << 8) + ptr[0] ) );
}
//inline double byte2Double( byte_t const* ptr ) {
//  return( (double)( (ptr[7] << 56) + (ptr[6] << 48) + (ptr[5] << 40) + (ptr[4] << 32) + (ptr[3] << 24) + (ptr[2] << 16) + (ptr[1] << 8) + ptr[0] ) );
//}
inline float byte2Float_SWAP( byte_t const* ptr ) {
  char c[4];
  c[0] = ptr[3];
  c[1] = ptr[2];
  c[2] = ptr[1];
  c[3] = ptr[0];
  float f;
  memcpy( &f, c, 4 );
  return f;
}

inline float byte_FloatIBM2Float_SWAP( byte_t const* ptr ) {
  unsigned char c[4];
  c[0] = ptr[3];
  c[1] = ptr[2];
  c[2] = ptr[1];
  c[3] = ptr[0];
  ibm2ieee( (unsigned char*)c, 1 );
  float f;
  memcpy( &f, c, 4 );
  return f;
}

inline float byte_FloatIBM2Float( byte_t const* ptr ) {
  unsigned char c[4];
  c[0] = ptr[0];
  c[1] = ptr[1];
  c[2] = ptr[2];
  c[3] = ptr[3];
  ibm2ieee( (unsigned char*)c, 1 );
  float f;
  memcpy( &f, c, 4 );
  return f;
}

/*
inline double byte2Double_SWAP( byte_t const* ptr ) {
  char c[8];
  c[0] = ptr[3];
  c[1] = ptr[2];
  c[2] = ptr[1];
  c[3] = ptr[0];
  c[4] = ptr[7];
  c[5] = ptr[6];
  c[6] = ptr[5];
  c[7] = ptr[4];
  double d;
  memcpy( &d, c, 8 );
  return d;
}
*/
inline void short2Byte_SWAP( short value, byte_t* outPtr ) {
  outPtr[1] = value & 0xff;
  outPtr[0] = (value & 0xff00) >> 8;
}
inline void short2Byte( short value, byte_t* outPtr ) {
  outPtr[0] = value & 0xff;
  outPtr[1] = (value & 0xff00) >> 8;
}
inline void int2Byte( int value, byte_t* outPtr ) {
  outPtr[0] = value & 0xff;
  outPtr[1] = (value & 0xff00) >> 8;
  outPtr[2] = (value & 0xff0000) >> 16;
  outPtr[3] = (value & 0xff000000) >> 24;
}
inline void int2Byte_SWAP( int value, byte_t* outPtr ) {
  outPtr[3] = value & 0xff;
  outPtr[2] = (value & 0xff00) >> 8;
  outPtr[1] = (value & 0xff0000) >> 16;
  outPtr[0] = (value & 0xff000000) >> 24;
}
inline void float2Byte_SWAP( float value, byte_t* outPtr ) {
  int valueInt;
  memcpy( &valueInt, &value, 4 );
  outPtr[3] = valueInt & 0xff;
  outPtr[2] = (valueInt & 0xff00) >> 8;
  outPtr[1] = (valueInt & 0xff0000) >> 16;
  outPtr[0] = (valueInt & 0xff000000) >> 24;
}
inline void float2Byte( float value, byte_t* outPtr ) {
  memcpy( outPtr, &value, 4 );
}
inline void double2Byte_SWAP( double value, byte_t* outPtr ) {
  unsigned char ptr[8];
  memcpy( ptr, &value, 8 );
  for( int i = 0; i < 8; i++ ) {
    memcpy( &outPtr[7-i], &ptr[i], 1 );
  }
}
inline void double2Byte( double value, byte_t* outPtr ) {
  memcpy( outPtr, &value, 8 );
}

} // end namespace
#endif

