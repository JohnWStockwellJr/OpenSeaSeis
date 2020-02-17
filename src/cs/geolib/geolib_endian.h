/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef GEOLIB_ENDIAN_H
#define GEOLIB_ENDIAN_H

namespace cseis_geolib {

/**
 * Endian related helper methods
 *
 * @author Bjorn Olofsson
 * @date 2007
 */

/**
*   Intel: LittleEndian
*   SPARC: BigEndian
* @return true if the current platform is a Little Endian machine, false for Big Endian
*/
bool isPlatformLittleEndian();

/**
* Endian-swap 4-byte words
* Swap every 4 bytes from (1234) to (4321)
* @param array  Data array
* @param size   Number of bytes in array (Must be multiple of 4)
*/
void swapEndian4( char* array, int size );

/**
* Endian-swap 2-byte words
* Swap every 2 bytes from (12) to (21)
* @param array  Data array
* @param size   Number of bytes in array (Must be multiple of 2)
*/
void swapEndian2( char* array, int size );

/**
* Endian-swap 8-byte words
* Swap every 8 bytes from (1234) to (4321)
* @param array  Data array
* @param size   Number of bytes in array (Must be multiple of 4)
*/
void swapEndian8( char* array, int size );

/**
* Endian-swap any-byte words
* Swap every for example 8 bytes from (12345678) to (87654321)
* @param array  Data array
* @param numBytesSwap      Number of bytes to swap (2,4,8...)
* @param numTotalNumBytes  Number of bytes in array (Must be multiple of numBytesSwap)
*/
void swapEndian( char* array, int numBytesSwap, int numTotalBytes );

} // end namespace

#endif

