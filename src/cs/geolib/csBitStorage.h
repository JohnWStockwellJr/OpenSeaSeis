/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_BIT_STORAGE_H
#define CS_BIT_STORAGE_H

#include "geolib_defines.h"
#include "csException.h"
#include <string>
#include <cstdio>

namespace cseis_geolib {

/**
 * Memory-efficient storage of information
 */
class csBitStorage {
  static int const NUM_BITS_PER_ELEMENT = 32;
public:
  csBitStorage() {
    myArray = NULL;
    myNumElements = 0;
    myMaxUsedBit  = 0;
    reallocate(1);
  }
  csBitStorage( int numBits ) {
    myArray = NULL;
    myNumElements = 0;
    myMaxUsedBit  = numBits;
    reallocate( (int)( (numBits-1)/NUM_BITS_PER_ELEMENT) + 1 );
  }
  ~csBitStorage() {
    delete [] myArray;
    myArray = NULL;
  }

  /**
   * Set one bit to 1
   */
  void setBit( int bit ) {
    int elemIndex = (int)( bit / NUM_BITS_PER_ELEMENT );
    if( elemIndex >= myNumElements ) {
      reallocate( elemIndex+1 );
    }
    int bitRed = bit - elemIndex * NUM_BITS_PER_ELEMENT;
    int num = 1 << bitRed;
    myArray[elemIndex] |= num;
    if( bit > myMaxUsedBit ) myMaxUsedBit = bit;
  }

  /**
   * Unset one bit (set to 0)
   */
  void unsetBit( int bit ) {
    int elemIndex = (int)( bit / NUM_BITS_PER_ELEMENT );
    if( elemIndex >= myNumElements ) {
      reallocate( elemIndex+1 );
    }
    int bitRed = bit - elemIndex * NUM_BITS_PER_ELEMENT;
    int num = 1 << bitRed;
    int dummy = 0;
    dummy |= num;
    myArray[elemIndex] |= num;   // Set bit to 1 (if not set yet)
    myArray[elemIndex] -= dummy; // Remove bit again, set to zero
    if( bit > myMaxUsedBit ) myMaxUsedBit = bit;
  }

  /**
   * @return True if indicated bit is set
   */
  bool isBitSet( int bit ) const {
    int elemIndex = (int)( bit / NUM_BITS_PER_ELEMENT );
    if( elemIndex >= myNumElements ) return false;
    int bitRed = bit - elemIndex * NUM_BITS_PER_ELEMENT;
    int num = 1 << bitRed;
    bool isSet = myArray[elemIndex] & num;
    return isSet;
  }

  /**
   * @return Number of internal storage elements
   */
  int numElements() const { return myNumElements; }

  /**
   * @return Number of used bits (set or not set)
   */
  int numBits() const {
    return myMaxUsedBit;
  }

  /**
   * @return Number of allocated bits (set or not set) in internal storage
   */
  int numAllocatedBits() const {
    return myNumElements * NUM_BITS_PER_ELEMENT;
  }

  /**
   * @return Number of bits that are set
   */
  int bitCount() const {
    int counter = 0;
    for( int ielem = 0; ielem < myNumElements; ielem++ ) {
      unsigned int num = myArray[ielem];
      for( ; num; counter++ ) {
        num &= num - 1; // clear the least significant bit set
      }
    }
    return counter;
  }

  void dump( FILE* stream ) {
    for( int i = 0; i < myNumElements; i++ ) {
      fprintf(stream,"#%d %-11d  ", i, myArray[i]);
    }
  }

private:
  void reallocate( int numElements ) {
    if( numElements <= myNumElements ) {
      myNumElements = numElements;
      return;
    }
    int* newArray = new int[numElements];
    for( int i = 0; i < numElements; i++ ) {
      newArray[i] = 0;
    }
    if( myArray != NULL ) {
      for( int i = 0; i < myNumElements; i++ ) {
        newArray[i] = myArray[i];
      }
      delete [] myArray;
    }
    myNumElements = numElements;
    myArray = newArray;
  }

  int myMaxUsedBit;
  int myNumElements;
  int* myArray;
};

} // end namespace
#endif

