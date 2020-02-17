/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_COMPARE_VECTOR_HH
#define CS_COMPARE_VECTOR_HH

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "csVector.h"

namespace cseis_geolib {

/**
* Vector implementation of general collection
* Vector where elements can be compared using '==' operator
* This means, this vector should only be used with types that allow comparison with '=='
*
* @author Bjorn Olofsson
* @date 2006
*/
template <typename T>
class csCompareVector : public csVector<T> {
public:
  csCompareVector( int initial_capacity );
  csCompareVector();
  csCompareVector( csCompareVector<T> const& obj );
  virtual ~csCompareVector();
  virtual inline bool contains( T const& value ) const;
};

//-------------------------------------------------
template<typename T>csCompareVector<T>::csCompareVector() : 
  csVector<T>() {
}
//-------------------------------------------------
template<typename T>csCompareVector<T>::csCompareVector( int initialCapacity ) :
  csVector<T>(initialCapacity) {
}
//-------------------------------------------------
template<typename T>csCompareVector<T>::csCompareVector( csCompareVector<T> const& obj ) : 
  csVector<T>( obj ) {
}
//-------------------------------------------------
template<typename T>csCompareVector<T>::~csCompareVector() { 
}
//-------------------------------------------------
template<typename T> inline bool csCompareVector<T>::contains( T const& value ) const {
  for( int i = 0; i < csCollection<T>::mySize; i++ ) {
    if( csCollection<T>::myArray[i] == value ) {
      return true;
    }
  }
  return false;
}

} // END namespace

#endif


