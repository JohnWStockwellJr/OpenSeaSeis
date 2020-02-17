/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_STACK_H
#define CS_STACK_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "csException.h"
#include "csCollection.h"

namespace cseis_geolib {

/**
* Stack implementation of general collection
*
* @author Bjorn Olofsson
* @date 2007
*/
template <typename T>
class csStack : public csCollection<T> {
public:
  csStack( int initial_capacity );
  csStack();
  csStack( csStack<T> const& obj );
  virtual ~csStack();
  inline void push( T const& value );
  inline T pop();
  inline T peek() const;
};

//-------------------------------------------------
template<typename T>csStack<T>::csStack() : 
  csCollection<T>() {
}
//-------------------------------------------------
template<typename T>csStack<T>::csStack( int initialCapacity ) :
  csCollection<T>(initialCapacity) {
}
//-------------------------------------------------
template<typename T>csStack<T>::csStack( csStack<T> const& obj ) : 
  csCollection<T>( obj ) {
}
//-------------------------------------------------
template<typename T>csStack<T>::~csStack() { 
}
//-------------------------------------------------
template<typename T> inline void csStack<T>::push( T const& value ) {
  if( csCollection<T>::mySize == csCollection<T>::myCapacity ) {
    this->reallocate();
  }
  csCollection<T>::myArray[csCollection<T>::mySize++] = value;
}
//-------------------------------------------------
template<typename T> inline T csStack<T>::peek() const {
  if( csCollection<T>::mySize > 0 ) {
    return csCollection<T>::myArray[csCollection<T>::mySize-1];
  }
  else {
    throw( csException("Peek() attempted on empty stack") );
  }
}
//-------------------------------------------------
template<typename T> inline T csStack<T>::pop() {
  if( csCollection<T>::mySize > 0 ) {
    return csCollection<T>::myArray[--csCollection<T>::mySize];
  }
  else {
    throw( csException("Pop() attempted on empty stack") );
  }
}

} // END namespace

#endif


