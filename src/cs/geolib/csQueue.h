/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_QUEUE_H
#define CS_QUEUE_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "csException.h"
#include "csCollection.h"

namespace cseis_geolib {

/**
* Queue implementation of general collection
*
* @author Bjorn Olofsson
* @date 2006
*/
template <typename T>
class csQueue : public csCollection<T> {
public:
  csQueue( int initial_capacity );
  csQueue();
  csQueue( csQueue<T> const& obj );
  virtual ~csQueue();
  virtual inline void push( T const& value );
  virtual inline T& peek() const;
  virtual inline T pop();

  virtual void dump() const;
  virtual void setMinCapacity( int minCapacity );

// Only for debugging:
  virtual inline T& at( int index ) const;
private:
  void init();  // Overrides base class init() function
  /// Index of first element. May be greater or equal to 0. All elements with smaller index are obsolete.
  int myFirstIndex;
  /// Minimum capacity of array
  int myMinCapacity;

  virtual inline void reallocate();
  virtual inline void reallocate( int capacity );
};

//-------------------------------------------------
template<typename T>csQueue<T>::csQueue() : 
  csCollection<T>() {
  init();
}
//-------------------------------------------------
template<typename T>csQueue<T>::csQueue( int initialCapacity ) :
  csCollection<T>(initialCapacity) {
  init();
}
//-------------------------------------------------
template<typename T>csQueue<T>::csQueue( csQueue<T> const& obj ) : 
  csCollection<T>( obj ) {
  init();
  if( csCollection<T>::mySize > 0 ) {
    for( int i = 0; i < csCollection<T>::mySize; i++ ) {
      csCollection<T>::myArray[i] = obj.myArray[i+obj.myFirstIndex];
    }
  }
}
//-------------------------------------------------
template<typename T>csQueue<T>::~csQueue() { 
}
//-------------------------------------------------
template<typename T> void csQueue<T>::init() {
  myFirstIndex = 0;
  myMinCapacity = std::max(10,csCollection<T>::myCapacity);
}
//-------------------------------------------------
template<typename T> void csQueue<T>::setMinCapacity( int minCapacity ) {
  myMinCapacity = minCapacity;
}
//-------------------------------------------------
template<typename T> inline void csQueue<T>::push( T const& value ) {
  if( csCollection<T>::mySize+myFirstIndex == csCollection<T>::myCapacity ) {
    if( csCollection<T>::mySize < myFirstIndex && csCollection<T>::myCapacity > myMinCapacity ) {
      // Reshuffle elements to beginning of array instead of reallocating
      int nElements = csCollection<T>::mySize;
      for( int i = 0; i < nElements; i++ ) {
        csCollection<T>::myArray[i] = csCollection<T>::myArray[i+myFirstIndex];
      }
      myFirstIndex = 0;
    }
    else {
      csCollection<T>::reallocate();
    }
  }
  csCollection<T>::myArray[myFirstIndex+csCollection<T>::mySize++] = value;
}
//-------------------------------------------------
template<typename T> inline T& csQueue<T>::peek() const {
  if( csCollection<T>::mySize > 0 ) {
    return csCollection<T>::myArray[myFirstIndex];
  }
  else {
    throw( csException("Attempted peek() on empty Queue object") );
  }
}
//-------------------------------------------------
template<typename T> inline T csQueue<T>::pop() {
  if( csCollection<T>::mySize > 0 ) {
    T obj = csCollection<T>::myArray[myFirstIndex];
    csCollection<T>::mySize--;
    if( csCollection<T>::mySize > 0 ) {
      myFirstIndex++;
    }
    else {
      myFirstIndex = 0;
    }
    return obj;
  }
  else {
    throw( csException("Attempted pop() on empty Queue object") );
  }
}
template<typename T> void csQueue<T>::reallocate() {
  csCollection<T>::reallocate();
}
template<typename T> void csQueue<T>::reallocate( int capacity ) {
  csCollection<T>::myCapacity = capacity;
  T* array_new = new T[csCollection<T>::myCapacity];
  int save = MIN( csCollection<T>::myCapacity, csCollection<T>::mySize );
  for( int i = 0; i < save; i++ ) {
    array_new[i] = csCollection<T>::myArray[i+myFirstIndex];   // Careful: May cause problems if T is a pointer type..?
  }
  if( csCollection<T>::myArray ) {
    delete [] csCollection<T>::myArray;   // Careful: May cause problems if T is a pointer type..?
  }
  myFirstIndex = 0;
  csCollection<T>::myArray = array_new;
  csCollection<T>::myNumReallocateCalls++;
}

template<typename T> void csQueue<T>::dump() const{
  fprintf(stdout,"*** Queue DUMP ****\n");
  fprintf(stdout,"myFirstIndex: %2d, mySize: %2d, myCapacity: %2d, myMinCapacity: %2d\n",
      myFirstIndex, csCollection<T>::mySize, csCollection<T>::myCapacity, myMinCapacity);
}

//-------------------------------------------------
template<typename T> inline T& csQueue<T>::at( int index ) const {
  if( index >= 0 && index < csCollection<T>::mySize )
    return csCollection<T>::myArray[index+myFirstIndex];
  else {
    throw( csException("Wrong index passed to vector object") );
  }
}

} // END namespace

#endif


