/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_VECTOR_H
#define CS_VECTOR_H

#include <cstdlib>
#include <cstdio>
#include "csException.h"
#include "csCollection.h"

namespace cseis_geolib {

/**
 * Vector implementation of general collection
 *
 * @author Bjorn Olofsson
 * @date 2006
 */
template <typename T>
class csVector : public csCollection<T> {
public:
  csVector( int initial_capacity );
  csVector();
  csVector( csVector<T> const& obj );
  virtual ~csVector();
  virtual inline void insertStart( T const& value );
  /// Returns index of inserted value
  virtual inline int insertEnd( T const& value );
  virtual inline int insert( T const& value );
  virtual inline void insert( T const& value, int pos );
  virtual inline void set( T const& value, int pos );
  virtual inline T& at( int index ) const;
  /// @return last element in list
  virtual inline T& last() const;
  virtual inline void remove( int index );
  virtual inline void remove( int firstIndex, int nElements );
  /// Allocate more elements and initialize them to the given value
  virtual inline void allocate( int nElements, T value );
};

//-------------------------------------------------
template<typename T>csVector<T>::csVector() : 
  csCollection<T>() {
}
//-------------------------------------------------
template<typename T>csVector<T>::csVector( int initialCapacity ) :
  csCollection<T>(initialCapacity) {
}
//-------------------------------------------------
template<typename T>csVector<T>::csVector( csVector<T> const& obj ) : 
  csCollection<T>( obj ) {
}
//-------------------------------------------------
template<typename T>csVector<T>::~csVector() { 
}
//-------------------------------------------------
template<typename T> inline void csVector<T>::set( T const& value, int pos ) {
  if( pos >= 0 && pos < csCollection<T>::mySize ) {
    csCollection<T>::myArray[pos] = value;
  }
}
//-------------------------------------------------
template<typename T> inline void csVector<T>::insertStart( T const& value ) {
  insert( value, 0 );
}
template<typename T> inline int csVector<T>::insertEnd( T const& value ) {
  if( csCollection<T>::mySize == csCollection<T>::myCapacity ) {
    this->reallocate();
  }
  csCollection<T>::myArray[csCollection<T>::mySize] = value;
  return csCollection<T>::mySize++;
}
template<typename T> inline int csVector<T>::insert( T const& value ) {
  return insertEnd( value );
}
template<typename T> inline void csVector<T>::insert( T const& value, int pos ) {
  if( pos > csCollection<T>::mySize ) pos = csCollection<T>::mySize;
  else if( pos < 0 ) pos = 0;
  if( csCollection<T>::mySize == csCollection<T>::myCapacity ) {
    this->reallocate();
  }
  csCollection<T>::mySize++;
  for( int i = csCollection<T>::mySize-1; i > pos; i-- ) {
    csCollection<T>::myArray[i] = csCollection<T>::myArray[i-1];
  }
  csCollection<T>::myArray[pos] = value;
}
//-------------------------------------------------
template<typename T> inline T& csVector<T>::at( int index ) const {
  //  if( index >= 0 && index < csCollection<T>::mySize )
    return csCollection<T>::myArray[index];
    // else {
    // throw( csException(1,"Wrong index passed to vector object") );
    // }
}
//-------------------------------------------------
template<typename T> inline T& csVector<T>::last() const {
  //  if( index >= 0 && index < csCollection<T>::mySize )
  if( csCollection<T>::mySize > 0 ) {
    return csCollection<T>::myArray[csCollection<T>::mySize-1];
  }
  else {
    throw( csException("csVector::last(): Called on empty vector object") );
  }
}
//-------------------------------------------------
template<typename T> inline void csVector<T>::remove( int index ) {
  if( index >= 0 && index < csCollection<T>::mySize ) {
    for( int i = index+1; i < csCollection<T>::mySize; i++ ) {
      csCollection<T>::myArray[i-1] = csCollection<T>::myArray[i];
    }
    csCollection<T>::mySize--;
  }
  else {
    throw( csException(1,"Wrong index passed to vector object, remove function") );
  }
}
template<typename T> inline void csVector<T>::remove( int firstIndex, int nElements ) {
  if( nElements <= 0 ) return;
  if( firstIndex >= 0 && firstIndex+nElements <= csCollection<T>::mySize ) {
    for( int i = firstIndex+nElements; i < csCollection<T>::mySize; i++ ) {
      csCollection<T>::myArray[i-nElements] = csCollection<T>::myArray[i];
    }
    csCollection<T>::mySize -= nElements;
  }
  else {
    throw( csException(1,"Wrong firstIndex or nElements passed to vector object, remove(firstIndex,nElements) function") );
  }
}
 template<typename T> inline void csVector<T>::allocate( int nElements, T value ) {
  if( nElements <= 0 ) return;
  int size_old = csCollection<T>::mySize;
  csCollection<T>::mySize += nElements;
  if( csCollection<T>::mySize >= csCollection<T>::myCapacity ) {
    this->reallocate( csCollection<T>::mySize );
  }
  for( int i = 0; i < nElements; i++ ) {
    csCollection<T>::myArray[size_old+i] = value;
  }
}

} // END namespace

#endif


