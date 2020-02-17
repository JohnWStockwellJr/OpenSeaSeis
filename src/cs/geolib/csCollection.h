/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_COLLECTION_H
#define CS_COLLECTION_H

namespace cseis_geolib {

/**
 * Base class for collections.
 * Is of no use as-is, as there is no possibility to add/access collected objects.
 * Derive from this class and add access methods according to needs.
 * May be replaced by C++ STL collection. The main reason for coding a new collection class suite is that
 * in my experience, the STL collection classes were (still are?) not always supported by different compilers on all platforms.
 *
 * @author Bjorn Olofsson
 * @date 2006
 */
template <typename T>
class csCollection {
public:
  /// Default constructor
  csCollection();
  /**
   * Pre-allocate at least the given capacity (=number of elements)
   * @param initial_capacity Number of elements to pre-allocate
   */
  csCollection( int initial_capacity );
  /// Copy constructor
  csCollection( csCollection const& obj );
  /// Assign operator
  virtual csCollection<T>& operator=( const csCollection<T>& obj );
  /// Destructor
  virtual ~csCollection();
  /// @return true if collection is empty
  virtual inline bool isEmpty() const;
  virtual inline int capacity() const;
  // @return true if collection contains 'value'. Does a simple '==' test only
  // virtual inline bool contains( T const& value ) const;
  virtual inline int size() const;
  /// Return const pointer to object array
  virtual inline T const* toArray() const;
  /// Clear collection. This may or may not free some memory
  virtual inline void clear();
  /// Dispose collection. This first clears the collection and then releases as much memory as possible
  virtual inline void dispose();
  /**
   * Reduce collection size (number of elements).
   * @param nElements Number fo elements
   */
  virtual inline void reduceSizeTo( int nElements );
  /// Dump contents to standard (or error) output
  virtual void dump() const;
protected:
  static const int DEFAULT_CAPACITY = 10;
  int myCapacity;
  int mySize;
  int myChunkSize;
  int myNumReallocateCalls;
  T* myArray;
  virtual inline void reallocate();
  virtual inline void reallocate( int capacity );
  /// Checks if capacity is filled up. If yes, reallocate another chunk
  //  virtual inline void checkCapacityReallocate();
  void init();  // Not a virtual method!
};

//int template<typename T>csCollection<T>::DEFAULT_CAPACITY = 10;

//-------------------------------------------------
template<typename T>csCollection<T>::csCollection() {
  myCapacity = DEFAULT_CAPACITY;
  init();
}
//-------------------------------------------------
template<typename T>csCollection<T>::csCollection( int initialCapacity ) {
  if( initialCapacity > 0 ) {
    myCapacity = initialCapacity;
  }
  else {
    myCapacity = DEFAULT_CAPACITY;
  }
  init();
}
//-------------------------------------------------
template<typename T>csCollection<T>::csCollection( csCollection const& obj ) {
  myCapacity = obj.myCapacity;  
  init();
  mySize = obj.size();
  for( int i = 0; i < mySize; i++ ) {
    myArray[i] = obj.myArray[i];
  }
}
//-------------------------------------------------
template <typename T> csCollection<T>& csCollection<T>::operator=( const csCollection<T>& obj ) { 
  //  printf("= constructor called %x...\n", this);
  mySize = 0;
  if( myCapacity < obj.mySize ) {
    myCapacity = obj.myCapacity;
    myChunkSize = obj.myChunkSize;
    // myArray = NULL;  // necessary in case this instance is just being allocated?
    reallocate( myCapacity );
  }
  myNumReallocateCalls = 0;
  mySize = obj.size();
  for( int i = 0; i < mySize; i++ ) {
    myArray[i] = obj.myArray[i];
  }
  //  printf("x = constructor left %x...\n", this);
  return *this;
}
//-------------------------------------------------
template<typename T>csCollection<T>::~csCollection() { 
  //  printf("Destructor called %x...\n", this);
  if( myArray != NULL ) {
    delete [] myArray;
    myArray = NULL;
  }
  mySize = 0;
  myCapacity = 0;
}
//-------------------------------------------------
template<typename T> void csCollection<T>::init() {
  mySize = 0;
  myChunkSize = myCapacity;
  myNumReallocateCalls = 0;
  myArray = new T[myCapacity];
}
//-------------------------------------------------
template<typename T> inline T const* csCollection<T>::toArray() const {
  return myArray;
}
//-------------------------------------------------
template<typename T> inline void csCollection<T>::clear() {
  mySize = 0;
}
//-------------------------------------------------
template<typename T> inline void csCollection<T>::dispose() {
  clear();
  reallocate( 1 );
}
//-------------------------------------------------
template<typename T> inline int csCollection<T>::capacity() const {
  return myCapacity;
}
//-------------------------------------------------
template<typename T> inline int csCollection<T>::size() const {
  return mySize;
}
//-------------------------------------------------
template<typename T> inline bool csCollection<T>::isEmpty() const {
  return( mySize == 0 );
}
//-------------------------------------------------
/*template<typename T> inline bool csCollection<T>::contains( T const& value ) const {
  for( int i = 0; i < mySize; i++ ) {
    if( myArray[i] == value ) return true;
  }
  return false;
}*/
//-------------------------------------------------
/*
template<typename T> void csCollection<T>::checkCapacityReallocate() {
  if( mySize == myCapacity ) {
    if( myNumReallocateCalls == 2 ) {
      myChunkSize *= 2;
      myNumReallocateCalls = 0;
    }
    reallocate( myCapacity + myChunkSize );
  }
}
*/
//-------------------------------------------------
template<typename T> void csCollection<T>::reallocate() {
  if( myNumReallocateCalls == 2 ) {
    myChunkSize *= 2;
    myNumReallocateCalls = 0;
  }
  reallocate( myCapacity + myChunkSize );
}
template<typename T> void csCollection<T>::reallocate( int capacity ) {
  myCapacity = capacity;
  T* array_new = new T[myCapacity];
  int save = std::min( myCapacity, mySize );
  for( int i = 0; i < save; i++ ) {
    array_new[i] = myArray[i];   // Careful: May cause problems if T is a pointer type..?
  }
  if( myArray ) {
    delete [] myArray;   // Careful: May cause problems if T is a pointer type..?
  }
  myArray = array_new;
  myNumReallocateCalls++;
}
//-------------------------------------------------
template<typename T> void csCollection<T>::dump() const {
//  for( int i = 0; i < mySize; i++ ) {
//    cout << "Value " << i << ": " << myArray[i] << endl;
//  }
}
//-------------------------------------------------
template<typename T> inline void csCollection<T>::reduceSizeTo( int nElements ) {
  if( nElements < 0 ) {
    mySize = 0;
  }
  else if( nElements < mySize ) {
    mySize = nElements;
  }
}

} // namespace

#endif

