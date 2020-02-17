/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SORT_H
#define CS_SORT_H

#include <cstdlib>

namespace cseis_geolib {

/**
 * Sort class
 *
 * Provides two sort algorithms
 * Types that can be used with this class must implement operators < and = 
 *
 * @author Bjorn Olofsson
 * @date 2006
 */
template <typename T>
class csSort {

 private:
  template <typename W>
    struct treeNode_t {
      W value;
      struct treeNode_t* smaller;
      struct treeNode_t* greater;
      struct treeNode_t* up;
    };

 public:
  csSort();
  ~csSort();

  /**
  * Simple sort algorithm
  * Use this sort routine for small arrays or arrays which are already partially sorted
  *
  * @param values (i/o) Array of values/objects that shall be sorted
  * @param nValues (i)  Number of elements in array
  */
  void simpleSort( T* values, int nValues );

  /**
  * Simple sort algorithm
  * Sorts index array as well, in order to keep track of change in sort order
  *
  * @param values (i/o) Array of values/objects that shall be sorted
  * @param nValues (i)  Number of elements in array
  * @param index  (i/o) Integer array. This will be sorted in the same way as the 'value' array
  *                     For example, this array can be set to the numbers 1..numValues beforehand in order to
  *                     obtain the sort order
  */
  void simpleSortIndex( T* values, int nValues, int* index );

  /**
  * Tree sort algorithm
  * Use this sort routine for big arrays that are randomly sorted
  *
  * @param values (i/o) Array of values/objects that shall be sorted
  * @param nValues (i)  Number of elements in array
  */
  void treeSort( T* values, int nValues );

 private:
  /// Tree nodes. Used for tree sort
  treeNode_t<T>* myTreeNodes;
  /// Number of currently allocated tree nodes
  int myNumAllocatedNodes;
};

/**
* Sort item
* Consists of sort item/value/object, and an index
* Can be used with sort algorithms, in order to keep track of item indices after sorting
*/
template <typename T> class csSortItem {
 public:
  csSortItem();
  csSortItem( T item, int index );
  csSortItem( csSortItem const& sortItem );
  ~csSortItem();

  csSortItem& operator=( const csSortItem& obj );
  inline bool operator==( csSortItem const& obj ) const;
  inline bool operator!=( csSortItem const& obj ) const;
  inline bool operator>=( csSortItem const& obj ) const;
  inline bool operator<=( csSortItem const& obj ) const;
  inline bool operator>( csSortItem const& obj ) const;
  inline bool operator<( csSortItem const& obj ) const;

  inline void set( T const& value, int index );
  inline T value() const;
  inline int index() const;
public:
  T myItem;
  int myIndex;
};

//================================================================================
//
// csSort implementation
//
//

 template<typename T>csSort<T>::csSort() {
   myNumAllocatedNodes = 0;
   myTreeNodes = NULL;
 }
 //--------------------------------------------------------------------------------
 template<typename T>csSort<T>::~csSort() {
   if( myTreeNodes != NULL ) {
     delete [] myTreeNodes;
     myTreeNodes = NULL;
   }
 }
 //--------------------------------------------------------------------------------
 template<typename T> void csSort<T>::simpleSort( T* values, int nValues ) {
   T tmp;
   bool hasMoved = false;
   do {
     hasMoved = false;
     for( int i = 1; i < nValues; i++ ) {
       if( values[i] < values[i-1] ) {
         hasMoved = true;
         tmp = values[i-1];
         values[i-1] = values[i];
         values[i] = tmp;
       }
     }
   } while( hasMoved );
 }
 //--------------------------------------------------------------------------------
 template<typename T> void csSort<T>::simpleSortIndex( T* values, int nValues, int* index ) {
   T tmp;
   int tmpIndex;
   bool hasMoved = false;
   do {
     hasMoved = false;
     for( int i = 1; i < nValues; i++ ) {
       if( values[i] < values[i-1] ) {
         hasMoved = true;
         tmp = values[i-1];
         values[i-1] = values[i];
         values[i] = tmp;
         
         tmpIndex   = index[i-1];
         index[i-1] = index[i];
         index[i]   = tmpIndex;
       }
     }
   } while( hasMoved );
 }
 //--------------------------------------------------------------------------------
 template<typename T> void csSort<T>::treeSort( T* values, int nValues ) {
   if( myNumAllocatedNodes < nValues ) {
     if( myTreeNodes != NULL ) {
       delete [] myTreeNodes;
     }
     myTreeNodes = new treeNode_t<T>[nValues];  // Create one tree node per value
     myNumAllocatedNodes = nValues;
   }
   treeNode_t<T>* ptr;
   treeNode_t<T>* ptrNew;
   treeNode_t<T>* firstNodePtr;
   // Initialize empty tree nodes
   for( int i = 0; i < nValues; i++ ) {
     myTreeNodes[i].value = values[i];
     myTreeNodes[i].smaller = NULL;
     myTreeNodes[i].greater = NULL;
     myTreeNodes[i].up      = NULL;
   }
   firstNodePtr = &myTreeNodes[0];
   
   // Insert all values into tree, at correct place. This is the loop that takes most time
   for( int i = 1; i < nValues; i++ ) {
     ptr    = firstNodePtr;
     ptrNew = &myTreeNodes[i];
     while( true ) {
       if( ptrNew->value < ptr->value ) {
         if( ptr->smaller != NULL ) {
           ptr = ptr->smaller;
         }
         else {  // Add smaller value
           ptr->smaller = ptrNew;
           ptrNew->up = ptr;
           break;
         }
       }
       else {
         if( ptr->greater != NULL ) {
           ptr = ptr->greater;
         }
         else {  // Add greater value
           ptr->greater = ptrNew;
           ptrNew->up = ptr;
           break;
         }
       }
     }
   }
   
   // Output sorted array
   ptr = firstNodePtr;
   int counter = 0;
   while( ptr != NULL ) {
     while( ptr->smaller != NULL ) {
       ptr = ptr->smaller;
     }
     values[counter++] = ptr->value;
     if( ptr->greater != NULL ) {
       ptr = ptr->greater;
     }
     else {
       do {
         while( ptr->up != NULL && ptr->up->smaller != ptr ) {
           ptr = ptr->up;
         }
         ptr = ptr->up;
         if( ptr == NULL ) break;
         values[counter++] = ptr->value;
       } while( ptr->greater == NULL );
       if( ptr != NULL ) ptr = ptr->greater;
     }
   }
 }

//================================================================================
//================================================================================
// csSortItem implementation
//


template<typename T>csSortItem<T>::csSortItem() : myItem(0), myIndex(0) {
}
template<typename T>csSortItem<T>::csSortItem( T item, int index ) {
  myItem  = item;
  myIndex = index;
}
template<typename T>csSortItem<T>::csSortItem( csSortItem const& sortItem ) {
  myItem  = sortItem.myItem;
  myIndex = sortItem.myIndex;
}
template<typename T>csSortItem<T>::~csSortItem() {
}

template<typename T> inline void csSortItem<T>::set( T const& value, int index ) {
  myItem  = value;
  myIndex = index;
}
template<typename T> inline T csSortItem<T>::value() const {
  return myItem;
}
template<typename T> inline int csSortItem<T>::index() const {
  return myIndex;
}

template<typename T> csSortItem<T>& csSortItem<T>::operator=( const csSortItem& obj ) {
  myItem  = obj.myItem;
  myIndex = obj.myIndex;
  return *this;
}
template<typename T> inline bool csSortItem<T>::operator==( const csSortItem& obj ) const {
  return( myItem == obj.myItem );
}
template<typename T> inline bool csSortItem<T>::operator!=( const csSortItem& obj ) const {
  return( myItem != obj.myItem );
}
template<typename T> inline bool csSortItem<T>::operator>=( const csSortItem& obj ) const {
  return( myItem >= obj.myItem );
}
template<typename T> inline bool csSortItem<T>::operator<=( const csSortItem& obj ) const {
  return( myItem <= obj.myItem );
}
template<typename T> inline bool csSortItem<T>::operator>( const csSortItem& obj ) const {
  return( myItem > obj.myItem );
}
template<typename T> inline bool csSortItem<T>::operator<( const csSortItem& obj ) const {
  return( myItem < obj.myItem );
}


} // end namespace
#endif

