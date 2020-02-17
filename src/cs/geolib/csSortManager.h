/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SORT_MANAGER_H
#define CS_SORT_MANAGER_H

#include "geolib_defines.h"


namespace cseis_geolib {

  template<typename T> class csSort;
  class csFlexNumber;

/**
 * Sort array item.
 *
 * This class is designed for use in conjunction with class csSort.
 * Use SortArrayItem to sort a number of values by a one or several (sort) keys. If the data to be sorted is
 * a 2D table with rows and columns, place all values for one row into one object of type
 * SortArrayItem, and all SortArrayItem objects into one array.
 * Before sorting the array, set the 'compare index' to the sort key index that shall be sorted,
 * then sort the array. Repeat sorting for all keys in case more than one shall be sorted.
 *
 * @author Bjorn Olofsson
 * @date 2007
 */
class SortArrayItem {
 public:
  SortArrayItem();
  /**
   * Constructor
   * @param numSortKeys  Sort data by 'numSortKeys' number of keys.
   */
  SortArrayItem( int numSortKeys );
  /**
   * Copy constructor
   */
  SortArrayItem( SortArrayItem const& obj );
  ~SortArrayItem();
  /**
   * Set number of keys to be soted on
   * @param numSortKeys  Sort data by 'numSortKeys' number of keys
   */
  void setNumSortKeys( int numSortKeys );
  /**
   * Set value for one key
   * @param keyIndex  Index of sort key
   * @param value     Value of sort key
   */
  void set( int keyIndex, cseis_geolib::csFlexNumber const& value );
  void setIndexToSort( int indexToSort );
  /**
   * Set index of key that shall be sorted next
   * @param keyIndex  Index of sort key
   */
  void setCompareIndex( int keyIndex );
  /**
   * Retrieve value
   * @param keyIndex  Sort key index
   */
  cseis_geolib::csFlexNumber const* get( int keyIndex ) const;
  /**
   * Dump all values. Use for debugging
   */
  void dump() const;
  bool operator<( const SortArrayItem& obj ) const;
  SortArrayItem& operator=( SortArrayItem const& obj );
  void flipPolarity( int keyIndex );
private:
  /// Index of sort key that shall be sorted this time.
  int myCompareIndex;
  /// Number of sort keys to be sorted by
  int myNumSortKeys;
  /// Values that shall be sorted, one for each sort key. Do not free memory pointed to by this array pointer.
  cseis_geolib::csFlexNumber*  myValues;
  /// Pointer to allocated value array. Free memory using this pointer.
  cseis_geolib::csFlexNumber*  myAllocatedPtr;
};

class csSortManager {
public:
  static int const SIMPLE_SORT = 101;
  static int const TREE_SORT   = 201;
public:
  /**
   * @param numSortKeys  Data shall be sorted by this number of keys.
   * @param sortMethod   Sort method: csSortManager::SIMPLE_SORT or csSortManager::TREE_SORT
   */
  csSortManager( int numSortKeys, int sortMethod = csSortManager::SIMPLE_SORT );
  ~csSortManager();
  /**
   * Reset number of values to be sorted next
   * @param numValues  Number of values to be sorted
   */
  void resetValues( int numValues );
  /**
   * Set one value to be sorted.
   * Call consecutively to set all values which shall be sorted.
   * @param valueIndex  Index of value to be set. Must be smaller than numValues set previously by 'resetValues(numValues)'.
   * @param keyIndex    Index of sort key that this value belongs to.
   * @param value       The value
   * @param indexToSort Additional index field which shall be sorted at the same time. -1: Do not use index field
   */
  void setValue( int valueIndex, int keyIndex, cseis_geolib::csFlexNumber const& value, int indexToSort = -1 );
  /**
   * Perform the sort.
   */
  void sort();
  /**
   * @param index 
   * @return
   */
  int sortedIndex( int index ) const;
  double sortedDoubleValue( int valIndex, int keyIndex = 0 ) const;
  float sortedFloatValue( int valIndex, int keyIndex = 0 ) const;
  csInt64_t sortedInt64Value( int valIndex, int keyIndex = 0 ) const;
  int sortedIntValue( int valIndex, int keyIndex = 0 ) const;
  void dump() const;
  int numValues() const;
  cseis_geolib::csFlexNumber const* sortedValue( int valIndex, int keyIndex = 0 ) const;
  void flipPolarity( int keyIndex = 0 );
private:

  /// Number of keys to be sorted by.
  int myNumSortKeys;
  /// Number of values that shall be sorted
  int myNumValues;
  /// Current number of values that shall be sorted
  int myNumCurrentValues;
  cseis_geolib::csSort<SortArrayItem>* mySort;
  SortArrayItem* mySortItems;
  /// Sort method
  int mySortMethod;
};

} // end namespace
#endif

