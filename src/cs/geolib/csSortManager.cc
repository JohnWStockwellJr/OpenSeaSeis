/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>

#include "csSortManager.h"
#include "csException.h"
#include "csFlexNumber.h"
#include "csSort.h"

using namespace cseis_geolib;

SortArrayItem::SortArrayItem() {
  myValues       = NULL;
  myNumSortKeys  = 0;
  myAllocatedPtr = NULL;
}
SortArrayItem::SortArrayItem( int numSortKeys ) {
  myValues       = NULL;
  myAllocatedPtr = NULL;
  setNumSortKeys( numSortKeys );
}
SortArrayItem::SortArrayItem( SortArrayItem const& obj ) {
  myValues       = obj.myValues;
  myNumSortKeys  = obj.myNumSortKeys;
  myCompareIndex = obj.myCompareIndex;
  myAllocatedPtr = NULL;
}
void SortArrayItem::setNumSortKeys( int numSortKeys ) {
  myNumSortKeys   = numSortKeys;
  myValues       = new csFlexNumber[myNumSortKeys];
  if( myAllocatedPtr != NULL ) {
    delete [] myAllocatedPtr;
  }
  myAllocatedPtr = myValues;
}
SortArrayItem::~SortArrayItem() {
  if( myAllocatedPtr != NULL ) {
    delete [] myAllocatedPtr;
    myAllocatedPtr = NULL;
  }
}
SortArrayItem& SortArrayItem::operator=( SortArrayItem const& obj ) {
  myNumSortKeys  = obj.myNumSortKeys;
  myValues       = obj.myValues;
  myCompareIndex = obj.myCompareIndex;
  return *this;
}

void SortArrayItem::set( int keyIndex, csFlexNumber const& value ) {
  myValues[keyIndex] = value;
}
void SortArrayItem::setIndexToSort( int indexToSort ) {
  myValues[myNumSortKeys-1] = indexToSort;
}

void SortArrayItem::setCompareIndex( int keyIndex ) {
  myCompareIndex = keyIndex;
}
bool SortArrayItem::operator<( const SortArrayItem& obj ) const {
  return( myValues[myCompareIndex] < obj.myValues[myCompareIndex] );
}
void SortArrayItem::dump() const {
  for( int ikey = 0; ikey < myNumSortKeys; ikey++ ) {
    if( myValues[ikey].type() == TYPE_INT ) {
      fprintf(stdout,"%6d ", myValues[ikey].intValue() );
    }
    else if( myValues[ikey].type() == TYPE_INT64 ) {
      fprintf(stdout,"%14lld ", myValues[ikey].int64Value() );
    }
    else {
      fprintf(stdout,"%6f ", myValues[ikey].doubleValue() );
    }
  }
  fprintf(stdout,"\n");
}
cseis_geolib::csFlexNumber const* SortArrayItem::get( int keyIndex ) const {
  return &myValues[keyIndex];
}
void SortArrayItem::flipPolarity( int keyIndex ) {
  myValues[keyIndex] = myValues[keyIndex] * -1;
}

//================================================================================
//================================================================================
//
//
csSortManager::csSortManager( int numSortKeys, int sortMethod ) {
  myNumValues   = 0;
  mySort        = new csSort<SortArrayItem>();
  mySortItems   = NULL;
  myNumSortKeys  = numSortKeys;
  mySortMethod    = sortMethod;
}
csSortManager::~csSortManager() {
  if( mySortItems != NULL ) {
    delete [] mySortItems;
    mySortItems = NULL;
  }
  if( mySort != NULL ) {
    delete mySort;
    mySort = NULL;
  }
}
int csSortManager::numValues() const {
  return myNumValues;
}
void csSortManager::resetValues( int numValues ) {
  if( numValues > myNumValues ) {
    if( mySortItems != NULL ) {
      delete [] mySortItems;
    }
    mySortItems = new SortArrayItem[numValues];
  }
  myNumValues = numValues;
  myNumCurrentValues = 0;
  for( int ival = 0; ival < myNumValues; ival++ ) {
    mySortItems[ival].setNumSortKeys( myNumSortKeys+1 );  // +1 to add index field to end of key list
    mySortItems[ival].set( myNumSortKeys, ival ); // Preset index field to consecutive value index
  }
}
void csSortManager::setValue( int valueIndex, int keyIndex, csFlexNumber const& value, int indexToSort ) {
  if( valueIndex >= myNumValues ) {
    throw( csException("csSortManager::setValue: Program bug found in calling function. Too many values added. Call resetValues() first") );
  }
  mySortItems[valueIndex].set( keyIndex, value );
  if( indexToSort >= 0 ) {
    mySortItems[valueIndex].setIndexToSort( indexToSort );
  }
}
void csSortManager::sort() {
  for( int ihdr = 0; ihdr < myNumSortKeys; ihdr++ ) {
    for( int ival = 0; ival < myNumValues; ival++ ) {
      mySortItems[ival].setCompareIndex( ihdr );
    }
    if( mySortMethod == csSortManager::TREE_SORT ) {
      mySort->treeSort( mySortItems, myNumValues );
    }
    else {
      mySort->simpleSort( mySortItems, myNumValues );
    }
  }
}
int csSortManager::sortedIndex( int index ) const {
  return mySortItems[index].get( myNumSortKeys )->intValue();
}
cseis_geolib::csFlexNumber const* csSortManager::sortedValue( int valIndex, int keyIndex ) const {
  return mySortItems[valIndex].get( keyIndex );
}
double csSortManager::sortedDoubleValue( int valIndex, int keyIndex ) const {
  return sortedValue( valIndex, keyIndex )->doubleValue();
}
float csSortManager::sortedFloatValue( int valIndex, int keyIndex ) const {
  return sortedValue( valIndex, keyIndex )->floatValue();
}
csInt64_t csSortManager::sortedInt64Value( int valIndex, int keyIndex ) const {
  return sortedValue( valIndex, keyIndex )->int64Value();
}
int csSortManager::sortedIntValue( int valIndex, int keyIndex ) const {
  return sortedValue( valIndex, keyIndex )->intValue();
}
void csSortManager::flipPolarity( int keyIndex ) {
  for( int ival = 0; ival < myNumValues; ival++ ) {
    mySortItems[ival].flipPolarity( keyIndex );
  }
}
void csSortManager::dump() const {
  for( int ival = 0; ival < myNumValues; ival++ ) {
    fprintf(stdout,"Value #%2d: ", ival);
    mySortItems[ival].dump();
  }
}
