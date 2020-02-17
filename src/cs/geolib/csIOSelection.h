/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_IO_SELECTION_H
#define CS_IO_SELECTION_H

#include <string>
#include <vector>
#include <algorithm>
#include "geolib_defines.h"

namespace cseis_geolib {
  class csIReader;
  class csFlexHeader;
  class csFlexNumber;
  template <typename T> class csVector;
  class csSortManager;
  class csSelection;

class csIOSortItem {
 public:
  csIOSortItem( int traceIndex_in, csFlexNumber* flexNum_in ) {
    traceIndex = traceIndex_in;
    flexNum    = flexNum_in;
  }
  ~csIOSortItem() {
  }
  csFlexNumber const* flexNum;
  int traceIndex;
 private:
};

/**
 * IOSelection.
 *
 * Helper class to select certain traces from input file.
 */
class csIOSelection {
 public:
  static const int SORT_NONE = 0;
  static const int SORT_INCREASING = 1;
  static const int SORT_DECREASING = 2;

 public:
  /**
   * Constructor
   *
   * @param headerName: Header name to select/sort on
   * @param sortOrder:  csIOSelection::SORT_NONE, SORT_INCREASING, or SORT_DECREASING
   * @param sortMethod: csSortManager::SIMPLE_SORT or csSortManager::TREE_SORT
   */
  csIOSelection( std::string const& headerName, int sortOrder, int sortMethod );
  ~csIOSelection();

  /**
   * Initialize header selection
   * Note that this method may take a long time to complete, since it will determine pointers to all selected traces in the input file, and perform the necessary sort
   *
   * @param reader The data reader object
   * @param hdrValueSelectionText  A string defining the selection of certain trace header values. Selection syntax is explained on www.seaseis.com
   */
  bool initialize( cseis_geolib::csIReader* reader, std::string const& hdrValueSelectionText );
  /**
   * step1, step2, step3: Same as initialize but in three separate steps
   * Purpose: Be able to scan partially through reader
   */
  void step1( cseis_geolib::csIReader* reader, std::string const& hdrValueSelectionText );
  bool step2( cseis_geolib::csIReader* reader, int& numTracesToRead );
  bool step3( cseis_geolib::csIReader* reader );
  /**
   * Call this method to reset header selection.
   * After the call, method getNextTraceIndex() will return the first selected trace
   */
  void resetHeaderSelection();
  /**
   * Call consecutively to retrieve all selected traces
   * @return trace index of next selected (and potentially sorted) trace
   */
  int getNextTraceIndex();
  /**
   * @return trace index of currently selected trace (by previously calling getNextTrace())
   */
  int getCurrentTraceIndex();
  /**
   * @return Number of selected traces
   */
  int getNumSelectedTraces() const;
  /**
   * @param traceIndex  Index of consecutively selected trace
   * @return Value of selected trace
   */
  cseis_geolib::csFlexNumber const* getSelectedValue( int traceIndex ) const;
  int getSelectedIndex( int traceIndex ) const;

 private:
  void flipSortOrder();

  std::string myHdrName;
  int myNumSelectedTraces;
  int mySortOrder;
  int myCurrentSelectedIndex;
  cseis_geolib::csVector<int>* mySelectedTraceIndexList;
  cseis_geolib::csVector<csFlexNumber*>* mySelectedValueList;
  std::vector<csIOSortItem*>* mySortList;

  int mySTEPCurrentTraceIndex;
  int mySTEPTraceIndex;
  cseis_geolib::csSelection* mySTEPSelection;
};

} // end namespace

#endif
