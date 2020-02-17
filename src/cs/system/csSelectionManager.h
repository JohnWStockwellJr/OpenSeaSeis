/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SELECTION_MANAGER_H
#define CS_SELECTION_MANAGER_H

#include <string>
#include "geolib_defines.h"

namespace cseis_geolib {
  // Forward declarations:
  class csFlexNumber;
  class csFlexHeader;
  class csSelection;
  template <typename T> class csVector;
}

namespace cseis_system {

class csTraceHeaderDef;
class csTraceHeader;

/**
 * Manages all aspects of a user header 'selection'
 *
 * Header value selections occur in many Cseis modules, for example in IF and ELSEIF statements.
 * A value selection consists of a list of trace header names, and a selection field according to the format
 * described in @class csSelection.
 *
 * This class serves as a helper class between a module and the selection field, with an easy interface.
 *
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
class csSelectionManager {
public:
  /**
   * Constructor
   */
  csSelectionManager();
  ~csSelectionManager();
  /**
   * Prepares the selection field
   * @param headerList     List of trace header names specified in input flow
   * @param selectionText  Selection text specified in input flow
   * @param hdef           Trace header definition object
   * @throws csException
   */
  void set( cseis_geolib::csVector<std::string> const* headerList, std::string const* selectionText,
            cseis_system::csTraceHeaderDef const* hdef );
  /**
   * Is the current trace contained in the selection?
   *
   * @param trcHeader  Trace header object
   * @return true if the header values passed to this method are contained in the trace header value selections
   */
  bool contains( cseis_system::csTraceHeader const* trcHeader );
  /**
   * Is the current trace contained in the selection?
   * All trace headers that this selection is based upon are checked against this value.
   * Use this method with caution.
   *
   * @param hdrValue  Trace header value
   * @return true if the header value passed to this method is contained in the trace header value selections
   */
  bool contains( cseis_geolib::csFlexHeader const* hdrValue );
  /**
   * @param index Index of trace header to return (>0 in case selection is based on more than one trace header)
   * @return Trace header name
   */
  std::string headerName( int index ) const;

  void dump();
private:
  csSelectionManager( csSelectionManager const& obj );
  /// Number of trace headers in this selection
  int  myNumHeaders;
  /// Trace header indexes 
  int* myHeaderIndex;
  /// Trace header formats (TYPE_INT, TYPE_FLOAT, TYPE_DOUBLE)
  cseis_geolib::type_t* myHeaderType;
  /// Trace header names
  std::string* myHeaderNames;
  /// The value selection field object that this instance manages
  cseis_geolib::csSelection* mySelection;
  /// Temporary field used to store the trace header values of the current trace
  cseis_geolib::csFlexNumber* myValues;
  ///
  std::string myStringSelectionText;
};

} // namespace

#endif

