/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TABLE_MANAGER_H
#define CS_TABLE_MANAGER_H

#include <string>
#include "csTable.h"

namespace cseis_geolib {
  // Forward declarations:
  class csFlexNumber;
  template <typename T> class csVector;
  template <typename T> class csTimeFunction;
}

namespace cseis_system {

  class csTraceHeaderDef;
  class csTraceHeader;

/**
* Table manager
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTableManager {
public:
  csTableManager( cseis_geolib::csTable const* table, cseis_system::csTraceHeaderDef const* hdef );
  ~csTableManager();
  void set( cseis_geolib::csTable const* table, cseis_system::csTraceHeaderDef const* hdef );
//  void set( cseis_geolib::csVector<std::string> const* headerList, std::string const* selectionText, cseis_system::csTraceHeaderDef const* hdef );
  /// @return true   if the header values passed to this method are equal to the trace header value selections
//  bool contains( cseis_system::csTraceHeader const* trcHeader );
//  void dump();
  int headerIndex( int keyIndex ) const;

  cseis_geolib::csTimeFunction<double> const* getFunction( csTraceHeader const* trcHdr );
  double getValue( csTraceHeader const* trcHdr, int valueIndex = 0 );
  double getValue( csTraceHeader const* trcHdr, double time );

  // Simple methods -> just redirect calls to table object
  inline int numLocations() const { return myTable->numLocations(); }
  inline int numKeys() const { return myTable->numKeys(); }
  inline std::string tableName() const { return myTable->tableName(); }
  inline std::string const keyName( int indexKey ) const { return myTable->keyName(indexKey); }
  inline std::string const valueName( int indexValue = 0 ) const { return myTable->valueName(indexValue); }
  inline int dimension() const { return myTable->dimension(); }

  inline cseis_geolib::csTable const* table() const { return myTable; }

private:
  void setKeyValueBuffer( csTraceHeader const* trcHdr );
  cseis_geolib::csTable const* myTable;

  int* myHeaderIndex;
  double* myKeyValueBuffer;
};

} // namespace

#endif

