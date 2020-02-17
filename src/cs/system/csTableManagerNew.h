/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TABLE_MANAGER_NEW_H
#define CS_TABLE_MANAGER_NEW_H

#include <string>

namespace cseis_geolib {
  // Forward declarations:
  class csFlexNumber;
  class csTableAll;
  template <typename T> class csVector;
  template <typename T> class csTimeFunction;
}

namespace cseis_system {

  class csTraceHeaderDef;
  class csTraceHeader;

/**
* Table manager.
* This class superseeds the older version of this class, i.e. csTableManager
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTableManagerNew {
public:
  csTableManagerNew( std::string table_filename, int table_type, cseis_system::csTraceHeaderDef const* hdef );
  ~csTableManagerNew();
  //  void set( cseis_geolib::csTable const* table, cseis_system::csTraceHeaderDef const* hdef );
  void dump();
  int headerIndex( int keyIndex ) const;

  cseis_geolib::csTimeFunction<double> const* getFunction( csTraceHeader const* trcHdr );
  double getValue( csTraceHeader const* trcHdr, int valueIndex = 0 );

  // Simple methods -> just redirect calls to table object
  int numLocations() const;
  int numKeys() const;
  std::string tableName() const;
  std::string const keyName( int indexKey ) const;
  std::string const valueName( int indexValue = 0 ) const;
  int type() const;
  
  inline cseis_geolib::csTableAll const* table() const { return myTable; }

private:
  void setKeyValueBuffer( csTraceHeader const* trcHdr );
  cseis_geolib::csTableAll* myTable;

  int* myHeaderIndex;
  double* myKeyValueBuffer;
};

} // namespace

#endif

