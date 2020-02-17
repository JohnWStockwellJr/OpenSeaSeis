/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstring>
#include "csTableManagerNew.h"
#include "csTableAll.h"
#include "csTraceHeaderDef.h"
#include "csTraceHeader.h"
#include "csException.h"

using namespace cseis_system;

csTableManagerNew::csTableManagerNew( std::string table_filename, int table_type, cseis_system::csTraceHeaderDef const* hdef )
{
  myTable = new cseis_geolib::csTableAll( table_type );
  myTable->initialize( table_filename );

  int nKeys = myTable->numKeys();
  myHeaderIndex    = new int[nKeys];
  myKeyValueBuffer = new double[nKeys];

  for( int ikey = 0; ikey < nKeys; ikey++ ) {
    if( !hdef->headerExists( myTable->keyName(ikey) ) ) {
      throw( cseis_geolib::csException("No matching trace header found for velocity table key '%s'", myTable->keyName(ikey).c_str()) );
    }
    myHeaderIndex[ikey] = hdef->headerIndex( myTable->keyName(ikey) );
  }

  myTable->readTableContents();
}
csTableManagerNew::~csTableManagerNew() {
  if( myKeyValueBuffer != NULL ) {
    delete [] myKeyValueBuffer;
    myKeyValueBuffer = NULL;
  }
  if( myHeaderIndex != NULL ) {
    delete [] myHeaderIndex;
    myHeaderIndex = NULL;
  }
  if( myTable != NULL ) {
    delete myTable;
    myTable = NULL;
  }
}
void csTableManagerNew::dump() {
  myTable->dump();
}
//-------------------------------------------------------------------------------
int csTableManagerNew::headerIndex( int keyIndex ) const {
  if( keyIndex >= 0 && keyIndex < numKeys() ) {
    return myHeaderIndex[keyIndex];
  }
  else {
    throw( cseis_geolib::csException("csTableManagerNew::headerIndex(): Wrong key index passed. This is probably due to a program bug in the calling method.") );
  }
}
//-------------------------------------------------------------------------------
cseis_geolib::csTimeFunction<double> const* csTableManagerNew::getFunction( csTraceHeader const* trcHdr ) {
  setKeyValueBuffer( trcHdr );
  return myTable->getFunction( myKeyValueBuffer );
}
//-------------------------------------------------------------------------------
double csTableManagerNew::getValue( csTraceHeader const* trcHdr, int valueIndex ) {
  setKeyValueBuffer( trcHdr );
  return myTable->getValue( myKeyValueBuffer, valueIndex );
}
//-------------------------------------------------------------------------------
void csTableManagerNew::setKeyValueBuffer( csTraceHeader const* trcHdr ) {
  for( int ikey = 0; ikey < numKeys(); ikey++ ) {
    myKeyValueBuffer[ikey] = trcHdr->doubleValue( myHeaderIndex[ikey] );
    //    printf("Table key %d %d: %f\n", ikey, myHeaderIndex[ikey], myKeyValueBuffer[ikey] );
  }
}
int csTableManagerNew::numLocations() const {
  return myTable->numLocations();
}
int csTableManagerNew::numKeys() const { return myTable->numKeys(); }
std::string csTableManagerNew::tableName() const { return myTable->tableName(); }
std::string const csTableManagerNew::keyName( int indexKey ) const { return myTable->keyName(indexKey); }
std::string const csTableManagerNew::valueName( int indexValue ) const { return myTable->valueName(indexValue); }
int csTableManagerNew::type() const { return myTable->type(); }

