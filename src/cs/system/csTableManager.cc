/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstring>
#include "csTableManager.h"
#include "csTable.h"
#include "csTraceHeaderDef.h"
#include "csTraceHeader.h"
#include "csException.h"

using namespace cseis_system;

csTableManager::csTableManager( cseis_geolib::csTable const* table, cseis_system::csTraceHeaderDef const* hdef )
{
  myTable = table;

  int nKeys = myTable->numKeys();
  myHeaderIndex    = new int[nKeys];
  myKeyValueBuffer = new double[nKeys];

  for( int ikey = 0; ikey < nKeys; ikey++ ) {
    if( !hdef->headerExists( myTable->keyName(ikey) ) ) {
      throw( cseis_geolib::csException("No matching trace header found for velocity table key '%s'", myTable->keyName(ikey).c_str()) );
    }
    myHeaderIndex[ikey] = hdef->headerIndex( myTable->keyName(ikey) );
  }
}

csTableManager::~csTableManager() {
  if( myKeyValueBuffer != NULL ) {
    delete [] myKeyValueBuffer;
    myKeyValueBuffer = NULL;
  }
  if( myHeaderIndex != NULL ) {
    delete [] myHeaderIndex;
    myHeaderIndex = NULL;
  }
}
//-------------------------------------------------------------------------------
int csTableManager::headerIndex( int keyIndex ) const {
  if( keyIndex >= 0 && keyIndex < numKeys() ) {
    return myHeaderIndex[keyIndex];
  }
  else {
    throw( cseis_geolib::csException("csTableManager::headerIndex(): Wrong key index passed. This is probably due to a program bug in the calling method.") );
  }
}
//-------------------------------------------------------------------------------
cseis_geolib::csTimeFunction<double> const* csTableManager::getFunction( csTraceHeader const* trcHdr ) {
  setKeyValueBuffer( trcHdr );
  return myTable->getFunction( myKeyValueBuffer );
}
//-------------------------------------------------------------------------------
double csTableManager::getValue( csTraceHeader const* trcHdr, int valueIndex ) {
  setKeyValueBuffer( trcHdr );
  return myTable->getValue( myKeyValueBuffer, valueIndex );
}
//-------------------------------------------------------------------------------
double csTableManager::getValue( csTraceHeader const* trcHdr, double time ) {
  setKeyValueBuffer( trcHdr );
  return myTable->getValue( myKeyValueBuffer, time );
}
//-------------------------------------------------------------------------------
void csTableManager::setKeyValueBuffer( csTraceHeader const* trcHdr ) {
  for( int ikey = 0; ikey < numKeys(); ikey++ ) {
    myKeyValueBuffer[ikey] = trcHdr->doubleValue( myHeaderIndex[ikey] );
//    printf("Table key %d %d: %f\n", ikey, myHeaderIndex[ikey], myKeyValueBuffer[ikey] );
  }
}

