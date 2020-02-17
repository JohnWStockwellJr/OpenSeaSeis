/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csInitExecEnv.h"
#include "csTable.h"
#include "csException.h"

using namespace cseis_system;

csExecPhaseEnv::csExecPhaseEnv(
    csTraceHeaderDef const* hDef,
    csExecPhaseDef*   eDef,
    csSuperHeader const*    sHdr ) {
  headerDef = hDef;
  execPhaseDef = eDef;
  superHeader = sHdr;
  errorCounter = 0;
}
csExecPhaseEnv::~csExecPhaseEnv() {
}

csInitPhaseEnv::csInitPhaseEnv(
    csTraceHeaderDef* hDef,
    csExecPhaseDef*   eDef,
    csSuperHeader*    sHdr,
    cseis_geolib::csTable const** tables,
    int numTables ) {
  headerDef = hDef;
  execPhaseDef = eDef;
  superHeader  = sHdr;
  errorCounter = 0;
  myNumTables  = numTables;
  myTables     = tables;
}
csInitPhaseEnv::~csInitPhaseEnv() {
}
cseis_geolib::csTable const* csInitPhaseEnv::getTable( std::string const& tableName ) const {
  try {
    for( int i = 0; i < myNumTables; i++ ) {
      if( !tableName.compare( myTables[i]->tableName()) ) {
        return myTables[i];
      }
    }
  }
  catch( cseis_geolib::csException& exc ) {
    printf("Error occurred: %s\n", exc.getMessage() );
    return NULL;
  }
  return NULL;
}


