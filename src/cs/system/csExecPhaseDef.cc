/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csExecPhaseDef.h"
#include "csTraceGather.h"
#include "csTraceHeaderDef.h"
#include "csSuperHeader.h"
#include "cseis_defines.h"
#include "csFlexNumber.h"
#include "csException.h"

using namespace cseis_system;

csExecPhaseDef::csExecPhaseDef( std::string& moduleName, int mpiProcID ) {
  myModuleName = moduleName;
  myMPIProcID  = mpiProcID;
  numTraces  = 1;
  traceMode  = TRCMODE_UNKNOWN;
  myExecType   = EXEC_TYPE_NORMAL;
  varPtr       = NULL;
  myIsCleanup  = false;
  myIsDebug    = false;
  myTracesAreWaiting = false;
  myIsLastCall  = false;
  myIsMPISupported = false;
}
csExecPhaseDef::~csExecPhaseDef() {
}
void csExecPhaseDef::setTraceSelectionMode( int mode ) {
  traceMode = mode;
}
void csExecPhaseDef::setTraceSelectionMode( int mode, int nTraces ) {
  traceMode = mode;
  numTraces = nTraces;
  if( nTraces <= 0 ) {
    throw( cseis_geolib::csException("csExecPhaseDef::setTraceSelectionMode: Incorrect number of traces specified: %d. This is most likely a program bug in the calling function", nTraces) );
  }
}
void csExecPhaseDef::setExecType( int theExecType ) {
  myExecType = theExecType;
}
void csExecPhaseDef::setTracesAreWaiting() {
  setTracesAreWaiting( true );
}
void csExecPhaseDef::setTracesAreWaiting( bool areWaiting ) {
  myTracesAreWaiting = areWaiting;
}
bool csExecPhaseDef::tracesAreWaiting() const {
  return myTracesAreWaiting;
}
void csExecPhaseDef::setTraceMode( int traceMode_in ) {
  traceMode = traceMode_in;
}
void csExecPhaseDef::setLastCall( bool isLastCall ) {
  myIsLastCall = isLastCall;
}
bool csExecPhaseDef::isLastCall() const {
  return myIsLastCall;
}
int csExecPhaseDef::getTraceMode() const {
  return traceMode;
}
void csExecPhaseDef::setCleanUp( bool isCleanUp ) {
  myIsCleanup = isCleanUp;
}
void csExecPhaseDef::setMPISupport( bool isSupported ) {
  myIsMPISupported = isSupported;
}
