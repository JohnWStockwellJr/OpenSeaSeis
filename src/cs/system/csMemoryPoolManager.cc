/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csMemoryPoolManager.h"
#include "geolib_defines.h"
#include "csException.h"
#include "csTrace.h"
#include "csTracePool.h"
#include "csTraceData.h"
#include "csTraceHeader.h"
#include "csTraceHeaderDef.h"
#include "csTraceHeaderInfo.h"
#include "csTraceHeaderInfoPool.h"

using namespace cseis_system;
using namespace cseis_geolib;

csMemoryPoolManager::csMemoryPoolManager( int policy ) {
  init(policy);
}
csMemoryPoolManager::csMemoryPoolManager() {
  init(csMemoryPoolManager::POLICY_SPEED);
}
void csMemoryPoolManager::init( int policy ) {
  myPolicy = policy;
  myTracePool = new csTracePool( myPolicy );
  myTraceHeaderInfoPool = new csTraceHeaderInfoPool();
  myMaxNumBytes = MAX_NUM_MEGABYTES * 1024L * 1024L;
  myMaxNumBytesAllocated = 0;
}
csMemoryPoolManager::~csMemoryPoolManager() {
  if( myTracePool ) {
    delete myTracePool;
    myTracePool = NULL;
  }
  if( myTraceHeaderInfoPool ) {
    delete myTraceHeaderInfoPool;
    myTraceHeaderInfoPool = NULL;
  }
}
csTrace* csMemoryPoolManager::getNewTrace() {
  checkMemory();
  //  static int counter = 0;
  //  static int const maxCount = 1;
  csTrace* trace = myTracePool->getNewTrace();
  return trace;
}
csTrace* csMemoryPoolManager::getNewTrace( csTraceHeaderDef const* hdefPtr, int numSamples ) {
  csTrace* traceNew = getNewTrace();
  traceNew->getTraceHeader()->setHeaders( hdefPtr );
  traceNew->getTraceDataObject()->set( numSamples );
  return traceNew;
}
csTrace* csMemoryPoolManager::getNewTrace( csTrace const* traceOld ) {
  csTrace* traceNew = getNewTrace();
  traceNew->getTraceHeader()->copyFrom( traceOld->getTraceHeader() );
  traceNew->getTraceDataObject()->setData( traceOld->getTraceDataObject() );
  return traceNew;
}
csTraceHeaderInfo const* csMemoryPoolManager::getNewTraceHeaderInfo( cseis_geolib::type_t type, std::string const& name, std::string const& description ) {
  csTraceHeaderInfo const* info = myTraceHeaderInfoPool->createTraceHeaderInfo( type, name, description );
  return info;
}
csTraceHeaderInfo const* csMemoryPoolManager::getNewTraceHeaderInfo( cseis_geolib::type_t type, int nElements, std::string const& name, std::string const& description ) {
  csTraceHeaderInfo const* info = myTraceHeaderInfoPool->createTraceHeaderInfo( type, name, description, nElements );
  return info;
}
void csMemoryPoolManager::dumpSummary( FILE* fout ) const {
  myTracePool->dumpSummary( fout );
  fprintf( fout," Total number of allocated (trace) memory:  %.2fkb  (= %.2fMb)\n",
           (double)myMaxNumBytesAllocated/(1024.0), (double)myMaxNumBytesAllocated/(1024.0*1024.0) );
}
bool csMemoryPoolManager::checkMemory() {
  if( myTracePool->numAvailableTraces() > 1 ) {
    return true;
  }
  csInt64_t numBytes = myTracePool->computeNumBytes();
  if( numBytes > myMaxNumBytesAllocated ) myMaxNumBytesAllocated = numBytes;
  //  fprintf(stderr,"Memory usage: %d traces (%d), %fkb, max allowed: %fkb\n",
  //        myTracePool->numAvailableTraces(), myTracePool->myNumAllocatedTraces, (double)numBytes/(1024.0), (double)myMaxNumBytes/(1024.0) );

  if( numBytes > myMaxNumBytes ) {
    throw( csException("\ncsMemoryPoolManager::checkMemory(): Bytes allocated: %fkb, max allowed: %fkb\n"\
                       "\nSeaSeis dynamically allocates traces for each module. Most trace-buffering usually occurs for multi-trace modules. " \
                       "Modules working on 'ensembles' are most prone to buffer a large amount of traces, such as sorting, stacking... " \
                       "In that respect, it is important that the user is aware of the correct setting of the 'ensemble trace header'. " \
                       "The ensemble header can be set using module 'ENS_DEFINE'. All consecutive traces for which the ensemble trace header "\
                       "value does not change constitute one ensemble.",
                       (double)numBytes/(1024.0), (double)myMaxNumBytes/(1024.0) ) );
  }
  return true;
}

