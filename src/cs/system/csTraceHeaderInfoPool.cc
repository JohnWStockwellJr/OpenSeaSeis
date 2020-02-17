/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csTraceHeaderInfoPool.h"
#include "csTraceHeaderInfo.h"
#include "csVector.h"

using namespace cseis_system;

csTraceHeaderInfoPool::csTraceHeaderInfoPool() {
  myInfoList = new cseis_geolib::csVector<csTraceHeaderInfo*>(50);
}
csTraceHeaderInfoPool::~csTraceHeaderInfoPool() {
  if( myInfoList != NULL ) {
    for( int i = 0; i < myInfoList->size(); i++ ) {
      delete myInfoList->at(i);
    }
    myInfoList->clear();
    delete myInfoList;
    myInfoList = NULL;
  }
}
csTraceHeaderInfo const* csTraceHeaderInfoPool::createTraceHeaderInfo(
  cseis_geolib::type_t type, std::string const& name, std::string const& description, int nElements ) {
  csTraceHeaderInfo* info = new csTraceHeaderInfo( type, nElements, name, description );
  myInfoList->insertEnd( info );
  return info;
}

