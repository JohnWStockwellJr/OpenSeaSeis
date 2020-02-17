/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csTraceGather.h"
#include "csTrace.h"
#include "csVector.h"
#include "csMemoryPoolManager.h"
#include "csTraceHeader.h"
#include "csTraceHeaderDef.h"
#include "csTraceData.h"
#include "csException.h"

using namespace cseis_system;

csTraceGather::csTraceGather() {
  myMemoryManager = NULL;
  myTraceList = new cseis_geolib::csVector<csTrace*>(10);
}
void csTraceGather::setMemoryManager( csMemoryPoolManager* memManager ) {
  myMemoryManager = memManager;
}
csTraceGather::csTraceGather( csMemoryPoolManager* memManager ) {
  myMemoryManager = memManager;
  myTraceList = new cseis_geolib::csVector<csTrace*>(10);
}
csTraceGather::csTraceGather( csTraceHeaderDef* hdef ) {
  myMemoryManager = hdef->getMemoryManager();
  myTraceList = new cseis_geolib::csVector<csTrace*>(10);
}
csTraceGather::~csTraceGather() {
  if( myTraceList != NULL ) { delete myTraceList; myTraceList = NULL; }
}
csTrace*& csTraceGather::operator []( int itrc ) {
  return myTraceList->at(itrc);  // May throw exception
}
csTrace* csTraceGather::trace( int itrc ) {
  return myTraceList->at(itrc);  // May throw exception
}
csTrace const* csTraceGather::trace( int itrc ) const {
  return myTraceList->at(itrc);  // May throw exception
}

//--------------------------------------------------------------
csTrace* csTraceGather::createTrace( int atTraceIndex, csTraceHeaderDef const* hdefPtr, int nSamples ) {
  if( myMemoryManager == NULL ) throw( cseis_geolib::csException("csTraceGather::createTrace: This trace gather object has no memory manager") );
  csTrace* trace = myMemoryManager->getNewTrace( hdefPtr, nSamples );
  myTraceList->insert( trace, atTraceIndex );
  return trace;
}
//--------------------------------------------------------------
void csTraceGather::createTraces( int atTraceIndex, int nTraces, csTraceHeaderDef const* hdefPtr, int nSamples ) {
  if( myMemoryManager == NULL ) throw( cseis_geolib::csException("csTraceGather::createTraces: This trace gather object has no memory manager") );
  if( atTraceIndex > myTraceList->size() ) atTraceIndex = myTraceList->size();
  int lastTraceIndex = atTraceIndex + nTraces - 1;
  for( int index = atTraceIndex; index <= lastTraceIndex; index++ ) {
    csTrace* trace = myMemoryManager->getNewTrace( hdefPtr, nSamples );
    myTraceList->insert( trace, atTraceIndex );
  }
}
//--------------------------------------------------------------
void csTraceGather::freeTraces( int firstTraceIndex, int nTraces ) {
  for( int i = 0; i < nTraces; i++ ) {
    myTraceList->at(firstTraceIndex+i)->free();
  }
  deleteTraces( firstTraceIndex, nTraces );
}
void csTraceGather::freeAllTraces() {
  reduceNumTracesTo(0);
}
//------------------------------------------------------------------
void csTraceGather::deleteTrace( int traceIndex ) {
  deleteTraces( traceIndex, 1 );
}
void csTraceGather::deleteTraces( int firstTraceIndex, int nTraces ) {
  myTraceList->remove( firstTraceIndex, nTraces );
}
//------------------------------------------------------------------
//
void csTraceGather::reduceNumTracesTo( int nTraces ) {
  if( nTraces < 0 ) {
    nTraces = 0;
  }
  else if( nTraces >= myTraceList->size() ) {
    return;
  }
  // !! Free traces before removing them from list:
  for( int i = nTraces; i < myTraceList->size(); i++ ) {
    myTraceList->at(i)->free();
  }
  myTraceList->reduceSizeTo( nTraces );
}
int csTraceGather::numTraces() const {
  return myTraceList->size();
}
//------------------------------------------------------------------
void csTraceGather::moveTraceTo( int traceIndex, csTraceGather* traceGather, int toTraceIndex ) {
  if( toTraceIndex < 0 || toTraceIndex > traceGather->numTraces() ) {
    toTraceIndex = traceGather->numTraces();
  }
  traceGather->addTrace( trace(traceIndex), toTraceIndex );
  deleteTrace( traceIndex );
}
void csTraceGather::moveTracesTo( int firstTraceIndex, int nTraces, csTraceGather* traceGather ) {
  for( int itrc = firstTraceIndex; itrc < firstTraceIndex+nTraces; itrc++ ) {
    traceGather->addTrace( trace(itrc) );
  }
  deleteTraces( firstTraceIndex, nTraces );
}
void csTraceGather::copyTraceTo( int traceIndex, csTraceGather* traceGather ) {
  if( myMemoryManager == NULL ) {
    throw( cseis_geolib::csException("csTraceGather::copyTraceTo(): This instance has no memory manager. To use this function, Construct the trace gather instance using one of the other constructors") );
  }
  csTrace* traceOld = myTraceList->at( traceIndex );
  csTrace* traceNew = myMemoryManager->getNewTrace( traceOld );
  traceGather->addTrace( traceNew );
}
//------------------------------------------------------------------
//
void csTraceGather::addTrace( csTrace* trace ) {
  myTraceList->insertEnd( trace );
}
void csTraceGather::addTrace( csTrace* trace, int atTraceIndex ) {
  myTraceList->insert( trace, atTraceIndex );
}

