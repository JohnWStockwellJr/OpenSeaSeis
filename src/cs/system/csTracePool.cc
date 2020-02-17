/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdlib>
#include <cstdio>
#include "csTracePool.h"
#include "csMemoryPoolManager.h"
#include "csException.h"
#include "csTrace.h"
#include <map>

using namespace cseis_system;

csTracePool::csTracePool( int policy ) {
  myNumAllocatedTraces = 0;
  myIndexNextFreeTrace = 0;
  myNumUsedTraces      = 0;
  myMaxNumUsedTraces   = 0;
  myPolicy             = policy;
  myTraceIndexMap = new std::map<int,int>();
  if( myPolicy == csMemoryPoolManager::POLICY_SPEED ) {
    myBlockSize = 4;
  }
  else {
    myBlockSize = 1;
  }
  reallocate( myBlockSize );
}
csTracePool::~csTracePool() {
  if( myTraces != NULL ) {
    for( int i = 0; i < myNumAllocatedTraces; i++ ) {
      if( myTraces[i] ) delete myTraces[i];
    }
    delete [] myTraces;
    myTraces = NULL;
  }
  if( myIsTraceFree != NULL ) {
    delete [] myIsTraceFree;
    myIsTraceFree = NULL;
  }
  if( myTraceIndexMap != NULL ) {
    delete myTraceIndexMap;
    myTraceIndexMap = NULL;
  }

}
//----------------------------------------------------
void csTracePool::reallocate( int numTracesNew ) {
  csTrace** trcs = new csTrace*[numTracesNew];
  bool* isFree = new bool[numTracesNew];
  if( myNumAllocatedTraces > 0 ) {
    for( int i = 0; i < myNumAllocatedTraces; i++ ) {
      trcs[i]   = myTraces[i];
      isFree[i] = myIsTraceFree[i];
    }
    delete [] myIsTraceFree;
    delete [] myTraces;
  }

// Only fully allocate one trace. Leave other traces to be allocated later
  isFree[myNumAllocatedTraces] = true;
  trcs[myNumAllocatedTraces]   = new csTrace(this);
  myTraceIndexMap->insert(std::pair<int,int>(trcs[myNumAllocatedTraces]->getIdentNumber(),myNumAllocatedTraces) );
  for( int i = myNumAllocatedTraces+1; i < numTracesNew; i++ ) {
    isFree[i] = true;  // is free yet not allocated yet
    trcs[i]   = NULL;
  }

  myTraces = trcs;
  myIsTraceFree  = isFree;
  myNumAllocatedTraces = numTracesNew;
}
//----------------------------------------------------
csInt64_t csTracePool::computeNumBytes() {
  csInt64_t numSamples = 0;
  for( int i = 0; i < myNumAllocatedTraces; i++ ) {
    if( !myIsTraceFree[i] && myTraces[i] != NULL ) {
      numSamples += myTraces[i]->numSamples();
    }
  }
  return( numSamples * (csInt64_t)sizeof(float) );
}
//----------------------------------------------------
//
void csTracePool::freeTrace( csTrace* trace ) {
  int identNumber = trace->getIdentNumber();
//  bool isOK = false;
  std::map<int,int>::iterator iter = myTraceIndexMap->find(identNumber);
  if( iter != myTraceIndexMap->end() ) {
    myIsTraceFree[iter->second] = true;
    myIndexNextFreeTrace = iter->second;
  }
  else {
    throw( cseis_geolib::csException("csTracePool::freeTrace: Error...") );
  }

//  if( !isOK ) {
//    printf("*** Free trace:  %3d %3d,  used: %3d   (NOT FOUND!!!!)\n", myIndexNextFreeTrace, myNumAllocatedTraces, myNumUsedTraces);
//  }
  myNumUsedTraces--;
  //  fprintf(stdout,"******* Free trace from trace pool *******\n");
  //  dump();
}
//----------------------------------------------------
//
csTrace* csTracePool::getNewTrace() {
  //  printf("*** Get new trace:  %3d %3d,   used: %3d", myIndexNextFreeTrace, myNumAllocatedTraces, myNumUsedTraces);
  if( myTraces[myIndexNextFreeTrace] == NULL ) {
    myTraces[myIndexNextFreeTrace] = new csTrace( this );
    myTraceIndexMap->insert(std::pair<int,int>(myTraces[myIndexNextFreeTrace]->getIdentNumber(),myIndexNextFreeTrace));
  }
  csTrace* trace = myTraces[myIndexNextFreeTrace];
//  printf(",   ID: %3d\n", myTraces[myIndexNextFreeTrace]->getIdentNumber() );
  myIsTraceFree[myIndexNextFreeTrace] = false;
  myNumUsedTraces++;
  myIndexNextFreeTrace++;
  if( myNumUsedTraces > myMaxNumUsedTraces ) myMaxNumUsedTraces= myNumUsedTraces;

  if( myNumUsedTraces >= myNumAllocatedTraces ) {
    reallocate( myNumAllocatedTraces+myBlockSize );
    if( myPolicy == csMemoryPoolManager::POLICY_SPEED ) myBlockSize += BLOCK_SIZE_ATOM;
  }

  if( myIndexNextFreeTrace >= myNumAllocatedTraces || !myIsTraceFree[myIndexNextFreeTrace] ) {
    for( int i = 0; i < myNumAllocatedTraces; i++ ) {
      if( myIsTraceFree[i] ) {
        myIndexNextFreeTrace = i;
        break;
      }
    }
    if( myIndexNextFreeTrace >= myNumAllocatedTraces || !myIsTraceFree[myIndexNextFreeTrace] ) {
      fprintf(stderr,"ERROR in trace pool. %d %d\n", myIndexNextFreeTrace, myNumAllocatedTraces );
    }
  }
  return trace;
}
//----------------------------------------------------
//
void csTracePool::dumpSummary( FILE* fout ) const {
  std::string speed  = "SPEED";
  std::string memory = "MEMORY";
  fprintf( fout," Total number of used/allocated traces:  %d/%d\n", myMaxNumUsedTraces, myNumAllocatedTraces );
}
void csTracePool::dump() {
  std::string speed  = "SPEED";
  std::string memory = "MEMORY";
  fprintf(stdout,"****** TracePool dump *******\n");
  fprintf(stdout," Number of traces allocated/used/max:   %d / %d / %d\n", myNumAllocatedTraces, myNumUsedTraces, myMaxNumUsedTraces );
  fprintf(stdout," Memory pool policy used:           %s\n", myPolicy == csMemoryPoolManager::POLICY_SPEED ? speed.c_str() : memory.c_str() );
  fprintf(stdout," Index of next free trace: %d\n", myIndexNextFreeTrace);
}

