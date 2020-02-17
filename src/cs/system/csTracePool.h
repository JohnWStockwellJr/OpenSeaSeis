/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TRACE_POOL_H
#define CS_TRACE_POOL_H

#include <cstdio>
#include <map>
#include "geolib_defines.h"

namespace cseis_system {

class csTrace;

/**
* Trace pool
*
* Buffers seismic traces. From this pool, new traces are retrieved and returned after use, and can
* thereafter be re-used again.
* The main idea behind this pool is to minimise the amount of memory allocation/free applications.
*  ...the pool is used instead of allocating/freeing traces on the fly
*
* However, this increase in speed is bought by a greater need for memory.
*
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTracePool {
public:
  csTracePool( int policy );
  ~csTracePool();
  /**
  * @return pointer to new trace object
  */
  csTrace* getNewTrace();
  friend class csTrace;
  friend class csMemoryPoolManager;
  void dumpSummary( FILE* fout ) const;
  /// For debugging purposes
  void dump();
  int numAvailableTraces() { return myNumAllocatedTraces-myNumUsedTraces;  };

protected:
  csInt64_t computeNumBytes();

private:

  /// Traces
  csTrace** myTraces;
  /// Number of traces that are currently in use
  int myNumUsedTraces;
  /// Maximum number of traces in use at one time
  int myMaxNumUsedTraces;
  /// Number of traces allocated in trace pool
  int myNumAllocatedTraces;
  /// Index field, for each trace. If 'yes' for any particular trace, it means this trace is still 'free' (=available to be used)
  bool* myIsTraceFree;
  /// Index in trace array of a free trace
  int myIndexNextFreeTrace;
  /**
   * Size of block that is allocated anew when new trace shall be allocated. Increases steadily when new traces are allocated repeatibly.
   * Reason is to allocate some additional traces in advance to save time. This is faster than allocating each trace separately, but uses more memory
   */
  int myBlockSize;
  /// Memory pool policy: Optimised for speed or memory usage
  int myPolicy;

  void reallocate( int newNumAllocatedTraces );
  /**
  * 'Free' the according trace from the buffer pool
  * This does not free any memory!
  * It makes the trace buffer available again, which means the trace is free'd to be used elsewhere, see 'getNewTrace()'.
  */
  void freeTrace( csTrace* trace );
  static int const BLOCK_SIZE_ATOM = 4;
  std::map<int,int>* myTraceIndexMap;
};

} // namespace
#endif

