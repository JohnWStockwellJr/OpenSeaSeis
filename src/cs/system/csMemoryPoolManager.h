/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_MEMORY_POOL_MANAGER_H
#define CS_MEMORY_POOL_MANAGER_H

#include <string>
#include <cstdio>
#include "geolib_defines.h"

namespace cseis_system {

  class csTracePool;
  class csTrace;
  class csTraceHeaderInfo;
  class csTraceHeaderInfoPool;
  class csTraceHeaderDef;

/**
 * Memory pool manager
 * Manages memory pools for seismic traces and trace header definitions.
 * ..could be widened to include more memory management, such as normal allocate requests.
 *
 * @author Bjorn Olofsson
 * @date   2007
 */

class csMemoryPoolManager {
 public:
  static int const POLICY_SPEED  = 651;
  static int const POLICY_MEMORY = 156;

  static csInt64_t const MAX_NUM_MEGABYTES = 12192;
public:
  csMemoryPoolManager();
  csMemoryPoolManager( int policy );
  ~csMemoryPoolManager();
  /**
  * @return new trace from memory pool
  */
  csTrace* getNewTrace();
  /**
   * Get new trace from memory pool. Initialize all necessary things: Trace headers, number of samples...
   * @return new trace from memory pool
   */
  csTrace* getNewTrace( csTraceHeaderDef const* hdefPtr, int numSamples );
  /**
   * Get new trace from memory pool. Initialize all necessary things: Trace headers, number of samples...
   * @return new trace from memory pool
   */
  csTrace* getNewTrace( csTrace const* trace_old );
  /**
  * @return new trace header info from memory pool
  */
  csTraceHeaderInfo const* getNewTraceHeaderInfo( cseis_geolib::type_t type, std::string const& name, std::string const& description );
  /**
  * @return new trace header info from memory pool
  */
  csTraceHeaderInfo const* getNewTraceHeaderInfo( cseis_geolib::type_t type, int nElements, std::string const& name, std::string const& description );
  /**
   * Dump summary
   * @param fout Output stream where dump shall be written to
   */
  void dumpSummary( FILE* fout ) const;
private:

  csMemoryPoolManager( csMemoryPoolManager const& obj );
  void init( int policy );
  bool checkMemory();
  /// Memory pool policy: Optimised for speed or memory usage
  int myPolicy;
  /// Trace pool buffer where all seismic traces are stored
  csTracePool* myTracePool;
  /// Trace header info pool buffer where all trace header info objects are stored
  csTraceHeaderInfoPool* myTraceHeaderInfoPool;
  /// Maximum number of bytes to allocate
  csInt64_t myMaxNumBytes;
  csInt64_t myMaxNumBytesAllocated;
};

} // namespace

#endif

