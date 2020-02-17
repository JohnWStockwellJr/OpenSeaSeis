/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TRACE_GATHER_H
#define CS_TRACE_GATHER_H

namespace cseis_geolib {
  template <typename T> class csVector;
}

namespace cseis_system {

  class csTrace;
  class csModule;
  class csMemoryPoolManager;
  class csTraceHeaderDef;

/**
* Collection of traces (pointers to the traces)
* This class should be optimized, maybe as a queue..?
* Also, access to trace generation (retrieval from pool) may be required
*
* This is the trace gather that is passed to the module exec methods
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTraceGather {
public:
  csTraceGather();
  csTraceGather( csMemoryPoolManager* memManager );
  csTraceGather( csTraceHeaderDef* hdef );
  ~csTraceGather();
  /**
  * Retrieve trace in gather
  * @param trace index (starting at 0)
  * @return pointer to trace at given trace index
  */
  csTrace*& operator [] ( int index );
  /**
  * Retrieve trace in gather
  * @param trace index (starting at 0)
  * @return pointer to trace at given trace index
  */
  csTrace* trace( int index );
  csTrace const* trace( int index ) const;
  /**
  * Create new trace in gather, add at end of gather.
  * Note that new trace sample values may not have been initialized.
  *
  * @param hdefPtr       Trace header definition that defines the trace header block of the new traces
  * @param nSamples      Number of samples in new traces
  */
  inline csTrace* createTrace( csTraceHeaderDef const* hdefPtr, int nSamples) { return createTrace( numTraces(), hdefPtr, nSamples ); }
  /**
  * Create new trace in gather, at specified trace index
  * Note that new trace sample values may not have been initialized.
  *
  * @param atTraceIndex  Create new trace, insert into gather at this trace index
  * @param hdefPtr       Trace header definition that defines the trace header block of the new traces
  * @param nSamples      Number of samples in new traces
  */
  csTrace* createTrace( int atTraceIndex, csTraceHeaderDef const* hdefPtr, int nSamples );
  /**
  * Create new traces in gather, at specified trace index
  * Note that new trace sample values may not have been initialized.
  *
  * @param atTraceIndex  Create new traces, insert into gather beginning at this trace index
  * @param nTraces       Number of traces to create
  * @param hdefPtr       Trace header definition that defines the trace header block of the new traces
  * @param nSamples      Number of samples in new traces
  */
  void createTraces( int atTraceIndex, int nTraces, csTraceHeaderDef const* hdefPtr, int nSamples );
  void reduceNumTracesTo( int nTraces );
  /**
  * Free specified trace
  * @param traceIndex  Index of trace to free
  */
  inline void freeTrace( int traceIndex ) {
    freeTraces( traceIndex, 1 );
  }
  /**
  * Free specified traces
  * @param firstTraceIndex  Index of first trace to free
  * @param nTraces          Number of traces to free
  */
  void freeTraces( int firstTraceIndex, int nTraces );
  /**
  * Free all traces
  */
  void freeAllTraces();
  /**
  * Move specified trace to another trace gather.
  * Insert trace at specified trace location, or at end of trace gather.
  */
  void moveTraceTo( int traceIndex, csTraceGather* traceGather, int toTraceIndex = -1 );
  void moveTracesTo( int firstTraceIndex, int nTraces, csTraceGather* traceGather );
  void copyTraceTo( int traceIndex, csTraceGather* traceGather );
  /**
  * @return number of traces currently in gather
  */
  int numTraces() const;
  /**
  * @return true if gather is empty
  */
  inline bool isEmpty() const { return numTraces() == 0; }
  friend class csModule;
private:
  void setMemoryManager( csMemoryPoolManager* memManager );
  csMemoryPoolManager* myMemoryManager;

  /**
  * Add trace to trace gather. Insert trace at end of gather.
  */
  void addTrace( csTrace* trace );
  void addTrace( csTrace* trace, int atTraceIndex );
  /**
  * Delete trace from trace gather
  * This removes the specified trace from this gather, without free-ing the trace
  */
  void deleteTrace( int traceIndex );
  /**
  * Delete traces from trace gather
  * This removes the specified traces from this gather, without free-ing the trace
  */
  void deleteTraces( int firstTraceIndex, int nTraces );
  csTraceGather( csTraceGather const& obj );
  cseis_geolib::csVector<csTrace*>* myTraceList;
};

} // namespace
#endif

