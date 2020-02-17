/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TRACE_H
#define CS_TRACE_H

namespace cseis_system {

class csTraceHeader;
class csTraceHeaderData;
class csTraceHeaderDef;
class csTraceData;
class csTracePool;

/**
* Seismic trace
* Contains all bits related to a seismic trace in Cseis:
*  - Trace samples/trace data
*  - Header values
*  - Header definitons (pointer only)
*
* Also:
*  - Unique trace identifier
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTrace {
public:
  csTrace();
  csTrace( csTracePool* const tracePool );
  ~csTrace();
  /// @return pointer to trace header object
  csTraceHeader* getTraceHeader();
  csTraceHeader const* getTraceHeader() const;
  /// @return pointer to trace data object
  csTraceData* getTraceDataObject();
  csTraceData const* getTraceDataObject() const;
  /// @return pointer to trace data object
  float* getTraceSamples();
  float const* getTraceSamples() const;
  /// @return number of samples in trace
  int numSamples() const;
  /// @return number of trace headers in trace
  int numHeaders() const;
  inline int getIdentNumber() {
    return myIdentNumber;
  }
  /**
  * 'Free' trace, so that trace can be used elsewhere
  * Does not necessarily free any memory (this depends on how memory management is set up)
  * NOTE: Do not use this trace anymore after making this call!
  */
  void free();
  /**
   * 'Trim' data
   * Free unused memory if possible
   */
  void trim();

  int mpi_getByteSize() const;
  int mpi_compress( char* data ) const;
  void mpi_decompress( char const* data );

 private:
  /// Object containing seismic trace header definitions (name,format,...)
  csTraceHeader* myTraceHeader;
  /// Object containing seismic data samples
  csTraceData*   myData;

//---------------------------------------
// Trace buffer mangement
public:
  csTracePool* const myTracePoolPtr;
  /// Unique trace identifier
  int const myIdentNumber;
  /// Sequential trace counter
  static int myIdentCounter;
};

} // namespace
#endif

