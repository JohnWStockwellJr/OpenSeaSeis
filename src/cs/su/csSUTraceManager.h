/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SU_TRACE_MANAGER_H
#define CS_SU_TRACE_MANAGER_H

#include <cstring>
#include <limits>
#include <cstdio>

#ifdef _WIN32
 extern "C" {
   #include <windows.h>
 }
#else
 extern "C" {
   #include <unistd.h>
 }
 #include "su.h"
 #include "segy.h"
#endif

#include "csVector.h"

namespace cseis_geolib {
  class csSegyTraceHeader;
  class csSegyHdrMap;
}

namespace cseis_su {

  static int const STATUS_NONE  = -1;
  static int const STATUS_TRACE_WAITING  = 1;
  static int const STATUS_EMPTY  = 2;
  static int const SU_FALSE = 0;
  static int const SU_TRUE  = 1;

  static int const HDRBYTES = 240;

/**
 * SU "trace manager"
 * - Manages the seismic traces exchanged between CSEIS and SU.
 * - Helps to pass on SU module self-documentation to CSEIS.
 * - Certain methods are designed to used only by CSEIS, and others by SU.
 * - Only one seismic trace is held at any one time
 *
 * Enables running of SU module within CSEIS flow, with minor modifications
 * to SU source code.
 */
class csSUTraceManager {
 public:
  csSUTraceManager();
  ~csSUTraceManager();

  /**
   * Set log file pointer
   * @param logFile Stream/file pointer to log file
   */
  void setLogFile( FILE* logFile );

  /**
   * @return Status: Either STATUS_TRACE_WAITING or STATUS_EMPTY (no trace available)
   */
  int getStatus() const { return myStatus; }

  /**
   * Retrieve seismic trace if available. Otherwise, return SU_FALSE.
   *
   * @return Success flag: SU_TRUE(1) or SU_FALSE(0)
   */
  int getTraceMaybe( int* trace );

  //-------------------- Methods to be used in SU program --------------------
  /**
   * Retrieve seismic trace (to be used by SU).
   * The program will wait within this method and only return when
   *  a) A trace is available, or
   *  b) The EOF flag was set.
   * @return Success flag: SU_TRUE(1) or SU_FALSE(0)
   */
  int getTrace( segy* trace );
  /**
   * Put a trace in the manager (to be used by SU)
   * 
   * @return Success flag: SU_TRUE(1) or SU_FALSE(0)
   */
  int putTrace( segy const* trace );

  //-------------------- Methods to be used in SU program --------------------
  /**
   * Retrieve pointer to trace (to be used by CSEIS)
   * The pointer is alive until method 'freeTrace' is called.
   * Make sure to call 'freeTrace()', otherwise the manager cannot accept any new traces
   *
   * @return Success flag: SU_TRUE(1) or SU_FALSE(0)
   */
  int getTracePtr( unsigned char const** bufferPtr );
  /**
   * Free seismic trace whose pointer was previously retrieved using method 'getTracePointer()'
   */
  void freeTrace();
  /**
   * Put a trace in the manager (to be used by CSEIS)
   *
   * @param suTrcHdr    SU trace header
   * @param samplesPtr  Pointer to data samples
   * @param numSamples  Number of samples in trace
   *
   * @return Success flag: SU_TRUE(1) or SU_FALSE(0)
   */
  int putTrace( cseis_geolib::csSegyTraceHeader const* suTrcHdr, float const* samplesPtr, int numSamples );

  /**
   * Set EOF flag
   * - No more traces will be placed in the manager. 
   * - Indicates that all trace retrieval processes can terminate.
   */
  void setEOF();
  /**
   * Set error message, error flag, and EOF flag
   *
   * @param text Error message, and optional argument list
   */
  void setError( char const* text, ... );
  /**
   * @return true if error flag has been set
   */
  bool isError() const;
  /**
   * @return  true if EOF flag has been set
   */
  bool isEOF() const { return myIsEOF; }
  /**
   * @return true if manager currently doesn't hold any traces
   */
  bool isEmpty() const;
  /**
   * Indicate that this call to the SU module is for retrieval of the self-doc only
   */
  void setDocRequestOnly();
  /**
   * @return Flag indicating whether this SU call is for retrieval of the self-doc only
   */
  bool isDocRequestOnly() const;

  /**
   * Helper method: Sleep specified number of milliseconds
   */
  static void sleep_ms( float milliseconds ) {
 #ifdef _WIN32
    Sleep( std::max(1,(unsigned int)(milliseconds/1000.0f)) );
 #else
    usleep( (unsigned int)(milliseconds * 1000.0f) ); // takes microseconds
 #endif
  }

  /**
   * Set SU module self-documentation
   */
  void setSUDoc( std::string& sdoc );
  /**
   * @return Pointer to SU module self-documentation
   */
  char const* getSUDoc() const;

  /**
   * @return Number of samples output by SU (if trace is available, -1 otherwise)
   */
  int numSamples() const;

  /**
   * @return Number of sample interval [ms] output by SU (if trace is available, -1 otherwise)
   */
  float sampleInt() const;

 private:
  void reallocateBuffer( int numSamples );

  /// Pointer to log file/stream: Do not allocate or deallocate
  std::FILE* myLogFilePtr;
  int myTrace;
  int myStatus;
  ///
  float mySleepMilliSec;
  /// EOF flag
  bool myIsEOF;
  /// Error flag
  bool myIsError;
  /// Flag indicating that this call to SU module is to retrieve self-doc only
  bool myIsDocRequestOnly;
  ///
  unsigned char* myBuffer;
  /// Numbwer of samples in trace
  int myNumSamples;
  std::string mySUDoc;
};


} // END namespace

#endif
