/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEISMIC_READER_VER_H
#define CS_SEISMIC_READER_VER_H

#include <cstdio>
#include <string>
#include <fstream>
#include "geolib_defines.h"

namespace cseis_geolib {
  class csHeaderInfo;
  class csFlexHeader;
}

namespace cseis_io {

class csSeismicIOConfig;

/**
 * Seismic file Reader, Cseis format
 *
 * Base class for version dependent seismic Reader
 * Note: As this is the first version, the versioning capabilities have not really been implemented yet.
 * Currently, this class handles everything, which should in part be moved to the version dependent class
 * 
 * @author Bjorn Olofsson
 * @date 2006
 */
class csSeismicReader_ver {
 public:
  static int const DEFAULT_BUFFERED_SAMPLES = 100000;
 public:
  csSeismicReader_ver( std::string filename, bool enableRandomAccess, int numTracesBuffer = 0 );
  virtual ~csSeismicReader_ver();
  
  virtual bool readFileHeader( csSeismicIOConfig* config ) = 0;     // pure virtual -> must be implemented in derived class
  
  virtual bool readTrace( float* samples, char* hdrValueBlock );
  virtual bool readTrace( float* samples, char* hdrValueBlock, int numSamples );
  virtual bool moveToTrace( int firstTraceIndex );
  virtual bool moveToTrace( int firstTraceIndex, int numTracesToRead );
  virtual bool moveToNextTrace();
  virtual void closeFile();
  static csSeismicReader_ver* createReaderObject( std::string filename, bool enableRandomAccess, int numTracesBuffer = 0 );
  static void extractVersionString( std::string filename, std::string& versionString );
  /**
   * Peek ahead in data file, at current position
   * @param byteOffset (i) Byte offset from current position
   * @param byteSize   (i) Byte size to peek
   * @param buffer     (o) Return buffer containing header value
   * @param traceIndex (i) Trace index at which header shall be peeked. Set to < 0 to peek header in current trace
   */
  virtual bool peek( int byteOffset, int byteSize, char* buffer, int traceIndex = -1 );
  /**
   *  Reset file pointer to start of current trace, from whereever it is at the moment
   */
  bool revertFromPeekPosition();

  int numTraces(void) { return myNumTraces; };
  int currentTraceIndex() const { return myCurrentTraceIndex; }

  int numTracesCapacity() const { return myBufferCapacityNumTraces; };

  int numSamples() const { return myNumSamples; }

protected:
  short int myVersionMinor;
  short int myVersionMajor;
  char myVersionText[4];
  std::string myFilename;
  /// C++ file
  std::ifstream* myFile;
  /// Size in bytes of trace header value block = All trace values of one trace
  int myByteSizeHdrValueBlock;
  /// Size in bytes of data samples of one trace
  int myByteSizeSamples;
  /// Size in bytes of one data sample
  int myByteSizeOneSample;
  /// Additional bytes needed for data compression
  int myByteSizeCompression;

  /// Buffer holding myNumTracesBuffer traces
  char* myDataBuffer;

  /// Number of samples in each trace
  int myNumSamples;
  
  /// Trace index where file pointer currently stands. Reading from input file the next time will read the trace with the given trace index.
  int myCurrentTraceIndex;
  /// Size in bytes of one full SEGY trace, including trace header and samples
  int myTraceByteSize;
  /// Capacity in number of traces of big buffer = Maximum number of traces that can be buffered at once
  int myBufferCapacityNumTraces;
  /// Total number of traces that are currently stored in the buffer
  int myBufferNumTraces;
  /// Index pointer to trace in buffer that is currently being read
  int myBufferCurrentTrace;
  ///
  int myLastTraceIndex;

protected:
  void open();
  bool initialize();
  void resizeDataBuffer( int newSize );
  bool readDataBuffer();
  bool readSingleTrace( float* samples, char* hdrValueBlock, int numSamples );
  bool seekg_relative( csInt64_t bytePosRelative );
  void decompressBuffer( float* samples, int numSamples );

  char* myTempBuffer;
  int myByteLoc;
  bool myIsReadFileHeader;
  /// File size in bytes
  csInt64_t myFileSize;
  /// Header size in bytes
  int myHeaderByteSize;
  /// Number of traces in input file
  int myNumTraces;
  /// true if random access is enabled
  bool myEnableRandomAccess;
  /// true if file pointer is currently elsewhere and has to be brought back before reading in a trace
  bool myPeekIsInProgress;
  int myCurrentPeekTraceIndex;
  int myCurrentPeekByteOffset;
  int myCurrentPeekByteSize;
  /// Index of first trace in current buffer
  int myBufferFirstTrace;

};

} // end namespace
#endif

