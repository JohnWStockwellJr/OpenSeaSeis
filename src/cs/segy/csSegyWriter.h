/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGY_WRITER_H
#define CS_SEGY_WRITER_H

#include <cstdio>
#include <string>
#include "geolib_defines.h"
#include "csHeaderInfo.h"

namespace cseis_geolib {

  template<typename T> class csVector;
  class csSegyBinHeader;
  class csSegyTraceHeader;
  class csFlexNumber;
  class csSegyHdrMap;

/**
 * SEGY Writer
 *
 * Segy files are stored in Big Endian format, independent of which platform they were created on.
 * --> Determine platform. If Little Endian, swap endian format on output.
 *
 * @author Bjorn Olofsson
 * @date 2006
 */

class csSegyWriter {
//---------------------------------------------------------------------------------------
public:
  csSegyWriter( std::string filename, int nTracesBuffer,  bool reverseByteOrder, bool autoscale_hdrs, bool isSUFormat );
  ~csSegyWriter();
  /// Initialize
  void initialize( csSegyHdrMap const* hdrMap, char const* newCharHdr );
  /// Open SEGY file
  void openFile();
  /// Close SEGY file
  void closeFile();

  /// @return number of bytes per sample
  inline int sampleByteSize() const { return mySampleByteSize; }
  inline int totalTraceSize() const { return myTotalTraceSize; }
  int numTraceHeaders() const;
  inline byte_t* binHdrBlock() const { return myBinHdrBlock; }
  inline char* charHdrBlock() const { return myCharHdrBlock; }
  cseis_geolib::csSegyBinHeader* binHdr() const { return myBinHdr; }

  /// Free memory used to store EBCDIC and binary headers
  void freeCharBinHdr();
  /// @return SEGY file name
  inline char const* filename() const { return myFilename.c_str(); }
  inline char const* getFilename() const { return myFilename.c_str(); }

//  void writeNextTrace( byte_t const* buffer, int nSamples = 0 );
  void writeNextTrace( byte_t const* buffer, int nSamples );
  void writeNextTrace( byte_t const* buffer, csSegyTraceHeader const* trcHdr, int nSamples );

  void setIntValue( int hdrIndex, int value );
  void setFloatValue( int hdrIndex, float value );
  void setDoubleValue( int hdrIndex, double value );
  void setStringValue( int hdrIndex, std::string const& value );

  csSegyTraceHeader* getTraceHeader() const { return myTrcHdr; }

private:
  /// Actual SEGY header map, after user has added all relevant headers to be written to output file
  csSegyHdrMap* myTrcHdrMap;
  cseis_geolib::csFlexNumber* myHdrValues;
  csSegyTraceHeader* myTrcHdr;

  int myNumTraceHeaders; // Kept up to date, = myInfoList->size()
  char* myCharHdrBlock;
  byte_t* myBinHdrBlock;
  csSegyBinHeader* myBinHdr;
  int   mySampleByteSize;
  int   myTotalTraceSize;
  bool  myIsAtEOF;
  char* myBigBuffer;
  /// Index pointer to current trace in buffer
  int myCurrentTrace;
  /// Number of traces that are currently stored in the buffer
  int myNumSavedTraces;
  /// Total trace counter of all traces accessed via the getNextTrace() method
  int myTraceCounter;
  /// true if endian swapping shall be performed.
  bool myDoSwapEndian;
  /// true if Writer has been initialized, i.e. char & bin headers have been read in, and trace header is sealed
  bool myHasBeenInitialized;
  bool myIsEBCDIC;
  bool myIsAutoScaleHeaders;
  bool myForceToWrite;
  bool myIsSUFormat;

  std::FILE*  myFile;
  std::string myFilename;

  float mySampleInt;
  int   myNumSamples;
  int myDataSampleFormat;
  int const NTRACES_BUFFER;

//-----------------------------------------------------------------------------------------
// Private access methods
//
private:
  void setCharHdr( char const* newCharHdr );
  void writeCharBinHdr();

  csSegyWriter();
  csSegyWriter( csSegyWriter const& obj );
  csSegyWriter& operator=( csSegyWriter const& obj );
};

} // namespace

#endif

