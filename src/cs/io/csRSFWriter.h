/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_RSF_WRITER_H
#define CS_RSF_WRITER_H

#include <cstdio>
#include <fstream>
#include <string>
#include "geolib_defines.h"

namespace cseis_io {

  template<typename T> class csVector;
  class csRSFHeader;

class csRSFWriter {
//---------------------------------------------------------------------------------------
public:
  csRSFWriter( std::string filename, int nTracesBuffer, bool reverseByteOrder, int numTraces_swapDim3, bool outputGrid, double tolerance = 1e-5 );
  ~csRSFWriter();
  //
  void initialize( csRSFHeader const* );
  void openFile();
  void closeFile();

  /// @return number of bytes per sample
  inline int sampleByteSize() const { return mySampleByteSize; }
  inline int totalTraceSize() const { return myTotalTraceSize; }
  inline char const* filename() const { return myFilename.c_str(); }
  inline char const* getFilename() const { return myFilename.c_str(); }

  void writeNextTrace( cseis_geolib::byte_t const* buffer, int nSamples, double valDim2, double valDim3 );
  void finalize();
private:
  void finalizeSwapDim3();
  static int const STAGE_1_INIT      = 1;
  static int const STAGE_2_SET_ORIG  = 2;
  static int const STAGE_3_DIM2_STEP = 3;
  static int const STAGE_4_DIM3_STEP = 4;
  static int const STAGE_5_WRITE_RSF = 5;
  static int const STAGE_6_COMPLETE  = 6;
  bool check( double valDim2, double valDim3 );

  csRSFHeader* myHdr;

  int   mySampleByteSize;
  int   myTotalTraceSize;
  bool  myIsAtEOF;
  char* myBigBuffer;
  /// Index pointer to current trace in buffer
  int myCurrentTrace;
  /// Number of traces that are currently stored in the buffer
  int myNumSavedTraces;
  /// Total trace counter of all traces written to disk
  int myTraceCounter;
  /// true if endian swapping shall be performed.
  bool myDoSwapEndian;
  /// Stage of writer object: Initialized, first trace, second trace... finalized
  int myStage;
  /// true if Writer has been initialized, 
  //  bool myHasBeenInitialized;
  /// true if Writer has been finalized, i.e. rsf header file has been written
  //  bool myHasBeenFinalized;
  double myPreviousValueDim2;
  double myPreviousValueDim3;
  int myCurrentCounterDim2;

  std::ofstream*  myFile;
  std::string myFilename;

  float mySampleInt;
  int   myNumSamples;
  int   myNumTracesBuffer;
  double const TOLERANCE;
  int myNumTraces_swapDim3; // true if dimension 3 and dimension 1 shall be swapped
  bool myOutputGrid; // true if grid definiton shall be written to rsf header file
//-----------------------------------------------------------------------------------------
// Private access methods
//
private:
  void writeRSFHdr();

  csRSFWriter();
  csRSFWriter( csRSFWriter const& obj );
  csRSFWriter& operator=( csRSFWriter const& obj );
};

} // namespace

#endif

