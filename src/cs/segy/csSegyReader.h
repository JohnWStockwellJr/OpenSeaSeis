/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGY_READER_H
#define CS_SEGY_READER_H

#include <cstdio>
#include <string>
#include <fstream>
#include "csSegyHdrMap.h"
#include "csSegyDefines.h"
#include "csIReader.h"


namespace cseis_geolib {

  template<typename T> class csVector;
  class csSegyBinHeader;
  class csHeaderInfo;
  class csSegyHeaderInfo;
  class csFlexNumber;
  class csFlexHeader;
  class csSegyTraceHeader;
  class csIOSelection;

/**
 * SEGY reader
 *
 * Segy files are stored in Big Endian format, independent of which platform they were created on.
 * a) Determine platform.
 * b) If Little Endian, swap endian format on input.
 *
 * @author Bjorn Olofsson
 * @date 2006
 */


class csSegyReader : public cseis_geolib::csIReader {
public:
  static int const DEFAULT_BUFFERED_SAMPLES = 250000;
  struct SegyReaderConfig {
    SegyReaderConfig() {
      numTracesBuffer   = 20;
      segyHeaderMapping = csSegyHdrMap::SEGY_STANDARD;
      reverseByteOrderData = false;
      reverseByteOrderHdr  = false;
      autoscaleHdrs     = true;
      overrideSampleFormat = cseis_segy::AUTO;
      enableRandomAccess = false;
      isSUFormat = false;
    }
    int  numTracesBuffer;
    int  segyHeaderMapping;
    bool reverseByteOrderData;
    bool reverseByteOrderHdr;
    bool autoscaleHdrs;
    int  overrideSampleFormat;
    bool enableRandomAccess;
    bool isSUFormat;
  };

//---------------------------------------------------------------------------------------
public:
  /**
  * Constructor.
  * @param filename (i) Filename
  * @param config   (i) Configuration parameters for SEGY reader
  */
  csSegyReader( std::string filename, SegyReaderConfig const& config, csSegyHdrMap const* hdrMap = NULL );
  ~csSegyReader();
  /**
  * This method initialises the reader object, by reading in the SEGY char and binary headers and setting
  * all internal buffers etc.
  *
  * Call this method before starting to read in data.
  */
  void initialize();
  /// Call this method if the number of samples in the SEGY file shall be overridden
  void initialize( int numSamples );

   /**
   * Set specific trace selection.
   * Only traces matching the selection will be read in and sorted, if requested.
   * 
   */
  bool setSelection( std::string const& hdrValueSelectionText, std::string const& headerName, int sortOrder, int sortMethod );
  void setSelectionStep1( std::string const& hdrValueSelectionText, std::string const& headerName, int sortOrder, int sortMethod );
  bool setSelectionStep2( int& numTracesToRead );
  bool setSelectionStep3();
  /**
   * @return Number of selected traces
   */
  int getNumSelectedTraces() const;
  /**
   * @param traceIndex  Index of consecutively selected trace
   * @return Value of selected trace
   */
  cseis_geolib::csFlexNumber const* getSelectedValue( int traceIndex ) const;
  int getSelectedIndex( int traceIndex ) const;

  /// @return Number of bytes per sample
  inline int sampleByteSize() const { return mySampleByteSize; }
  /// @return Number of bytes in one trace (including trace header and samples)
  inline int traceByteSize() const { return myTraceByteSize; }
  /// @return Number of defined trace headers
  int numTraceHeaders() const;
  /// @return Number of traces in file
  int numTraces() const { return myNumTraces; }
  /// @return SEGY file name
  inline char const* filename() const { return myFilename.c_str(); }
  inline byte_t const* binHdrBlock() const { return myBinHdrBlock; }
  inline char const* charHdrBlock() const { return myCharHdrBlock; }
  void dumpTrcHdr( FILE* fout ) const;
  cseis_geolib::csSegyBinHeader const* binHdr() const { return myBinHdr; }
  /// Set char header format. true for EBCDIC, false for ASCII
  void setCharHdrFormat( bool isEBCDIC );


  /// Free memory used to store EBCDIC and binary headers
  void freeCharBinHdr();
  /// Close SEGY file
  void closeFile();
  /**
  * Same as moveToTrace( int traceIndex ), with one additional parameter:
  *
  * @param numTracesToRead (i) Expected(!) number of traces that will be read at once by repeated calls to getNextTrace().
  *                            This parameter makes sure that no more than numTracesToRead will be read into a temporary buffer. This
  *                            will speed up the read process somewhat, depending on the number of traces to pre-buffer.
  *                            If unsure, use method call without this parameter
  */
  bool moveToTrace( int traceIndex, int numTracesToRead );
  bool moveToTrace( int traceIndex );
  /**
  * Retrieve next trace. Read in next trace, return sample values.
  * Trace header values must be retrieved by calls to methods hdrIntValue() etc...
  *
  * @param buffer   (o) Output buffer (sample values only, no trace header). Buffer must be fully allocated to hold one full trace.
  * @param nSamples (i) Number of samples to read
  * @return true if operation was successful, false if no more trace was found (end of file reached)
  */
  bool getNextTrace( byte_t* sampleBufferOut, int nSamples );
  bool getNextTrace( byte_t* sampleBufferOut );
  float const* getNextTracePointer( int nSamples );
  float const* getNextTracePointer();

  bool setHeaderToPeek( std::string const& headerName );
  bool setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType );
  bool peekHeaderValue( cseis_geolib::csFlexHeader* hdrValue, int traceIndex = -1 );   
  bool peek( int byteOffset, int byteSize, char* buffer, int traceIndex = -1 );
  bool revertFromPeekPosition();
  int getPeekByteOffset() const;

  csSegyTraceHeader const* getTraceHeader() const { return myTrcHdr; }

  csSegyHeaderInfo const* header( int hdrIndex ) const;

  bool addHeader( int byteLoc, int byteSize, type_t inType, cseis_geolib::csHeaderInfo const& hdr );
  bool addHeader( cseis_geolib::csSegyHeaderInfo const& hdr );
  bool addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theType );
  bool addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theType, std::string const& theDesc );

  int numSamples() const { return myNumSamples; }
  float sampleIntMS() const { return mySampleInt; }
  int dataSampleFormat() const { return myDataSampleFormat; }

  int numTracesCapacity() { return myBufferCapacityNumTraces; }
  int getCurrentTraceIndex() const { return myCurrentTraceInFile; }

  void setVariableTraceLength();
private:
  csSegyHdrMap*      myTrcHdrMap;
  csSegyTraceHeader* myTrcHdr;
  
  char* myCharHdrBlock;
  byte_t* myBinHdrBlock;
  csSegyBinHeader* myBinHdr;
  // Size in byte of one sample. Depends on data sample format, e.g. 4-byte floating point.
  int mySampleByteSize;
  /// Size in bytes of one full SEGY trace, including trace header and samples
  int myTraceByteSize;

  ///
  int myLastTraceIndex;

  /// Big buffer holding all pre-read traces
  char* myBigBuffer;
  /// Capacity in number of traces of big buffer = Maximum number of traces that can be buffered at once
  int myBufferCapacityNumTraces;
  /// Total number of traces that are currently stored in the buffer
  int myBufferNumTraces;
  /// Index pointer to trace in buffer that is currently being read
  int myBufferCurrentTrace;
  /// Trace index where file pointer currently stands. Reading from SEGY file the next time will read the trace with the given trace index.
  int myCurrentTraceInFile;
  /// Trace index of trace that is going to be retrieved next (either from buffer or input file)
  int myCurrentTraceRead;

  /// 'Copy' buffer. Holds data samples of one trace. Is only allocated on request (typically call from JNI interface)
  byte_t* myCopyBuffer;

//  Total trace counter of all traces accessed via the getNextTrace() method

  /// C++ file
  std::ifstream* myFile;
  /// File name
  std::string myFilename;
  /// Size of file in bytes
  csInt64_t myFileSize;

  /// true if endian swapping shall be performed, for data samples.
  bool myDoSwapEndianData;
  /// true if endian swapping shall be performed, for binary and trace headers.
  bool myDoSwapEndianHdr;
  /// true if reader has been initialized, i.e. char & bin headers have been read in, and trace header is sealed
  bool myHasBeenInitialized;
  /// true if SEGY char header is formatted in EBCDIC code
  bool myIsEBCDIC;
  /// true if some particular trace headers shall be scaled automatically (geometry headers)
  bool myIsAutoScaleHeaders;
  /// true if sample format code in SEGY binary header shall be overriden by manual format code
  int myOverrideSampleFormat;
  /// true if random access shall be enabled
  bool myEnableRandomAccess;
  /// Override number of samples found in in header
  bool myIsOverrideNumSamples;

  /// Sample interval [ms]
  float mySampleInt;
  /// Number of samples
  int myNumSamples;
  /// Preset definition for SEGY trace header mapping, e.g. cseis_segy::HDR_MAPPING_STANDARD
  //int mySegyHeaderMapping;
  /// SEGY data sample format, e.g. cseis_segy::DATA_FORMAT_IEEE
  int myDataSampleFormat;
  /// Number of traces in file as determined from file size
  int myNumTraces;
  /// Residual number of bytes at end of SEGY file (typically 0). This may be the case for corrupted SEGY files, or when SEGY file is being written at the same time.
  int myResidualNumBytesAtEnd;

  /// Header to check/peek
  int  myHdrCheckByteOffset;
  cseis_geolib::type_t myHdrCheckInType;
  cseis_geolib::type_t myHdrCheckOutType;
  int   myHdrCheckByteSize;
  char* myHdrCheckBuffer;
  bool myPeekIsInProgress;
  int myCurrentPeekTraceIndex;
  int myCurrentPeekByteOffset;
  int myCurrentPeekByteSize;
  bool myIsSUFormat;
  cseis_geolib::csIOSelection* myIOSelection;
  bool myIsVariableTraceLength;
public:
  void resetTrcHdrMap( csSegyHdrMap* map );
  csSegyHdrMap const* getTrcHdrMap() const { return myTrcHdrMap; }
//-----------------------------------------------------------------------------------------
// Private access methods
//
private:
  void readCharBinHdr();
  void openFile();

  csSegyReader();
  csSegyReader( csSegyReader const& obj );
  csSegyReader& operator=( csSegyReader const& obj );
};

} // namespace


#endif
