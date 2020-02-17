/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGD_READER_H
#define CS_SEGD_READER_H

#include <fstream>
#include <cstdio>
#include <string>
#include "csSegdDefines.h"

namespace cseis_segd {

class csSegdBuffer;
class csBaseHeader;
class csGeneralHeader1;
class csGeneralHeader2;
class csGeneralHeaderN;
class csChanSetHeader;
class csSampleSkew;
class csTraceHeader;
class csTraceHeaderExtension;
class csGeneralTrailer;
class csExtendedHeader;
class csExternalHeader;
class csSegdHdrValues;

/**
* SEGD file reader
*
* Example usage:
*  csSegdReader::configuration config;
*  config.isDebug           = true;
*  config.navSystemID       = vars->navSystemID;
*  config.navInterfaceID    = cseis_segd::UNKNOWN;
*  config.recordingSystemID = cseis_segd::UNKNOWN;
*  config.thisIsRev0        = thisIsRev0;
*
*  try {
*    vars->segdReader = new csSegdReader( filename );
*    vars->segdReader->setConfiguration( config );
*    vars->segdReader->open();
*    vars->segdReader->readNewRecord()
*  }
*  catch( csException& e ) {
*  }
*
*  float* samples = ...;
*  cseis_segd::commonTraceHeaderStruct  comTrcHdr;
*  while( !vars->segdReader->getNextTrace( samples, comTrcHdr ) ) {  // Retrieve one trace at a time
*    if( vars->segdReader->readNextRecord( comRecordHdr ) ) {   // read in next SEGD record for multi-record files
*    }
*  }
*
* @author Bjorn Olofsson
*/

class csSegdReader {
public:
  struct configuration {
    int   recordingSystemID;
    int   navSystemID;
    int   navInterfaceID;
    bool  isDebug;
    bool  thisIsRev0;   // True if this is a revision 0 SEGD file
    bool  readAuxTraces;
    bool  numSamplesAddOne;
    bool  overrideTraceHdrExtensions;
    int numTraceHdrExtensions;
    int numBytes_genhdr2_extHdrBlocks;
  };
  struct bytePosition {
    int generalHdr1;
    int generalHdr2;
    int generalHdrN;
    int chanSetHdr;
    int extendedHdr;
    int externalHdr;
    int firstTrace;
    int firstTraceExt;
    int firstTraceData;
    int currentTraceData;
  };

public:
  csSegdReader();
  ~csSegdReader();
  void setChanSetToRead( int chanSetIndexToRead ) { myChanSetIndexToRead = chanSetIndexToRead; }
  /**
  * Set parameters that determine which additional tasks the SEGD reader shall perform
  */
  void setConfiguration( csSegdReader::configuration& config );
  /**
  * Open SEGD file.
  */
  bool open( std::string filename );
  /**
  * Close SEGD file.
  */
  void closeFile();
  /**
   * Read in all file headers, check consistency etc.
   */
  bool readNewRecordHeaders();
  /**
  * Read in next record.
  * Call only when next SEGD record shall be read in.
  * First, read in all traces of previous record.
  *
  * @param comRecHdr
  * @return false if no further record is present in input file
  */
  bool readNextRecord( commonRecordHeaderStruct& comRecHdr );
  /**
  * Retrieve next consecutive trace
  * Repeated calls to this method return traces one by one
  * @param trace      (o)
  * @param comTrcHdr  (o)
  * @return false if no further traces are found in this record, or if current record has not been fully read in yet
  */
  bool getNextTrace( float* trace, commonTraceHeaderStruct& comTrcHdr );
  float const* getNextTracePointer( commonTraceHeaderStruct& comTrcHdr );
  /**
  * @return common file headers
  */
  commonFileHeaderStruct const* getCommonFileHeaders() const {
    return &myComFileHdr;
  }
  /**
  * Return byte size of one full SEGD record
  */
  inline int recordByteSize() const { return myRecordByteSize; }
  /**
  * @return byte size of one seismic trace (data only)
  */
  int traceDataByteSize() const;
  int traceDataByteSize( int chanSetIndex ) const;
  int numTraces() const;
  int numHeaders() const;
  /**
  * Dump file header information for current record
  * @param dumpFlag  flag that determines which headers will be dumped
  */
  void dumpFileHeaders( int dumpFlag, std::ofstream* outStream );
  /**
  * Dump all trace header & trace header extensions
  */
  void dumpAllTraces( std::ofstream* outStream );
  /**
   * Dump last read trace
   */
  void dumpTrace( std::ofstream* outStream );

  /**
  * Dump essential information
  */
  void dumpEssentialFileInfo( FILE* dumpFile );
  /**
  * Dump record headers
  */
  void dumpRecordHeaders( FILE* dumpFile );
  
  csExtendedHeader const* extendedHdr() const {
    return myExtendedHdr;
  }
  csTraceHeaderExtension const* traceHdrExtension() const {
    return myTraceHdrExtension;
  }
  csSegdHdrValues const* hdrValues() const {
    return mySegdHdrValues;
  }
  csSegdHdrValues* hdrValues() {
    return mySegdHdrValues;
  }
  int numChanSets() const { return myNumChanSetsPerScanType; }
  void retrieveChanSetInfo( int chanSetIndex, commonChanSetStruct& info ) const;
  csGeneralHeader1 const* generalHdr1() const { return myGeneralHdr1; }
  csGeneralHeader2 const* generalHdr2() const { return myGeneralHdr2; }
  bool readCharHdr( char* charHdr, int numBytes );
private:
  /// 
  void setDefaultConfiguration();
  /// Read numBytes into buffer
  /// @return false if problem occurred
  bool readBuffer( byte* buffer, int numBytes );
  /// Extract common set of headers for current record
  void extractCommonRecordHeaders( commonRecordHeaderStruct& comRecHdr );
  /// Extract common set of headers for current trace
  void extractCommonTraceHeaders( byte const* bufferPtr, commonTraceHeaderStruct& commonTrcHdr );

  std::string myFileName;
  FILE* myFile;
  /// File where header information may be dumped
  FILE* myDumpFile;
  int myDumpFlag;

  //--------------------------------------------------------------
  // SEGD Headers
  
  csGeneralHeader1* myGeneralHdr1;  // Single pointer
  csGeneralHeader2* myGeneralHdr2;  // Single pointer
  csGeneralHeaderN* myGeneralHdrN;
  csChanSetHeader*  myChanSetHdr;
  int* myChanSetNumSamples;
  int* myChanSetSampleInt_us;
  int myChanSetIndexToRead;
  csSampleSkew*     mySampleSkew;
  csTraceHeader*    myTraceHdr;   // Single trace header
  csTraceHeaderExtension*  myTraceHdrExtension;   // Trace header extension for single trace
  csExtendedHeader* myExtendedHdr;
  csExternalHeader* myExternalHdr;

  /// Start byte position of all SEGD headers
  bytePosition myBytePos;
  
  //---------------------------------------------------------------
  // Base parameters extracted from SEGD headers, needed to read in SEGD file (needed to compute sizes..)

  int myNumScanTypes;
  int myNumChanSetsPerScanType;  // Maximum number of chan sets per scan types (if less, dummy chan sets are present)
  int myNumExtendedHdrBlocks;
  int myNumExternalHdrBlocks;
  int myNumGeneralTrailerBlocks;
  int myNumGeneralHdrBlocks;  // All general header blocks
  int myNumExtraGeneralHdrBlocks;  // Additional blocks to general header 1 & 2
  int myNumSampleSkewBlocks;
  int myNumTraceHdrExtensions;

  //---------------------------------------------------------------
  //
  // Configuration parameters for SEGD Reader
  configuration myConfig;
  /// Recording system identifier. Depends on manufacturer
  int myRecordingSystemID;
  // Byte size of one full trace, including trace header & extensions
  int* myFullTraceByteSize;
  // Common parameters for whole file
  commonFileHeaderStruct myComFileHdr;
  /// Chan type ID for all traces
  int* myChanTypeID;
  /// Byte size of one full SEGD record
  int myRecordByteSize;
  /// Buffer holding one full SEGD record, including all headers and trace samples
  byte* myBuffer_oneRecord;
  bool myHasJustBeenInitialized;
  /// Number of bytes read during initialization
  int myInitializeNumBytesRead;
  bool myIsSupported;
  int mySequentialTraceCounter;
  int mySequentialChanSetIndexCounter;
  int mySequentialChanSetAccTraces;
  //
  double* myMPDescaleOperator;
  // Segd header values. Collects (most) available headers into one object
  csSegdHdrValues* mySegdHdrValues;
  /// 'Copy' buffer. Holds data samples of one trace. Is only allocated on request (typically call from JNI interface)
  float* myCopyBuffer;

  byte* myBuffer_gen1;
  byte* myBuffer_gen2;

  csSegdBuffer* myBuffer_traceHdr;
  csSegdBuffer* myBuffer_traceHdrExtension;

  csSegdBuffer* buffer_genN;
  csSegdBuffer* buffer_scantypes;
  csSegdBuffer* buffer_extendedHdr;
  csSegdBuffer* buffer_externalHdr;

  // TEMP
  //  FILE* myTEMPFile;
};

} // end namespace
#endif

