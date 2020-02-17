/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_BASE_HEADER_H
#define CS_BASE_HEADER_H

#include <iostream>
#include "csSegdDefines.h"
#include "geolib_defines.h"

namespace cseis_segd {

/**
 * Base class for various SEGD headers
 *
 *
 * @author Bjorn Olofsson
 * @date 2006
 */
class csBaseHeader {
public:
  csBaseHeader();
  csBaseHeader( int numBytes );
  virtual ~csBaseHeader();
  virtual void extractHeaders( byte const* buffer ) = 0;
  virtual void dump( std::ostream& cs );
  static int const BLOCK_SIZE = 32;
private:
  csBaseHeader( csBaseHeader const& obj );
};


//---------------- General header block #1 -------------------------
class csGeneralHeader1 : public csBaseHeader {
public:
  csGeneralHeader1();
  /// Get sample interval [us]
  int getSampleInt_us() {
    return (int)((1000.0*(float)(baseScanInterval))/16.0);
  }
  virtual ~csGeneralHeader1();
  virtual void extractHeaders( byte const* buffer );
  virtual void dump( std::ostream& cs );
public:
  int fileNum;                     // 01-02 BCD   File number (0-9999)
  int formatCode;                  // 03-04 BCD   Format code BCD
  int genConsts[12];               // 05-10 BCD   General constants  (1 for 408UL, 2 for SEAL)
  int year;                        // 11    BCD   Year (0-99)
  int numGeneralHdrBlocks;         // 12    UINT4 Number of additional general header blocks
  int julianDay;                   // 12-13 BCD   Day of year (xxx-)
  int hour;                        // 14    BCD   Hour of day
  int min;                         // 15    BCD   Minute of hour
  int sec;                         // 16    BCD   Second of minute
  int manufactCode;                // 17    BCD   Manufacturer's code
  int manufactSerialNum;           // 18-19 BCD   and serial number
  int numBytesPerScan;             // 20-22 BCD   Bytes per scan (multiplexed only)
  int baseScanInterval;            // 23    UINT8 Base scan interval
                                   // (0.25ms: 4, 0.5ms: 8, 1ms: 10, 2ms: 20, 4ms: 40)
  int polarityCode;                // 24    UINT4 Polarity
  int numScansPerBlock;            // 24-25 UINT4 Number of scans per block
  int recordTypeCode;              // 26    UINT4 Record type (8: normal, 2: test record)
  int recordLen;                   // 26-27 BCD   Record length (if FFFx, extended record length is used), in increments of 0.5x1024ms, from 0.5
  int recordLenMS;                 //             --> recordLen*0.5*1024ms
  int numScanTypesPerRecord;       // 28    BCD   Scan types per record
  int numChanSetsPerScanType;      // 29    BCD   Number of channels sets per scan type
  int numSkewBlocks;               // 30    BCD   Skew blocks
  int numExtendedHdrBlocks;        // 31    BCD   Extended header length
  int numExternalHdrBlocks;        // 32    BCD   External header length
};

//---------------- General header block #2 -------------------------
class csGeneralHeader2 : public csBaseHeader {
public:
  csGeneralHeader2();
  virtual ~csGeneralHeader2();
  virtual void extractHeaders( byte const* buffer );
  virtual void dump( std::ostream& cs );
  void set( int numBytes_externalHdrBlocks );
public:
  int expandedFileNum;              // 01-03 UINT24 Expanded file number
  int numExtendedChanSets;       // 04-05 UIN16  Extended channel sets and scan types
  int numExtendedHdrBlocks;         // 06-07 UINT16 Extended header blocks
  int numExternalHdrBlocks;         // 08-09 UINT16 external header blocks
  int revisionNum[2];               // 11-12 UINT8  SEG-D revision number
  int numGeneralTrailerBlocks;      // 13-14 UINT16 General trailer number
  int extendedRecordLen;            // 15-17 UINT24 Extended record length (0-128000ms)
  int generalHeaderBlockNum;        // 19    UINT8  General header block number
 private:
  int myNumBytes_externalHdrBlocks; // Number of bytes of field where 'number of external header blocks' is stored. In standard SEGD rev2, only 2 bytes are read in.
};

//---------------- General header block #n -------------------------
class csGeneralHeaderN : public csBaseHeader {
public:
  csGeneralHeaderN();
  virtual ~csGeneralHeaderN();
  virtual void extractHeaders( byte const* buffer );
  virtual void dump( std::ostream& cs );
//  static int const BLOCK_SIZE;
public:
  int expandedFileNum;            // 01-03 UINT24 Expanded file number
  double sourceLineNum;           // 04-08 INT24/INT16  Source line number
  double sourcePointNum;          // 09-13 INT24/INT16  Source point number
  int sourcePointIndex;           // 14    UINT8  Source point index
  int phaseControl;               // 15    UINT8  Phase control
  int vibratorType;               // 16    UINT8  Type vibrator
  int phaseAngle;                 // 17-18 INT16  Phase angle
  int generalHeaderBlockNum;      // 19    UINT8  General header block number
  int sourceSetNum;               // 20    UINT8  Source set number
};

//---------------- Channel set descriptor -------------------------
class csChanSetHeader : public csBaseHeader {
public:
  csChanSetHeader();
  virtual ~csChanSetHeader();
  virtual void extractHeaders( byte const* buffer );
  virtual void dump( std::ostream& cs );
public:
  int scanTypeNum;                // 01    BCD Scan type
  int chanSetNum;                 // 02    BCD Channel set number
  int chanSetStartTime;           // 03-04 UINT16 Channel set start time
  int chanSetEndTime;             // 05-06 UINT16 Channel set end time
  double mpFactor;                // 07-08 UINT16 Descaling exponent, converted to 2^mp
  int numChannels;                // 09-10 BCD Channels in this channel set
  int chanTypeID;                 // 11    UINT4 Channel type identification
                                  //  0 (0000) Unused
                                  //  1 (0001) Seismic
                                  //  2 (0010) Time break
                                  //  3 (0011) Up hole
                                  //  4 (0100) Water break
                                  //  5 (0101) Time counter
                                  //  6 (0110) External Data
                                  //  7 (0111) Other
                                  //  8 (1000) Signature/unfiltered
                                  //  9 (1001) Signature/filtered
                                  // 12 (1100) Auxiliary trailer
  int subScanPerBaseScanBinExp;   // 12    BCD Sample/channel gain
  int gainControlCode;            // 12    UINT4 - Gain control method
  int aliasFilterFreq;            // 13-14 BCD Alias filter frequency
  int aliasFilterSlope;           // 15-16 BCD Alias filter slope
  int lowCutFilterFreq;           // 17-18 BCD Low cut filter frequency
  int lowCutFilterSlope;          // 19-20 BCD Low cut filter slope
  int firstNotchFilterFreq;       // 21-22 BCD 1st notch filter frequency
  int secondNotchFilterFreq;      // 23-24 BCD 2nd notch filter frequency
  int thirdNotchFilterFreq;       // 25-26 BCD 3rd notch filter frequency
  int extendedChanSetNum;         // 27-28 UIN16 Extended channel set number
  int extendedHeaderFlag;         // 29    UINT4 Extended header flag
  int numTraceHeaderExtensions_rev2; // 29    UINT4 Number of trace header extensions (rev 2 only)
  int verticalStack;              // 30    UINT8 Vertical stack
  int streamerCableNum;           // 31    UINT8 Streamer number
  int arrayForming;               // 32    UINT8 Array forming
};

//---------------- Sample skew -------------------------
class csSampleSkew : public csBaseHeader {
public:
  csSampleSkew();
  virtual ~csSampleSkew();
  virtual void extractHeaders( byte const* buffer );
  virtual void dump( std::ostream& cs );
public:
  int* skewFactor;
};

//--------------- Demultiplexed trace header -----------------------

class csTraceHeader : public csBaseHeader {
public:
  csTraceHeader();
  virtual ~csTraceHeader();
  virtual void extractHeaders( byte const* buffer );
  virtual void dump( std::ostream& cs );
public:
  static int const BLOCK_SIZE = 20;
  int fileNumber;                 // 01-02 BCD file number
  int scanType;                   // 03    BCD scan type (1-2)
  int chanSetNumber;              // 04    BCD channel set number
  int traceNumber;                // 05-06 BCD trace number
  int firstTimingWord;            // 07-09 UINT24 timing word of the first sample if
                                  //       the data were written in the multiplexed
                                  //       format
  int numTraceHeaderExtensions;   // 10    UINT4 number of trace header extensions
  int sampleSkew;                 // 11    UINT8 sample skew of the first sample
                                  //       of the trace. It is a part of
                                  //       the fractional part of the base scan interval.
                                  //       Res : 1/256 scan interval
  int traceEdit;                  // 12    UINT8 trace edit, 00: no edit, 01,02: Trace intentionally muted/zeroed, 03: Edited by acquisition system
  int timeBreakWindow;            // 13-15 UINT16 time from time break to the end of
                                  //       the internal time break window (binary
                                  //       number, inc 1ms)
  int extendedChanSetNumber;      // 16-17 UINT16 extended channel set number
  int extendedFileNumber;         // 18-20 UINT24 expanded file number
};

//---------------- Trace header extension -------------------------
class csTraceHeaderExtension : public csBaseHeader {
public:
  csTraceHeaderExtension( int numBlocks );
  virtual ~csTraceHeaderExtension();
  virtual void extractHeaders( byte const* buffer );
  virtual void extractHeaders( byte const* buffer, commonTraceHeaderStruct* comTrcHdr );
  virtual void dump( std::ostream& cs );
  virtual int getNumSamples() const;
  virtual int numBlocks() const;
  virtual int nanoSeconds() const;
  virtual csInt64_t timeSamp1_us() const;

public:
  struct block0Struct {
    int rcvLineNumber;              // 01-03 INT24 Receiver line number
    int rcvPointNumber;             // 04-06 INT24 Receiver point number
    int rcvPointIndex;              // 07    INT8 Receiver point index
    int numSamples;                 // 08-10 UINT24 Number of samples per traces
    int extendedRcvLineNumber;      // Rev 2: 11-15, only valid if rcvLineNumber set to FFFFFF
    int extendedRcvPointNumber;     // Rev 2: 16-20, only valid if rcvPointNumber set to FFFFFF
    int sensorType;                 // Rev 2: 21
  } block0;
 protected:
  int myNumBlocks;
};

//---------------- General trailer ---------------------------------

class csGeneralTrailer : public csBaseHeader {
public:
  csGeneralTrailer();
  virtual ~csGeneralTrailer();
  virtual void extractHeaders( byte const* buffer );
  virtual void dump( std::ostream& cs );
public:
  int generalTrailerNumber;       // 01-02 UINT16 General trailer number
  int chanTypeId;                 // 11    UINT8 Channel type identification
};

//---------------- Extended header ---------------------------------

class csExtendedHeader : public csBaseHeader {
public:
  csExtendedHeader();
  virtual ~csExtendedHeader();
  virtual void extractHeaders( byte const* buffer );
  virtual void extractHeaders( byte const* buffer, int totalNumBytes = 0 );
  virtual void dump( std::ostream& cs );
  virtual int nanoSeconds() const;
  virtual int numTraces() const;
};


class csBaseExternalHeader : public csBaseHeader {
public:
  csBaseExternalHeader( int numBlocks );
  virtual ~csBaseExternalHeader();
  virtual void extractHeaders( byte const* buffer ) = 0;
  virtual void extractHeaders( byte const* buffer, int totalNumBytes ) = 0;
  virtual void dump( std::ostream& cs );
private:
  int myNumBlocks;
};

} // end namespace

#endif

