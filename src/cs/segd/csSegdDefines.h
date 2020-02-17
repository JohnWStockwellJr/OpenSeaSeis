/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef SEGD_DEFINES_H
#define SEGD_DEFINES_H

#include <string>
#include <iostream>
#include <cstdio>

/* char[] to different integer types conversion */

#define UINT4L(x)  ((*(x) & 0xf))
#define UINT4H(x)  ((*(x))>>4 & 0xf)
#define  INT4L(x)  ((*(x) & 0x8) == 0x8 ? UINT4L(x)-0x10 : UINT4L(x))
#define  INT4H(x)  ((*(x) & 0x80) == 0x80 ? UINT4H(x)-0x10 : UINT4H(x))

#define UINT12L(x) ((UINT4L(x)<<8) + *(x+1))
#define UINT12H(x) ((*(x)<<4) + UINT4H(x+1))
#define  INT12L(x) ((*(x) & 0x8) == 0x8 ? UINT12L(x)-0x1000 : UINT12L(x))
#define  INT12H(x) ((*(x) & 0x80) == 0x80 ? UINT12H(x)-0x1000 : UINT12H(x))

// Not sure if the following two methods are correct
#define UINT8(x)  (*(x))
#define  INT8(x)  ((*(x) & 0x8) == 0x8 ? UINT8(x)-0x10000 : UINT8(x))

#define UINT16(x)  ((*(x)<<8) + *(x+1))
#define  INT16(x)  ((*(x) & 0x80) == 0x80 ? UINT16(x)-0x10000 : UINT16(x))

#define UINT16_swap(x)  ((*(x+1)<<8) + *(x))

#define UINT24(x)  ((((*(x)<<8) + *(x+1))<<8) + *(x+2))
#define  INT24(x)  ((*(x) & 0x80) == 0x80 ? UINT24(x)-0x1000000 : UINT24(x))

//#define UINT32(x)  ((((((*(x)<<8) + *(x+1))<<8) + *(x+2))<<8) + *(x+3))
#define  INT32(x)  ((*(x) & 0x80) == 0x80 ? UINT32(x)-0x1000000 : UINT32(x))
#define UINT32(x)  (*(x+3) + (*(x+2)<<8) + (*(x+1)<<16) + (*(x)<<24))

// Does not work for some reason...use csSegdFunctions-->convert2uint64 instead
//#define UINT64(x) ( ((csInt64_t)*(buffer+7)) + (((csInt64_t)*(buffer+6)<<8UL)) + (((csInt64_t)*(buffer+5)<<16UL)) + (((csInt64_t)*(buffer+4)<<24UL)) + (((csInt64_t)*(buffer+3)<<32UL)) + (((csInt64_t)*(buffer+2)<<40UL)) + (((csInt64_t)*(buffer+1)<<48UL)) + (((csInt64_t)*(buffer)<<56UL)) )

namespace cseis_segd {
  typedef unsigned char byte;

  static int const UNKNOWN                  = -1;

  /// Manufacturer/recording system IDs
  static int const RECORDING_SYSTEM_SEAL    = 121;
  static int const RECORDING_SYSTEM_GEORES  = 122;
  static int const RECORDING_SYSTEM_GECO    = 123;
  static int const RECORDING_SYSTEM_DIGISTREAMER  = 124;
  static int const RECORDING_SYSTEM_FAIRFIELD  = 125;

  /// External header/nav system IDs
  static int const NAV_HEADER_NONE     = 200;
  static int const NAV_HEADER_PSI      = 201;
  static int const NAV_HEADER_LABO     = 202;
  static int const NAV_HEADER_HYDRONAV_VER1 = 203;
  static int const NAV_HEADER_HYDRONAV_VER6 = 204;

  /// Nav interface IDs
  static int const CM_DIGI_COMP_A = 101;
  static int const CM_DIGI_PSIB   = 102;
  static int const CM_DIGI_TS     = 103;
  static int const CM_DIGI_BIRD_DEPTHS = 104;

  static int const DUMP_NONE     = 0;
  static int const DUMP_GENERAL  = 2;
  static int const DUMP_CHANSET  = 4;
  static int const DUMP_EXTENDED = 8;
  static int const DUMP_EXTERNAL = 16;
  static int const DUMP_TRACE    = 32;
  static int const DUMP_EXTERNAL_BIRDS = 64;
  static int const DUMP_RAW      = 128;
  static int const DUMP_ALL      = DUMP_GENERAL + DUMP_CHANSET + DUMP_EXTENDED + DUMP_EXTERNAL + DUMP_TRACE;

  struct shotTimeStruct {
    shotTimeStruct() {
      year        = 0;
      julianDay   = 0;
      hour        = 0;
      minute      = 0;
      second      = 0;
      microSecond = 0;
    }
    shotTimeStruct( shotTimeStruct const& s ) {
      year        = s.year;
      julianDay   = s.julianDay;
      hour        = s.hour;
      minute      = s.minute;
      second      = s.second;
      microSecond = s.microSecond;
    }
    int year;
    int julianDay;
    int hour;
    int minute;
    int second;
    int microSecond;
  };

// general parameters that need to be requested
//  int seqNumber;

  /// Most common SEGD file headers that stay the same for all records in one sail line
  class commonFileHeaderStruct {
  public:
    commonFileHeaderStruct() {
      formatCode     = 0;
      formatCode     = 0;
      revisionNum[0] = 0;
      revisionNum[1] = 0;
      manufactCode   = 0;
      manufactSerialNum = 0;
      for( int i = 0; i < 14; i++ ) {
        srcLineName[i]    = ' ';
      }
      srcLineName[13]    = '\0';
      recordLength_ms = 0;
      sampleInt_us  = 0;
      numSamples    = 0;
      sampleBitSize = 0;
      totalNumChan  = 0;
      numSeismicChan= 0;
      numAuxChan    = 0;
    }
    /// SEGD file format code
    int formatCode;
    /// SEGD revision number
    int revisionNum[2];
    /// Manufacturer code
    int manufactCode;
    /// Manufacturer serial number
    int manufactSerialNum;
    /// Source line name
    char srcLineName[16];
    /// Record length [ms]  (= time of last sample..?)
    int recordLength_ms;
    /// Sample interval [us]
    int sampleInt_us;
    /// Number of samples
    int numSamples;
    /// Sample bit size
    int sampleBitSize;
    /// Total number of channels
    int totalNumChan;
    /// Number of seismic channels
    int numSeismicChan; 
    /// Number of auxiliary channels
    int numAuxChan;
  };

  /// Most common SEGD record headers that change from shot record to shot record
  class commonRecordHeaderStruct {
  public:
  commonRecordHeaderStruct() : shotTime()
      {
        fileNum = 0;
        shotNum = 0;
        srcEasting = 0.0;
        srcNorthing = 0.0;
        srcElev = 0.0;
      }
    commonRecordHeaderStruct( commonRecordHeaderStruct const& hdr )
      {
        fileNum     = hdr.fileNum;
        shotNum     = hdr.shotNum;
        srcEasting  = hdr.srcEasting;
        srcNorthing = hdr.srcNorthing;
        srcElev     = hdr.srcElev;
        shotTime    = hdr.shotTime;
      }
    /// Shot record field file number (FFID)
    int fileNum;
    /// Shot time
    shotTimeStruct shotTime;

    /// Shot point number
    int shotNum;
    /// Source easting [m]
    double srcEasting;
    /// Source northing [m]
    double srcNorthing;
    /// Source elevation [m]
    float  srcElev;
    void dump( FILE* fout ) {
      fprintf(fout,"fileNum: %d, shotTime: %d %d %d %d %d\n", fileNum, shotTime.year, shotTime.julianDay,
        shotTime.minute, shotTime.second, shotTime.microSecond );
      fprintf(fout,"shotNum: %d, srcEast: %f, srcNorth: %f, srcElev: %f\n", shotNum, srcEasting, srcNorthing, srcElev );
    }
  };

  /// Most common SEGD trace headers that change from trace to trace
  class commonTraceHeaderStruct {
  public:
    commonTraceHeaderStruct() {
      chanNum = 0;
      chanTypeID = 0;
      traceEdit = 0;
      rcvLineNumber = 0;
      rcvPointNumber = 0;
      rcvPointIndex = 0;
      rcvEasting   = 0.0;
      rcvNorthing  = 0.0;
      rcvElevation = 0.0;
      serialNumber = 0;
      sensor = 0;
      incl_i = 0;
      incl_c = 0;
      incl_v = 0;
      numSamples = 0;
      sampleInt_us = 0;
      chanSet = 0;
    }
    /// Channel number
    int chanNum;
    /// Channel type identifier
    int chanTypeID;
    /// Trace edit flag...
    int traceEdit;
    /// Receiver line number
    int rcvLineNumber;
    /// Receiver point number
    int rcvPointNumber;
    /// Receiver point index
    int rcvPointIndex;
    /// Receiver easting [m]
    double rcvEasting;
    /// Receiver northing [m]
    double rcvNorthing;
    /// Receiver elevation [m]
    float  rcvElevation;
    /// Sensor serial number
    int serialNumber;
    /// Sensor ID
    int sensor;
    /// Inclinometer 'ICV' values
    int incl_i;
    int incl_c;
    int incl_v;
    int numSamples;
    int sampleInt_us;
    int chanSet;
    void dump( FILE* fout ) {
      fprintf(fout,"chanNum: %d, chanTypeID: %d, traceEdit: %d\n", chanNum, chanTypeID, traceEdit);
      fprintf(fout,"rcvLineNum: %d, rcvPointNum: %d, rcvPointIndex: %d\n", rcvLineNumber, rcvPointNumber, rcvPointIndex );
      fprintf(fout,"sensor: %d, incl_i: %d, incl_c: %d, incl_v: %d\n", sensor, incl_i, incl_c, incl_v );
      fprintf(fout,"chanSet: %d, numSamples: %d, sampleInt_us: %d\n", chanSet, numSamples, sampleInt_us );
      //      fprintf(fout,"serialNum: %d, rcvEast: %f, rcvNorth: %f, rcvElev: %f, serialNum: %d\n", serialNumber, rcvEasting, rcvNorthing, rcvElevation );
    }
  };
//    int navStatus;             // 0: no status available, 1: Offline, 2: Approach, 3: Seabed, 4: Runout
//    int numDeadSeismicTraces;
//    int numLiveSeismicTraces;
//    float  waterDepth;
//    float  masterGyro;
//    float  masterCMG;
//    float  masterSpeed;
//    double masterEasting;
//    double masterNorthing;
//    double sourceDeltaEasting;
//    double sourceDeltaNorthing;
//    float  lineBearing;
//    int    lineStatus;
    //int testRecordType;       // 1 Field noise, 2 Field Impulse, 3....14

  /// Information in channel set header
  class commonChanSetStruct {
  public:
    commonChanSetStruct() {
      numSamples   = 0;
      sampleInt_us = 0;
      numChannels  = 0;
    }
    int numSamples;
    int sampleInt_us;
    int numChannels;
  };

int bcd( byte const* bytes, int startNibble, int numNibbles );
}

struct header
{
  std::string s;
  header(const std::string _s) : s(_s)
  {}
  friend std::ostream& operator<<(std::ostream &os, const header &h)
  {
    return os << std::string(16,'*') << "  " << h.s << "  " << std::string(16,'*') << std::endl;
  }
};

#endif

