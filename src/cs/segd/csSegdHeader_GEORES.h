/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGD_HEADERS_GEORES_H
#define CS_SEGD_HEADERS_GEORES_H

#include "csSegdDefines.h"
#include "csSegdHeader.h"

namespace cseis_segd {

/**
 * Extended Header - OYO Geospace GeoRes system
 *
 *
 */

class csExtendedHeader_GEORES : public csExtendedHeader {
public:
  csExtendedHeader_GEORES();
  virtual ~csExtendedHeader_GEORES();
  virtual void extractHeaders( byte const* buffer, int numBytes=0 );
  virtual void dump( std::ostream& os );
  virtual int nanoSeconds() const {
    return fracAcqNanoSeconds;
  }
public:
  int segMfgCode;       // 1:   SEG manufacturer code = 9
  int sponsorshipCode;  // 2-3: 0-255: Geospace, 256-...: other
  int contentType;      // 4:   Unique code identifying inf in this header module
  int contentVersion;   // 5:   Version of header content type
  int contentCompatibilityVersion; // 6: ...
  int totalCount32Blocks;// 7-8: Number of consecutive 32-byte blocks in this header module
  
  int fracAcqNanoSeconds; // 9-12:  Additional number of nanoseconds past second mark (in general header 1)
  int shotSerialNumber;   // 13-16: Serial number assigned to each shot
  int recordIndex;        // 17:    0-255. 
  int recordSerialNumber; // 18-19: 1-4095. Unique shot record index for each shot serial number
  int triggerMethod;      // 19
  int sampleClockSynch;   // 20 high nibble
  int sampleClockSynchStatus; // 20 low nibble
  int sampleClockSynchErr;// 21-23: Synchronization phase error [us]
  int numAuxChan;         // 24:    0-255
  int numReceiverChan;    // 25-26: 0-65535
  int recordTimeOffset;   // 27-28: Amount of milliseconds from trigger event [ms]
  int sourceType;         // 29
  int extTraceHeaderVersion;// 30
  int extTraceHeaderCompVersion; // 31
  int numExtTraceHeaderBlocks; // 32:   Number of extended trace header blocks per trace

// Block 2:
  int segMfgCode2;       // 1:   SEG manufacturer code = 9
  int sponsorshipCode2;  // 2-3: 0-255: Geospace, 256-...: other
  int contentType2;      // 4:   Unique code identifying inf in this header module
  int contentVersion2;   // 5:   Version of header content type
  int contentCompatibilityVersion2; // 6: ...
  int totalCount32Blocks2;// 7-8: Number of consecutive 32-byte blocks in this header module

  double srcEasting;  // 9-16  (64)
  double srcNorthing; // 17-24 (64)
  float srcElev;     // 25-29  (32)

// Block 3:
  int segMfgCode3;       // 1:   SEG manufacturer code = 9
  int sponsorshipCode3;  // 2-3: 0-255: Geospace, 256-...: other
  int contentType3;      // 4:   Unique code identifying inf in this header module
  int contentVersion3;   // 5:   Version of header content type
  int contentCompatibilityVersion3; // 6: ...
  int totalCount32Blocks3;// 7-8: Number of consecutive 32-byte blocks in this header module

  int fireMethod;
  int seqStartMethod;
  int timebreakMethod;
  int armDelay;
  int fireDelay;
  int timebreakWindowDelay;
  int timebreakWindowDuration;
  int startToTimebreak;
  int startToTimebreakFrac;
  int timebreakToSampleClock;

private:
  bool myIsLittleEndian;
};

/**
 * Demultiplexed trace header extension,all blocks
 *
 *
 *
 */
class csTraceHeaderExtension_GEORES : public csTraceHeaderExtension {
public:
  csTraceHeaderExtension_GEORES( int numBlocks );
  virtual ~csTraceHeaderExtension_GEORES();
  virtual void extractHeaders( byte const* buffer );
  virtual void extractHeaders( byte const* buffer, commonTraceHeaderStruct* comTrcHdr );
  virtual int getNumSamples() const {
    return block1.numSamples;
  }
  virtual void dump( std::ostream& os );
public:
  struct block1Struct {
    int rcvLineNumber;              // 01-03 INT24 Receiver line number
    int rcvPointNumber;             // 04-06 INT24 Receiver point number
    int rcvPointIndex;              // 07    INT8 Receiver point index
    int numSamples;                 // 08-10 UINT24 Number of samples per traces
    int extendedRcvLineNumber;      // 11-13 INT24 Extended rcv line Integer
    int extendedRcvLineFrac;        // 14-15 UINT24 Extended rcv line Fraction
    int extendedRcvPointNumber;     // 16-18 INT24 Extended rcv point Integer
    int extendedRcvPointFrac;       // 19-20 UINT24 Extended rcv point Fraction
    int sensorType;                 // 21    UNIT8 Sensor type
  } block1;

  struct block2Struct {
    double rcvEasting;        // 01-08 IEEE64
    double rcvNorthing;       // 09-16 IEEE64
    float rcvElevation;       // 17-20 IEEE32
    int recordIndex;          // 21    UINT8
    int traceDataProcessDesc; // 22    UINT8
    int secLowCutFilter;      // 23    UINT8
    int analogGainApplied;    // 24-25H INT12
    int sensorOrientDesc;     // 25L   UINT4
    // If sensorOrientDesc = 0
    //    --> No further headers
    // If sensorOrientDesc = 1:
    int sensorAzim;           // 27-28H INT12
    int sensorIncl;           // 28L-29 INT12
    int sensorAzimDev;        // 30-31H INT12
    int sensorInclDev;        // 31L-32 INT12
    // If sensorOrientDesc = 2:
    int sineOfSensorIncl;     // 27-28 INT16
    // If sensorOrientDesc = 3:
    int sensorVertContrib;    // 27-28 INT16
    int sensorInlineContrib;  // 29-30 INT16
    int sensorXlineContrib;   // 31-32 INT16
    // If sensorOrientDesc = 4:
    //Already defined: int sensorVertContrib;    // 27-28 INT16
    int sensorNorthContrib;   // 29-30 INT16
    int sensorEastContrib;    // 31-32 INT16
  } block2;

  struct block3Struct {
    int moduleSerialNumber;   // 01-03 UINT24
    int appliedSensorCalib;   // 04L-05 INT12
  } block3;
private:
  bool myIsLittleEndian;
};


} // end namespace
#endif

