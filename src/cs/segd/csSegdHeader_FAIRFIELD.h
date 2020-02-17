/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGD_HEADER_FAIRFIELD_H
#define CS_SEGD_HEADER_FAIRFIELD_H

#include "csSegdDefines.h"
#include "csSegdHeader.h"
#include "geolib_defines.h"

namespace cseis_segd {

/**
 * Extended Header - Fairfield
 *
 */

class csExtendedHeader_FAIRFIELD : public csExtendedHeader {
public:
  csExtendedHeader_FAIRFIELD();
  virtual ~csExtendedHeader_FAIRFIELD();
  virtual void extractHeaders( byte const* buffer, int numBytes=0 );
  virtual void dump( std::ostream& os );
  virtual int nanoSeconds() const;
  virtual int numTraces() const;

public:
// Block 1:
  csInt64_t unitID;               //   1-8 ID number of remote unit  
  csInt64_t timeDeploy_epoch_us;  //  9-16 Epoch deployment time
  csInt64_t timeRetrieve_epoch_us;// 17-24 Epoch pickup time
  csInt64_t timeStart_epoch_us;   // 25-32 Remote unit Epoch start time

// Block 2:
  float driftWindow_us;        // 1-4  IEEE Acquisition drift window [us]
  csInt64_t clockDrift_ns;   // 5-12 Clock drift [ns]
  int clockStopMethod;       //   13 Clock stop method (0:normal,1:storage full,2:power loss,3:reboot)
  int freqDriftFlag;         //   14 Frequency drift flag (0:not within specs, 1:within specs)
  int oscillatorType;        //   15 Oscillator type (0:control based,1:atomic,2:ovenized,3:double overnized,4:disciplined)
  int dataCollectionMethod;  //   16 Data collection method (0:normal shots,1:continuous,2:shot sliced with guard bands)
  int numRecordsInFile;      // 17-20 Number of records (shots or time slices) in this file
  int numFilesInUnit;        // 21-24 Total number of files to be acquired for this remote unit
  int fileNum;               // 25-28 File number X
  int dataDecimationFlag;    //   29  Data decimation flag (0:not decimated,1:decimated)
  int baseScanInterval;      //   30  Original base scan interval
  int numDecimationFiltCoef; // 31-32 Number of decimation flter coefficients

// Block 3:
  int rec_line;    //   1-4 Receiver line number (replicated from other places)
  int rec_point;   //   5-8 Receiver point (replicated from other places)
  int rec_index;   //     9 Receiver point index (replicated)
  int sou_line1;   // 10-13 First shot line written to this file
  int sou_point1;  // 14-17 First shot point written to this file
  int sou_index1;  //    18 First shot point index written to this file
  int sou_line2;   // 19-22 Last shot line written to this file
  int sou_point2;  // 23-26 Last shot point written to this file
  int sou_index2;  //    27 Last shot point index written to this file
//  int reserved;   // 28-32 (set to 0)

private:
  bool myIsLittleEndian;
};

/**
 * Demultiplexed trace header extension, all blocks
 *
 *
 *
 */
class csTraceHeaderExtension_FAIRFIELD : public csTraceHeaderExtension {
public:
  csTraceHeaderExtension_FAIRFIELD( int numBlocks );
  virtual ~csTraceHeaderExtension_FAIRFIELD();
  virtual void extractHeaders( byte const* buffer );
  virtual void extractHeaders( byte const* buffer, commonTraceHeaderStruct* comTrcHdr );
  virtual int getNumSamples() const {
    return block1.numSamples;
  }
  virtual void dump( std::ostream& os );
  virtual int nanoSeconds() const;
  virtual csInt64_t timeSamp1_us() const;

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
    // 00 - n/a
    // 01 - hydrophone
    // 02 - Vertical geophone
    // 03 - Inline geophone
    // 04 - Crossline geophone
    // 05 - Other horizontal geophone
    // 10 - Geophone X
    // 11 - Geophone Y
    // 12 - Geophone Z
    int traceCount; // 22-25 Trace count within file
  } block1;

  struct block2Struct { // "Shot" block
    int sou_line;  // 1-4
    int sou_point; // 5-8
    int sou_index; // 9
    int sou_x_preplot; // 10-13
    int sou_y_preplot; // 14-17
    int sou_x_fin; // 18-21 x10
    int sou_y_fin; // 22-25 x10
    int sou_z_fin;  // 26-29 x10
    int souInfoSrcFlag;    // 30 1:preplan, 2:as-shot, 3:post-processed
    int sou_type; // 31 Energy source type, 0:undefined,1:vibroseis,2:dynamite,3:airgun
    int sou_id; // 32 Source identifier
  } block2;

  struct block3Struct { // "Time" block
    csInt64_t shotTime_epoch_us;  // 01-08 Int64
    csInt64_t shotSkewTime_us;    // 09-16 Int64 Shot skew time from sample boundary
    csInt64_t clockShiftApplied_ns; // 17-24 Applied clock correction time shift [ns]
    csInt64_t clockShiftRemainder_ns; // 25-32 Remaining (not applied) clock correction time shift [ns]
  } block3;

  struct block4Struct { // "Data" block
    int timePreGuard_ms; // 1-4 Pre Shot Guard Band [ms]
    int timePostGuard_ms;// 5-8 Post Shot Guard Band [ms]
    int gain_db;         // 9 Preamp gain [dB]
    int traceClipFlag;   // 10 0:not clipped, 1:digital clip, 2:analog clip
    int recordTypeCode;  // 11
    int shotStatusFlag;  // 12
    int external_sou_id; // 13-16
    int dataSampleFlag;  // 17
    // 0: First sample time is one sample interval after the Epoch time of start of the record
    // 1: First sample time is at Epoch time of the start of the record
    int fbptime;  // 25-28
    int rmsNoise; // 29-32
  } block4;

  struct block5Struct { // "Receiver" block
    int rec_line;  // 1-4  Receiver line number
    int rec_point; // 5-8  Receiver point
    int rec_index; // 9    Receiver point index
    int rec_x_preplot;       // 10-13 x10
    int rec_y_preplot;       // 14-17 x10
    int rec_x_fin;           // 18-21 x10
    int rec_y_fin;           // 22-25 x10
    float rec_z_fin;         // 26-29 x10
    int recInfoSourceFlag;   // 30 1:preplan, 2:as-laid, 3:...
  } block5;

  struct block6Struct { // "Orientation" block
    float tilt_h1x;           // 1-4 IEEE Float
    float tilt_h2x;           // 5-8 IEEE Float
    float tilt_vx;            // 9-12 IEEE Float
    float tilt_h1y;           // 13-16 IEEE Float
    float tilt_h2y;           // 17-20 IEEE Float
    float tilt_vy;            // 21-24 IEEE Float
    float tilt_h1z;           // 25-28 IEEE Float
    float tilt_h2z;           // 29-32 IEEE Float
  } block6;

  struct block7Struct { // "Orientation" block
    float tilt_vz;            // 1-4 IEEE Float
    float azim_deg;           // 5-8 IEEE Float
    float pitch_deg;          // 9-12 IEEE Float
    float roll_deg;           // 13-16 IEEE Float
    float temperature;        // 17-20 IEEE Float
    float humidity;           // 21-24 IEEE Float
    int orientMatrixVer;      // 25-28 Binary
    int gimbalCorrection;     // 29 Binary
    int headingRotation;      // 30-31 Binary Heading rotation angle (tenths of deg)
  } block7;

  struct block8Struct { // "Instrument Test" block
    int testCode; // 1-4 Magseis Fairfield Test Analysis code
    int test1_attn; // 5-8 First Test Oscillator Attenuation
    int test2_attn; // 9-12 Second Test Oscillator Attenuation
    int startDelay_us; // 13-16 Start Delay (in uSec)
    int dcFilterFlag; // 17-20 DC Filter Flag, 0 = No Filter, 1 = Apply Filter
    float dcFilterFreq; // 21-24 IEEE DC Filter frequency
    int preampPath; // 25-28 Preamp Path
    // 0 = external input selected (default)
    // 1 = simulated data selected
    // 2 = pre-amp input shorted to ground
    // 3 = test oscillator with sensors
    // 4 = test oscillator without sensors
    // 5 = common mode test oscillator with sensors
    // 6 = common mode test oscillator without sensors
    // 7 = test oscillator on positive sensors with neg sensor grounded
    // 8 = test oscillator on negative sensors with pos sensor grounded
    // 9 = test oscillator on positive PA input, with neg PA input ground
    // 10 = test oscillator on negative PA input, with pos PA input ground
    // 11 = test oscillator on positive PA input, with neg PA input ground, no sensors
    // 12 = test oscillator on negative PA input, with pos PA input ground, no sensors
    int testSignalType; // 29-32 Test Oscillator Signal Type
    // 0 = test oscillator path open
    // 1 = test signal selected
    // 2 = DC reference selected
    // 3 = test oscillator path grounded
    // 4 = DC reference toggle selected
  } block8;

  struct block9Struct { // "Instrument Test" block
  } block9;
  struct block10Struct { // "Instrument Test" block
  } block10;

private:
  bool myIsLittleEndian;
};


} // end namespace
#endif

