/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGD_HEADERS_SEAL_H
#define CS_SEGD_HEADERS_SEAL_H

#include "csSegdDefines.h"
#include "csSegdHeader.h"

namespace cseis_segd {

/**
 * Extended Header - Sercel SEAL system
 *
 *
 */

class csExtendedHeader_SEAL : public csExtendedHeader {
public:
  csExtendedHeader_SEAL();
  virtual ~csExtendedHeader_SEAL();
  virtual void extractHeaders( byte const* buffer, int numBytes = 0 );
  virtual void dump( std::ostream& os );
public:
  int acqLength;  // 1000-128000ms
  int sampleInt; // 250,500,1000,2000,4000 us
  int totalNumTraces;    // 1-100000
  int numAuxTraces;      // 1-60
  int numSeismicTraces; 
  int numDeadSeismicTraces;
  int numLiveSeismicTraces;
  int sourceType;           // 0: no source, 1: Impulsive, 2: Vibro
  int numSamples;           // 1-128000
  int shotNumber;           // 1-9999
  float TBwindow;           // 0-64s
  int testRecordType;       // 1 Field noise, 2 Field Impulse, 3....14
  int swathFirstLine;    // 
  int swathFirstNum;
  int spreadNumber;
  int spreadType;        // 1: Generic, 2: Absolute
  int timeBreak;         // 0-9999us
  int upholeTime;        // us
  int blasterID;
  int blasterStatus;
  int refractionDelay;   // ms
  int TBtoT0time;        // us
  int internalTimeBreak;      // 0: no, 1: yes
  int prestackWithinFieldUnits;// 0: no, 1: yes
  int noiseElimType;     // 1: None, 2: Diversity stack, 3: Historic, 4: Off
  int lowTracePerc;      // 0-100%
  int lowTraceValue;     // 0-132dB
  int numWindows;
  int historicEditType;
  int noisyTracePerc;
  // .....

  int historicRange;     // 0-36dB
  int historicTaperLengthExponent;
  int thresholdVar;
  int historicThresholdInitValue;
  int historicZeroingLength;
  int typeOfProcess;             // 1: No operation (raw data), 2: Stack, 3: Correlation after stack, 4: Correlation before stack
  int acquisitionTypeTables[32];
  int thresholdTypeTables[32];
  int stackFold;   // 1 to 32
  char auxCorrelationDesc[80];
  int investigationLength;    // 100-128000ms
  int autocorrelationPeakTime;
  int fftExponent;
  int sourceAuxNB;
  int pilotLength;
  int sweepLength;
  int acquisitionNum;
  float maxOfMaxAux;          // IEEE
  float maxOfMaxSeis;         // IEEE
  int dumpStackFold;

  char tapeLabel[16];
  int tapeNum;
  char softwareRevision[16];
  char date[12];
  double srcEasting;
  double srcNorthing;
  float srcElev;
  int filesPerTape;
  int fileCount;
  char acquisitionErrorDesc[160];
  int filterType;            // 1: 0.8 Minimum, 2: 0.8 Linear, 3: SU6R-like
  int profileNum;            
  int t0mode;                // 1: Internal closure, 2: External closure
  char srcLineName[16];
  int seqNumber;
  int navStatus;             // 0: no status available, 1: Offline, 2: Approach, 3: Seabed, 4: Runout
  int shotsToGo;

  int shotsBeforeSeabed;
  int totalNumEditions;
  float auxDigitalLCfilter;  // 0: OFF
  float seisDigitalLCfilter; // 0: OFF
  float seisLC3dbDown;
  int numCharsFromNav;
private:
  bool myIsLittleEndian;

};


/**
 * Demultiplexed trace header extension,all blocks
 *
 *
 *
 */
class csTraceHeaderExtension_SEAL : public csTraceHeaderExtension {
public:
  csTraceHeaderExtension_SEAL( int numBlocks );
  virtual ~csTraceHeaderExtension_SEAL();
  void extractHeaders( byte const* buffer, commonTraceHeaderStruct* comTrcHdr );
  virtual void extractHeaders( byte const* buffer );
  virtual int getNumSamples() const {
    return block1.numSamples;
  }
  virtual void dump( std::ostream& os );
public:
  struct block1Struct {
    int rcvLineNumber;
    int rcvPointNumber;
    int rcvPointIndex;
    int numSamples;
    int extendedRcvLineNumber;
    int extendedRcvPointNumber;
    int sensorType;
  } block1;

  struct block2Struct {
    double rcvEasting;
    double rcvNorthing;
    float rcvElevation;
    int sensorTypeNB;
  } block2;

  struct block3Struct {
    float resistanceLimitLow;
    float resistanceLimitHigh;
    float resistanceValue;
    float tiltLimit;
    float tiltValue;
    int resistanceError;
    int tiltError;
  } block3;

  struct block4Struct {
    float capacitanceLimitLow;
    float capacitanceLimitHigh;
    float capacitanceValue;
    float cutOffLimitLow;
    float cutOffLimitHigh;
    float cutOffValue;
    int capacitanceError;
    int cutOffError;
  } block4;

  struct block5Struct {
    float leakageLimit;
    float leakageValue;
    int leakageError;
  } block5;

  struct block6Struct {
    int unitType; // 0: not identified, 1: FDU, 2: SU6R
    int unitSerialNumber;
    int chanNumber;
    int fduAssemblyType;
    int fduAssemblySerialNumber;
    int fduAssemblyLocation;
    int fduUnitType;
    int chanType;
    float sensorSensitivity;
  } block6;

  struct block7Struct {
    int controlUnitType;
    int controlUnitSerialNumber;
    int chanGainScale;
    int chanFilter;
    int chanDataError;
    int chanEditStatus; // 1: dead, 2: acquisition/retrieve error, 3: noise edition
    float sampleToMVconversionFactor;
    int numStacksLow;
    int numStacksNoisy;
    int chanTypeID;
    int chanProcess;
    float analogLowcutFilter;
    float digitalLowcutFilter;
    int numEditions;
    float compoundLowcutFilter;
  } block7;
private:
  bool myIsLittleEndian;

};



} // end namespace
#endif

