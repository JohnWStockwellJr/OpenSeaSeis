/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_GCS90_HEADER_H
#define CS_GCS90_HEADER_H

#include "csSegdDefines.h"
#include <iostream>

namespace cseis_segd {

static int const SIZE_GCS90_FRONT   = 90;
static int const SIZE_GCS90_ARRAYS  = 4;
static int const SIZE_GCS90_GUNDATA = 22;

  struct gunDataStruct {
    int  portNumber;
    char gunMode;  // A[uto], M[anual], S[pare], O[ff]
    char detectMode;
    int seqNumber;
    char autoFire;
    float staticOffset; // [m]
    float gunDelay;     // [ms]
    float fireTime;     // [ms]
    float delta;        // [ms]
    float depth;        // [m]
  };
/**
*
* @author Bjorn Olofsson
*/
class csGCS90Header {
public:
  csGCS90Header();
  ~csGCS90Header();
  int totalSize();
  bool isOK() { return myHeaderSeemsOK; }
  void extractHeaders( byte const* buffer, int numBytes );
  void dump( std::ostream& os );
public:
  bool myHeaderSeemsOK;
  
  char gcs90HeaderFlag[7];   // *GCS90 Header Flag Bytes
  int  numBytesInBlock;
  char lineNumber[9];
  long shotPoint;
  char activeArrayMask[3];
  char triggerMode;        // I: Internal, E: External
  char date[9];
  char time[10];
  int currentSeqNumber;
  int numSubarrays;
  int numGunsInArray;
  int numActiveGuns;
  int numDeltaErrors;
  int numAutofires;
  int numMisfires;
  float deltaSpread;  // in [ms]
  int volumeFired;
  float averageDelta;
  float averageDeltaDeviation;
  int manifoldPressure;
  int deepTow;

// For each subarray  
  int* subarrayStringPressure;   // Strin pressure of each subarray

// For each gun
  gunDataStruct* gunData;
private:
  bool myHeadersAreExtracted;
  int myTotalByteSize;
};

} // end namespace

#endif

