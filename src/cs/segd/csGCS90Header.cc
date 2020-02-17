/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csGCS90Header.h"
#include "csSegdFunctions.h"
#include "csSegdDefines.h"

#include <cstring>
#include <string>

using namespace std;
using namespace cseis_segd;

csGCS90Header::csGCS90Header() {
  myHeadersAreExtracted = false;
  myTotalByteSize = 0;

  myHeaderSeemsOK = false;
  for( int i = 0; i < 6; i++ ) {
    gcs90HeaderFlag[i] = ' ';
  }
  gcs90HeaderFlag[6] = '\0';
  numBytesInBlock = 0;
  lineNumber[8] = '\0';
  shotPoint = 0;
  activeArrayMask[2] = '\0';
  triggerMode = '0';
  date[8] = '\0';
  time[9] = '\0';
  currentSeqNumber = 0;
  numSubarrays = 0;
  numGunsInArray = 0;
  numActiveGuns = 0;
  numDeltaErrors = 0;
  numAutofires = 0;
  numMisfires = 0;
  deltaSpread = 0.0f;
  volumeFired = 0;
  averageDelta = 0.0f;
  averageDeltaDeviation = 0.0f;
  manifoldPressure = 0;
  deepTow = 0;
  subarrayStringPressure = NULL;
  gunData = NULL;
}
csGCS90Header::~csGCS90Header() {
  if( subarrayStringPressure ) {
    delete [] subarrayStringPressure;
    subarrayStringPressure = NULL;
  }
  if( gunData ) {
    delete [] gunData;
    gunData = NULL;
  }
}
int csGCS90Header::totalSize() {
  if( !myHeadersAreExtracted ) return 0;
  else return myTotalByteSize; //numBytesInBlock;
}
void csGCS90Header::extractHeaders( byte const* buffer, int numBytes ) {
  if( numBytes < 90 ) return;
  
  myHeadersAreExtracted = true;
  myHeaderSeemsOK = true;

  memcpy( gcs90HeaderFlag, &buffer[0], 6 );
  if( strcmp( gcs90HeaderFlag, "*GCS90" ) ) {
    myHeaderSeemsOK = false;
    return;
  }
  numBytesInBlock = string2int( &buffer[6], 4 );
  memcpy( lineNumber, &buffer[10], 8 );
  shotPoint = string2long( &buffer[18], 10 );
  memcpy( activeArrayMask, &buffer[28], 2 );
  triggerMode = buffer[30];
  memcpy( date, &buffer[31], 8 );
  memcpy( time, &buffer[39], 9 );
  currentSeqNumber = string2int( &buffer[48], 1 );
  numSubarrays     = string2int( &buffer[49], 1 );
  numGunsInArray   = string2int( &buffer[50], 2 );
  numActiveGuns    = string2int( &buffer[52], 2 );
  numDeltaErrors   = string2int( &buffer[54], 2 );
  numAutofires     = string2int( &buffer[56], 2 );
  numMisfires      = string2int( &buffer[58], 2 );
  
  deltaSpread = (float)(string2int( &buffer[60], 2 )) * 0.1;
  volumeFired      = string2int( &buffer[63], 5 );
  averageDelta     = string2float( &buffer[68], 5 );
  averageDeltaDeviation = string2float( &buffer[73], 2 );
  
  manifoldPressure = string2int( &buffer[82], 4 );
  deepTow          = string2int( &buffer[86], 4 );

// For each subarray  
  int bytePos = 90;
  if( subarrayStringPressure ) delete [] subarrayStringPressure;
  subarrayStringPressure = NULL;
  if( numSubarrays <= 0 ) return;
  if( numBytes < bytePos+4*numSubarrays ) return;
  subarrayStringPressure = new int[numSubarrays];
  for( int i = 0; i < numSubarrays; i++ ) {
    subarrayStringPressure[i] = string2int( &buffer[bytePos], 4 );
    bytePos += 4;
  }

  if( gunData ) delete [] gunData;
  gunData = NULL;
  if( numGunsInArray <= 0 ) return;
  if( numBytes < bytePos+22*numGunsInArray+2 ) return;
  gunData = new gunDataStruct[numGunsInArray];
  for( int i = 0; i < numGunsInArray; i++ ) {
    gunDataStruct& gun = gunData[i];
    gun.portNumber = string2int( &buffer[bytePos], 2);
    gun.gunMode    = buffer[bytePos+2];
    gun.detectMode = buffer[bytePos+3];
    gun.seqNumber  = string2int( &buffer[bytePos+4], 1 );
    gun.autoFire   = buffer[bytePos+5];
    gun.staticOffset = string2float( &buffer[bytePos+7], 3 ) * 0.1;
    gun.gunDelay   = string2float( &buffer[bytePos+10], 3 )  * 0.1;
    gun.fireTime   = string2float( &buffer[bytePos+13], 3 )  * 0.1;
    gun.delta      = string2float( &buffer[bytePos+16], 3 )  * 0.1;
    gun.depth      = string2float( &buffer[bytePos+19], 3 )  * 0.1;
    bytePos += 22;
  }  
  bytePos += 2; // 2 extra bytes for CRLF
  if( numBytesInBlock != 0 ) {
    myTotalByteSize = numBytesInBlock;
  }
  else {
    myTotalByteSize = 92;
  }
}

void csGCS90Header::dump( std::ostream& os )
{
  os <<
  header("GCS90 header start") << '\n' << "GCS90 Header flag : " << gcs90HeaderFlag << '\n';
  if( !myHeaderSeemsOK ) {
    os << "...the flag is not recognized as a valid GCS90 header.. decoding aborted." << endl <<
    header("GCS90 header end") << endl;
    return;
  }
  os <<
  "numBytesInBlock   : " << numBytesInBlock << endl <<
  "lineNumber        : " << lineNumber << endl <<
  "shotPoint         : " << shotPoint << endl <<
  "activeArrayMask   : " << activeArrayMask << endl <<
  " triggerMode      : " << triggerMode << endl <<
  " date : " << date << endl <<
  " time : " << time << endl <<
  "currentSeqNumber  : " << currentSeqNumber << endl <<
  " numSubarrays     : " << numSubarrays << endl <<
  " numGunsInArray   : " << numGunsInArray << endl <<
  " numActiveGuns    : " << numActiveGuns << endl <<
  " numDeltaErrors   : " << numDeltaErrors << endl <<
  " numAutofires     : " << numAutofires << endl <<
  " numMisfires      : " << numMisfires << endl <<
  " deltaSpread      : " << deltaSpread << endl <<
  " volumeFired      : " << volumeFired << endl <<
  " averageDelta     : " << averageDelta << endl <<
  " averageDeltaDeviation : " << averageDeltaDeviation << endl <<
  " manifoldPressure : " << manifoldPressure << endl <<
  " deepTow          : " << deepTow << endl;
  
  for( int i = 0; i < numSubarrays; i++ ) {
    os << " String pressure subarray #" << i << ": " << subarrayStringPressure[i] << endl;
  }
  os << "Listing: Gun #num: portNum gunMode detectMode seqNum autoFire   statOff delay time delta depth\0" << endl;
  char line[256];
  for( int i = 0; i < numGunsInArray; i++ ) {
    gunDataStruct& gunTmp = gunData[i];
    sprintf( line, "Gun #%3d: %2d  %c  %c  %1d  %c  %5.1f %5.1f %5.1f  %5.1f  %5.1f",
      i+1, gunTmp.portNumber,
      gunTmp.gunMode,gunTmp.detectMode,
      gunTmp.seqNumber, gunTmp.autoFire,
      gunTmp.staticOffset, gunTmp.gunDelay,
      gunTmp.fireTime, gunTmp.delta, gunTmp.depth );
    os << line << endl;
  }
  
  os << header("GCS90 header end  ") << endl;
}

