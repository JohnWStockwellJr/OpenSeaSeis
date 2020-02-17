/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegdHeader_FAIRFIELD.h"
#include "csSegdFunctions.h"
#include "geolib_endian.h"
#include <ostream>
#include <iostream>
//#include <bitset>
#include <cmath>
#include <cstring>

using std::endl;
using namespace cseis_segd;
using namespace std;

csExtendedHeader_FAIRFIELD::csExtendedHeader_FAIRFIELD() : csExtendedHeader() {
  myIsLittleEndian = cseis_geolib::isPlatformLittleEndian();
}
csExtendedHeader_FAIRFIELD::~csExtendedHeader_FAIRFIELD() {
}
void csExtendedHeader_FAIRFIELD::extractHeaders( byte const* buffer, int numBytes ) {
// Block 1:
  unitID                = cseis_segd::convert2uint64(&buffer[0]);   //   1-8 ID number of remote UINT  
  timeDeploy_epoch_us   = cseis_segd::convert2uint64(&buffer[8]);   //  9-16 Epoch deployment time
  timeRetrieve_epoch_us = cseis_segd::convert2uint64(&buffer[16]);  // 17-24 Epoch pickup time
  timeStart_epoch_us    = cseis_segd::convert2uint64(&buffer[24]);  // 25-32 Remote UINT Epoch start time

  /*
  for( int i = 0; i < 4; i+=1 ) {
    csInt64_t c1 = (csInt64_t)*(buffer+7+i*8);
    csInt64_t c2 = (csInt64_t)*(buffer+6+i*8)<<8UL;
    csInt64_t c3 = (csInt64_t)*(buffer+5+i*8)<<16UL;
    csInt64_t c4 = (csInt64_t)*(buffer+4+i*8)<<24UL;
    csInt64_t c5 = (csInt64_t)*(buffer+3+i*8)<<32UL;
    csInt64_t c6 = (csInt64_t)*(buffer+2+i*8)<<40UL;
    csInt64_t c7 = (csInt64_t)*(buffer+1+i*8)<<48UL;
    csInt64_t c8 = (csInt64_t)*(buffer+i*8)<<56UL;

    csInt64_t c10 = c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8; // THIS WORKED!!!!
    csInt64_t c11 =
      ((csInt64_t)*(buffer+7+i*8)) +
      ((csInt64_t)*(buffer+6+i*8)<<8UL) +
      ((csInt64_t)*(buffer+5+i*8)<<16UL) +
      ((csInt64_t)*(buffer+4+i*8)<<24UL) +
      ((csInt64_t)*(buffer+3+i*8)<<32UL) +
      ((csInt64_t)*(buffer+2+i*8)<<40UL) +
      ((csInt64_t)*(buffer+1+i*8)<<48UL) +
      ((csInt64_t)*(buffer+i*8)<<56UL );
    csInt64_t c12 = UINT64(&buffer[i*8]);

    fprintf(stdout,"  UNIT64   c10 c11 c12:  %lld  %lld  %lld  (%lld)\n", c10, c11, c12, c12-c11 );
  }
*/

// Block 2:
  int bytePos = 32;

// 1-4  IEEE Acquisition drift window [us]
  if( myIsLittleEndian ) {
    char buffer4[4];
    memcpy( buffer4, &buffer[bytePos+0], 4 );
    cseis_geolib::swapEndian4( buffer4, 4 );
    memcpy( &driftWindow_us, &buffer4, 4 );
  }
  else {
    memcpy( &driftWindow_us, &buffer[bytePos+0], 4 );
  }
  clockDrift_ns        = cseis_segd::convert2uint64(&buffer[bytePos+4]);   // 5-12 Clock drift [ns]
  clockStopMethod      = UINT8(&buffer[bytePos+12]);   //   13 Clock stop method (0:normal,1:storage full,2:power loss,3:reboot)
  freqDriftFlag        = UINT8(&buffer[bytePos+13]);   //   14 Frequency drift flag (0:not within specs, 1:within specs)
  oscillatorType       = UINT8(&buffer[bytePos+14]);   //   15 Oscillator type (0:control based,1:atomic,2:ovenized,3:double overnized,4:disciplined)
  dataCollectionMethod = UINT8(&buffer[bytePos+15]);   //   16 Data collection method (0:normal shots,1:continuous,2:shot sliced with guard bands)
  numRecordsInFile     = UINT32(&buffer[bytePos+16]);  // 17-20 Number of records (shots or time slices) in this file
  numFilesInUnit       = UINT32(&buffer[bytePos+20]);  // 21-24 Total number of files to be acquired for this remote UINT
  fileNum              = UINT32(&buffer[bytePos+24]);  // 25-28 File number (out of total above)
  dataDecimationFlag   = UINT8(&buffer[bytePos+28]);   //   29  Data decimation flag (0:not decimated,1:decimated)
  baseScanInterval     = UINT8(&buffer[bytePos+29]);   //   30  Original base scan interval
  numDecimationFiltCoef= UINT16(&buffer[bytePos+30]);  // 31-32 Number of decimation flter coefficients

// Block 3:
  bytePos += 32;
  rec_line   = UINT32(&buffer[bytePos+0]);   //   1-4 Receiver line number (replicated from other places)
  rec_point  = UINT32(&buffer[bytePos+4]);   //   5-8 Receiver point (replicated from other places)
  rec_index  = UINT8(&buffer[bytePos+8]);    //     9 Receiver point index (replicated)
  sou_line1  = UINT32(&buffer[bytePos+9]);   // 10-13 First shot line written to this file
  sou_point1 = UINT32(&buffer[bytePos+13]);  // 14-17 First shot point written to this file
  sou_index1 = UINT8(&buffer[bytePos+17]);   //    18 First shot point index written to this file
  sou_line2  = UINT32(&buffer[bytePos+18]);  // 19-22 Last shot line written to this file
  sou_point2 = UINT32(&buffer[bytePos+22]);  // 23-26 Last shot point written to this file
  sou_index2 = UINT8(&buffer[bytePos+26]);   //    27 Last shot point index written to this file
}

int csExtendedHeader_FAIRFIELD::nanoSeconds() const {
  return( 1000 * (int)( timeStart_epoch_us - 1000000UL * (csInt64_t)(timeStart_epoch_us/1000000UL) ) );
}
int csExtendedHeader_FAIRFIELD::numTraces() const {
  return( numRecordsInFile );
}

//----------------------------------------------------------------------

csTraceHeaderExtension_FAIRFIELD::csTraceHeaderExtension_FAIRFIELD( int numBlocks ) : csTraceHeaderExtension( numBlocks ) {
  myIsLittleEndian = cseis_geolib::isPlatformLittleEndian();
}
csTraceHeaderExtension_FAIRFIELD::~csTraceHeaderExtension_FAIRFIELD() {
}
void csTraceHeaderExtension_FAIRFIELD::extractHeaders( byte const* buffer, commonTraceHeaderStruct* comTrcHdr ) {
  extractHeaders( buffer );

  comTrcHdr->rcvLineNumber  = block1.rcvLineNumber;
  comTrcHdr->rcvPointNumber = block1.rcvPointNumber;
  comTrcHdr->rcvPointIndex  = block1.rcvPointIndex;
  comTrcHdr->sensor         = block1.sensorType;

  if( myNumBlocks < 5 ) return;
  comTrcHdr->rcvEasting   = (double)block5.rec_x_fin / 10.0;
  comTrcHdr->rcvNorthing  = (double)block5.rec_y_fin / 10.0;
  comTrcHdr->rcvElevation = (float)block5.rec_z_fin / 10.0F;

  if( myNumBlocks < 3 ) return;
  /*  
  comTrcHdr->serialNumber = (int)unitID;
  comTrcHdr->incl_i = block2.sensorInlineContrib;
  comTrcHdr->incl_c = block2.sensorXlineContrib;
  comTrcHdr->incl_v = block2.sensorVertContrib;
  */

  //  fprintf(stdout,"Serial number: %d %d %d\n", comTrcHdr->serialNumber, comTrcHdr->sensor,  comTrcHdr->incl_i );
}
int csTraceHeaderExtension_FAIRFIELD::nanoSeconds() const {
  fprintf(stdout,"Shottime:  %lld %lld\n", block3.shotTime_epoch_us, block3.shotSkewTime_us );
  csInt64_t time_us = block3.shotSkewTime_us + block3.shotTime_epoch_us;
  return( 1000 * (int)( time_us - 1000000UL * (csInt64_t)(time_us/1000000UL) ) );
}
csInt64_t csTraceHeaderExtension_FAIRFIELD::timeSamp1_us() const {
  return( block3.shotSkewTime_us + block3.shotTime_epoch_us );
}
//---------------------------------------
void csTraceHeaderExtension_FAIRFIELD::extractHeaders( byte const* buffer ) {
  memset( &block1, 0, sizeof(block1Struct) );
  memset( &block2, 0, sizeof(block2Struct) );
  memset( &block3, 0, sizeof(block3Struct) );
  memset( &block4, 0, sizeof(block4Struct) );
  memset( &block5, 0, sizeof(block5Struct) );
  memset( &block6, 0, sizeof(block6Struct) );
  memset( &block7, 0, sizeof(block7Struct) );
  memset( &block8, 0, sizeof(block8Struct) );
  memset( &block9, 0, sizeof(block9Struct) );
  memset( &block10, 0, sizeof(block10Struct) );

  // Block 1: Standard SEGD rev2 block
  block1.rcvLineNumber  = UINT24(&buffer[0]);
  block1.rcvPointNumber = UINT24(&buffer[3]);
  block1.rcvPointIndex  = UINT8(&buffer[6]);
  block1.numSamples     = UINT24(&buffer[7]);
  block1.extendedRcvLineNumber  = INT24(&buffer[10]);
  block1.extendedRcvLineFrac    = UINT24(&buffer[13]);
  block1.extendedRcvPointNumber = INT24(&buffer[15]);
  block1.extendedRcvPointFrac   = UINT24(&buffer[18]);
  block1.sensorType = UINT8(&buffer[20]);
  block1.traceCount = UINT32(&buffer[21]);

  if( myNumBlocks < 2 ) return;
  int bytePos = 32;

  // Block 2: "Shot" block
  block2.sou_line      = UINT32(&buffer[bytePos+0]);   // 1-4
  block2.sou_point     = UINT32(&buffer[bytePos+4]);   // 5-8
  block2.sou_index     = UINT8(&buffer[bytePos+8]);    // 9
  block2.sou_x_preplot = UINT32(&buffer[bytePos+9]);   // 10-13
  block2.sou_y_preplot = UINT32(&buffer[bytePos+13]);  // 14-17
  block2.sou_x_fin     = UINT32(&buffer[bytePos+17]);  // 18-21 x10
  block2.sou_y_fin     = UINT32(&buffer[bytePos+21]);  // 22-25 x10
  block2.sou_z_fin     = UINT32(&buffer[bytePos+25]);  // 26-29 x10
  block2.souInfoSrcFlag = UINT8(&buffer[bytePos+29]);  // 30 1:preplan, 2:as-shot, 3:post-processed
  block2.sou_type      = UINT8(&buffer[bytePos+30]);   // 31 Energy source type, 0:undefined,1:vibroseis,2:dynamite,3:airgun
  block2.sou_id        = UINT8(&buffer[bytePos+31]);   // 32 Source identifier

  if( myNumBlocks < 3 ) return;
  bytePos += 32;

  // Block 3: "Time" block
  block3.shotTime_epoch_us      = cseis_segd::convert2uint64(&buffer[bytePos+0]);  // 01-08 Int64
  block3.shotSkewTime_us        = cseis_segd::convert2uint64(&buffer[bytePos+8]);  // 09-16 Int64 Shot skew time from sample boundary
  block3.clockShiftApplied_ns   = cseis_segd::convert2uint64(&buffer[bytePos+16]); // 17-24 Applied clock correction time shift [ns]
  block3.clockShiftRemainder_ns = cseis_segd::convert2uint64(&buffer[bytePos+24]); // 25-32 Remaining (not applied) clock correction time shift [ns]

  if( myNumBlocks < 4 ) return;
  bytePos += 32;

// Block 4: "Data" block
  block4.timePreGuard_ms  = UINT32(&buffer[bytePos+0]);   // 1-4 Pre Shot Guard Band [ms]
  block4.timePostGuard_ms = UINT32(&buffer[bytePos+4]);   // 5-8 Post Shot Guard Band [ms]
  block4.gain_db          = UINT8(&buffer[bytePos+8]);    // 9 Preamp gain [dB]
  block4.traceClipFlag    = UINT8(&buffer[bytePos+9]);    // 10 0:not clipped, 1:digital clip, 2:analog clip
  block4.recordTypeCode   = UINT8(&buffer[bytePos+10]);   // 11
  block4.shotStatusFlag   = UINT8(&buffer[bytePos+11]);   // 12
  block4.external_sou_id  = UINT32(&buffer[bytePos+12]);  // 13-16
  block4.dataSampleFlag   = UINT8(&buffer[bytePos+16]);   // 17
  block4.fbptime          = UINT32(&buffer[bytePos+24]);  // 25-28
  block4.rmsNoise         = UINT32(&buffer[bytePos+28]);  // 29-32

  if( myNumBlocks < 5 ) return;
  bytePos += 32;

  // Block 5: "Receiver" block
  block5.rec_line      = UINT32(&buffer[bytePos+0]);   // 1-4  Receiver line number
  block5.rec_point     = UINT32(&buffer[bytePos+4]);   // 5-8  Receiver point
  block5.rec_index     = UINT8(&buffer[bytePos+8]);    // 9    Receiver poblock5.index
  block5.rec_x_preplot = UINT32(&buffer[bytePos+9]);   // 10-13 x10
  block5.rec_y_preplot = UINT32(&buffer[bytePos+13]);  // 14-17 x10
  block5.rec_x_fin     = UINT32(&buffer[bytePos+17]);  // 18-21 x10
  block5.rec_y_fin     = UINT32(&buffer[bytePos+21]);  // 22-25 x10
  block5.rec_z_fin     = UINT32(&buffer[bytePos+25]);  // 26-29 x10
//  float rec_z_fin;              // 26-29 x10
  block5.recInfoSourceFlag = UINT32(&buffer[bytePos+29]); // 30 1:preplan, 2:as-laid, 3:...
  /*
  char buffer4[4];
  memcpy( buffer4, &buffer[bytePos+25], 4 );
  if( myIsLittleEndian ) {
    cseis_geolib::swapEndian( buffer4, 4, 4 );
  }
  memcpy( &block5.rec_z_fin, buffer4, 4 );  // 1-4 IEEE Float
  */
  if( myNumBlocks < 6 ) return;
  bytePos += 32;

  // Block 6: "Orientation" block
  char buffer32[32];
  memcpy( buffer32, &buffer[bytePos+0], 32 );
  if( myIsLittleEndian ) {
    cseis_geolib::swapEndian4( buffer32, 32 );
  }
  memcpy( &block6.tilt_h1x, &buffer32[0], 4 );  // 1-4 IEEE Float
  memcpy( &block6.tilt_h2x, &buffer32[4], 4 );  // 5-8 IEEE Float
  memcpy( &block6.tilt_vx,  &buffer32[8], 4 );  // 9-12 IEEE Float
  memcpy( &block6.tilt_h1y, &buffer32[12], 4 ); // 13-16 IEEE Float
  memcpy( &block6.tilt_h2y, &buffer32[16], 4 ); // 17-20 IEEE Float
  memcpy( &block6.tilt_vy,  &buffer32[20], 4 ); // 21-24 IEEE Float
  memcpy( &block6.tilt_h1z, &buffer32[24], 4 ); // 25-28 IEEE Float
  memcpy( &block6.tilt_h2z, &buffer32[28], 4 ); // 29-32 IEEE Float

  if( myNumBlocks < 7 ) return;
  bytePos += 32;

  memcpy( buffer32, &buffer[bytePos+0], 24 );
  if( myIsLittleEndian ) {
    cseis_geolib::swapEndian4( buffer32, 24 );
  }
  memcpy( &block7.tilt_vz,     &buffer32[0], 4 );  // 1-4 IEEE Float
  memcpy( &block7.azim_deg,    &buffer32[4], 4 );  // 5-8 IEEE Float
  memcpy( &block7.pitch_deg,   &buffer32[8], 4 );  // 9-12 IEEE Float
  memcpy( &block7.roll_deg,    &buffer32[12], 4 ); // 13-16 IEEE Float
  memcpy( &block7.temperature, &buffer32[16], 4 ); // 17-20 IEEE Float
  memcpy( &block7.humidity,    &buffer32[20], 4 ); // 21-24 IEEE Float
  block7.orientMatrixVer  = UINT32(&buffer[bytePos+24]); // 25-28 Binary
  block7.gimbalCorrection = UINT8(&buffer[bytePos+28]);  // 29 Binary
  block7.headingRotation  = UINT16(&buffer[bytePos+29]); // 30-31 Binary Heading rotation angle (tenths of deg)

  if( myNumBlocks < 8 ) return;
  bytePos += 32;
  // Block 8 not implemented yet

  if( myNumBlocks < 9 ) return;
  bytePos += 32;
  // Block 9 not implemented yet

  if( myNumBlocks < 10 ) return;
  bytePos += 32;
  // Block 10 not implemented yet
}

namespace cseis_segd {
void csExtendedHeader_FAIRFIELD::dump( std::ostream& os )
{
  os <<
    header("ExtendedHeader start") << '\n' <<
// Block 1:
  "ID number of remote unit        : " << unitID << endl << 
  "Epoch deployment time [us]      : " << timeDeploy_epoch_us << endl <<
  "Epoch pickup time [us]          : " << timeRetrieve_epoch_us << endl <<
  "Remote unit Epoch start time[us]: " << timeStart_epoch_us << endl <<
  endl <<

// Block 2:
  "Acquisition drift window [us]                       : " << driftWindow_us << endl <<
  "Clock drift [ns]                                    : " << clockDrift_ns << endl <<
  "Clock stop method (0:normal,1:storage full,2:power loss,3:reboot)           : " << clockStopMethod << endl <<       //   13 
  "Frequency drift flag (0:not within specs, 1:within specs)                   : " << freqDriftFlag << endl <<
  "Oscillator type (0:control based,1:atomic,2:oven,3:dbl oven,4:disciplined)  : " << oscillatorType << endl <<
  "Data collection method (0:shots,1:continuous,2:shot sliced with guard bands): " << dataCollectionMethod << endl <<
  "Number of records (shots or time slices) in file    : " << numRecordsInFile << endl <<
  "Total number of files to be acquired for remote unit: " << numFilesInUnit << endl <<
  "File number (out of total above)                    : " << fileNum << endl <<
  "Data decimation flag (0:not decimated,1:decimated)  : " << dataDecimationFlag << endl <<
  "Original base scan intervalbaseScanInterval         : " << endl <<
  "Number of decimation flter coefficients             : " << numDecimationFiltCoef << endl <<
  endl <<

// Block 3:
  "Receiver line number (replicated)     : " << rec_line << endl <<
  "Receiver point (replicated)           : " << rec_point << endl <<
  "Receiver point index (replicated)     : " << rec_index << endl <<
  "First shot line written to file       : " << sou_line1 << endl <<
  "First shot point written to file      : " << sou_point1 << endl << 
  "First shot point index written to file: " << sou_index1 << endl <<
  "Last shot line written to file        : " << sou_line2 << endl <<
  "Last shot point written to file       : " << sou_point2 << endl <<
  "Last shot point index written to file : " << sou_index2 << endl <<
  header("ExtendedHeader end  ");
}

void csTraceHeaderExtension_FAIRFIELD::dump( std::ostream& os )
{
  os <<
    header("TraceHeaderExtension start") << endl;

  os << " Block 1 - Standard block: " << endl <<
    "rcvLineNumber          : " << block1.rcvLineNumber << endl <<
    "rcvPointNumber         : " << block1.rcvPointNumber << endl <<
    "rcvPointIndex          : " << block1.rcvPointIndex << endl <<
    "numSamples             : " << block1.numSamples << endl <<
    "extendedRcvLineNumber  : " << block1.extendedRcvLineNumber << endl <<
    "extendedRcvLineFrac    : " << block1.extendedRcvLineFrac << endl <<
    "extendedRcvPointNumber : " << block1.extendedRcvPointNumber << endl <<
    "extendedRcvPointFrac   : " << block1.extendedRcvPointFrac << endl <<
    "sensorType             : " << block1.sensorType << endl <<
    "traceCount             : " << block1.traceCount << endl;

  if( myNumBlocks < 2 ) goto endofdump;
  os << " Block 2 - Shot block: " << endl <<
    "sou_line               : " << block2.sou_line << endl <<
    "sou_point              : " << block2.sou_point << endl <<
    "sou_index              : " << block2.sou_index << endl <<
    "sou_x_preplot          : " << block2.sou_x_preplot << endl <<
    "sou_y_preplot          : " << block2.sou_y_preplot << endl <<
    "sou_x_fin              : " << block2.sou_x_fin << endl <<
    "sou_y_fin              : " << block2.sou_y_fin << endl <<
    "sou_z_fin              : " << block2.sou_z_fin << endl <<
    "souInfoSrcFlag         : " << block2.souInfoSrcFlag << endl <<
    "sou_type               : " << block2.sou_type << endl <<
    "sou_id                 : " << block2.sou_id << endl;

  if( myNumBlocks < 3 ) goto endofdump;
  os << " Block 3 - Time block: " << endl <<
    "shotTime_epoch_us      : " << block3.shotTime_epoch_us << endl <<
    "shotSkewTime_us        : " << block3.shotSkewTime_us << endl <<
    "clockShiftApplied_ns   : " << block3.clockShiftApplied_ns << endl <<
    "clockShiftRemainder_ns : " << block3.clockShiftRemainder_ns << endl;

  if( myNumBlocks < 4 ) goto endofdump;
  os << " Block 4 - Data block: " << endl <<
    "timePreGuard_ms       : " << block4.timePreGuard_ms << endl <<
    "timePostGuard_ms      : " << block4.timePostGuard_ms << endl <<
    "gain_db               : " << block4.gain_db << endl <<
    "traceClipFlag         : " << block4.traceClipFlag << endl <<
    "recordTypeCode        : " << block4.recordTypeCode << endl <<
    "shotStatusFlag        : " << block4.shotStatusFlag << endl <<
    "external_sou_id       : " << block4.external_sou_id << endl <<
    "dataSampleFlag        : " << block4.dataSampleFlag << endl <<
    "fbptime               : " << block4.fbptime << endl <<
    "rmsNoise              : " << block4.rmsNoise << endl;

  if( myNumBlocks < 5 ) goto endofdump;
  os << " Block 5 - Receiver block: " << endl <<
    "rec_line              : " << block5.rec_line << endl <<
    "rec_point             : " << block5.rec_point      << endl <<
    "rec_index             : " << block5.rec_index      << endl <<
    "rec_x_preplot         : " << block5.rec_x_preplot  << endl <<
    "rec_y_preplot         : " << block5.rec_y_preplot  << endl <<
    "rec_x_fin             : " << block5.rec_x_fin      << endl <<
    "rec_y_fin             : " << block5.rec_y_fin      << endl <<
    "rec_z_fin             : " << block5.rec_z_fin      << endl <<
    "recInfoSourceFlag     : " << block5.recInfoSourceFlag << endl;

  if( myNumBlocks < 6 ) goto endofdump;
  os << " Block 6 - Orientation block: " << endl <<
    "tilt_h1x         : " << block6.tilt_h1x << endl <<
    "tilt_h2x         : " << block6.tilt_h2x << endl <<
    "tilt_vx          : " << block6.tilt_vx << endl <<
    "tilt_h1y         : " << block6.tilt_h1y << endl <<
    "tilt_h2y         : " << block6.tilt_h2y << endl <<
    "tilt_vy          : " << block6.tilt_vy << endl <<
    "tilt_h1z         : " << block6.tilt_h1z << endl <<
    "tilt_h2z         : " << block6.tilt_h2z << endl;
  
  if( myNumBlocks < 7 ) goto endofdump;
  os << " Block 7 - Orientation block: " << endl <<
    "tilt_vz          : " << block7.tilt_vz << endl <<
    "azim_deg         : " << block7.azim_deg << endl <<
    "pitch_deg        : " << block7.pitch_deg << endl <<
    "roll_deg         : " << block7.roll_deg << endl <<
    "temperature      : " << block7.temperature << endl <<
    "humidity         : " << block7.humidity << endl <<
    "orientMatrixVer  : " << block7.orientMatrixVer << endl <<
    "gimbalCorrection : " << block7.gimbalCorrection << endl <<
    "headingRotation  : " << block7.headingRotation << endl;

endofdump:
    os << header("TraceHeaderExtension end  ") << endl;

}
} // end namespace


