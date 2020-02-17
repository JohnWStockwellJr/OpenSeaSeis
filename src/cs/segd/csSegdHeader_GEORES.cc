/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegdHeader_GEORES.h"
#include "geolib_endian.h"
#include <ostream>
#include <cmath>
#include <cstring>

using std::endl;
using namespace cseis_segd;
using namespace std;

csExtendedHeader_GEORES::csExtendedHeader_GEORES() : csExtendedHeader() {
  myIsLittleEndian = cseis_geolib::isPlatformLittleEndian();
}
csExtendedHeader_GEORES::~csExtendedHeader_GEORES() {
}

void csExtendedHeader_GEORES::extractHeaders( byte const* buffer, int numBytes ) {
  segMfgCode = UINT8(&buffer[0]);          // 1:   SEG manufacturer code = 9
  sponsorshipCode = UINT16(&buffer[1]);    // 2-3: 0-255: Geospace, 256-...: other
  contentType = UINT8(&buffer[3]);         // 4:   Unique code identifying inf in this header module
  contentVersion = UINT8(&buffer[4]);      // 5:   Version of header content type
  contentCompatibilityVersion = UINT8(&buffer[5]); // 6: ...
  totalCount32Blocks = UINT16(&buffer[6]);  // 7-8: Number of consecutive 32-byte blocks in this header module
  
  fracAcqNanoSeconds = UINT32(&buffer[8]);   // 9-12:  Additional number of nanoseconds past second mark (in general header 1)
  shotSerialNumber = UINT32(&buffer[12]);    // 13-16: Serial number assigned to each shot
  recordIndex = UINT8(&buffer[16]);          // 17:    0-255. 
  recordSerialNumber = UINT12H(&buffer[17]); // 18-19 high: 1-4095. Unique shot record index for each shot serial number
  triggerMethod = UINT4L(&buffer[18]);       // 19 low nibble, 0-15
  sampleClockSynch = UINT4H(&buffer[19]);    // 20 high nibble, 0-15
  sampleClockSynchStatus = UINT4L(&buffer[19]); // 20 low nibble, 0-15
  sampleClockSynchErr = UINT24(&buffer[20]); // 21-23: Synchronization phase error [us]
  numAuxChan = UINT8(&buffer[23]);           // 24:    0-255
  numReceiverChan = UINT16(&buffer[24]);     // 25-26: 0-65535
  recordTimeOffset = INT16(&buffer[26]);     // 27-28: Amount of milliseconds from trigger event [ms]
  sourceType = UINT8(&buffer[28]);           // 29
  extTraceHeaderVersion = UINT8(&buffer[29]);// 30
  extTraceHeaderCompVersion = UINT8(&buffer[30]); // 31
  numExtTraceHeaderBlocks = UINT8(&buffer[31]); // 32:   Number of extended trace header blocks per trace

// Block 2:
  int bytePos = 32;
  segMfgCode2 = UINT8(&buffer[bytePos+0]);       // 1:   SEG manufacturer code = 9
  sponsorshipCode2 = UINT16(&buffer[bytePos+1]); // 2-3: 0-255: Geospace, 256-...: other
  contentType2 = UINT8(&buffer[bytePos+3]);      // 4:   Unique code identifying inf in this header module
  contentVersion2 = UINT8(&buffer[bytePos+4]);   // 5:   Version of header content type
  contentCompatibilityVersion2 = UINT8(&buffer[bytePos+5]); // 6: ...
  totalCount32Blocks2 = UINT16(&buffer[bytePos+6]);// 7-8: Number of consecutive 32-byte blocks in this header module

  if( myIsLittleEndian ) {
    char buffer8[8];
    memcpy( buffer8, &buffer[bytePos+8], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &srcEasting, buffer8, 8 );

    memcpy( buffer8, &buffer[bytePos+16], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &srcNorthing, buffer8, 8 );

    memcpy( buffer8, &buffer[bytePos+24], 4 );
    cseis_geolib::swapEndian4( buffer8, 4 );
    memcpy( &srcElev, &buffer8, 4 );
  }
  else {
    memcpy( &srcEasting, &buffer[bytePos+8], 8 );
    memcpy( &srcNorthing, &buffer[bytePos+16], 8 );
    memcpy( &srcElev, &buffer[bytePos+24], 4 );
  }

// Block 3:
  bytePos = 32+32;
  segMfgCode3 = UINT8(&buffer[bytePos+0]);       // 1:   SEG manufacturer code = 9
  sponsorshipCode3 = UINT16(&buffer[bytePos+1]);  // 2-3: 0-255: Geospace, 256-...: other
  contentType3 = UINT8(&buffer[bytePos+3]);      // 4:   Unique code identifying inf in this header module
  contentVersion3 = UINT8(&buffer[bytePos+4]);   // 5:   Version of header content type
  contentCompatibilityVersion3 = UINT8(&buffer[bytePos+5]); // 6: ...
  totalCount32Blocks3 = UINT16(&buffer[bytePos+6]);// 7-8: Number of consecutive 32-byte blocks in this header module

  fireMethod = UINT4H(&buffer[bytePos+8]);
  seqStartMethod = UINT4L(&buffer[bytePos+8]);
  timebreakMethod = UINT4H(&buffer[bytePos+9]);
  armDelay = UINT16(&buffer[bytePos+10]);
  fireDelay = UINT16(&buffer[bytePos+12]);
  timebreakWindowDelay = UINT16(&buffer[bytePos+14]);
  timebreakWindowDuration = UINT16(&buffer[bytePos+16]);
  startToTimebreak = UINT8(&buffer[bytePos+23]);
  startToTimebreakFrac = UINT32(&buffer[bytePos+24]);
  timebreakToSampleClock = UINT32(&buffer[bytePos+28]);
}
//----------------------------------------------------------------------

csTraceHeaderExtension_GEORES::csTraceHeaderExtension_GEORES( int numBlocks ) : csTraceHeaderExtension( numBlocks ) {
  myIsLittleEndian = cseis_geolib::isPlatformLittleEndian();
}
csTraceHeaderExtension_GEORES::~csTraceHeaderExtension_GEORES() {
}
void csTraceHeaderExtension_GEORES::extractHeaders( byte const* buffer, commonTraceHeaderStruct* comTrcHdr ) {
  int blockSize = 32;
  //    fprintf(stdout,"Extract trace header extension 1...\n");

  comTrcHdr->rcvLineNumber  = UINT24(&buffer[0]);
  comTrcHdr->rcvPointNumber = UINT24(&buffer[3]);
  comTrcHdr->rcvPointIndex  = UINT8(&buffer[6]);
  comTrcHdr->sensor         = block1.sensorType;

  if( myNumBlocks < 2 ) return;
  if( myIsLittleEndian ) {
    char buffer8[8];
    memcpy( buffer8, &buffer[0+blockSize], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &comTrcHdr->rcvEasting, buffer8, 8 );

    memcpy( buffer8, &buffer[8+blockSize], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &comTrcHdr->rcvNorthing, buffer8, 8 );

    memcpy( &buffer8, &buffer[16+blockSize], 4 );
    cseis_geolib::swapEndian4( buffer8, 4 );
    memcpy( &comTrcHdr->rcvElevation, buffer8, 4 );
  }
  else {
    memcpy( &comTrcHdr->rcvEasting, &buffer[0+blockSize], 8 );
    memcpy( &comTrcHdr->rcvNorthing, &buffer[8+blockSize], 8 );
    memcpy( &comTrcHdr->rcvElevation, &buffer[16+blockSize], 4 );
  }

  if( myNumBlocks < 3 ) return;
  comTrcHdr->serialNumber = UINT24(&buffer[0+blockSize*2]);
  comTrcHdr->incl_i = block2.sensorInlineContrib;
  comTrcHdr->incl_c = block2.sensorXlineContrib;
  comTrcHdr->incl_v = block2.sensorVertContrib;

  //  fprintf(stdout,"Serial number: %d %d %d\n", comTrcHdr->serialNumber, comTrcHdr->sensor,  comTrcHdr->incl_i );
}
//---------------------------------------
void csTraceHeaderExtension_GEORES::extractHeaders( byte const* buffer ) {
  int blockSize = 32;
  //  fprintf(stdout,"Extract trace header extension 2...\n");

  memset( &block2, 0, sizeof(block2Struct) );
  memset( &block3, 0, sizeof(block3Struct) );

  block1.rcvLineNumber  = UINT24(&buffer[0]);
  block1.rcvPointNumber = UINT24(&buffer[3]);
  block1.rcvPointIndex  = UINT8(&buffer[6]);
  block1.numSamples     = UINT24(&buffer[7]);
  block1.extendedRcvLineNumber  = INT24(&buffer[10]);
  block1.extendedRcvLineFrac    = UINT24(&buffer[13]);
  block1.extendedRcvPointNumber = INT24(&buffer[15]);
  block1.extendedRcvPointFrac   = UINT24(&buffer[18]);
  block1.sensorType = UINT8(&buffer[20]);

  if( myNumBlocks < 2 ) return;


  if( myIsLittleEndian ) {
    char buffer8[8];
    memcpy( buffer8, &buffer[0+blockSize], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &block2.rcvEasting, buffer8, 8 );

    memcpy( buffer8, &buffer[8+blockSize], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &block2.rcvNorthing, buffer8, 8 );

    memcpy( &buffer8, &buffer[16+blockSize], 4 );
    cseis_geolib::swapEndian4( buffer8, 4 );
    memcpy( &block2.rcvElevation, buffer8, 4 );
  }
  else {
    memcpy( &block2.rcvEasting, &buffer[0+blockSize], 8 );
    memcpy( &block2.rcvNorthing, &buffer[8+blockSize], 8 );
    memcpy( &block2.rcvElevation, &buffer[16+blockSize], 4 );
  }

  block2.recordIndex       = UINT8(&buffer[20+blockSize]);
  block2.traceDataProcessDesc = UINT8(&buffer[21+blockSize]);
  block2.secLowCutFilter   = UINT8(&buffer[22+blockSize]);
  block2.analogGainApplied = INT12H(&buffer[23+blockSize]);
  block2.sensorOrientDesc  = UINT4L(&buffer[24+blockSize]);

  //  fprintf(stdout,"Sensor orient desc: %d\n", block2.sensorOrientDesc);

  if( block2.sensorOrientDesc == 0 ) {
    // Nothing...
  }
  else if( block2.sensorOrientDesc == 1 ) {
    block2.sensorAzim       = INT12H(&buffer[26+blockSize]);
    block2.sensorIncl       = INT12L(&buffer[27+blockSize]);
    block2.sensorAzimDev    = INT12H(&buffer[29+blockSize]);
    block2.sensorInclDev    = INT12L(&buffer[30+blockSize]);
  }
  else if( block2.sensorOrientDesc == 2 ) {
    block2.sineOfSensorIncl = INT16(&buffer[26+blockSize]);
  }
  else if( block2.sensorOrientDesc == 3 ) {
    block2.sensorVertContrib   = INT16(&buffer[26+blockSize]);
    block2.sensorInlineContrib = INT16(&buffer[28+blockSize]);
    block2.sensorXlineContrib  = INT16(&buffer[30+blockSize]);
  }
  else if( block2.sensorOrientDesc == 4 ) {
    block2.sensorVertContrib   = INT16(&buffer[26+blockSize]);
    block2.sensorNorthContrib = INT16(&buffer[28+blockSize]);
    block2.sensorEastContrib  = INT16(&buffer[30+blockSize]);
  }

  if( myNumBlocks < 3 ) return;
  block3.moduleSerialNumber = UINT24(&buffer[0+blockSize*2]);
  block3.appliedSensorCalib = UINT12L(&buffer[3+blockSize*2]);

  //  fprintf(stdout,"Module serial number: %d\n", block3.moduleSerialNumber);

}

namespace cseis_segd {
void csExtendedHeader_GEORES::dump( std::ostream& os )
{
  os <<
    header("ExtendedHeader start") << '\n' <<
  "SEG manufacturer code = 9            : " << segMfgCode << endl <<
  "sponsorshipCode                      : " << sponsorshipCode << endl <<
  "Unique code identifying header       : " << contentType << endl <<
  "Version of header content type       : " << contentVersion << endl <<
  "contentCompatibilityVersion          : " << contentCompatibilityVersion << endl <<
  "Number of consecutive 32-byte blocks : " << totalCount32Blocks << endl <<
  "Additional nanoseconds past second   : " << fracAcqNanoSeconds << endl <<
  "Serial number assigned to each shot  : " << shotSerialNumber << endl <<
  "recordIndex                          : " << recordIndex << endl <<
  "Unique shot record serial number     : " << recordSerialNumber << endl <<
  "triggerMethod                        : " << triggerMethod << endl <<
  "sampleClockSynch                     : " << sampleClockSynch << endl <<
  "sampleClockSynchStatus               : " << sampleClockSynchStatus << endl <<
  "Synchronization phase error [us]     : " << sampleClockSynchErr << endl <<
  "numAuxChan                           : " << numAuxChan << endl <<
  "numReceiverChan                      : " << numReceiverChan << endl <<
  "Milliseconds from trigger event [ms] : " << recordTimeOffset << endl <<
  "sourceType                           : " << sourceType << endl <<
  "extTraceHeaderVersion                : " << extTraceHeaderVersion << endl <<
  "extTraceHeaderCompVersion            : " << extTraceHeaderCompVersion << endl <<
  "Number of extended trace header blocks: " << numExtTraceHeaderBlocks << endl <<
  endl <<

  "SEG manufacturer code = 9            : " << segMfgCode2 << endl <<
  "sponsorshipCode                      : " << sponsorshipCode2 << endl <<
  "Unique code identifying header       : " << contentType2 << endl <<
  "Version of header content type       : " << contentVersion2 << endl <<
  "contentCompatibilityVersion          : " << contentCompatibilityVersion2 << endl <<
  "Number of consecutive 32-byte blocks : " << totalCount32Blocks2 << endl <<
  "srcEasting                           : " << srcEasting << endl <<
  "srcNorthing                          : " << srcNorthing << endl <<
  "srcElev                              : " << srcElev << endl <<
  endl <<

  "SEG manufacturer code = 9            : " << segMfgCode3 << endl <<
  "sponsorshipCode                      : " << sponsorshipCode3 << endl <<
  "Unique code identifying header       : " << contentType3 << endl <<
  "Version of header content type       : " << contentVersion3 << endl <<
  "contentCompatibilityVersion          : " << contentCompatibilityVersion3 << endl <<
  "Number of consecutive 32-byte blocks : " << totalCount32Blocks3 << endl <<
  
  "fireMethod                           : " << fireMethod << endl <<
  "seqStartMethod                       : " << seqStartMethod << endl <<
  "timebreakMethod                      : " << timebreakMethod << endl <<
  "armDelay                             : " << armDelay << endl <<
  "fireDelay                            : " << fireDelay << endl <<
  "timebreakWindowDelay                 : " << timebreakWindowDelay << endl <<
  "timebreakWindowDuration              : " << timebreakWindowDuration << endl <<
  "startToTimebreak                     : " << startToTimebreak << endl <<
  "startToTimebreakFrac                 : " << startToTimebreakFrac << endl <<
  "timebreakToSampleClock               : " << timebreakToSampleClock << endl <<

    header("ExtendedHeader end  ");
}

void csTraceHeaderExtension_GEORES::dump( std::ostream& os )
{
  os <<
    header("TraceHeaderExtension start") << endl;

    os <<
      "rcvLineNumber          : " << block1.rcvLineNumber << endl <<
      "rcvPointNumber         : " << block1.rcvPointNumber << endl <<
      "rcvPointIndex          : " << block1.rcvPointIndex << endl <<
      "numSamples             : " << block1.numSamples << endl <<
      "extendedRcvLineNumber  : " << block1.extendedRcvLineNumber << endl <<
      "extendedRcvLineFrac    : " << block1.extendedRcvLineFrac << endl <<
      "extendedRcvPointNumber : " << block1.extendedRcvPointNumber << endl <<
      "extendedRcvPointFrac   : " << block1.extendedRcvPointFrac << endl <<
      "sensorType             : " << block1.sensorType << endl;
      if( myNumBlocks < 2 ) goto endofdump;

    os <<
      "rcvEasting          : " << block2.rcvEasting << endl <<
      "rcvNorthing         : " << block2.rcvNorthing << endl <<
      "rcvElevation        : " << block2.rcvElevation << endl <<
      "recordIndex         : " << block2.recordIndex << endl <<
      "traceDataProcessDesc: " << block2.traceDataProcessDesc << endl <<
      "secLowCutFilter     : " << block2.secLowCutFilter << endl <<
      "analogGainApplied   : " << block2.analogGainApplied << endl <<
      "sensorOrientDesc    : " << block2.sensorOrientDesc << endl;

    if( block2.sensorOrientDesc == 1 ) {
      os <<
      "sensorAzim          : " << block2.sensorAzim << endl <<
      "sensorIncl          : " << block2.sensorIncl << endl <<
      "sensorAzimDev       : " << block2.sensorAzimDev << endl <<
      "sensorInclDev       : " << block2.sensorInclDev << endl;
    }
    else if( block2.sensorOrientDesc == 2 ) {
      os <<
      "sineOfSensorIncl    : " << block2.sineOfSensorIncl << endl;
    }
    else if( block2.sensorOrientDesc == 3 ) {
      os <<
      "sensorVertContrib   : " << block2.sensorVertContrib << endl <<
      "sensorInlineContrib : " << block2.sensorInlineContrib << endl <<
      "sensorXlineContrib  : " << block2.sensorXlineContrib << endl;
    }
    else if( block2.sensorOrientDesc == 4 ) {    
      os <<
      "sensorVertContrib   : " << block2.sensorVertContrib << endl <<
      "sensorNorthContrib  : " << block2.sensorNorthContrib << endl <<
      "sensorEastContrib   : " << block2.sensorEastContrib << endl;
    }

    if( myNumBlocks < 3 ) goto endofdump;
    os <<
      "moduleSerialNumber  : " << block3.moduleSerialNumber << endl <<
      "appliedSensorCalib  : " << block3.appliedSensorCalib << endl;

endofdump:
    os << header("TraceHeaderExtension end  ") << endl;

}
} // end namespace


