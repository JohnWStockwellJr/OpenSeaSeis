/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegdHeader_SEAL.h"
#include "geolib_endian.h"
#include <ostream>
#include <cstring>

using namespace std;
using namespace cseis_segd;

csExtendedHeader_SEAL::csExtendedHeader_SEAL() : csExtendedHeader() {
  myIsLittleEndian = cseis_geolib::isPlatformLittleEndian();
}
csExtendedHeader_SEAL::~csExtendedHeader_SEAL() {
}
void csExtendedHeader_SEAL::extractHeaders( byte const* buffer, int numBytes ) {
  acqLength = UINT32(&buffer[0]);
  sampleInt = UINT32(&buffer[4]);
  totalNumTraces = UINT32(&buffer[8]);
  numAuxTraces = UINT32(&buffer[12]);
  numSeismicTraces = UINT32(&buffer[16]);
  numDeadSeismicTraces = UINT32(&buffer[20]);
  numLiveSeismicTraces = UINT32(&buffer[24]);
  sourceType = UINT32(&buffer[28]);
  numSamples = UINT32(&buffer[32]);
  shotNumber = UINT32(&buffer[36]);
  memcpy( &TBwindow, &buffer[40], 4 );
  testRecordType = UINT32(&buffer[44]);
  swathFirstLine = UINT32(&buffer[48]);
  swathFirstNum = UINT32(&buffer[52]);
  spreadNumber = UINT32(&buffer[56]);
  spreadType = UINT32(&buffer[60]);
  timeBreak = UINT32(&buffer[64]);
  upholeTime = UINT32(&buffer[68]);
  blasterID = UINT32(&buffer[72]);
  blasterStatus = UINT32(&buffer[76]);
  refractionDelay = UINT32(&buffer[80]);

  TBtoT0time = INT32(&buffer[84]);

  internalTimeBreak = UINT32(&buffer[88]);
  prestackWithinFieldUnits = UINT32(&buffer[92]);
  noiseElimType = UINT32(&buffer[96]);
  lowTracePerc = UINT32(&buffer[100]);
  lowTraceValue = UINT32(&buffer[104]);
  numWindows = UINT32(&buffer[108]);
  historicEditType = UINT32(&buffer[112]);
  noisyTracePerc = UINT32(&buffer[116]);

  historicRange               = UINT32(&buffer[120]);
  historicTaperLengthExponent = UINT32(&buffer[124]);
  thresholdVar                = UINT32(&buffer[128]);
  historicThresholdInitValue  = UINT32(&buffer[132]);
  historicZeroingLength       = UINT32(&buffer[136]);
  typeOfProcess               = UINT32(&buffer[140]);
  for( int i = 0; i < 32; i++ ) {
    acquisitionTypeTables[i]  = UINT32(&buffer[144+i]);
  }
  for( int i = 0; i < 32; i++ ) {
    thresholdTypeTables[i]    = UINT32(&buffer[272+i]);
  }
  stackFold                   = UINT32(&buffer[400]);   // 1 to 32
  memcpy( auxCorrelationDesc, &buffer[404], 80 );
  investigationLength         = UINT32(&buffer[484]);    // 100-128000ms
  autocorrelationPeakTime     = UINT32(&buffer[488]);
  fftExponent                 = UINT32(&buffer[492]);
  sourceAuxNB                 = UINT32(&buffer[497]);
  pilotLength                 = UINT32(&buffer[500]);
  sweepLength                 = UINT32(&buffer[504]);
  acquisitionNum              = UINT32(&buffer[508]);
  memcpy( &maxOfMaxAux, &buffer[512], 4 );
  memcpy( &maxOfMaxSeis, &buffer[516], 4 );
  dumpStackFold               = UINT32(&buffer[520]);

  memcpy( tapeLabel, &buffer[524], 16 );
  tapeNum = UINT32(&buffer[540]);
  memcpy( softwareRevision, &buffer[544], 16 );
  memcpy( date, &buffer[560], 16 );

  if( myIsLittleEndian ) {
    char buffer8[8];
    memcpy( buffer8, &buffer[572], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &srcEasting, buffer8, 8 );

    memcpy( buffer8, &buffer[580], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &srcNorthing, buffer8, 8 );

    memcpy( buffer8, &buffer[589], 4 );
    cseis_geolib::swapEndian4( buffer8, 4 );
    memcpy( &srcElev, &buffer8, 4 );
  }
  else {
    memcpy( &srcEasting, &buffer[572], 8 );
    memcpy( &srcNorthing, &buffer[580], 8 );
    memcpy( &srcElev, &buffer[589], 4 );
  }


  filesPerTape = UINT32(&buffer[596]);
  fileCount = UINT32(&buffer[600]);
  memcpy( &acquisitionErrorDesc, &buffer[604], 160 );
  filterType = UINT32(&buffer[764]);
  profileNum = UINT32(&buffer[768]);
  t0mode = UINT32(&buffer[780]);
  memcpy( &srcLineName, &buffer[784], 16 );
  seqNumber = UINT32(&buffer[800]);
  navStatus = UINT32(&buffer[804]);
  shotsToGo = UINT32(&buffer[808]);
}
//----------------------------------------------------------------------

csTraceHeaderExtension_SEAL::csTraceHeaderExtension_SEAL( int numBlocks ) : csTraceHeaderExtension( numBlocks ) {
  myIsLittleEndian = cseis_geolib::isPlatformLittleEndian();
}
csTraceHeaderExtension_SEAL::~csTraceHeaderExtension_SEAL() {
}
void csTraceHeaderExtension_SEAL::extractHeaders( byte const* buffer, commonTraceHeaderStruct* comTrcHdr ) {
  int blockSize = 32;

  comTrcHdr->rcvLineNumber  = UINT24(&buffer[0]);
  comTrcHdr->rcvPointNumber = UINT24(&buffer[3]);
  comTrcHdr->rcvPointIndex  = UINT8(&buffer[6]);

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

  if( myNumBlocks < 6 ) return;
  comTrcHdr->serialNumber = UINT24(&buffer[1+blockSize*5]);
}
//--------------------------------------------------------
//
void csTraceHeaderExtension_SEAL::extractHeaders( byte const* buffer ) {
  int blockSize = 32;
  memset( &block2, 0, sizeof(block2Struct) );
  memset( &block3, 0, sizeof(block3Struct) );
  memset( &block4, 0, sizeof(block4Struct) );
  memset( &block5, 0, sizeof(block5Struct) );
  memset( &block6, 0, sizeof(block6Struct) );
  memset( &block7, 0, sizeof(block7Struct) );

  block1.rcvLineNumber  = UINT24(&buffer[0]);
  block1.rcvPointNumber = UINT24(&buffer[3]);
  block1.rcvPointIndex  = UINT8(&buffer[6]);
  block1.numSamples     = UINT24(&buffer[7]);
  block1.extendedRcvLineNumber  = UINT32(&buffer[10]);
  block1.extendedRcvPointNumber = UINT32(&buffer[15]);
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
  block2.sensorTypeNB = UINT8(&buffer[20+blockSize]);

  if( myNumBlocks < 3 ) return;
  char buffer128[128];
  memcpy( buffer128, &buffer[0+blockSize*2], 20 );
  if( myIsLittleEndian ) {
    cseis_geolib::swapEndian( buffer128, 4, 20 );
  }
  memcpy( &block3.resistanceLimitLow, &buffer128[0], 4 );
  memcpy( &block3.resistanceLimitHigh, &buffer128[4], 4 );
  memcpy( &block3.resistanceValue, &buffer128[8], 4 );
  memcpy( &block3.tiltLimit, &buffer128[12], 4 );
  memcpy( &block3.tiltValue, &buffer128[16], 4 );
  block3.resistanceError = UINT8(&buffer[20+blockSize*2]);
  block3.tiltError = UINT8(&buffer[21+blockSize*2]);

  if( myNumBlocks < 4 ) return;
  
  memcpy( buffer128, &buffer[0+blockSize*3], 24 );
  if( myIsLittleEndian ) {
    cseis_geolib::swapEndian( buffer128, 4, 24 );
  }
  memcpy( &block4.capacitanceLimitLow, &buffer128[0], 4 );
  memcpy( &block4.capacitanceLimitHigh, &buffer128[4], 4 );
  memcpy( &block4.capacitanceValue, &buffer128[8], 4 );
  memcpy( &block4.cutOffLimitLow, &buffer128[12], 4 );
  memcpy( &block4.cutOffLimitHigh, &buffer128[16], 4 );
  memcpy( &block4.cutOffValue, &buffer128[20], 4 );

  block4.capacitanceError = UINT8(&buffer[24+blockSize*3]);
  block4.cutOffError = UINT8(&buffer[25+blockSize*3]);

  if( myNumBlocks < 5 ) return;
  memcpy( buffer128, &buffer[0+blockSize*4], 8 );
  if( myIsLittleEndian ) {
    cseis_geolib::swapEndian( buffer128, 4, 8 );
  }
  memcpy( &block5.leakageLimit, &buffer128[0], 4 );
  memcpy( &block5.leakageValue, &buffer128[4], 4 );
  block5.leakageError = UINT8(&buffer[24+blockSize*4]);

  if( myNumBlocks < 6 ) return;
  block6.unitType   = UINT8(&buffer[0+blockSize*5]); // 0: not identified, 1: FDU, 2: SU6R
  block6.unitSerialNumber   = UINT24(&buffer[1+blockSize*5]);
  block6.chanNumber   = UINT8(&buffer[4+blockSize*5]);
  block6.fduAssemblyType   = UINT8(&buffer[8+blockSize*5]);
  block6.fduAssemblySerialNumber   = UINT24(&buffer[9+blockSize*5]);
  block6.fduAssemblyLocation   = UINT8(&buffer[12+blockSize*5]);
  block6.fduUnitType   = UINT8(&buffer[16+blockSize*5]);
  block6.chanType   = UINT8(&buffer[17+blockSize*5]);
  //  block6.sensorSensitivity   = (float)UINT32(&buffer[20+blockSize*5]);
  if( myIsLittleEndian ) {
    char buffer4[4];
    memcpy( buffer4, &buffer[20+blockSize*5], 4 );
    cseis_geolib::swapEndian4( buffer4, 4 );
    memcpy( &block6.sensorSensitivity, buffer4, 4 );
  }
  else {
    memcpy( &block6.sensorSensitivity, &buffer[20+blockSize*5], 4 );
  }

  if( myNumBlocks < 7 ) return;
  block7.controlUnitType = UINT8(&buffer[0+blockSize*6]);
  block7.controlUnitSerialNumber = UINT8(&buffer[1+blockSize*6]);
  block7.chanGainScale = UINT8(&buffer[4+blockSize*6]);
  block7.chanFilter = UINT8(&buffer[5+blockSize*6]);
  block7.chanDataError = UINT8(&buffer[6+blockSize*6]);
  block7.chanEditStatus = UINT8(&buffer[7+blockSize*6]);

  if( myIsLittleEndian ) {
    char buffer4[4];
    memcpy( buffer4, &buffer[8+blockSize*6], 4 );
    cseis_geolib::swapEndian4( buffer4, 4 );
    memcpy( &block7.sampleToMVconversionFactor, buffer4, 4 );
  }
  else {
    memcpy( &block7.sampleToMVconversionFactor, &buffer[8+blockSize*6], 4 );
  }
  //  memcpy( &block7.sampleToMVconversionFactor, &buffer[8+blockSize*6], 4 );
  block7.numStacksLow = UINT8(&buffer[12+blockSize*6]);
  block7.numStacksNoisy = UINT8(&buffer[13+blockSize*6]);
  block7.chanTypeID = UINT8(&buffer[14+blockSize*6]);
  block7.chanProcess = UINT8(&buffer[15+blockSize*6]);

  memcpy( buffer128, &buffer[16+blockSize*6], 8 );
  if( myIsLittleEndian ) {
    cseis_geolib::swapEndian( buffer128, 4, 8 );
  }
  memcpy( &block7.analogLowcutFilter, &buffer128[16], 4 );
  memcpy( &block7.digitalLowcutFilter, &buffer128[20], 4 );
  block7.numEditions = UINT32(&buffer[24+blockSize*6]);
  memcpy( buffer128, &buffer[28+blockSize*6], 4 );
  if( myIsLittleEndian ) {
    cseis_geolib::swapEndian( buffer128, 4, 4 );
  }
  memcpy( &block7.compoundLowcutFilter, &buffer128[28], 4 );
}

namespace cseis_segd {
void csExtendedHeader_SEAL::dump( std::ostream& os )
{
  os <<
    header("ExtendedHeader start") << '\n' <<
    "acqLength (1000-128000ms)     : " << acqLength << endl <<
    "sampleInt (250-4000us)        : " << sampleInt << endl <<
    "totalNumTraces                : " << totalNumTraces << endl <<
    "numAuxTraces                  : " << numAuxTraces << endl <<
    "numSeismicTraces              : " << numSeismicTraces << endl <<
    "numDeadSeismicTraces          : " << numDeadSeismicTraces << endl <<
    "numLiveSeismicTraces          : " << numLiveSeismicTraces << endl <<
    "sourceType                    : " << sourceType << endl <<
    "numSamples                    : " << numSamples << endl <<
    "shotNumber                    : " << shotNumber << endl <<
    "TBwindow                      : " << TBwindow  << endl <<
    "testRecordType                : " << testRecordType << endl <<
    "swathFirstLine                : " << swathFirstLine << endl <<
    "swathFirstNum                 : " << swathFirstNum << endl <<
    "spreadNumber                  : " << spreadNumber << endl <<
    "spreadType                    : " << spreadType << endl <<
    "timeBreak                     : " << timeBreak << endl <<
    "upholeTime                    : " << upholeTime << endl <<
    "blasterID                     : " << blasterID  << endl <<
    "blasterStatus                 : " << blasterStatus << endl <<
    "refractionDelay               : " << refractionDelay << endl <<
    "TBtoT0time                    : " << TBtoT0time << endl <<
    "internalTimeBreak             : " << internalTimeBreak << endl <<
    "prestackWithinFieldUnits      : " << prestackWithinFieldUnits << endl <<
    "noiseElimType                 : " << noiseElimType<< endl <<
    "lowTracePerc                  : " << lowTracePerc << endl <<
    "lowTraceValue                 : " << lowTraceValue << endl <<
    "numWindows                    : " << numWindows << endl <<
    "historicEditType              : " << historicEditType << endl <<
    "noisyTracePerc                : " << noisyTracePerc << endl <<
    "historicRange                 : " << historicRange   << endl <<
    "historicTaperLengthExponent   : " << historicTaperLengthExponent   << endl <<
    "thresholdVar                  : " << thresholdVar   << endl <<
    "historicThresholdInitValue    : " << historicThresholdInitValue   << endl <<
    "historicZeroingLength         : " << historicZeroingLength   << endl <<
    "typeOfProcess                 : " << typeOfProcess   << endl <<
    "acquisitionTypeTables         : ";
  for( int i = 0; i < 32; i++ ) {
    os << acquisitionTypeTables[i] << " ";
  }
  os << endl <<
    "thresholdTypeTables           : ";
  for( int i = 0; i < 32; i++ ) {
    os << thresholdTypeTables[i] << " ";
  }
  os << endl <<
    "stackFold                     : " << stackFold   << endl <<   // 1 to 32
    "auxCorrelationDesc            : " << auxCorrelationDesc   << endl <<
    "investigationLength           : " << investigationLength   << endl <<    // 100-128000ms
    "autocorrelationPeakTime       : " << autocorrelationPeakTime   << endl <<
    "fftExponent                   : " << fftExponent   << endl <<
    "sourceAuxNB                   : " << sourceAuxNB   << endl <<
    "pilotLength                   : " << pilotLength   << endl <<
    "sweepLength                   : " << sweepLength   << endl <<
    "acquisitionNum                : " << acquisitionNum   << endl <<
    "maxOfMaxAux                   : " << maxOfMaxAux   << endl <<
    "maxOfMaxSeis                  : " << maxOfMaxSeis   << endl <<
    "dumpStackFold                 : " << dumpStackFold   << endl <<

    "tapeLabel                     : " << tapeLabel << endl <<
    "tapeNum                       : " << tapeNum << endl <<
    "softwareRevision              : " << softwareRevision << endl <<
    "date                          : " << date  << endl <<
    "srcEasting                    : " << srcEasting<< endl <<
    "srcNorthing                   : " << srcNorthing << endl <<
    "srcElev                       : " << srcElev << endl <<
    "filesPerTape                  : " << filesPerTape << endl <<
    "fileCount                     : " << fileCount<< endl <<
    "acquisitionErrorDesc          : " << acquisitionErrorDesc << endl <<
    "filterType                    : " << filterType  << endl <<
    "profileNum                    : " << profileNum  << endl <<
    "t0mode                        : " << t0mode  << endl <<
    "srcLineName                   : " << srcLineName<< endl <<
    "seqNumber                     : " << seqNumber << endl <<
    "navStatus                     : " << navStatus<< endl <<
    "shotsToGo                     : " << shotsToGo << endl <<
    header("ExtendedHeader end  ");
}

void csTraceHeaderExtension_SEAL::dump( std::ostream& os ) {
  os <<
    header("TraceHeaderExtension start") << endl;
  {
    os <<
      "rcvLineNumber          : " << block1.rcvLineNumber << endl <<
      "rcvPointNumber         : " << block1.rcvPointNumber << endl <<
      "rcvPointIndex          : " << block1.rcvPointIndex << endl <<
      "numSamples             : " << block1.numSamples << endl <<
      "extendedRcvLineNumber  : " << block1.extendedRcvLineNumber << endl <<
      "extendedRcvPointNumber : " << block1.extendedRcvPointNumber << endl <<
      "sensorType             : " << block1.sensorType << endl;
      if( myNumBlocks < 2 ) goto endofdump;
    os <<
      "rcvEasting          : " << block2.rcvEasting << endl <<
      "rcvNorthing         : " << block2.rcvNorthing << endl <<
      "rcvElevation        : " << block2.rcvElevation << endl <<
      "sensorTypeNB        : " << block2.sensorTypeNB << endl;
      if( myNumBlocks < 3 )  goto endofdump;
    os <<
      "resistanceLimitLow  : " << block3.resistanceLimitLow << endl <<
      "resistanceLimitHigh : " << block3.resistanceLimitHigh << endl <<
      "resistanceValue     : " << block3.resistanceValue << endl <<
      "tiltLimit           : " << block3.tiltLimit << endl <<
      "tiltValue           : " << block3.tiltValue << endl <<
      "resistanceError     : " << block3.resistanceError << endl <<
      "tiltError           : " << block3.tiltError << endl;
      if( myNumBlocks < 4 )  goto endofdump;
    os <<
      "capacitanceLimitLow  : " << block4.capacitanceLimitLow << endl <<
      "capacitanceLimitHigh : " << block4.capacitanceLimitHigh << endl <<
      "capacitanceValue     : " << block4.capacitanceValue << endl <<
      "cutOffLimitLow       : " << block4.cutOffLimitLow << endl <<
      "cutOffLimitHigh      : " << block4.cutOffLimitHigh << endl <<
      "cutOffValue          : " << block4.cutOffValue << endl <<
      "capacitanceError     : " << block4.capacitanceError << endl <<
      "cutOffError          : " << block4.cutOffError << endl;
      if( myNumBlocks < 5 )  goto endofdump;
    os <<
      "leakageLimit         : " << block5.leakageLimit << endl <<
      "leakageValue         : " << block5.leakageValue << endl <<
      "leakageError         : " << block5.leakageError << endl;
      if( myNumBlocks < 6 )  goto endofdump;
    os <<
      "unitType                : " << block6.unitType << endl <<
      "unitSerialNumber        : " << block6.unitSerialNumber << endl <<
      "chanNumber              : " << block6.chanNumber << endl <<
      "fduAssemblyType         : " << block6.fduAssemblyType << endl <<
      "fduAssemblySerialNumber : " << block6.fduAssemblySerialNumber << endl <<
      "fduAssemblyLocation     : " << block6.fduAssemblyLocation << endl <<
      "fduUnitType             : " << block6.fduUnitType << endl <<
      "chanType                : " << block6.chanType << endl <<
      "sensorSensitivity       : " << block6.sensorSensitivity << endl;
      if( myNumBlocks < 7 )  goto endofdump;
    os <<
      "controlUnitType            : " << block7.controlUnitType << endl <<
      "controlUnitSerialNumber    : " << block7.controlUnitSerialNumber << endl <<
      "chanGainScale              : " << block7.chanGainScale << endl <<
      "chanFilter                 : " << block7.chanFilter << endl <<
      "chanDataError              : " << block7.chanDataError << endl <<
      "chanEditStatus             : " << block7.chanEditStatus << endl <<
      "sampleToMVconversionFactor : " << block7.sampleToMVconversionFactor << endl <<
      "numStacksLow               : " << block7.numStacksLow << endl <<
      "numStacksNoisy             : " << block7.numStacksNoisy << endl <<
      "chanTypeID                 : " << block7.chanTypeID << endl <<
      "chanProcess                : " << block7.chanProcess << endl <<
      "analogLowcutFilter         : " << block7.analogLowcutFilter << endl <<
      "digitalLowcutFilter        : " << block7.digitalLowcutFilter << endl <<
      "numEditions                : " << block7.numEditions << endl <<
      "compoundLowcutFilter       : " << block7.compoundLowcutFilter << endl;
    }
endofdump:
    os << header("TraceHeaderExtension end  ") << endl;

}

}

