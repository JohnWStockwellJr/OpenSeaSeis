/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegdHeader.h"
#include "csException.h"
#include "csSegdDefines.h"
#include <cmath>

using cseis_geolib::csException;
using std::endl;
using namespace cseis_segd;

//----------------------------------------------------------
csBaseHeader::csBaseHeader() {}
csBaseHeader::~csBaseHeader() {}
void csBaseHeader::extractHeaders( byte const* buffer ) {
}

//----------------------------------------------------------
csGeneralHeader1::csGeneralHeader1() : csBaseHeader() {}
csGeneralHeader1::~csGeneralHeader1() {}
void csGeneralHeader1::extractHeaders( byte const* buffer ) {
  fileNum                 = bcd(&buffer[0], 0, 4);
  formatCode              = bcd(&buffer[2], 0, 4);
  for (int i=0, j=0; i<6; i++) {
    genConsts[j++]        = bcd(&buffer[4+i], 0, 2);
    genConsts[j++]        = bcd(&buffer[4+i], 1, 2);
  }
  year                    = bcd(&buffer[10], 0, 2);
  numGeneralHdrBlocks     = UINT4H(&buffer[11]);
  julianDay               = bcd(&buffer[11], 1, 3);
  hour                    = bcd(&buffer[13], 0, 2);
  min                     = bcd(&buffer[14], 0, 2);
  sec                     = bcd(&buffer[15], 0, 2);
  manufactCode            = bcd(&buffer[16], 0, 2);
  manufactSerialNum       = bcd(&buffer[17], 0, 4);
  numBytesPerScan         = bcd(&buffer[19], 0, 6);
  baseScanInterval        = UINT8(&buffer[22]); //  (buffer[22]*1000) >> 4; gives ms
  polarityCode            = UINT4H(&buffer[23]);
  numScansPerBlock        = bcd(&buffer[23], 0, 3);
  recordTypeCode          = UINT4H(&buffer[25]);
  recordLen               = bcd(&buffer[25], 1, 3);
  //  recordLenMS             = ((recordLen * 1024) / 2) / 10;  // If this rounds down, this is not supported
  recordLenMS             = (recordLen * 1024) / 10;  // Fix 12/04/16: Multiple by 1024ms instead of 512ms
  numScanTypesPerRecord    = bcd(&buffer[27], 0, 2);
  numChanSetsPerScanType  = bcd(&buffer[28], 0, 2);
  numSkewBlocks           = bcd(&buffer[29], 0, 2);
  numExtendedHdrBlocks    = bcd(&buffer[30], 0, 2);
  numExternalHdrBlocks    = bcd(&buffer[31], 0, 2);

}

//----------------------------------------------------------
csGeneralHeader2::csGeneralHeader2() : csBaseHeader() {
  myNumBytes_externalHdrBlocks = 2;
}
csGeneralHeader2::~csGeneralHeader2() {}
void csGeneralHeader2::set( int numBytes_externalHdrBlocks ) {
  myNumBytes_externalHdrBlocks = numBytes_externalHdrBlocks;
}
void csGeneralHeader2::extractHeaders( byte const* buffer )
{

  expandedFileNum                = UINT24(&buffer[0]);
  numExtendedChanSets            = UINT16(&buffer[3]);
  numExtendedHdrBlocks           = UINT16(&buffer[5]);
  numExternalHdrBlocks           = UINT16(&buffer[7]);
  if( myNumBytes_externalHdrBlocks == 3 ) {
    int numExternalHdrBlocks_NEW = UINT24(&buffer[7]);
    // fprintf(stderr,"External header blocks: %d --> %d\n", numExternalHdrBlocks, numExternalHdrBlocks_NEW );
    numExternalHdrBlocks = numExternalHdrBlocks_NEW;
  }
  revisionNum[0]                 = buffer[10];
  revisionNum[1]                 = buffer[11];
  numGeneralTrailerBlocks        = UINT16(&buffer[12]);
  extendedRecordLen              = UINT24(&buffer[14]);
  generalHeaderBlockNum          = buffer[18];
}

//----------------------------------------------------------
csGeneralHeaderN::csGeneralHeaderN() : csBaseHeader() {}
csGeneralHeaderN::~csGeneralHeaderN() {}
void csGeneralHeaderN::extractHeaders( byte const* buffer )
{
  expandedFileNum       = UINT24(&buffer[0]);

  sourceLineNum = INT24(&buffer[3]);
  double tmp = (double)UINT16(&buffer[6]);
  double tmp2 = tmp;
  if( tmp2 < 1.0 ) tmp2 = 1.0;
  sourceLineNum += tmp / std::pow( 10.0, int(std::log(tmp2/10.0)) );

  sourcePointNum = INT24(&buffer[8]);
  tmp = (double)UINT16(&buffer[11]);
  tmp2 = tmp;
  if( tmp2 < 1.0 ) tmp2 = 1.0;
  sourcePointNum += tmp / std::pow( 10.0, int(std::log(tmp2/10.0)) );

  sourcePointIndex      = buffer[13];
  phaseControl          = buffer[14];
  vibratorType          = buffer[15];
  phaseAngle            = UINT16(&buffer[16]);
  generalHeaderBlockNum = buffer[18];
  sourceSetNum          = buffer[19];
}

void showbits( byte a ) {
  int i  , k , mask;

  for( i =7 ; i >= 0 ; i--) {
    mask = 1 << i;
    k = a & mask;
    if( k == 0)
      fprintf(stderr,"0 ");
    else
      fprintf(stderr,"1 ");
  }
  fprintf(stderr,"\n");
}

//----------------------------------------------------------
csChanSetHeader::csChanSetHeader() : csBaseHeader() {}
csChanSetHeader::~csChanSetHeader() {}
void csChanSetHeader::extractHeaders( byte const* buffer )
{
  scanTypeNum              = bcd(&buffer[0], 0, 2);
  chanSetNum               = bcd(&buffer[1], 0, 2);
  chanSetStartTime         = UINT16(&buffer[2]) * 2;  // Fix 12/04/16: multiply by 2ms
  chanSetEndTime           = UINT16(&buffer[4]) * 2;  // Fix 12/04/16: multiply by 2ms

  // Computation of mp descale factor. From SU
  unsigned char mp1 = (unsigned char)buffer[6];
  unsigned char mp2 = (unsigned char)buffer[7];
  double mpNew = ((double) ((mp2 & 0x7f) << 8 | mp1)) / 1024.;
  if( mp2 >> 7 ) mpNew *= -1.;
  mpFactor = pow ((double) 2., (double) mpNew);
  
  numChannels              = bcd(&buffer[8], 0, 4);
  chanTypeID               = UINT4H(&buffer[10]);
  subScanPerBaseScanBinExp = bcd(&buffer[11], 0, 1);  // Fixed 11/08/24
  gainControlCode          = UINT4L(&buffer[11]);
  aliasFilterFreq          = bcd(&buffer[12], 0, 4);
  aliasFilterSlope         = bcd(&buffer[14], 0, 4);
  lowCutFilterFreq         = bcd(&buffer[16], 0, 4);
  lowCutFilterSlope        = bcd(&buffer[18], 0, 4);
  firstNotchFilterFreq     = bcd(&buffer[20], 0, 4);
  secondNotchFilterFreq    = bcd(&buffer[22], 0, 4);
  thirdNotchFilterFreq     = bcd(&buffer[24], 0, 4);
  extendedChanSetNum       = UINT16(&buffer[26]);
  extendedHeaderFlag       = UINT4H(&buffer[28]);
  numTraceHeaderExtensions_rev2 = UINT4L(&buffer[28]);  // Only valid for rev 2
  verticalStack            = buffer[29];
  streamerCableNum         = buffer[30];
  arrayForming             = buffer[31];
}

//----------------------------------------------------------
csSampleSkew::csSampleSkew() : csBaseHeader() {
  skewFactor = new int[csBaseHeader::BLOCK_SIZE];
}
csSampleSkew::~csSampleSkew() {
  if( skewFactor != NULL ) {
    delete [] skewFactor;
    skewFactor = NULL;
  }
}
void csSampleSkew::extractHeaders( byte const* buffer ) {
  //...
}

//----------------------------------------------------------
csExtendedHeader::csExtendedHeader() : csBaseHeader() {}
csExtendedHeader::~csExtendedHeader() {}
void csExtendedHeader::extractHeaders( byte const* buffer) {
}
void csExtendedHeader::extractHeaders( byte const* buffer, int totalNumBytes ) {
}
int csExtendedHeader::nanoSeconds() const {
  return 0;
}
int csExtendedHeader::numTraces() const {
  return 0;
}

//----------------------------------------------------------
csBaseExternalHeader::csBaseExternalHeader( int numBlocks ) : csBaseHeader() {
  myNumBlocks = numBlocks;
}
csBaseExternalHeader::~csBaseExternalHeader() {}
//void csBaseExternalHeader::extractHeaders( byte const* buffer, int totalNumBytes ) {
//}
//----------------------------------------------------------
csTraceHeader::csTraceHeader() : csBaseHeader() {}
csTraceHeader::~csTraceHeader() {}
void csTraceHeader::extractHeaders( byte const* buffer )
{
  fileNumber               = bcd(&buffer[0], 0, 4);
  scanType                 = bcd(&buffer[2], 0, 2);
  chanSetNumber            = bcd(&buffer[3], 0, 2);
  traceNumber              = bcd(&buffer[4], 0, 4);
  firstTimingWord          = bcd(&buffer[6], 0, 3);

  numTraceHeaderExtensions = UINT4L(&buffer[9]);
  sampleSkew               = UINT8(&buffer[10]);
  traceEdit                = UINT8(&buffer[11]);
  timeBreakWindow          = UINT24(&buffer[12]);
  extendedChanSetNumber    = UINT16(&buffer[15]);
  extendedFileNumber       = UINT24(&buffer[17]);
}

//----------------------------------------------------------
csGeneralTrailer::csGeneralTrailer() : csBaseHeader() {}
csGeneralTrailer::~csGeneralTrailer() {}
void csGeneralTrailer::extractHeaders( byte const* buffer ) {
  generalTrailerNumber = UINT16(&buffer[0]);
  chanTypeId           = buffer[10];
// My interpretation:
//  chanTypeId           = UINT16(&buffer[10]);
}

//----------------------------------------------------------
csTraceHeaderExtension::csTraceHeaderExtension( int numBlocks ) : csBaseHeader() {
  myNumBlocks = numBlocks;
}
csTraceHeaderExtension::~csTraceHeaderExtension() {}
int csTraceHeaderExtension::getNumSamples() const {
  return 0;
}
int csTraceHeaderExtension::numBlocks() const {
  return myNumBlocks;
}
int csTraceHeaderExtension::nanoSeconds() const {
  return 0;
}
csInt64_t csTraceHeaderExtension::timeSamp1_us() const {
  return 0;
}
void csTraceHeaderExtension::extractHeaders( byte const* buffer, commonTraceHeaderStruct* comTrcHdr ) {
  comTrcHdr->rcvLineNumber  = UINT24(&buffer[0]);
  comTrcHdr->rcvPointNumber = UINT24(&buffer[3]);
  comTrcHdr->rcvPointIndex  = UINT8(&buffer[6]);

  if( comTrcHdr->rcvLineNumber == 16777215 ) {  // == FF FF FF
    comTrcHdr->rcvLineNumber  = UINT32(&buffer[10]);
  }
  if( comTrcHdr->rcvPointNumber == 16777215 ) {  // == FF FF FF
    comTrcHdr->rcvPointNumber  = UINT32(&buffer[15]);
  }
}
void csTraceHeaderExtension::extractHeaders( byte const* buffer ) {
// My interpretation of the format:
  block0.rcvLineNumber      = UINT24(&buffer[0]);
  block0.rcvPointNumber     = UINT24(&buffer[3]);
  block0.rcvPointIndex      = UINT24(&buffer[6]);
  block0.numSamples         = UINT24(&buffer[9]);

  block0.extendedRcvLineNumber  = INT24(&buffer[11]) + UINT16(&buffer[14])/100;
  block0.extendedRcvPointNumber = INT24(&buffer[16]) + UINT16(&buffer[19])/100;
  block0.sensorType = UINT8(&buffer[21]);
}

//-----------------------------------------------------------
void csBaseHeader::dump( std::ostream& cs ) {}
void csGeneralHeader1::dump( std::ostream& cs )
{
  cs <<
  header("generalHeader1 start") << '\n' <<
  "file number                               : " << fileNum << '\n' <<
  "format code                               : " << formatCode << '\n' <<
  "general constants                         : ";
  for (int i=0; i<12; i++)
    cs << genConsts[i] << ' ';
  cs << '\n' <<
  "year                                      : " << year << '\n' <<
  "number of blks in general header extension: " << numGeneralHdrBlocks << '\n' <<
  "julian day                                : " << julianDay << '\n' <<
  "minute                                    : " << min << '\n' <<
  "hour                                      : " << hour << '\n' <<
  "second                                    : " << sec << '\n' <<
  "manufacturer's code                       : " << manufactCode << '\n' <<
  "manufacturer's serial number              : " << manufactSerialNum << '\n' <<
  "bytes per scan (multiplexed formats)      : " << numBytesPerScan << '\n' <<
  "base scan interval                        : " << baseScanInterval << " 1/16ms  (= "<< ((float)baseScanInterval/16.0) <<"ms)\n" <<
  "polarity                                  : " << polarityCode << '\n' <<
  "scan/block                                : " << numScansPerBlock << '\n' <<
  "record type (8=normal record)             : " << recordTypeCode << '\n' <<
  "record length (raw)                       : " << recordLen << '\n' <<
  "   --> record length (ms)                 : " << recordLenMS << '\n' <<
  "scan types / record                       : " << numScanTypesPerRecord << '\n' <<
  "channel sets per scan type                : " << numChanSetsPerScanType << '\n' <<
  "number of 32-bytes fields for sample skew : " << numSkewBlocks << '\n' <<
  "extended header length                    : " << numExtendedHdrBlocks << '\n' <<
  "external header length                    : " << numExternalHdrBlocks << '\n' <<
  header("generalHeader1 end  ");

}

void csGeneralHeader2::dump( std::ostream& cs )
{
  cs <<
  header("generalHeader2 start") << '\n' <<
  "expanded file number            : " << expandedFileNum << '\n' <<
  "extended channel sets/scantype  : " << numExtendedChanSets << '\n' <<
  "extended header blocks          : " << numExtendedHdrBlocks << '\n' <<
  "external header blocks          : " << numExternalHdrBlocks << '\n' <<
  "seg-d revision number           : " << revisionNum[0] <<'.'<< revisionNum[1] << '\n' <<
  "general trailer number of blocks: " << numGeneralTrailerBlocks << '\n' <<
  "extended record length          : " << extendedRecordLen << " msec\n" <<
  "general header block number     : " << generalHeaderBlockNum << '\n' <<
  header("generalHeader2 end  ");
}
void csGeneralHeaderN::dump( std::ostream& cs )
{
  cs << header("generalHeaderN") << '\n' <<
  "expanded file number       : " << expandedFileNum << '\n' <<
  "source line number         : " << sourceLineNum << '\n' <<
  "source point number        : " << sourcePointNum << '\n' <<
  "source point index         : " << sourcePointIndex << '\n' <<
  "phase control              : " << phaseControl << '\n' <<
  "vibrator type              : " << vibratorType << '\n' <<
  "phase angle                : " << phaseAngle << '\n' <<
  "general header block number: " << generalHeaderBlockNum << '\n' <<
  "source set number          : " << sourceSetNum << endl <<
  header("generalHeaderN end  ") << std::endl;
}
void csChanSetHeader::dump( std::ostream& cs )
{
  cs << header("channelSetHeader") << '\n' <<
  "scan type number                 : " << scanTypeNum << '\n' <<
  "channel set number               : " << chanSetNum << '\n' <<
  "channel set start time           : " << chanSetStartTime << '\n' <<
  "channel set end time             : " << chanSetEndTime << '\n' <<
  "2^(MP factor) = descale operator : " << mpFactor << '\n' <<
  "number of channels in channel set: " << numChannels << '\n' <<
  "channel type code                : " << chanTypeID << '\n' <<
  "subscans per base scan binary exp: " << subScanPerBaseScanBinExp << '\n' <<
  "gain control code                : " << gainControlCode << '\n' <<
  "alias filter frequency           : " << aliasFilterFreq << " Hz\n" <<
  "alias filter slope               : " << aliasFilterSlope << " dB/Octave\n" <<
  "low cut filter frequency         : " << lowCutFilterFreq << " Hz\n" <<
  "low cut filter slope             : " << lowCutFilterSlope << " dB/Octave\n" <<
  "first notch filter frequency     : " << firstNotchFilterFreq << " 1/10 Hz\n" <<
  "second notch filter frequency    : " << secondNotchFilterFreq << " 1/10 Hz\n" <<
  "third notch filter frequency     : " << thirdNotchFilterFreq << " 1/10 Hz\n" <<
  "extended channel set number      : " << extendedChanSetNum << '\n' <<
  "extended header flag             : " << extendedHeaderFlag << '\n' <<
  "number of trc hdr ext (rev >=1.5): " << numTraceHeaderExtensions_rev2 << '\n' <<
  "vertical stack                   : " << verticalStack << '\n' <<
  "streamer cable number            : " << streamerCableNum << '\n' <<
  "array forming                    : " << arrayForming << endl <<
  header("chanSetHeader end  ") << std::endl;
}

void csSampleSkew::dump( std::ostream& cs )
{
  cs << header("Sample skew factors...") << '\n' <<
  "...not implemented yet." << '\n';
}

void csTraceHeader::dump( std::ostream& cs )
{
  cs <<
    header("TraceHeader start") << '\n' <<
    "fileNumber            : " << fileNumber << std::endl <<
    "scanTypeNumber        : " << scanType << std::endl <<
    "chanSetNumber         : " << chanSetNumber << std::endl <<
    "traceNumber           : " << traceNumber << std::endl <<
    "firstTimingWord       : " << firstTimingWord << std::endl <<
    "numTraceHeaderExtensions : " << numTraceHeaderExtensions << std::endl <<
    "sampleSkew            : " << sampleSkew << std::endl <<
    "traceEdit             : " << traceEdit << std::endl <<
    "timeBreakWindow       : " << timeBreakWindow << std::endl <<
    "extendedChanSetNumber : " << extendedChanSetNumber << std::endl <<
    "extendedFileNumber    : " << extendedFileNumber << std::endl <<
    header("TraceHeader end  ") << '\n';
}

void csGeneralTrailer::dump( std::ostream& cs )
{
  cs <<
    header("GeneralTrailer start") << '\n' <<
    "generalTrailerNumber : " << generalTrailerNumber << std::endl <<
    "chanTypeId           : " << chanTypeId << std::endl <<
    header("GeneralTrailer end  ") << '\n';
}

void csTraceHeaderExtension::dump( std::ostream& cs )
{
  cs <<
    header("TraceHeaderExtension start") << std::endl <<
    "rcvLineNumber          : " << block0.rcvLineNumber << std::endl <<
    "rcvPointNumber         : " << block0.rcvPointNumber << std::endl <<
    "rcvPointIndex          : " << block0.rcvPointIndex << std::endl <<
    "numSamples             : " << block0.numSamples << std::endl <<
    "extendedRcvLineNumber  : " << block0.rcvLineNumber << std::endl <<
    "extendedRcvPointNumber : " << block0.rcvPointNumber << std::endl <<
    "sensorType             : " << block0.sensorType << std::endl <<
    header("TraceHeaderExtension end  ") << std::endl;
}

void csExtendedHeader::dump( std::ostream& cs )
{
  cs <<
    header("ExtendedHeader start") << std::endl <<
    "...manufacturer/recording system not supported..." << std::endl <<
    header("ExtendedHeader end  ") << std::endl;
}

void csBaseExternalHeader::dump( std::ostream& cs )
{
  cs <<
    header("ExternalHeader start") << std::endl <<
    "...external header format not supported..." << std::endl <<
    header("ExternalHeader end  ") << std::endl;
}


