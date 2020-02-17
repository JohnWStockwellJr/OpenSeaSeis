/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegyBinHeader.h"
#include "geolib_endian.h"
#include "csByteConversions.h"
#include <cstring>

using namespace cseis_geolib;

csSegyBinHeader::csSegyBinHeader( bool doSwapEndian ) {
  myDoSwapEndian = doSwapEndian;

  jobID = 0;
  lineNum = 0;
  reelNum = 0;
  numTraces = 0;
  numAuxTraces = 0;
  sampleIntUS = 0;
  sampleIntOrigUS = 0;
  numSamples = 0;
  numSamplesOrig = 0;
  dataSampleFormat = 0;
  fold = 0;
  sortCode = 0;
  vertSumCode = 0;
  sweepFreqStart = 0; // Hz
  sweepFreqEnd = 0;   // Hz
  sweepCode = 0;
  taperType = 0;
  correlatedTraces = 0;
  gainRecovered = 0;
  ampRecoveryMethod = 0;
  unitSystem = 0;
  polarity = 0;
  vibPolarityCode = 0;
  revisionNum = 0;
  fixedTraceLengthFlag = 0;
  numExtendedBlocks = 0;
}
csSegyBinHeader::~csSegyBinHeader() {
  // nothing...
}
void csSegyBinHeader::extractHeaders( byte_t const* buffer ) {

  if( myDoSwapEndian ) {
    jobID            = byte2Int_SWAP( &buffer[0] );
    lineNum          = byte2Int_SWAP( &buffer[4] );
    reelNum          = byte2Int_SWAP( &buffer[8] );

    numTraces        = byte2UShort_SWAP( &buffer[12] );
    numAuxTraces     = byte2UShort_SWAP( &buffer[14] );
    sampleIntUS      = byte2UShort_SWAP( &buffer[16] );
    sampleIntOrigUS  = byte2UShort_SWAP( &buffer[18] );
    numSamples       = byte2UShort_SWAP( &buffer[20] );
    numSamplesOrig   = byte2UShort_SWAP( &buffer[22] );
    dataSampleFormat = byte2UShort_SWAP( &buffer[24] );

    fold             = byte2UShort_SWAP( &buffer[26] );
    sortCode         = byte2UShort_SWAP( &buffer[28] );
    vertSumCode      = byte2UShort_SWAP( &buffer[30] );
    sweepFreqStart   = byte2UShort_SWAP( &buffer[32] );
    sweepFreqEnd     = byte2UShort_SWAP( &buffer[34] );
    sweepCode        = byte2UShort_SWAP( &buffer[36] );
    taperType        = byte2UShort_SWAP( &buffer[46] );
    correlatedTraces = byte2UShort_SWAP( &buffer[48] );
    gainRecovered    = byte2UShort_SWAP( &buffer[50] );
    ampRecoveryMethod= byte2UShort_SWAP( &buffer[52] );
    unitSystem       = byte2UShort_SWAP( &buffer[54] );

    polarity         = byte2UShort_SWAP( &buffer[56] );
    vibPolarityCode  = byte2UShort_SWAP( &buffer[58] );

    revisionNum      = byte2UShort_SWAP( &buffer[300] );
    fixedTraceLengthFlag  = byte2UShort_SWAP( &buffer[302] );
    numExtendedBlocks = byte2UShort_SWAP( &buffer[304] );
  }
  else {
    jobID            = byte2Int( &buffer[0] );
    lineNum          = byte2Int( &buffer[4] );
    reelNum          = byte2Int( &buffer[8] );

    numTraces        = byte2UShort( &buffer[12] );
    numAuxTraces     = byte2UShort( &buffer[14] );
    sampleIntUS      = byte2UShort( &buffer[16] );
    sampleIntOrigUS  = byte2UShort( &buffer[18] );
    numSamples       = byte2UShort( &buffer[20] );
    numSamplesOrig   = byte2UShort( &buffer[22] );
    dataSampleFormat = byte2UShort( &buffer[24] );

    fold             = byte2UShort( &buffer[26] );
    sortCode         = byte2UShort( &buffer[28] );
    vertSumCode      = byte2UShort( &buffer[30] );
    sweepFreqStart   = byte2UShort( &buffer[32] );
    sweepFreqEnd     = byte2UShort( &buffer[34] );
    sweepCode        = byte2UShort( &buffer[36] );
    taperType        = byte2UShort( &buffer[46] );
    correlatedTraces = byte2UShort( &buffer[48] );
    gainRecovered    = byte2UShort( &buffer[50] );
    ampRecoveryMethod= byte2UShort( &buffer[52] );
    unitSystem       = byte2UShort( &buffer[54] );
    polarity         = byte2UShort( &buffer[56] );
    vibPolarityCode  = byte2UShort( &buffer[58] );

    revisionNum      = byte2UShort( &buffer[300] );
    fixedTraceLengthFlag  = byte2UShort( &buffer[302] );
    numExtendedBlocks = byte2UShort( &buffer[304] );
  }
}
void csSegyBinHeader::encodeHeaders( byte_t* buffer ) const {
  if( myDoSwapEndian ) {
    int2Byte_SWAP( jobID, &buffer[0] );
    int2Byte_SWAP( lineNum, &buffer[4] );
    int2Byte_SWAP( reelNum, &buffer[8] );

    short2Byte_SWAP( (short)numTraces, &buffer[12] );
    short2Byte_SWAP( (short)numAuxTraces, &buffer[14] );
    short2Byte_SWAP( (short)sampleIntUS, &buffer[16] );
    short2Byte_SWAP( (short)sampleIntOrigUS, &buffer[18] );
    short2Byte_SWAP( (short)numSamples, &buffer[20] );
    short2Byte_SWAP( (short)numSamplesOrig, &buffer[22] );
    short2Byte_SWAP( (short)dataSampleFormat, &buffer[24] );

    short2Byte_SWAP( (short)fold, &buffer[26] );
    short2Byte_SWAP( (short)sortCode, &buffer[28] );
    short2Byte_SWAP( (short)vertSumCode, &buffer[30] );
    short2Byte_SWAP( (short)sweepFreqStart, &buffer[32] );
    short2Byte_SWAP( (short)sweepFreqEnd, &buffer[34] );
    short2Byte_SWAP( (short)sweepCode, &buffer[36] );
    short2Byte_SWAP( (short)taperType, &buffer[46] );
    short2Byte_SWAP( (short)correlatedTraces, &buffer[48] );
    short2Byte_SWAP( (short)gainRecovered, &buffer[50] );
    short2Byte_SWAP( (short)ampRecoveryMethod, &buffer[52] );
    short2Byte_SWAP( (short)unitSystem, &buffer[54] );

    short2Byte_SWAP( (short)polarity, &buffer[56] );
    short2Byte_SWAP( (short)vibPolarityCode, &buffer[58] );
    
    short2Byte_SWAP( (short)revisionNum, &buffer[300] );
    short2Byte_SWAP( (short)fixedTraceLengthFlag, &buffer[302] );
    short2Byte_SWAP( (short)numExtendedBlocks, &buffer[304] );
  }
  else {
    int2Byte( jobID, &buffer[0] );
    int2Byte( lineNum, &buffer[4] );
    int2Byte( reelNum, &buffer[8] );

    short2Byte( (short)numTraces, &buffer[12] );
    short2Byte( (short)numAuxTraces, &buffer[14] );
    short2Byte( (short)sampleIntUS, &buffer[16] );
    short2Byte( (short)sampleIntOrigUS, &buffer[18] );
    short2Byte( (short)numSamples, &buffer[20] );
    short2Byte( (short)numSamplesOrig, &buffer[22] );
    short2Byte( (short)dataSampleFormat, &buffer[24] );

    short2Byte( (short)fold, &buffer[26] );
    short2Byte( (short)sortCode, &buffer[28] );
    short2Byte( (short)vertSumCode, &buffer[30] );
    short2Byte( (short)sweepFreqStart, &buffer[32] );
    short2Byte( (short)sweepFreqEnd, &buffer[34] );
    short2Byte( (short)sweepCode, &buffer[36] );
    short2Byte( (short)taperType, &buffer[46] );
    short2Byte( (short)correlatedTraces, &buffer[48] );
    short2Byte( (short)gainRecovered, &buffer[50] );
    short2Byte( (short)ampRecoveryMethod, &buffer[52] );
    short2Byte( (short)unitSystem, &buffer[54] );

    short2Byte( (short)polarity, &buffer[56] );
    short2Byte( (short)vibPolarityCode, &buffer[58] );

    short2Byte( (short)revisionNum, &buffer[300] );
    short2Byte( (short)fixedTraceLengthFlag, &buffer[302] );
    short2Byte( (short)numExtendedBlocks, &buffer[304] );
  }
}
void csSegyBinHeader::dump( FILE* dumpFile ) const {
  fprintf( dumpFile,"Job identification number     : %u\n", jobID );
  fprintf( dumpFile,"Line number                   : %u\n", lineNum );
  fprintf( dumpFile,"Reel number                   : %u\n", reelNum );
  fprintf( dumpFile,"Number of traces/ensemble     : %u\n", numTraces );
  fprintf( dumpFile,"Number of aux traces/ensemble : %u\n", numAuxTraces );
  fprintf( dumpFile,"Sample interval [us]          : %u\n", sampleIntUS );
  fprintf( dumpFile,"Sample interval(field tape)[us]:%u\n", sampleIntOrigUS );
  fprintf( dumpFile,"Number of samples             : %u\n", numSamples );
  fprintf( dumpFile,"Number of samples (field tape): %u\n", numSamplesOrig );
  fprintf( dumpFile,"Data sample format            : %u\n", dataSampleFormat );
  fprintf( dumpFile,"Ensemble fold                 : %u\n", fold );
  fprintf( dumpFile,"Trace sort code               : %u\n", sortCode );
  fprintf( dumpFile,"Vertical sum code             : %u\n", vertSumCode );
  fprintf( dumpFile,"Sweep freq at start           : %u\n", sweepFreqStart );
  fprintf( dumpFile,"Sweep freq at end             : %u\n", sweepFreqEnd );
  fprintf( dumpFile,"Sweep type code               : %u\n", sweepCode );
  fprintf( dumpFile,"Taper type                    : %u\n", taperType );
  fprintf( dumpFile,"Correlated data traces        : %u\n", correlatedTraces );
  fprintf( dumpFile,"Binary gain recovered         : %u\n", gainRecovered );
  fprintf( dumpFile,"Amplitude recovery method     : %u\n", ampRecoveryMethod );
  fprintf( dumpFile,"Measurement system            : %u\n", unitSystem );
  fprintf( dumpFile,"Impulse signal polarity       : %u\n", polarity );
  fprintf( dumpFile,"Vibratory polarity code       : %u\n", vibPolarityCode );
  fprintf( dumpFile,"SEG Y Format Revision number  : %u\n", revisionNum );
  fprintf( dumpFile,"Fixed length trace flag       : %u\n", fixedTraceLengthFlag );
  fprintf( dumpFile,"Extended header blocks        : %u\n", numExtendedBlocks );
}

