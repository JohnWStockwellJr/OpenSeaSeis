/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegdHeader.h"
#include "csSegdReader.h"
#include "csSegdDefines.h"
#include "csSegdFunctions.h"
#include "csSegdBuffer.h"
#include "csSegdHdrValues.h"

#include "csExternalHeader.h"
#include "csSegdHeader_GEORES.h"
#include "csSegdHeader_SEAL.h"
#include "csSegdHeader_DIGISTREAMER.h"
#include "csSegdHeader_FAIRFIELD.h"

#include "csFileUtils.h"
#include "csTimer.h"
#include "csException.h"
#include "geolib_endian.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

using namespace cseis_segd;
using cseis_geolib::csException;
using std::cout;
using std::endl;
using std::string;
using std::memcpy;

void convertToFloat_20bit( short const* in, float* out, int numSamples );
void convertToFloat_24bit( unsigned char const* in, float* out, int numSamples );

#define SIGNBIT (1UL << 23)

csSegdReader::csSegdReader() {
  myRecordingSystemID = UNKNOWN;
  myFile = NULL;
  myCopyBuffer = NULL;
  myRecordByteSize = 0;
  myGeneralHdr1 = new csGeneralHeader1();
  myGeneralHdr2 = new csGeneralHeader2();
  myGeneralHdrN = NULL;
  myChanSetHdr  = NULL;
  myChanSetNumSamples = NULL;
  myChanSetSampleInt_us  = NULL;
  myMPDescaleOperator = NULL;
  mySampleSkew  = NULL;
  myExtendedHdr = NULL;
  myExternalHdr = NULL;
  myTraceHdr    = new csTraceHeader();
  myTraceHdrExtension = NULL;
  myBuffer_oneRecord = NULL;
  myIsSupported = true;
  myFullTraceByteSize = NULL;
  mySegdHdrValues = new csSegdHdrValues();
  myChanSetIndexToRead = -1;

  myBuffer_gen1 = new byte[csBaseHeader::BLOCK_SIZE];
  myBuffer_gen2 = new byte[csBaseHeader::BLOCK_SIZE];
  buffer_genN          = new csSegdBuffer();
  buffer_scantypes     = new csSegdBuffer();
  buffer_extendedHdr   = new csSegdBuffer();
  buffer_externalHdr   = new csSegdBuffer();

  myBuffer_traceHdr      = NULL;
  myBuffer_traceHdr      = new csSegdBuffer();

  myBuffer_traceHdrExtension = NULL;
  myBuffer_traceHdrExtension = new csSegdBuffer();

  myBytePos.chanSetHdr  = 0;
  myBytePos.extendedHdr = 0;
  myBytePos.externalHdr = 0;
  myBytePos.firstTrace  = 0;
  myBytePos.generalHdr1 = 0;
  myBytePos.generalHdr2 = 0;
  myBytePos.generalHdrN = 0;
  myBytePos.currentTraceData = 0;
  myBytePos.firstTraceData   = 0;

  myNumScanTypes             = 0;
  myNumChanSetsPerScanType   = 0;  // Maximum number of chan sets per scan types (if less, dummy chan sets are present)
  myNumExtendedHdrBlocks     = 0;
  myNumExternalHdrBlocks     = 0;
  myNumGeneralTrailerBlocks  = 0;
  myNumGeneralHdrBlocks      = 0;  // All general header blocks
  myNumExtraGeneralHdrBlocks = 0;  // Additional blocks to general header 1 & 2
  myNumSampleSkewBlocks      = 0;
  myNumTraceHdrExtensions    = -1;

  myChanTypeID = NULL;
  myDumpFile = NULL;
  mySequentialTraceCounter = 0;

  setDefaultConfiguration();
}
//---------------------------------------------------------
csSegdReader::~csSegdReader() {
  if( myFile != NULL ) {
    fclose( myFile );
    myFile = NULL;
  }
  if( myCopyBuffer != NULL ) delete [] myCopyBuffer;
  if( mySegdHdrValues != NULL ) delete mySegdHdrValues;
  if( myBuffer_gen1 ) delete [] myBuffer_gen1;
  if( myBuffer_gen2 ) delete [] myBuffer_gen2;
  myBuffer_gen1 = NULL; myBuffer_gen2 = NULL;
  if( buffer_genN ) delete buffer_genN;
  if( buffer_scantypes ) delete buffer_scantypes;
  if( buffer_extendedHdr ) delete buffer_extendedHdr;
  if( buffer_externalHdr ) delete buffer_externalHdr;
  if( myBuffer_traceHdr != NULL ) delete myBuffer_traceHdr;
  if( myBuffer_traceHdrExtension != NULL ) delete myBuffer_traceHdrExtension;
  buffer_genN          = NULL;
  buffer_scantypes     = NULL;
  buffer_extendedHdr   = NULL;
  buffer_externalHdr   = NULL;
  myBuffer_traceHdr      = NULL;
  myBuffer_traceHdrExtension = NULL;
    
  if( myBuffer_oneRecord != NULL ) {
    delete [] myBuffer_oneRecord;
    myBuffer_oneRecord = NULL;
  }
  if( myGeneralHdr1 != NULL ) {
    delete myGeneralHdr1;
    myGeneralHdr1 = NULL;
  }
  if( myGeneralHdr2 != NULL ) {
    delete myGeneralHdr2;
    myGeneralHdr2 = NULL;
  }
  if( myGeneralHdrN != NULL ) {
    delete [] myGeneralHdrN;
    myGeneralHdrN = NULL;
  }
  if( myChanSetHdr != NULL ) {
    delete [] myChanSetHdr;
    myChanSetHdr = NULL;
  }
  if( myChanSetNumSamples != NULL ) {
    delete [] myChanSetNumSamples;
    myChanSetNumSamples = NULL;
  }
  if( myChanSetSampleInt_us != NULL ) {
    delete [] myChanSetSampleInt_us;
    myChanSetSampleInt_us = NULL;
  }
  if( myFullTraceByteSize != NULL ) {
    delete [] myFullTraceByteSize;
    myFullTraceByteSize= NULL;
  }
  if( myMPDescaleOperator != NULL ) {
    delete [] myMPDescaleOperator;
    myMPDescaleOperator = NULL;
  }
  if( mySampleSkew != NULL ) {
    delete [] mySampleSkew;
    mySampleSkew = NULL;
  }
  if( myExtendedHdr != NULL ) {
    delete myExtendedHdr;
    myExtendedHdr = NULL;
  }
  if( myExternalHdr != NULL ) {
    delete myExternalHdr;
    myExternalHdr = NULL;
  }
  if( myTraceHdr != NULL ) {
    delete myTraceHdr;
    myTraceHdr = NULL;
  }
  if( myTraceHdrExtension != NULL ) {
    delete myTraceHdrExtension;
    myTraceHdrExtension = NULL;
  }
  if( myChanTypeID != NULL ) {
    delete [] myChanTypeID; myChanTypeID = NULL;
  }
}
void csSegdReader::closeFile() {
  if( myFile != NULL ) {
    fclose( myFile );
    myFile = NULL;
  }
}
//---------------------------------------------------------
void csSegdReader::setDefaultConfiguration() {
  myConfig.recordingSystemID = UNKNOWN;
  myConfig.navSystemID       = UNKNOWN;
  myConfig.navInterfaceID    = UNKNOWN;
  myConfig.isDebug           = false;
  myConfig.numSamplesAddOne  = false;
  myConfig.thisIsRev0        = false;
  myConfig.readAuxTraces     = false;
  myConfig.overrideTraceHdrExtensions = false;
  myConfig.numTraceHdrExtensions = 0;
}
//---------------------------------------------------------
void csSegdReader::setConfiguration( csSegdReader::configuration& config ) {
  if( &config != &myConfig ) memcpy( &myConfig, &config, sizeof(csSegdReader::configuration) );
  if( config.isDebug ) {
    fprintf(stdout,"Configuration:\n");
    fprintf(stdout,"recordingSystemID    : %d\n", config.recordingSystemID );
    fprintf(stdout,"navSystemID          : %d\n", config.navSystemID );
    fprintf(stdout,"navInterfaceID       : %d\n", config.navInterfaceID );
    fprintf(stdout,"thisIsRev0           : %s\n", config.thisIsRev0 ? "yes" : "no" );
    fprintf(stdout,"readAuxTraces        : %s\n", config.readAuxTraces ? "yes" : "no" );
    fprintf(stdout,"numSamplesAddOne     : %s\n", config.numSamplesAddOne ? "yes" : "no" );
    fprintf(stdout,"isDebug              : %s\n", config.isDebug ? "yes" : "no" );
  }
}
//---------------------------------------------------------
// Read in general headers and chan set headers of first file
// Check consistency of SEGD file
// Throw exception if problem occurs
// 
//
bool csSegdReader::open( std::string filename ) {
  myFileName = filename;
  myHasJustBeenInitialized = true; 
  if( myFile != NULL ) {
    fclose( myFile );
    myFile = NULL;
  }
  myFile = fopen( myFileName.c_str(), "rb" );
  if( myFile == NULL ) {
    throw( cseis_geolib::csException("Could not open SEGD file") );
  }
  return true;
}
bool csSegdReader::readNewRecordHeaders() {
  // 1. Read general header 1
  // 2. Read general header 2
  // 3. Read general header N...
  // 4. Read chan set headers for all scan types

  myHasJustBeenInitialized = true;

  try {
    if( !readBuffer( myBuffer_gen1, csBaseHeader::BLOCK_SIZE ) ) {
      return false;
    }
  }
  catch( csException& exc ) {
    return false;
  }
  try {
    // Extract headers for general header 1
    myGeneralHdr1->extractHeaders( myBuffer_gen1 );
    myBytePos.generalHdr1 = 0;

    if( !myConfig.thisIsRev0 ) {
      readBuffer( myBuffer_gen2, csBaseHeader::BLOCK_SIZE );
      // Extract headers for general header 2
      myGeneralHdr2->set( myConfig.numBytes_genhdr2_extHdrBlocks );
      myGeneralHdr2->extractHeaders( myBuffer_gen2 );
      myBytePos.generalHdr2 = csBaseHeader::BLOCK_SIZE;
      myBytePos.generalHdrN = myBytePos.generalHdr2 + csBaseHeader::BLOCK_SIZE;
    }
    else {
      myBytePos.generalHdr2 = myBytePos.generalHdr1;
      myBytePos.generalHdrN = myBytePos.generalHdr2;
    }

  // Set/compute header values + Check header values
    myNumScanTypes              = myGeneralHdr1->numScanTypesPerRecord;
    myNumGeneralHdrBlocks       = myGeneralHdr1->numGeneralHdrBlocks+1;
    myNumExtraGeneralHdrBlocks  = myGeneralHdr1->numGeneralHdrBlocks-1;
    myNumSampleSkewBlocks       = myGeneralHdr1->numSkewBlocks;

    myComFileHdr.sampleInt_us     = myGeneralHdr1->getSampleInt_us();
    myComFileHdr.sampleBitSize     = sampleBitSize( myGeneralHdr1->formatCode );
    myComFileHdr.formatCode        = myGeneralHdr1->formatCode;
    myComFileHdr.manufactCode      = myGeneralHdr1->manufactCode;
    myComFileHdr.manufactSerialNum = myGeneralHdr1->manufactSerialNum;

    if( !myConfig.thisIsRev0 ) {
      myNumChanSetsPerScanType    = myGeneralHdr1->numChanSetsPerScanType <= 99 ? myGeneralHdr1->numChanSetsPerScanType : myGeneralHdr2->numExtendedChanSets;
      myNumExtendedHdrBlocks      = myGeneralHdr1->numExtendedHdrBlocks <= 99 ? myGeneralHdr1->numExtendedHdrBlocks : myGeneralHdr2->numExtendedHdrBlocks;
      myNumExternalHdrBlocks      = myGeneralHdr1->numExternalHdrBlocks <= 99 ? myGeneralHdr1->numExternalHdrBlocks : myGeneralHdr2->numExternalHdrBlocks;
      myNumGeneralTrailerBlocks   = myGeneralHdr2->numGeneralTrailerBlocks;

      myComFileHdr.recordLength_ms   = myGeneralHdr1->recordLen <= 999 ? myGeneralHdr1->recordLenMS : myGeneralHdr2->extendedRecordLen;
      myComFileHdr.revisionNum[0]    = myGeneralHdr2->revisionNum[0];
      myComFileHdr.revisionNum[1]    = myGeneralHdr2->revisionNum[1];
      if( myComFileHdr.sampleInt_us == 0 ) throw( csException("Trace sample interval is 0") );
      myComFileHdr.numSamples        = (1000*myComFileHdr.recordLength_ms) / myComFileHdr.sampleInt_us;
      // +1 or not, this is inconsistent in different example SEGD files!
      if( myConfig.numSamplesAddOne ) {
        myComFileHdr.numSamples += 1;
      }
    }
    else {  // SEG-D revision 0 special
      myNumChanSetsPerScanType    = myGeneralHdr1->numChanSetsPerScanType;
      myNumExtendedHdrBlocks      = myGeneralHdr1->numExtendedHdrBlocks;
      myNumExternalHdrBlocks      = myGeneralHdr1->numExternalHdrBlocks;
      myNumGeneralTrailerBlocks   = 0;

      myComFileHdr.recordLength_ms   = myGeneralHdr1->recordLenMS;
      myComFileHdr.revisionNum[0]    = 0;
      myComFileHdr.revisionNum[1]    = 0;

      myComFileHdr.numSamples        = (1000*myComFileHdr.recordLength_ms) / myComFileHdr.sampleInt_us + 1;
      //      if( (int)(myComFileHdr.numSamples/4)*4 != myComFileHdr.numSamples ) {
      //  myComFileHdr.numSamples = ((int)(myComFileHdr.numSamples/4)+1)*4;
      // }
    }

    myIsSupported =
      isFormatCodeSupported(myGeneralHdr1->formatCode ) &&
      isRevisionSupported(myComFileHdr.revisionNum[0],myComFileHdr.revisionNum[1]);
    //      isManufacturerSupported(myGeneralHdr1->manufactCode) &&
    if( !isFormatCodeSupported(myGeneralHdr1->formatCode ) ) {
      throw( csException("SEG-D format code %d is not supported.", myGeneralHdr1->formatCode) );
    }
    if( !isRevisionSupported(myComFileHdr.revisionNum[0],myComFileHdr.revisionNum[1]) ) {
      throw( csException("SEG-D revision %d.%d is not supported.", myComFileHdr.revisionNum[0],myComFileHdr.revisionNum[1]) );
    }
    // Read & extract headers for general header N...
    if( !myConfig.thisIsRev0 ) {
      myBytePos.chanSetHdr = myBytePos.generalHdrN + myNumExtraGeneralHdrBlocks*csBaseHeader::BLOCK_SIZE;
      if( myNumExtraGeneralHdrBlocks > 0 ) {
        if( myGeneralHdrN ) {
          delete[] myGeneralHdrN; myGeneralHdrN = NULL;
        }
        myGeneralHdrN = new csGeneralHeaderN[myNumExtraGeneralHdrBlocks];
        buffer_genN->setNumBytes( csBaseHeader::BLOCK_SIZE*myNumExtraGeneralHdrBlocks );
        readBuffer( buffer_genN->buffer(), myNumExtraGeneralHdrBlocks*csBaseHeader::BLOCK_SIZE );
        for( int iblock = 0; iblock < myNumExtraGeneralHdrBlocks; iblock++ ) {
          myGeneralHdrN[iblock].extractHeaders( &(buffer_genN->buffer()[iblock*csBaseHeader::BLOCK_SIZE]) );
        }
      }
    }
    else {
      myBytePos.chanSetHdr = myBytePos.generalHdr1 + myNumGeneralHdrBlocks*csBaseHeader::BLOCK_SIZE;
    }
    
    // Read & extract headers for chan set headers for each scan type and sample skew...
    if( myChanSetHdr ) {
      delete [] myChanSetHdr;
      myChanSetHdr = NULL;
    }
    if( myMPDescaleOperator ) {
      delete [] myMPDescaleOperator;
      myMPDescaleOperator = NULL;
    }
    if( myChanSetNumSamples == NULL ) {
      delete [] myChanSetNumSamples;
      myChanSetNumSamples = NULL;
    }
    if( myChanSetSampleInt_us != NULL ) {
      delete [] myChanSetSampleInt_us;
      myChanSetSampleInt_us = NULL;
    }
    if( myFullTraceByteSize != NULL ) {
      delete [] myFullTraceByteSize;
      myFullTraceByteSize = NULL;
    }
    myChanSetHdr = new csChanSetHeader[myNumScanTypes*myNumChanSetsPerScanType];
    myChanSetSampleInt_us = new int[myNumScanTypes*myNumChanSetsPerScanType];
    myFullTraceByteSize   = new int[myNumScanTypes*myNumChanSetsPerScanType];
    myMPDescaleOperator   = new double[myNumScanTypes*myNumChanSetsPerScanType];
    myChanSetNumSamples   = new int[myNumScanTypes*myNumChanSetsPerScanType];
    for( int i = 0; i < myNumScanTypes*myNumChanSetsPerScanType; i++ ) {
      myChanSetNumSamples[i] = 0;
    }

    if( myNumSampleSkewBlocks == 0 ) {
      mySampleSkew = NULL;
    }
    else {
      if( mySampleSkew ) {
        delete [] mySampleSkew;
        mySampleSkew = NULL;
      }
      mySampleSkew = new csSampleSkew[myNumScanTypes*myNumSampleSkewBlocks];
    }

    int totalNumTraces   = 0;
    int numSeismicTraces = 0;
    int numBytes = csBaseHeader::BLOCK_SIZE * (myNumScanTypes * (myNumChanSetsPerScanType + myNumSampleSkewBlocks));
    buffer_scantypes->setNumBytes( numBytes );
    readBuffer( buffer_scantypes->buffer(), numBytes );
    myBytePos.extendedHdr = myBytePos.chanSetHdr + numBytes;
    for( int iscan = 0; iscan < myNumScanTypes; iscan++ ) {
      for( int ichanset = 0; ichanset < myNumChanSetsPerScanType; ichanset++ ) {
        int index = iscan*myNumChanSetsPerScanType + ichanset;
        myChanSetHdr[index].extractHeaders( &(buffer_scantypes->buffer()[index*csBaseHeader::BLOCK_SIZE]) );
        totalNumTraces += myChanSetHdr[index].numChannels;
        if( myChanSetHdr[index].chanTypeID == 1 ) numSeismicTraces += myChanSetHdr[index].numChannels;
        myMPDescaleOperator[index] = myChanSetHdr[index].mpFactor;
        myNumTraceHdrExtensions   = std::max( myNumTraceHdrExtensions, myChanSetHdr[index].numTraceHeaderExtensions_rev2 );
        myChanSetSampleInt_us[index] = (int)( myComFileHdr.sampleInt_us / pow(2,myChanSetHdr[index].subScanPerBaseScanBinExp) );
        myChanSetNumSamples[index] = (myChanSetHdr[index].chanSetEndTime - myChanSetHdr[index].chanSetStartTime) * 1000 / myChanSetSampleInt_us[index];
        if( myConfig.numSamplesAddOne ) {
          myChanSetNumSamples[index] += 1;
        }
        if( myConfig.isDebug ) {
          fprintf(stdout,"Chan set #%d: Sample interval: %dus, number of samples: %d, number of traces: %d   (comFileHdr.numSamples: %d), num trace header extensions: %d\n",
                  index, myChanSetSampleInt_us[index], myChanSetNumSamples[index], myChanSetHdr[index].numChannels, myComFileHdr.numSamples, myNumTraceHdrExtensions );
        }
      }
      for( int iskew = 0; iskew < myNumSampleSkewBlocks; iskew++ ) {
        mySampleSkew[iscan*myNumSampleSkewBlocks + iskew].extractHeaders( buffer_scantypes->buffer() );
      }
    }

    // Read extended header, extract if supported
    buffer_extendedHdr->setNumBytes( myNumExtendedHdrBlocks*csBaseHeader::BLOCK_SIZE );
    if( myConfig.isDebug ) {
      fprintf(stdout,"Reading in %d bytes extended header blocks (#%d)\n", myNumExtendedHdrBlocks*csBaseHeader::BLOCK_SIZE, myNumExtendedHdrBlocks );
    }

    readBuffer( buffer_extendedHdr->buffer(), myNumExtendedHdrBlocks*csBaseHeader::BLOCK_SIZE );
    myBytePos.externalHdr = myBytePos.extendedHdr + myNumExtendedHdrBlocks*csBaseHeader::BLOCK_SIZE;
    myRecordingSystemID = manufacturerRecordingSystem( myGeneralHdr1->manufactCode );
    if( myExtendedHdr ) {
      delete myExtendedHdr;
      myExtendedHdr = NULL;
    }
    switch( myRecordingSystemID ) {
      case RECORDING_SYSTEM_GEORES:
        myExtendedHdr = new csExtendedHeader_GEORES();
        break;
      case RECORDING_SYSTEM_SEAL:
        myExtendedHdr = new csExtendedHeader_SEAL();
        break;
      case RECORDING_SYSTEM_FAIRFIELD:
        myExtendedHdr = new csExtendedHeader_FAIRFIELD();
        break;
      case RECORDING_SYSTEM_DIGISTREAMER:
        myExtendedHdr = new csExtendedHeader_DIGI();
        break;
      default:
        myExtendedHdr = new csExtendedHeader();
        break;
    }
    myExtendedHdr->extractHeaders( buffer_extendedHdr->buffer(), myNumExtendedHdrBlocks*csBaseHeader::BLOCK_SIZE );

    if( myExtendedHdr->numTraces() > 0 ) {
      totalNumTraces   = 0;
      numSeismicTraces = 0;
      if( myConfig.isDebug ) {
        fprintf(stdout,"WARNING:  Changed numChannels in each chanset set to %d (from extended header)\n", myExtendedHdr->numTraces() );
      }
      for( int iscan = 0; iscan < myNumScanTypes; iscan++ ) {
        for( int ichanset = 0; ichanset < myNumChanSetsPerScanType; ichanset++ ) {
          int index = iscan*myNumChanSetsPerScanType + ichanset;
          myChanSetHdr[index].numChannels = myExtendedHdr->numTraces();
          totalNumTraces += myChanSetHdr[index].numChannels;
          if( myChanSetHdr[index].chanTypeID == 1 ) numSeismicTraces += myChanSetHdr[index].numChannels;
        }
      }
    }

    myComFileHdr.totalNumChan   = totalNumTraces;
    myComFileHdr.numSeismicChan = numSeismicTraces;
    myComFileHdr.numAuxChan     = totalNumTraces - numSeismicTraces;
    if( myChanTypeID ) {
      delete [] myChanTypeID;
      myChanTypeID = NULL;
    }
    myChanTypeID = new int[totalNumTraces];
    totalNumTraces = 0;
    for( int iscan = 0; iscan < myNumScanTypes; iscan++ ) {
      for( int ichanset = 0; ichanset < myNumChanSetsPerScanType; ichanset++ ) {
        int index = iscan*myNumChanSetsPerScanType + ichanset;
        for( int itrc = 0; itrc < myChanSetHdr[index].numChannels; itrc++ ) {
          myChanTypeID[totalNumTraces++] = myChanSetHdr[index].chanTypeID;
        }
      }
    }



  // Read in external header, do not extract headers yet
    numBytes = myNumExternalHdrBlocks*csBaseHeader::BLOCK_SIZE;

    if( myConfig.thisIsRev0 && myRecordingSystemID == RECORDING_SYSTEM_GECO ) {
      numBytes += 8;  // !CHANGE! Why 8 bytes not accounted for???
    }

    if( myExternalHdr ) {
      delete myExternalHdr;
      myExternalHdr = NULL;
    }

    if( myConfig.isDebug ) {
      fprintf(stdout,"Reading in %d bytes external header blocks (#%d)\n", numBytes, myNumExternalHdrBlocks);
    }

    if( myNumExternalHdrBlocks > 0 ) {
      buffer_externalHdr->setNumBytes( numBytes );
      readBuffer( buffer_externalHdr->buffer(), numBytes );
      myExternalHdr = new csExternalHeader( myNumExternalHdrBlocks, myConfig.navSystemID, myConfig.navInterfaceID );
      myExternalHdr->extractHeaders( buffer_externalHdr->buffer(), numBytes );
    }
    myBytePos.firstTrace = myBytePos.externalHdr + numBytes;
    
  // Read in first trace header, to get number of trace header extensions, and check some other details...
    myBuffer_traceHdr->setNumBytes( csTraceHeader::BLOCK_SIZE );
    readBuffer( myBuffer_traceHdr->buffer(), csTraceHeader::BLOCK_SIZE );
    myTraceHdr->extractHeaders( myBuffer_traceHdr->buffer() );

    if( myComFileHdr.revisionNum[0] >= 2 || (myComFileHdr.revisionNum[0] == 1 && myComFileHdr.revisionNum[1] >= 5) ) {
      if( myTraceHdr->numTraceHeaderExtensions != myNumTraceHdrExtensions ) {
        if( !myConfig.overrideTraceHdrExtensions ) throw csException("Unequal number of trace header extensions specified in chan set header and trace header: %d != %d",
                                                                     myNumTraceHdrExtensions, myTraceHdr->numTraceHeaderExtensions );
        myTraceHdr->numTraceHeaderExtensions = myConfig.numTraceHdrExtensions;
      }
    }
    myNumTraceHdrExtensions = myTraceHdr->numTraceHeaderExtensions;

    // Read in trace header extensions for first trace
    if( myTraceHdrExtension ) {
      delete myTraceHdrExtension;
      myTraceHdrExtension = NULL;
    }
    if( myNumTraceHdrExtensions > 0 ) {
      switch( myRecordingSystemID ) {
        case RECORDING_SYSTEM_GEORES:
          myTraceHdrExtension = new csTraceHeaderExtension_GEORES( myNumTraceHdrExtensions );
          break;
        case RECORDING_SYSTEM_SEAL:
          myTraceHdrExtension = new csTraceHeaderExtension_SEAL( myNumTraceHdrExtensions );
          break;
        case RECORDING_SYSTEM_FAIRFIELD:
          myTraceHdrExtension = new csTraceHeaderExtension_FAIRFIELD( myNumTraceHdrExtensions );
          break;
        case RECORDING_SYSTEM_DIGISTREAMER:
          myTraceHdrExtension = new csTraceHeaderExtension_DIGI( myNumTraceHdrExtensions );
          break;
        default:
          myTraceHdrExtension = new csTraceHeaderExtension( myNumTraceHdrExtensions );
          break;
      }
      myBuffer_traceHdrExtension->setNumBytes( csBaseHeader::BLOCK_SIZE*myNumTraceHdrExtensions );
      readBuffer( myBuffer_traceHdrExtension->buffer(), csBaseHeader::BLOCK_SIZE*myNumTraceHdrExtensions );
      myBytePos.firstTraceExt = myBytePos.firstTrace + csTraceHeader::BLOCK_SIZE;
      myTraceHdrExtension->extractHeaders( myBuffer_traceHdrExtension->buffer() );
      if( myTraceHdrExtension->getNumSamples() != myComFileHdr.numSamples ) {
        if( myConfig.isDebug) fprintf(stdout,"WARNING: Unequal number of samples in General Header 1 and trace hdr extension: %d != %d\n",
          myComFileHdr.numSamples, myTraceHdrExtension->getNumSamples() );
      }
      if( myTraceHdrExtension->getNumSamples() != 0 ) myComFileHdr.numSamples = myTraceHdrExtension->getNumSamples();
    }
    myBytePos.firstTraceData = myBytePos.firstTrace +
        csTraceHeader::BLOCK_SIZE + csBaseHeader::BLOCK_SIZE*myNumTraceHdrExtensions;
    myBytePos.currentTraceData = myBytePos.firstTraceData;

    if( myComFileHdr.revisionNum[0] >= 2 ) {
      for( int ichanset = 0; ichanset < myNumChanSetsPerScanType; ichanset++ ) {
        myFullTraceByteSize[ichanset] =
          traceDataByteSize(ichanset) + csTraceHeader::BLOCK_SIZE + myChanSetHdr[ichanset].numTraceHeaderExtensions_rev2 * csBaseHeader::BLOCK_SIZE;
      }
    }
    else {
      for( int ichanset = 0; ichanset < myNumChanSetsPerScanType; ichanset++ ) {
        myFullTraceByteSize[ichanset] =
          traceDataByteSize(ichanset) + csTraceHeader::BLOCK_SIZE + myNumTraceHdrExtensions * csBaseHeader::BLOCK_SIZE;
      }
    }
    if( !myConfig.thisIsRev0 ) {
      myInitializeNumBytesRead = myBytePos.firstTraceData;
      /*      myInitializeNumBytesRead = 2*csBaseHeader::BLOCK_SIZE +  // general hdr 1 & 2
        myNumExtraGeneralHdrBlocks*csBaseHeader::BLOCK_SIZE +  // general hdrs N
        (myNumScanTypes*(myNumChanSetsPerScanType + myNumSampleSkewBlocks))*csBaseHeader::BLOCK_SIZE +   // chanset hdr
        myNumExtendedHdrBlocks*csBaseHeader::BLOCK_SIZE +   // extended hdr
        myNumExternalHdrBlocks*csBaseHeader::BLOCK_SIZE +   // external hdr
        csTraceHeader::BLOCK_SIZE +   // trace hdr
        csBaseHeader::BLOCK_SIZE*myNumTraceHdrExtensions;   // trace hdr extension */
    }
    else {
      myInitializeNumBytesRead = myBytePos.firstTraceData;
    }

  // Compute full SEGD record byte size
    int newRecordByteSize = csBaseHeader::BLOCK_SIZE * (
      (myGeneralHdr1->numGeneralHdrBlocks+1) +  // General header blocks
      myNumScanTypes*(myNumChanSetsPerScanType+myNumSampleSkewBlocks) +
      myNumExtendedHdrBlocks +
      myNumExternalHdrBlocks );

    int t1 = newRecordByteSize;

    newRecordByteSize +=
      myComFileHdr.totalNumChan*(csTraceHeader::BLOCK_SIZE) +
      (myComFileHdr.totalNumChan*myNumTraceHdrExtensions +
       myNumGeneralTrailerBlocks) * csBaseHeader::BLOCK_SIZE;

    for( int ichanset = 0; ichanset < myNumChanSetsPerScanType; ichanset++ ) {
      newRecordByteSize += traceDataByteSize(ichanset) * myChanSetHdr[ichanset].numChannels;
    }

    if( myConfig.isDebug ) {
      fprintf(stdout,"myRecordSize: %d (%d %d %d traceDataByteSize(1): %d)\n", newRecordByteSize, myComFileHdr.totalNumChan, myComFileHdr.numSamples, myComFileHdr.sampleBitSize, traceDataByteSize(0));
      int t2 = t1 + myComFileHdr.totalNumChan*traceDataByteSize() + myComFileHdr.totalNumChan*(csTraceHeader::BLOCK_SIZE);
      int t3 = myComFileHdr.totalNumChan*traceDataByteSize();
      int t4 = (myComFileHdr.totalNumChan*myNumTraceHdrExtensions + myNumGeneralTrailerBlocks) * csBaseHeader::BLOCK_SIZE;
      fprintf(stdout,"  %d %d %d  %d   size: %d\n", t1, t2, t3, t4, traceDataByteSize() );
      //      csInt64_t fileSize = cseis_geolib::csFileUtils::retrieveFileSize( myFilename );
    }
  // Check & dump parameters
    if( myIsSupported ) {
      if( newRecordByteSize != myRecordByteSize ) {
        myRecordByteSize = newRecordByteSize;
        if( myBuffer_oneRecord ) {
          delete [] myBuffer_oneRecord;
          myBuffer_oneRecord = NULL;
        }
        myBuffer_oneRecord = new byte[myRecordByteSize*sizeof(byte)];
      }
  
      memcpy( &myBuffer_oneRecord[0], myBuffer_gen1, csBaseHeader::BLOCK_SIZE );
      if( !myConfig.thisIsRev0 ) {
        memcpy( &myBuffer_oneRecord[myBytePos.generalHdr2], myBuffer_gen2, csBaseHeader::BLOCK_SIZE );
      }
      if( myNumExtraGeneralHdrBlocks > 0 ) {
        memcpy( &myBuffer_oneRecord[myBytePos.generalHdrN], buffer_genN->buffer(), (myBytePos.chanSetHdr-myBytePos.generalHdrN) );
      }  
      memcpy( &myBuffer_oneRecord[myBytePos.chanSetHdr], buffer_scantypes->buffer(), (myBytePos.extendedHdr-myBytePos.chanSetHdr) );
  
      if( myNumExtendedHdrBlocks > 0 ) {
        memcpy( &myBuffer_oneRecord[myBytePos.extendedHdr], buffer_extendedHdr->buffer(), (myBytePos.externalHdr-myBytePos.extendedHdr) );
      }
      if( myNumExternalHdrBlocks > 0 ) {
        memcpy( &myBuffer_oneRecord[myBytePos.externalHdr], buffer_externalHdr->buffer(), (myBytePos.firstTrace-myBytePos.externalHdr) );
      }
      memcpy( &myBuffer_oneRecord[myBytePos.firstTrace], myBuffer_traceHdr->buffer(), csTraceHeader::BLOCK_SIZE );
      if( myNumTraceHdrExtensions > 0 ) {
        memcpy( &myBuffer_oneRecord[myBytePos.firstTraceExt], myBuffer_traceHdrExtension->buffer(), csBaseHeader::BLOCK_SIZE*myNumTraceHdrExtensions );
      }
    }
  }
  catch( csException& exc ) {
    throw exc;
//    printf("Exception occurred: %s\n", exc.getMessage() );
  }

  /*  if( buffer_genN ) delete [] buffer_genN;
  if( buffer_scantypes ) delete [] buffer_scantypes;
  if( buffer_extendedHdr ) delete [] buffer_extendedHdr;
  if( buffer_externalHdr ) delete [] buffer_externalHdr;
  if( myBuffer_traceHdr ) delete [] myBuffer_traceHdr;
  if( myBuffer_traceHdrExtension ) delete [] myBuffer_traceHdrExtension;
  */
  
  //  mySequentialTraceCounter = myComFileHdr.totalNumChan;  // Disable retrieval of data traces until data has been read in
  mySequentialTraceCounter = 0;
  mySequentialChanSetIndexCounter = 0;
  mySequentialChanSetAccTraces = 0;
  mySegdHdrValues->setFileHeaderValues( myComFileHdr );
  return true;
  // In Next method, read in whole SEGD record
}  // END: bool csSegdReader::readNewRecordHeaders() 
//---------------------------------------------------------

//
bool csSegdReader::readCharHdr( char* charHdr, int numBytes ) {
  return( readBuffer( (byte*)charHdr, numBytes ) );
}
bool csSegdReader::readBuffer( byte* buffer, int numBytes ) {
  int sizeRead = (int)fread( buffer, 1, numBytes, myFile );
  return( sizeRead != 0 );
}
//---------------------------------------------------------
//
bool csSegdReader::readNextRecord( commonRecordHeaderStruct& comRecHdr ) {
  if( myConfig.isDebug ) fprintf(stdout,"Read next record. Current trace #%d / %d   (%d)\n", mySequentialTraceCounter, myComFileHdr.totalNumChan, myHasJustBeenInitialized );

  int startByte = 0;
  if( !myHasJustBeenInitialized ) {
    try {
      if( !readNewRecordHeaders() ) {
        return false;
      }
    }
    catch( csException& e ) {
      // ..reading in header failed. Return false without printing error message.
      return false;
    }
  }
  if( myHasJustBeenInitialized ) {
    startByte = myInitializeNumBytesRead;
  }
  int numBytes = myRecordByteSize - startByte;
  if( myConfig.isDebug ) {
    fprintf(stdout,"Reading %d bytes (=full record of %d bytes + initial %d bytes)\n", numBytes, myRecordByteSize, startByte );
    csInt64_t fileSize = cseis_geolib::csFileUtils::retrieveFileSize( myFileName );
    csInt64_t remainBytes = fileSize - myBytePos.firstTrace; // startByte; // myInitializeNumBytesRead
    int numTraces = (int)( remainBytes / (csInt64_t)traceDataByteSize(0) );
    fprintf(stdout,"File size: %lld  , remaining bytes: %lld, num traces:  %d  %d  %d\n", fileSize, remainBytes, numTraces, numTraces/4, numTraces/(4*35) );
  }
  if( readBuffer( &myBuffer_oneRecord[startByte], numBytes ) ) {
    extractCommonRecordHeaders( comRecHdr );
    mySegdHdrValues->setRecordHeaderValues( comRecHdr );
    mySequentialTraceCounter = 0;   // Reset sequential trace counter for trace retrieval
    mySequentialChanSetIndexCounter = 0;
    mySequentialChanSetAccTraces = 0;
    myBytePos.currentTraceData = myBytePos.firstTraceData;
    if( myHasJustBeenInitialized ) myHasJustBeenInitialized = false;
    return true;
  }
  else if( myHasJustBeenInitialized ) {
    throw csException("Error occurred when reading SEGD trace data. SEGD file corrupt, or unsupported SEGD format..?");
  }
  else {
    return false;
  }
}
//---------------------------------------------------------
bool csSegdReader::getNextTrace( float* trace, commonTraceHeaderStruct& comTrcHdr ) {
  if( myConfig.isDebug ) fprintf(stdout,"Read next trace: #%d / %d\n", mySequentialTraceCounter, myComFileHdr.totalNumChan );

  if( myChanSetIndexToRead >= 0 ) {
    // Restrict to reading only specified channel set, if requested
    while( mySequentialChanSetIndexCounter < myChanSetIndexToRead ) {
      myBytePos.currentTraceData   += myChanSetHdr[mySequentialChanSetIndexCounter].numChannels * myFullTraceByteSize[mySequentialChanSetIndexCounter];
      mySequentialTraceCounter     += myChanSetHdr[mySequentialChanSetIndexCounter].numChannels;
      mySequentialChanSetAccTraces += myChanSetHdr[mySequentialChanSetIndexCounter].numChannels;
      mySequentialChanSetIndexCounter += 1;
    }
    if( mySequentialTraceCounter >= myChanSetHdr[mySequentialChanSetIndexCounter].numChannels+mySequentialChanSetAccTraces ) return false;
    if( mySequentialChanSetIndexCounter > myChanSetIndexToRead ) return false;
  }
  else {
    if( mySequentialTraceCounter >= myChanSetHdr[mySequentialChanSetIndexCounter].numChannels+mySequentialChanSetAccTraces ) {
      mySequentialChanSetAccTraces    += myChanSetHdr[mySequentialChanSetIndexCounter].numChannels;
      mySequentialChanSetIndexCounter += 1;
      if( mySequentialChanSetIndexCounter >= (myNumScanTypes * myNumChanSetsPerScanType) ) return false;
    }
  }
  if( mySequentialTraceCounter == myComFileHdr.totalNumChan ) return false;
  if( myChanSetIndexToRead >= 0 && mySequentialChanSetIndexCounter != myChanSetIndexToRead ) return false;
  int bytePosData = myBytePos.currentTraceData;

  int bytePosHdr  = bytePosData - csTraceHeader::BLOCK_SIZE; // - 32;
  if( myTraceHdrExtension ) bytePosHdr -= csBaseHeader::BLOCK_SIZE * myTraceHdrExtension->numBlocks();

  myBytePos.currentTraceData += myFullTraceByteSize[mySequentialChanSetIndexCounter];
  int dataByteSize = traceDataByteSize(mySequentialChanSetIndexCounter);
  int numSamples = myChanSetNumSamples[mySequentialChanSetIndexCounter];
  extractCommonTraceHeaders( &myBuffer_oneRecord[bytePosHdr], comTrcHdr );
  comTrcHdr.chanSet      = mySequentialChanSetIndexCounter + 1; // +1 to convert to 'user-domain' number starting at 1
  comTrcHdr.numSamples   = myChanSetNumSamples[mySequentialChanSetIndexCounter];
  comTrcHdr.sampleInt_us = myChanSetSampleInt_us[mySequentialChanSetIndexCounter];
  bool retValue = true;

  if( myConfig.isDebug ) { 
    fprintf(stdout,"Reading #%2d in from byte %d / %d   %d %d,  chanset : %d  first: %d, fullTraceSize: %d, numchan: %d\n", mySequentialTraceCounter,
            bytePosData, myRecordByteSize, dataByteSize,
            bytePosData+dataByteSize,
            mySequentialChanSetIndexCounter,
            myBytePos.firstTraceData, myFullTraceByteSize[mySequentialChanSetIndexCounter],
            myChanSetHdr[mySequentialChanSetIndexCounter].numChannels+mySequentialChanSetAccTraces );
  }
  mySequentialTraceCounter += 1;

  if( comTrcHdr.chanTypeID == 1 || myConfig.readAuxTraces ) {  // Only read in if this is a seismic trace, or if aux traces shall be read in as well
    switch( myComFileHdr.formatCode ) {
    case 8058:
      cseis_geolib::swapEndian4( reinterpret_cast<char*>(&myBuffer_oneRecord[bytePosData]), dataByteSize );
      memcpy( trace, &myBuffer_oneRecord[bytePosData], dataByteSize );
      break;
    case 8015:
      convertToFloat_20bit( reinterpret_cast<short*>(&myBuffer_oneRecord[bytePosData]), trace, numSamples );
      break;
    case 8036:
      fprintf(stdout,"Byte pos %d  numSamples %d  byteSize %d\n", bytePosData, numSamples, dataByteSize);
      convertToFloat_24bit( reinterpret_cast<unsigned char*>(&myBuffer_oneRecord[bytePosData]), trace, numSamples );
      break;
    }
  
    for( int i = 0; i < numSamples; i++ ) {
      trace[i] *= (float)myMPDescaleOperator[mySequentialChanSetIndexCounter];
    }
  }
  else { // try reading in next trace
    retValue = getNextTrace( trace, comTrcHdr );
    if( retValue ) mySegdHdrValues->setTraceHeaderValues( comTrcHdr );
  }
  
  return retValue;
}
//---------------------------------------------------------
//
void csSegdReader::extractCommonRecordHeaders( commonRecordHeaderStruct& comRecHdr ) {
  myGeneralHdr1->extractHeaders( myBuffer_oneRecord );

  comRecHdr.fileNum            = myGeneralHdr1->fileNum != 16665 ? myGeneralHdr1->fileNum : myGeneralHdr2->expandedFileNum;  // Fix 12/04/17
  comRecHdr.shotTime.year      = myGeneralHdr1->year;
  comRecHdr.shotTime.julianDay = myGeneralHdr1->julianDay;
  comRecHdr.shotTime.hour      = myGeneralHdr1->hour;
  comRecHdr.shotTime.minute    = myGeneralHdr1->min;
  comRecHdr.shotTime.second    = myGeneralHdr1->sec;

  comRecHdr.shotTime.microSecond = 0;
  comRecHdr.shotNum    = 0;
  comRecHdr.srcEasting = 0;
  comRecHdr.srcEasting = 0;
  comRecHdr.srcElev    = 0;

  if( myExternalHdr ) {
    myExternalHdr->setCommonRecordHeaders( comRecHdr );
  }
}
//---------------------------------------------------------
//
void csSegdReader::extractCommonTraceHeaders( byte const* bufferPtr, commonTraceHeaderStruct& comTrcHdr ) {
  myTraceHdr->extractHeaders( bufferPtr );
  comTrcHdr.chanNum    = myTraceHdr->traceNumber;
  comTrcHdr.traceEdit  = myTraceHdr->traceEdit;
  comTrcHdr.chanTypeID = myChanTypeID[mySequentialTraceCounter];

  //  myTraceHdr->dump( cout );

  if( myTraceHdrExtension ) {
    myTraceHdrExtension->extractHeaders( &bufferPtr[csTraceHeader::BLOCK_SIZE] );
    myTraceHdrExtension->extractHeaders( &bufferPtr[csTraceHeader::BLOCK_SIZE], &comTrcHdr );
  }
}
//---------------------------------------------------------
//
void csSegdReader::dumpAllTraces( std::ofstream* outStream ) {

  int traceCounter = 0;
  int bytePos = myBytePos.firstTrace;

  //  for( int iscan = 0; iscan < myNumScanTypes; iscan++ ) {
  for( int ichanset = 0; ichanset < myNumChanSetsPerScanType; ichanset++ ) {
    for( int itrc = 0; itrc < myChanSetHdr[ichanset].numChannels; itrc++ ) {
      bytePos += myFullTraceByteSize[ichanset];
      myTraceHdr->extractHeaders( &myBuffer_oneRecord[bytePos] );
      myTraceHdr->dump( *outStream );
      if( myNumTraceHdrExtensions > 0 ) {
        myTraceHdrExtension->extractHeaders( &myBuffer_oneRecord[bytePos+csTraceHeader::BLOCK_SIZE] );
        myTraceHdrExtension->dump( *outStream );
      }
      traceCounter += 1;
    }
  }
}
void csSegdReader::dumpTrace( std::ofstream* outStream ) {
  myTraceHdr->dump( *outStream );
  if( myNumTraceHdrExtensions > 0 ) {
    myTraceHdrExtension->dump( *outStream );
  }
}
//---------------------------------------------------------
//
void csSegdReader::dumpRecordHeaders( FILE* dumpFile ) {
  fprintf( dumpFile,"\n\n***** SEGD record headers *****\n" );
  fprintf( dumpFile,"FFID number                       : %d  (actually %d if 16665)\n", myGeneralHdr1->fileNum, myGeneralHdr2->expandedFileNum );
  fprintf( dumpFile,"Shot time (year day hour:min:sec) : %4d %3d %2d:%2d:%2d\n",
    myGeneralHdr1->year, myGeneralHdr1->julianDay, myGeneralHdr1->hour, myGeneralHdr1->min, myGeneralHdr1->sec );
}
//---------------------------------------------------------
//
void csSegdReader::dumpEssentialFileInfo( FILE* dumpFile ) {
  string manufacturerName;
  getManufacturerName( myGeneralHdr1->manufactCode, &manufacturerName );
  fprintf( dumpFile,"\n\n***** SEGD General information *****\n" );
  fprintf( dumpFile,"Format code                      : %d %s\n", myGeneralHdr1->formatCode,
           isFormatCodeSupported(myGeneralHdr1->formatCode ) ? "" : " (NOT SUPPORTED)" );
  fprintf( dumpFile,"Manufacturer code                : %d %s\n", myGeneralHdr1->manufactCode,
           isManufacturerSupported(myGeneralHdr1->manufactCode) ? "" : "(NOT SUPPORTED)" );
  fprintf( dumpFile,"                                   %s\n", manufacturerName.c_str() );
  fprintf( dumpFile,"SEGD revision                    : %d.%d %s\n", myComFileHdr.revisionNum[0], myComFileHdr.revisionNum[1],
           isRevisionSupported(myComFileHdr.revisionNum[0],myComFileHdr.revisionNum[1]) ? "" : "(NOT SUPPORTED)" );
  if( !myConfig.thisIsRev0 ) {
    fprintf( dumpFile,"Number of general hdr blocks     : 2 + %d\n", myNumExtraGeneralHdrBlocks);
  }
  else {
    fprintf( dumpFile,"Number of general hdr blocks     : %d\n", myNumGeneralHdrBlocks);
  }
  fprintf( dumpFile,"Number of scan types             : %d\n", myNumScanTypes);
  fprintf( dumpFile,"Number of channel sets/scan type : %d\n", myNumChanSetsPerScanType);
  fprintf( dumpFile,"Number of extended hdr blocks    : %d\n", myNumExtendedHdrBlocks);
  fprintf( dumpFile,"Number of external hdr blocks    : %d\n", myNumExternalHdrBlocks);
  fprintf( dumpFile,"Number of general hdr trailers   : %d\n", myNumGeneralTrailerBlocks);
  fprintf( dumpFile,"Number of sample skew blocks     : %d\n", myNumSampleSkewBlocks);
  fprintf( dumpFile,"Number of trace hdr extensions   : %d\n", myNumTraceHdrExtensions );
  fprintf( dumpFile,"Total number of traces           : %d\n", myComFileHdr.totalNumChan );
  fprintf( dumpFile,"\n" );
  fprintf( dumpFile,"Record length       : %d [ms]\n", myComFileHdr.recordLength_ms );
  fprintf( dumpFile,"Sample interval     : %d [us]\n", myComFileHdr.sampleInt_us );
  fprintf( dumpFile,"Number of samples   : %d\n", myComFileHdr.numSamples );
  fprintf( dumpFile,"Sample bit size     : %d\n", myComFileHdr.sampleBitSize );
  fprintf( dumpFile,"\n" );
  fprintf( dumpFile,"Number of bytes read in so far...: %d\n", myInitializeNumBytesRead );
  if( myIsSupported ) {
    int traceByteSize = traceDataByteSize() + csTraceHeader::BLOCK_SIZE +
      myNumTraceHdrExtensions * csBaseHeader::BLOCK_SIZE;
    fprintf( dumpFile,"Number of bytes in one full record : %d  (each trace: %d bytes = %d bits ?= %d bits)\n",
             myRecordByteSize, traceByteSize, traceByteSize*8, myChanSetNumSamples[0] * myComFileHdr.sampleBitSize );
  }
  else {
    fprintf( dumpFile,"...this SEGD file format is not supported\n" );
  }
  fprintf( dumpFile,"\n" );
}
//---------------------------------------------------------
//
void csSegdReader::dumpFileHeaders( int dumpFlag, std::ofstream* outStream ) {
  if( (dumpFlag & DUMP_GENERAL) && myGeneralHdr1 ) {
    *outStream << endl << "   General header 1   " << endl;
    myGeneralHdr1->dump( *outStream );
  }
  if( (dumpFlag & DUMP_GENERAL) && myGeneralHdr2 ) {
    *outStream << endl << "   General header 2   " << endl;
    myGeneralHdr2->dump( *outStream );
  }
  if( (dumpFlag & DUMP_GENERAL) && myGeneralHdrN ) {
    for( int iblock = 0; iblock < myNumExtraGeneralHdrBlocks; iblock++ ) {
      *outStream << endl << "   General header N (" << (iblock+1) << ")  " << endl;
      myGeneralHdrN[iblock].dump( *outStream );
    }
  }
  if( (dumpFlag & DUMP_CHANSET) && myChanSetHdr ) {
    for( int iscan = 0; iscan < myNumScanTypes; iscan++ ) {
      *outStream << endl << "   Scan type  #" << (iscan+1) << "  " << endl;
      for( int ichanset = 0; ichanset < myNumChanSetsPerScanType; ichanset++ ) {
        int index = iscan*myNumChanSetsPerScanType + ichanset;
        *outStream << endl << "   Chan set  #" << (ichanset+1) << "  " << endl;
        myChanSetHdr[index].dump( *outStream );
      }
    }
  }
  if( dumpFlag & DUMP_EXTENDED && myExtendedHdr ) {
    *outStream << endl << "   Extended header   " << endl;
    if( myExtendedHdr ) {
      myExtendedHdr->dump( *outStream );
    }
    else if( myBytePos.extendedHdr != 0 ) {
      if( dumpFlag & DUMP_RAW ) dumpRawHex( *outStream, &myBuffer_oneRecord[myBytePos.extendedHdr], myNumExtendedHdrBlocks*csBaseHeader::BLOCK_SIZE );
    }
    if( dumpFlag & DUMP_RAW ) dumpRawHex( *outStream, &myBuffer_oneRecord[myBytePos.extendedHdr], myNumExtendedHdrBlocks*csBaseHeader::BLOCK_SIZE );
  }
  if( dumpFlag & DUMP_EXTERNAL && myExternalHdr ) {
    *outStream << endl << "   External header   " << endl;
    if( myExternalHdr ) {
      myExternalHdr->dump( *outStream );
      if( dumpFlag & DUMP_RAW ) dumpRawASCII( *outStream, &myBuffer_oneRecord[myBytePos.externalHdr], myNumExternalHdrBlocks*csBaseHeader::BLOCK_SIZE );
    }
    else if( myBytePos.externalHdr != 0 ) {
      if( dumpFlag & DUMP_RAW ) dumpRawASCII( *outStream, &myBuffer_oneRecord[myBytePos.externalHdr], myNumExternalHdrBlocks*csBaseHeader::BLOCK_SIZE );
    }
  }
  if( (dumpFlag & DUMP_TRACE) && myTraceHdr ) {
    *outStream << endl << "   Trace header from first trace   " << endl;
    myTraceHdr->dump( *outStream );
  }
  if( (dumpFlag & DUMP_TRACE) && myTraceHdrExtension ) {
    *outStream << endl << "   Trace header extension from first trace   " << endl;
    myTraceHdrExtension->dump( *outStream );
  }
  if( myNumTraceHdrExtensions > 0 ) {
    *outStream << "-------- Trace header extension " << (myBytePos.firstTraceExt) << " Size: " << csBaseHeader::BLOCK_SIZE*myNumTraceHdrExtensions << endl;
    // TEMP COMMENTED OUT    dumpRawHex( *outStream, &myBuffer_oneRecord[myBytePos.firstTraceExt], csBaseHeader::BLOCK_SIZE*myNumTraceHdrExtensions );
  }
  //  *outStream << "-------- First trace data  " << myBytePos.firstTraceData << endl;
  // dumpRawHex( *outStream, &myBuffer_oneRecord[myBytePos.firstTraceData], 40 );

}
int csSegdReader::traceDataByteSize() const {
  return traceDataByteSize(0);
}

int csSegdReader::traceDataByteSize( int chanSetIndex ) const {
  if( !myConfig.thisIsRev0 ) {
    int byteSize = (myChanSetNumSamples[chanSetIndex] * myComFileHdr.sampleBitSize)/8;
    // byteSize *= 2;
    return byteSize;
  }
  else {
    int byteSize = (myComFileHdr.numSamples * myComFileHdr.sampleBitSize);
    byteSize *= 2; // Why double???
    byteSize /= 8; // convert to bytes
    byteSize += 3; // Why 3 bytes extra?
    return byteSize;
  }
}

void convertToFloat_20bit( short const* in, float* out, int numSamples ) {
  int counter_2bytes = 0;

  for( int isamp = 0; isamp < numSamples; isamp +=4 ) {
    short allExponents = (int)in[counter_2bytes++];
    int expo1 = (allExponents & 15);
    int expo2 = ((allExponents >> 4) & 15);
    int expo3 = ((allExponents >> 8) & 15);
    int expo4 = ((allExponents >> 12) & 15);
    short frac1 = in[counter_2bytes++];
    short frac2 = in[counter_2bytes++];
    short frac3 = in[counter_2bytes++];
    short frac4 = in[counter_2bytes++];

    if( frac1 < 0 ) frac1 = (short)-(~frac1);
    if( frac2 < 0 ) frac2 = (short)-(~frac2);
    if( frac3 < 0 ) frac3 = (short)-(~frac3);
    if( frac4 < 0 ) frac4 = (short)-(~frac4);

    out[isamp]   = ldexp( (double)frac1, expo1 );
    out[isamp+1] = ldexp( (double)frac2, expo2 );
    out[isamp+2] = ldexp( (double)frac3, expo3 );
    out[isamp+3] = ldexp( (double)frac4, expo4 );

    // fprintf(stdout,"Expo: %4x -- %4d %4d %4d %4d -- %6d %6d %6d %6d -- %f %f %f %f\n",
    // allExponents, expo1, expo2, expo3, expo4, frac1, frac2, frac3, frac4, out[isamp], out[isamp+1], out[isamp+2], out[isamp+3]);

  }
}

void convertToFloat_24bit( unsigned char const* in, float* out, int numSamples ) {
  int counter = 0;
  for( int isamp = 0; isamp < numSamples; isamp++ ) {
    unsigned long accum;
    accum = (unsigned long)in[counter+0] << 16; // Top most byte
    accum |= (unsigned)in[counter+1] << 8;
    accum |= in[counter+2];
    long val = (long)(accum ^ SIGNBIT) - (long)SIGNBIT;
    out[isamp] = (float)val;
    counter += 3;
  }
}

int csSegdReader::numTraces() const {
  if( myConfig.readAuxTraces ) {
    return myComFileHdr.totalNumChan;
  }
  else {
    return myComFileHdr.numSeismicChan;
  }
}
float const* csSegdReader::getNextTracePointer( commonTraceHeaderStruct& comTrcHdr ) {
  if( myCopyBuffer == NULL ) {
    myCopyBuffer = new float[myComFileHdr.numSamples];
  }
  if( getNextTrace( myCopyBuffer, comTrcHdr ) ) return myCopyBuffer;
  else return NULL;
}
void csSegdReader::retrieveChanSetInfo( int chanSetIndex, commonChanSetStruct& info ) const {
  if( chanSetIndex >= numChanSets() ) throw( csException("csSegdReader::retrieveChanSetInfo: Wrong chan set index passed: %d\n", chanSetIndex) );
  info.numSamples   = myChanSetNumSamples[chanSetIndex];
  info.sampleInt_us = myChanSetSampleInt_us[chanSetIndex];
  info.numChannels  = myChanSetHdr[chanSetIndex].numChannels;
}

