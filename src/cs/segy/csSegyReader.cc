/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <string>
#include <cstdio>
#include <cstring>
#include "csSegyReader.h"
#include "csSegyHdrMap.h"
#include "csSegyTraceHeader.h"
#include "csSegyDefines.h"
#include "csSegyHeaderInfo.h"
#include "csSegyBinHeader.h"
#include "csHeaderInfo.h"
#include "csVector.h"
#include "csException.h"
#include "csFlexNumber.h"
#include "csFlexHeader.h"
#include "csGeolibUtils.h"
#include "csByteConversions.h"
#include "geolib_endian.h"
#include "methods_number_conversions.h"
#include "csFileUtils.h"
#include "geolib_defines.h"
#include "csSegyHdrMap.h"
#include "csIOSelection.h"

using namespace cseis_geolib;

csSegyReader::csSegyReader( std::string filename, SegyReaderConfig const& config, csSegyHdrMap const* hdrMap ) {
  myFilename     = filename;
  myFileSize     = csFileUtils::FILESIZE_UNKNOWN;
  myCharHdrBlock = NULL;
  myBinHdrBlock  = NULL;
  myIsOverrideNumSamples = false;

  myTrcHdr       = NULL;
  myTrcHdrMap    = NULL;
  myBinHdr       = NULL;
  myBigBuffer    = NULL;
  myCopyBuffer   = NULL;
  myFile         = NULL;
  myIOSelection     = NULL;
  myNumTraces       = 0;
  myLastTraceIndex  = 0;
  myTraceByteSize   = 0;
  myCurrentTraceInFile = 0;  // Pointer to current trace in input file
  myCurrentTraceRead   = 0;  // Current trace index (read from either buffer or input file)
  myBufferCurrentTrace = 0;  // Pointer to current trace in buffer
  myBufferNumTraces    = 0;  // Number of traces currently stored in buffer


  myResidualNumBytesAtEnd   = 0;
  myIsEBCDIC                = true;
  myHasBeenInitialized      = false;
  myIsAutoScaleHeaders      = config.autoscaleHdrs;
  myBufferCapacityNumTraces = config.numTracesBuffer;
  myOverrideSampleFormat    = config.overrideSampleFormat;
  myEnableRandomAccess      = config.enableRandomAccess;
  myIsSUFormat              = config.isSUFormat;

  myHdrCheckByteOffset   = 0;
  myHdrCheckInType       = cseis_geolib::TYPE_UNKNOWN;
  myHdrCheckOutType       = cseis_geolib::TYPE_UNKNOWN;
  myHdrCheckByteSize   = 0;
  myHdrCheckBuffer     = NULL;
  myPeekIsInProgress = false;
  myCurrentPeekTraceIndex = 0;
  myCurrentPeekByteOffset = 0;
  myCurrentPeekByteSize   = 0;
  myIsVariableTraceLength = false;

  if( myOverrideSampleFormat != cseis_segy::AUTO &&
     (myOverrideSampleFormat != cseis_segy::DATA_FORMAT_IEEE &&
      myOverrideSampleFormat != cseis_segy::DATA_FORMAT_IBM  &&
      myOverrideSampleFormat != cseis_segy::DATA_FORMAT_INT32  &&
      myOverrideSampleFormat != cseis_segy::DATA_FORMAT_INT16) ) {
    throw csException( "Overriden SEGY data sample format not supported: %d  (Supported formats are: 1:IBM, 2: 32bit INT, 3: 16bit INT, 5:IEEE)\nFile name: %s", myOverrideSampleFormat, myFilename.c_str() );
  }
  else if( myOverrideSampleFormat == cseis_segy::DATA_FORMAT_INT16 ) {
    // Override number of traces to read in at once. For INT16, it is required that only one trace is read in at once
    myBufferCapacityNumTraces = 1;
  }
  if( hdrMap == NULL ) {
    myTrcHdrMap     = new csSegyHdrMap(config.segyHeaderMapping,false);
  }
  else {
    myTrcHdrMap = new csSegyHdrMap( hdrMap );
  }

  myDoSwapEndianData = isPlatformLittleEndian();
  myDoSwapEndianHdr  = isPlatformLittleEndian();
  if( config.reverseByteOrderHdr ) myDoSwapEndianHdr = !myDoSwapEndianHdr;
  if( config.reverseByteOrderData ) myDoSwapEndianData = !myDoSwapEndianData;

  openFile();

  myCharHdrBlock = new char[cseis_segy::SIZE_CHARHDR];
  myBinHdrBlock  = new byte_t[cseis_segy::SIZE_BINHDR];
  myBinHdr       = new csSegyBinHeader( myDoSwapEndianHdr );
  for( int i = 0; i < cseis_segy::SIZE_CHARHDR; i++ ) {
    myCharHdrBlock[i] = 0;
  }
  for( int i = 0; i < cseis_segy::SIZE_BINHDR; i++ ) {
    myBinHdrBlock[i] = 0;
  }
}
//-----------------------------------------------------------------------------------------
csSegyReader::~csSegyReader() {
  if( myIOSelection != NULL ) {
    delete myIOSelection;
    myIOSelection = NULL;
  }
  freeCharBinHdr();
  if( myBigBuffer ) {
    delete [] myBigBuffer;
    myBigBuffer = NULL;
  }
  if( myFile != NULL ) {
    closeFile();
    delete myFile;
    myFile = NULL;
  }
  if( myTrcHdrMap ) {
    delete myTrcHdrMap;
    myTrcHdrMap = NULL;
  }
  if( myTrcHdr ) {
    delete myTrcHdr;
    myTrcHdr = NULL;
  }
  if( myBinHdr ) {
    delete myBinHdr;
    myBinHdr = NULL;
  }
  if( myCopyBuffer ) {
    delete [] myCopyBuffer;
    myCopyBuffer = NULL;
  }
  if( myHdrCheckBuffer != NULL ) {
    delete [] myHdrCheckBuffer;
    myHdrCheckBuffer = NULL;
  }    
}
//-----------------------------------------------------------------------------------------
void csSegyReader::setVariableTraceLength() {
  myIsVariableTraceLength = true;
}
//-----------------------------------------------------------------------------------------
void csSegyReader::openFile() {
  myFile = new std::ifstream();
  myFile->open( myFilename.c_str(), std::ios::in | std::ios::binary );
  if( myFile->fail() ) {
    throw csException("Could not open SEGY file: %s", myFilename.c_str());
  }

  // Determine file size
  if( myEnableRandomAccess ) {
    try {
      myFileSize = csFileUtils::retrieveFileSize( myFilename );
    }
    catch( csException& e ) {
      throw csException("Error occurred while determining file size. System message: %s\n", e.getMessage() );
    }
  }
}
//-----------------------------------------------------------------------------------------
void csSegyReader::closeFile() {
  if( myFile != NULL ) {
    myFile->close();
    // Do not clear ios flags..
  }
}
//-----------------------------------------------------------------------------------------
void csSegyReader::initialize( int numSamples ) {
  myIsOverrideNumSamples = true;
  myNumSamples = numSamples;
  initialize();
}
void csSegyReader::initialize() {

  if( myHasBeenInitialized ) return;

  // First, read char & bin hdrs. This is necessary to set total trace size etc..
  readCharBinHdr();

  if( myBufferCapacityNumTraces <= 0 ) {  // Set number of buffered traces if it wasn't set explicitely before
    myBufferCapacityNumTraces = DEFAULT_BUFFERED_SAMPLES / myNumSamples;
    if( myBufferCapacityNumTraces <= 0 ) myBufferCapacityNumTraces = 1;
    else if( myBufferCapacityNumTraces > 20 ) myBufferCapacityNumTraces = 20;
  }

  myTraceByteSize  = myNumSamples*mySampleByteSize+cseis_segy::SIZE_TRCHDR;
  if( myFileSize != csFileUtils::FILESIZE_UNKNOWN ) {
    long long sizeTraces  = myFileSize - ( (long long)cseis_segy::SIZE_CHARHDR + (long long)cseis_segy::SIZE_BINHDR );
    if( myIsSUFormat ) sizeTraces = myFileSize;  // Special for SU format: No bin & char header
    myNumTraces      = (int)(sizeTraces/(long long)myTraceByteSize);  // Total number of full traces in input file
    myLastTraceIndex = myNumTraces-1;
    myResidualNumBytesAtEnd = (int)(sizeTraces - (long long)myNumTraces*(long long)myTraceByteSize);  // Spurious bytes at end of file. File corrupted?

    // Do not allocate buffer that is larger than input file!
    if( myNumTraces < myBufferCapacityNumTraces ) {
      myBufferCapacityNumTraces = myNumTraces;
    }
    if( myNumTraces == 0 ) {
      throw( csException("SEGY file does not contain any traces. File corrupted..?") );
    }
  }

  myBigBuffer = new char[ myBufferCapacityNumTraces * myTraceByteSize ];
  if( !myBigBuffer ) {
    throw( csException("csSegyReader::initialize(): Not enough memory...") );
  }
  
  myTrcHdrMap->initScalars();
  myTrcHdr = new csSegyTraceHeader( myTrcHdrMap );

  myHasBeenInitialized = true;
}

//*******************************************************************
//
// Char & bin headers
//
//*******************************************************************
void csSegyReader::setCharHdrFormat( bool isEBCDIC ) {
  myIsEBCDIC = isEBCDIC;
}

void csSegyReader::readCharBinHdr() {
  if( myTrcHdrMap->mapID() != csSegyHdrMap::SEGY_PSEGY && !myIsSUFormat ) {
    myFile->read( myCharHdrBlock, cseis_segy::SIZE_CHARHDR );
    if( myFile->fail() ) {
      throw csException("Unexpected error occurred when reading SEGY char header");
    }

    if( myIsEBCDIC ) {
      for( int i = 0; i < cseis_segy::SIZE_CHARHDR; i++) {
        myCharHdrBlock[i] = ebcdic2char( (short)myCharHdrBlock[i] );
      }
    }

    myFile->read( (char*)myBinHdrBlock, cseis_segy::SIZE_BINHDR );
    if( myFile->fail() ) {
      throw csException("Unexpected error occurred when reading SEGY binary header");
    }

    myBinHdr->extractHeaders( myBinHdrBlock );

    myDataSampleFormat = myBinHdr->dataSampleFormat;
    mySampleInt        = (float)myBinHdr->sampleIntUS * 0.001f;
    if( !myIsOverrideNumSamples ) {
      myNumSamples       = myBinHdr->numSamples;
    }

  }
  //---------------------------------------------------------------------------------
  // SU format does not have bin & char headers. Instead, number of samples and sample interval is read from first trace header
  else if( myIsSUFormat ) {
    // Determine file size
    if( !myEnableRandomAccess ) {
      // Random access needed to peek number of samples & sample interval from first trace
      myEnableRandomAccess = true;
      try {
        myFileSize = csFileUtils::retrieveFileSize( myFilename );
      }
      catch( csException& e ) {
        throw csException("Error occurred while determining file size. System message: %s\n", e.getMessage() );
      }
    }
    myHasBeenInitialized = true; // Trick to make peek function work. Reset afterwards

    myDataSampleFormat = cseis_segy::DATA_FORMAT_IEEE;
    cseis_geolib::csFlexHeader hdrValue;
    bool success;
    if( !myIsOverrideNumSamples ) {
      if( myTrcHdrMap->mapID() != csSegyHdrMap::SEGY_SU_ONLY ) {      
        setHeaderToPeek( "nsamp" );
      }
      else {
        setHeaderToPeek( "ns" );
      }
      success = peekHeaderValue( &hdrValue, 0 );
      if( !success ) throw csException("Could not read number of samples from SU file: %s", myFilename.c_str());
      myNumSamples = hdrValue.intValue();
    }
    if( myTrcHdrMap->mapID() != csSegyHdrMap::SEGY_SU_ONLY ) {
      setHeaderToPeek( "sampint_us" );
    }
    else {
      setHeaderToPeek( "dt" );
    }
    success = peekHeaderValue( &hdrValue, 0 );
    if( !success ) throw csException("Could not read sample interval from SU file: %s", myFilename.c_str());
    mySampleInt = (float)hdrValue.intValue() * 0.001f;

    myHasBeenInitialized = false;

    // Set fake SU bin header information:
    myBinHdr->dataSampleFormat = (short)myDataSampleFormat;
    myBinHdr->numSamples       = myNumSamples;
    myBinHdr->sampleIntUS      = (int)(mySampleInt*1000);
  }
  //---------------------------------------------------------------------------------
  else {  // PSEGY 'special' format
    for( int i = 0; i < cseis_segy::SIZE_CHARHDR; i++ ) {
      myCharHdrBlock[i] = ' ';
    }
    byte_t* trcHdrBlock = new byte_t[cseis_segy::SIZE_TRCHDR];
    
    myFile->read( (char*)trcHdrBlock, cseis_segy::SIZE_TRCHDR );
    if( myFile->fail() ) {
      throw csException("Unexpected error occurred when reading P-SEGY trace header block");
    }

    int numSamples_2byte = byte2Short( &trcHdrBlock[114] );
    int sampleInt_2byte  = byte2Short( &trcHdrBlock[116] );
    int numSamples_4byte = byte2Int( &trcHdrBlock[228] );
    int sampleInt_4byte  = byte2Int( &trcHdrBlock[200] );
    int dataFormat       = byte2Short( &trcHdrBlock[204] );

    if( myDoSwapEndianHdr ) {
      numSamples_2byte = byte2Short_SWAP( &trcHdrBlock[114] );
      sampleInt_2byte  = byte2Short_SWAP( &trcHdrBlock[116] );
      numSamples_4byte = byte2Int_SWAP( &trcHdrBlock[228] );
      sampleInt_4byte  = byte2Int_SWAP( &trcHdrBlock[200] );
      dataFormat       = byte2Short_SWAP( &trcHdrBlock[204] );
    }

    if( sampleInt_2byte != 1 ) { 
      mySampleInt = (float)sampleInt_2byte * 0.001f;
    }
    else {
      mySampleInt = (float)sampleInt_4byte * 0.001f;
    }
    if( myIsOverrideNumSamples ) {
      // Nothing
    }
    else if( numSamples_2byte < 32767 && numSamples_2byte > 0 ) { 
      myNumSamples = numSamples_2byte;
    }
    else {
      myNumSamples = numSamples_4byte;
    }
    if( dataFormat == 0 ) {  // 16bit int
      myDataSampleFormat = cseis_segy::DATA_FORMAT_INT16;
//      if( myDataSampleFormat != myBinHdr->dataSampleFormat ) {  // comment out: binary header of PSEGY contains only garbage
//        throw( csException("PSEGY file: Inconsistent data sample format in trace header (=%d) and binary header (=%d)", myDataSampleFormat, myBinHdr->dataSampleFormat) );
//      }
    }
    else if( dataFormat == 1 ) {  // 32bit int
      myDataSampleFormat = cseis_segy::DATA_FORMAT_INT32;
    }
    else {
      myDataSampleFormat = -1;
    }
    // Endian format of input data file seems to be wrong...

    /*
    fprintf(stdout,"PSEGY: Standard 2 byte headers --> PSEGY 4 byte headers\n" );
    fprintf(stdout," Sample interval   : %d  -->  %d\n", sampleInt_2byte, sampleInt_4byte );
    fprintf(stdout," Number of samples : %d  -->  %d\n", numSamples_2byte, numSamples_4byte );
    fprintf(stdout," Data sample format: %d  -->  %d\n", dataFormat, myDataSampleFormat );
    */
    delete [] trcHdrBlock;
    myFile->clear();
    myFile->seekg(0, std::ios::beg);

    myBinHdr->dataSampleFormat = (short)myDataSampleFormat;
    myBinHdr->numSamples       = myNumSamples;
    myBinHdr->sampleIntUS      = (int)(mySampleInt*1000);
  }

  //  if( mySampleInt <= 0 ) {
  // Do not throw exception. Sample interval is not really needed in segy reader... Calling program should check if it is correct
  // }
  if( myOverrideSampleFormat == cseis_segy::AUTO &&
     (myDataSampleFormat != cseis_segy::DATA_FORMAT_IEEE &&
      myDataSampleFormat != cseis_segy::DATA_FORMAT_IBM  &&
      myDataSampleFormat != cseis_segy::DATA_FORMAT_INT32  &&
      myDataSampleFormat != cseis_segy::DATA_FORMAT_INT16) ) {
    throw csException( "SEGY data sample format not supported: %d  (Supported formats are: 1:IBM, 2: 32bit INT, 3: 16bit INT, 5:IEEE)\nFile name: %s", myDataSampleFormat, myFilename.c_str() );
  }
  else if( myNumSamples <= 0 && !myIsOverrideNumSamples ) {
    throw csException( "SEGY file corrupted: non-physical number of samples: %d. Try reversing the byte order.", myNumSamples );
  }
  if( myOverrideSampleFormat != cseis_segy::AUTO ) {
    myDataSampleFormat = myOverrideSampleFormat;
  }

  if( myDataSampleFormat != cseis_segy::DATA_FORMAT_INT16 ) {
    mySampleByteSize = 4;
  }
  else {
    mySampleByteSize = 2;
  }
}
//-----------------------------------------------------------------------------------------
void csSegyReader::freeCharBinHdr() {
  if( myCharHdrBlock ) {
    delete [] myCharHdrBlock;
    myCharHdrBlock = NULL;
  }
  if( myBinHdrBlock ) {
    delete [] myBinHdrBlock;
    myBinHdrBlock = NULL;
  }
}

//*******************************************************************
//
// Reading & retrieving traces
// Random access/moving to specified trace index
//
//*******************************************************************

//-----------------------------------------------------
//
bool csSegyReader::moveToTrace( int traceIndex ) {
  return moveToTrace( traceIndex, myNumTraces-traceIndex );
}
bool csSegyReader::moveToTrace( int traceIndex, int numTracesToRead ) {
  //  fprintf(stdout,"SEGY movetotrace IN: %d %d %d\n", myCurrentTraceInFile, myBufferNumTraces, traceIndex );
  //  fflush(stdout);

  if( myEnableRandomAccess == false ) {
    throw( csException("csSegyReader::moveToTrace: Random access not enabled. Set enableRandomAccess to true. This is a program bug in the calling function." ) );
  }
  else if( myFileSize == csFileUtils::FILESIZE_UNKNOWN ) {
    throw( csException("csSegyReader::moveToTrace: File size unknown. This may be due to a compatibility problem of this compiled version of the program on the current platform." ) );
  }
  else if( traceIndex < 0 || traceIndex >= myNumTraces ) {
    throw( csException("csSegyReader::moveToTrace: Incorrect trace index: %d (number of traces in input file: %d). This is a program bug in the calling method", traceIndex, myNumTraces) );
  }
  else if( !myHasBeenInitialized ) {
    initialize();
//    throw( csException("csSegyReader::moveToTrace: Called method before initializing SEGY reader. This is a program bug in the calling method") );
  }
  if (myPeekIsInProgress ) revertFromPeekPosition();

  //  fprintf(stdout,"SEGY movetotrace: %d %d %d\n", myCurrentTraceInFile, myBufferNumTraces, traceIndex );
  //  fflush(stdout);
  myBufferNumTraces = 0;   // New line since v1.80: Not sure if this line belongs here or if it was only for testing
  if( myCurrentTraceInFile-myBufferNumTraces != traceIndex ) {
    csInt64_t bytePosRelative = (csInt64_t)(traceIndex-myCurrentTraceInFile) * (csInt64_t)myTraceByteSize;
    if( (traceIndex > myCurrentTraceInFile && bytePosRelative < 0) ||
        (traceIndex < myCurrentTraceInFile && bytePosRelative > 0) ) {
      return false;
    }
    myFile->clear(); // Clear all flags
    if( !csFileUtils::seekg_relative( bytePosRelative, myFile ) ) return false;
    if( myFile->fail() ) return false;
    //    fprintf(stdout,"SEGY movetotrace: RESET myBufferNumTraces\n");
    //    fflush(stdout);
    myBufferNumTraces = 0;
  }
  myCurrentTraceInFile  = traceIndex;
  myBufferCurrentTrace = 0;
  myLastTraceIndex = std::min( traceIndex + numTracesToRead - 1, myNumTraces-1 );  // Index of last trace that will be read in one go by consecutive calls to getNextTrace()
  return !myFile->fail();
}

//-----------------------------------------------------
// This method expects a fully allocated buffer as input, which will be overwritten
//
bool csSegyReader::getNextTrace( byte_t* sampleBufferOut ) {
  return getNextTrace( sampleBufferOut, myNumSamples );
}

//--------------------------------------------------------------------------------
bool csSegyReader::getNextTrace( byte_t* sampleBufferOut, int nSamples ) {
  if( nSamples <= 0 || nSamples > myNumSamples ) {
    throw( csException("csSegyReader::getNextTrace: Incorrect number of samples: %d (true number of samples: %d). This is a program bug in the calling method", nSamples, myNumSamples) );
  }
  else if( !myHasBeenInitialized ) {
    initialize();
  }
  if (myPeekIsInProgress ) revertFromPeekPosition();

  if( myIOSelection ) {
    int traceIndex = myIOSelection->getNextTraceIndex();
    if( traceIndex < 0 ) return false;
    bool success = moveToTrace( traceIndex );
    if( !success ) return false;
  }

  //  fprintf(stdout,"SEGY %d %d %d %d   %d\n", myBufferCurrentTrace, myBufferNumTraces, myCurrentTraceInFile, myNumTraces, myFileSize );
  //========================================================
  // START: Read in new traces into big buffer
  // Read in as many traces as possible at once, to avoid extensive I/O
  //
  //  fprintf(stdout,"SEGY %d %d:  %d  %d\n", myBufferCurrentTrace, myCurrentTraceInFile, myBufferNumTraces, myResidualNumBytesAtEnd );
  //  fflush(stdout);

  if( myBufferCurrentTrace == myBufferNumTraces ) {
    if( myFile->eof() ) {
      return false;
    }

    if( myFileSize == csFileUtils::FILESIZE_UNKNOWN ) {

      //      fprintf(stderr,"Segyreader: Read new traces %d\n", myBufferNumTraces);
      myBufferNumTraces = myBufferCapacityNumTraces;

      if( !myIsVariableTraceLength ) {
        myFile->read( myBigBuffer, myTraceByteSize*myBufferNumTraces );
      }
      else {
        int numBytes = 118;
        myFile->read( myBigBuffer, numBytes );
        int nsamp   = 0;
        if( myDoSwapEndianHdr ) {
          nsamp = (int)(unsigned short)byte2Short_SWAP( (byte_t*)(&myBigBuffer[114]) );
        }
        else {
          nsamp = (int)(unsigned short)byte2Short( (byte_t*)(&myBigBuffer[114]) );
        }
        myFile->read( &myBigBuffer[numBytes], (240-numBytes)+4*nsamp );
      }
      /*      else {  // Read in corrupted file
        int numBytes = 118;
        myFile->read( myBigBuffer, numBytes );
        int nsamp   = 0;
        int sampint = 0;
        bool found = false;
        do {        
          if( myDoSwapEndianHdr ) {
            nsamp = (int)(unsigned short)byte2Short_SWAP( (byte_t*)(&myBigBuffer[114]) );
            sampint = (int)(unsigned short)byte2Short_SWAP( (byte_t*)(&myBigBuffer[116]) );
          }
          else {
            nsamp = (int)(unsigned short)byte2Short( (byte_t*)(&myBigBuffer[114]) );
            sampint = (int)(unsigned short)byte2Short( (byte_t*)(&myBigBuffer[116]) );
          }
          fprintf(stdout,"%4d  NSAMP = %d, SAMPINT = %d\n", myTEMPNumTracesRead, nsamp, sampint );
          if( nsamp == 2277 && sampint == 4000 ) {
            break;
          }
          for( int i = 0; i < 116; i++ ) {
            myBigBuffer[i] = myBigBuffer[i+2];
          }
          myFile->read( &myBigBuffer[116], 2 );
        } while( !found );
        myFile->read( &myBigBuffer[118], (240-numBytes) + 2277*4 );
        myTEMPNumTracesRead += 1;
        } */
      //      fprintf(stdout,"SEGY %d %d  read: %d in\n", myBufferCurrentTrace, myCurrentTraceInFile, myCurrentTraceRead );
      //      fflush(stdout);
      if( myFile->fail() ) {
        //        fprintf(stdout,"SEGY %d %d read: %d fail\n", myBufferCurrentTrace, myCurrentTraceInFile, myCurrentTraceRead );
        //        fflush(stdout);
        if( myFile->eof() ) {
          int numBytesRead        = (int)myFile->gcount();
          myBufferNumTraces       = numBytesRead / myTraceByteSize;
          myResidualNumBytesAtEnd = numBytesRead % myTraceByteSize;
          //          fprintf(stdout,"SEGY %d %d eof:  %d  %d\n", myBufferCurrentTrace, myCurrentTraceInFile, myBufferNumTraces, myResidualNumBytesAtEnd );
          //          fflush(stdout);
          if( myBufferNumTraces == 0 ) return false;
        }
        else {
          //          fprintf(stdout,"SEGY %d %d read: %d NOT eof\n", myBufferCurrentTrace, myCurrentTraceInFile, myCurrentTraceRead );
          //          fflush(stdout);
          closeFile();
          throw( csException("csSegyReader::getNextTrace: Unexpected error occurred when reading in data from input file '%s'", myFilename.c_str()) );
        }
      }
    }
    else { // file size is known == random access enabled
      if( myCurrentTraceInFile == myNumTraces ) return false;

      // Set myBufferNumTraces: Number of traces to be read into buffer
      if( myCurrentTraceInFile < myLastTraceIndex ) {
        // Make sure only as many traces are read in as necessary:
        myBufferNumTraces = std::min( myBufferCapacityNumTraces, myLastTraceIndex-myCurrentTraceInFile+1 );
      }
      else {
        myBufferNumTraces = myBufferCapacityNumTraces;
      }
      if( myCurrentTraceInFile+myBufferNumTraces > myNumTraces ) myBufferNumTraces = myNumTraces - myCurrentTraceInFile;
      //    fprintf(stdout,"Reading in traces: %d (saved: %d)\n", myBufferCurrentTrace, myBufferNumTraces );

      // The following read statement only makes sense if the user specified argument 'nSamples' is not too different from myNumSamples
      // Otherwise (when nSamples << myNumSamples), the read process may be sped up significantly by skipping over the bytes not needed.
      myFile->clear(); // Clear all flags
      myFile->read( myBigBuffer, myTraceByteSize*myBufferNumTraces );

      if( myFile->fail() ) {
        //        fprintf(stdout,"SEGY %d %d fail2: %d\n", myBufferCurrentTrace, myCurrentTraceInFile, myCurrentTraceRead );
        //        fflush(stdout);
        closeFile();
        throw( csException("csSegyReader::getNextTrace: Unexpected error occurred when reading in data from input file '%s'", myFilename.c_str()) );
      }
      else if( myFile->eof() ) {
        fprintf(stdout,"End of file reached...\n" );
        fflush(stdout);
      }
    }
    myBufferCurrentTrace = 0;
    myCurrentTraceInFile += myBufferNumTraces;

    // Perform ENDIAN processing if necessary
    // Convert samples from endian format in input file to internal endian format
    if( myDoSwapEndianData ) {
      if( myDataSampleFormat != cseis_segy::DATA_FORMAT_INT16 ) {
        for( int itrc = 0; itrc < myBufferNumTraces; itrc++ ) {
          swapEndian4( myBigBuffer+myTraceByteSize*itrc+cseis_segy::SIZE_TRCHDR, nSamples*mySampleByteSize );
        }
      }
      else {
        for( int itrc = 0; itrc < myBufferNumTraces; itrc++ ) {
          swapEndian2( myBigBuffer+myTraceByteSize*itrc+cseis_segy::SIZE_TRCHDR, nSamples*mySampleByteSize );
        }
      }
    }

    // Convert sample value format if necessary
    if( myDataSampleFormat == cseis_segy::DATA_FORMAT_IBM ) {
      for( int itrc = 0; itrc < myBufferNumTraces; itrc++ ) {
        ibm2ieee( (unsigned char*)(myBigBuffer+myTraceByteSize*itrc+cseis_segy::SIZE_TRCHDR), nSamples );
      }
    }
    else if( myDataSampleFormat == cseis_segy::DATA_FORMAT_IEEE ) {
      // Nothing to be done here...
    }
    else if( myDataSampleFormat== cseis_segy::DATA_FORMAT_INT32 ) {
      for( int itrc = 0; itrc < myBufferNumTraces; itrc++ ) {
        convertInt2Float( (int*)(myBigBuffer+myTraceByteSize*itrc+cseis_segy::SIZE_TRCHDR), nSamples );
      }
    }
    else if( myDataSampleFormat== cseis_segy::DATA_FORMAT_INT16 ) {
      //For 16bit files, only single trace is read in at once
      convertShort2Float( (short*)(myBigBuffer+cseis_segy::SIZE_TRCHDR), (float*)sampleBufferOut, nSamples );
    }
  }
  //
  // END: Read in new traces into big buffer
  //========================================================

  // Set output buffer
  if( myDataSampleFormat != cseis_segy::DATA_FORMAT_INT16 ) {
    memcpy( sampleBufferOut, &myBigBuffer[myBufferCurrentTrace*myTraceByteSize+cseis_segy::SIZE_TRCHDR], nSamples*mySampleByteSize );
  }

  // Extract header values
  byte_t* trcHdrPtr = reinterpret_cast<byte_t*>( &myBigBuffer[myBufferCurrentTrace*myTraceByteSize] );
  myTrcHdr->readHeaderValues( trcHdrPtr, myDoSwapEndianHdr, myIsAutoScaleHeaders );

  myBufferCurrentTrace += 1;
  return true;
}
// This method returns a constant pointer to the trace buffer
//
float const* csSegyReader::getNextTracePointer() {
  return getNextTracePointer( myNumSamples );
}
float const* csSegyReader::getNextTracePointer( int nSamples ) {
  if( myCopyBuffer == NULL ) {
    myCopyBuffer = new byte_t[myNumSamples*mySampleByteSize];
  }
  if( getNextTrace( myCopyBuffer, nSamples ) ) {
    return (float*)myCopyBuffer;
  }
  else {
    return NULL;
  }
}
//-------------------------------------------------------------
//
void csSegyReader::dumpTrcHdr( FILE* fout ) const {
  if( !myBigBuffer ) return;
  byte_t* trcHdrPtr = reinterpret_cast<byte_t*>( &myBigBuffer[myBufferCurrentTrace*myTraceByteSize] );
  csSegyTraceHeader::dump( trcHdrPtr, myDoSwapEndianHdr, fout );
}
//--------------------------------------------------------------
//

bool csSegyReader::setHeaderToPeek(std::string const& headerName) {
  int idx = myTrcHdrMap->headerIndex(headerName);    // returns -1 on failure

  if (idx < 0) {
    return false;
  }
  else {
    revertFromPeekPosition();
    cseis_geolib::csSegyHeaderInfo const* info = myTrcHdrMap->header(idx);
    myHdrCheckInType     = info->inType;
    myHdrCheckOutType    = info->outType;
    myHdrCheckByteOffset = info->byteLoc;
    
    if( info->byteSize != myHdrCheckByteSize && myHdrCheckBuffer != NULL ) {
      delete [] myHdrCheckBuffer;
      myHdrCheckBuffer = NULL;
    }
    if( myHdrCheckBuffer == NULL ) {
      myHdrCheckByteSize = info->byteSize;
      myHdrCheckBuffer = new char[myHdrCheckByteSize];
    }
    return true;
  }
}
bool csSegyReader::setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType ) {
  bool success = setHeaderToPeek( headerName );
  headerType = myHdrCheckOutType;
  return success;
}
int csSegyReader::getPeekByteOffset() const {
  return myHdrCheckByteOffset;
}
bool csSegyReader::peekHeaderValue(cseis_geolib::csFlexHeader* hdrValue, int traceIndex) {
  if (myHdrCheckBuffer == NULL) {
    throw(cseis_geolib::csException("csSegyReader::peekHeaderValue: No header has been set for checking. This is a program bug in the calling function"));
  }
  
  bool success = true;
  if (traceIndex < 0) {
    success = peek(myHdrCheckByteOffset, myHdrCheckByteSize, myHdrCheckBuffer);
  }
  else {
    success = peek(myHdrCheckByteOffset, myHdrCheckByteSize, myHdrCheckBuffer, traceIndex);
  }
  
  if (success) {
    if( myHdrCheckInType == cseis_geolib::TYPE_INT ) {
      if( myDoSwapEndianHdr ) {
        if( myHdrCheckByteSize == 2 ) {
          hdrValue->setIntValue( (int)byte2Short_SWAP( (byte_t*)myHdrCheckBuffer ) );
        }
        else if( myHdrCheckByteSize == 4 ) {
          hdrValue->setIntValue( byte2Int_SWAP( (byte_t*)myHdrCheckBuffer ) );
        }
      }
      else {
        if( myHdrCheckByteSize == 2 ) {
          hdrValue->setIntValue( (int)byte2Short( (byte_t*)myHdrCheckBuffer ) );
        }
        else if( myHdrCheckByteSize == 4 ) {
          hdrValue->setIntValue( byte2Int( (byte_t*)myHdrCheckBuffer ) );
        }
      }
    }
    else if( myHdrCheckInType == cseis_geolib::TYPE_USHORT ) {
      if( myDoSwapEndianHdr ) {
        hdrValue->setIntValue( (int)(unsigned short)byte2Short_SWAP( (byte_t*)myHdrCheckBuffer ) );
      }
      else {
        hdrValue->setIntValue( (int)(unsigned short)byte2Short( (byte_t*)myHdrCheckBuffer ) );
      }
    }
    else if( myHdrCheckInType == cseis_geolib::TYPE_FLOAT ) {
      if( myDoSwapEndianHdr ) {
        hdrValue->setFloatValue( byte2Float_SWAP( (byte_t*)myHdrCheckBuffer ) );
      }
      else {
        hdrValue->setFloatValue( byte2Float( (byte_t*)myHdrCheckBuffer ) );
      }
    }
    else if( myHdrCheckInType == cseis_geolib::TYPE_DOUBLE ) {
      if( myDoSwapEndianHdr ) {
        hdrValue->setDoubleValue( (double)byte2Int_SWAP( (byte_t*)myHdrCheckBuffer ) );
      }
      else {
        hdrValue->setDoubleValue( (double)byte2Int( (byte_t*)myHdrCheckBuffer ) );
      }
    }
    else if( myHdrCheckInType == cseis_geolib::TYPE_STRING ) {
      char* text1 = new char[myHdrCheckByteSize+1];
      memcpy(text1,myHdrCheckBuffer,myHdrCheckByteSize);
      text1[myHdrCheckByteSize] = '\0';
      std::string text = text1;
      hdrValue->setStringValue( text );
      delete [] text1;
    }
    //    fprintf(stderr,"peekHeaderValue: %d\n", hdrValue->intValue());
  }
  
  return success;
}

// Reset file pointer to start of current trace, from whereever it is at the moment
bool csSegyReader::revertFromPeekPosition() {
  if (myPeekIsInProgress ) {
    csInt64_t bytePosRelative( (csInt64_t)(myCurrentTraceInFile - myCurrentPeekTraceIndex ) * (csInt64_t)myTraceByteSize - (csInt64_t)(myCurrentPeekByteOffset+myCurrentPeekByteSize) ); 
    //    fprintf(stderr,"revertFrom...: %lld\n", bytePosRelative);
    if( !csFileUtils::seekg_relative(bytePosRelative,myFile) ) return false;
  
    myPeekIsInProgress      = false;
    myCurrentPeekTraceIndex = 0;
  }    
  
  return true;
}

//--------------------------------------------------------------------
bool csSegyReader::peek( int byteOffset, int byteSize, char* buffer, int traceIndex ) {
  if( !myHasBeenInitialized ) {
    throw( cseis_geolib::csException("csSegyReader::peek: reader object has not been initialized. This is a program bug in the calling function") );
  }
  if( myEnableRandomAccess == false ) {
    throw( cseis_geolib::csException("csSegyReader::peek: Random access not enabled. Set enableRandomAccess to true. This is a program bug in the calling function." ) );
  }
  if( myFileSize == cseis_geolib::csFileUtils::FILESIZE_UNKNOWN ) {
    throw( cseis_geolib::csException("csSegyReader::peek: File size unknown. This may be due to a compatibility problem of this compiled version of the program on the current platform." ) );
  }

  //  myCurrentPeekByteOffset = byteOffset;
  //  myCurrentPeekByteSize   = byteSize;

  if( traceIndex < 0 && myCurrentTraceInFile >= myNumTraces && myBufferCurrentTrace == myBufferNumTraces ) {
    //    fprintf(stderr,"Seismic_ver: PEEK LAST TRACE  %d %d %d\n", myNumTraces, myBufferCurrentTrace, myBufferNumTraces );
    return false;  // Last trace has been reached. Cannot peek ahead.
  }

  bool success = true;
  // No traces are buffered: Read ahead in input file:
  if( traceIndex >= 0 ) {
    csInt64_t bytePosRelative = 0;
  
    if( myPeekIsInProgress ) {      // subsequent peeks, jump from header block offset     
      bytePosRelative = (csInt64_t)(traceIndex - myCurrentPeekTraceIndex) * (csInt64_t)myTraceByteSize - myCurrentPeekByteSize + (byteOffset-myCurrentPeekByteOffset);
      // OLD:                bytePosRelative = (csInt64_t)(traceIndex - myCurrentPeekTraceIndex) * (csInt64_t)myTraceByteSize - byteSize;
    }
    else {                          // 1st peek, sets file pointer to offset in header block    
      bytePosRelative = (csInt64_t)(traceIndex - myCurrentTraceInFile) * (csInt64_t)myTraceByteSize + (csInt64_t)byteOffset;
      myPeekIsInProgress = true;
    }  

    if( !csFileUtils::seekg_relative(bytePosRelative,myFile) ) return false;  
   
    myCurrentPeekTraceIndex = traceIndex;    
    //    fprintf(stderr,"peek(1): bytePosRelative: %lld, Current peek pos: %d\n", bytePosRelative, traceIndex );

    myFile->read( buffer, byteSize );
    success = !myFile->fail();
  }
  else if( myBufferCapacityNumTraces == 1 || myBufferCurrentTrace == myBufferNumTraces ) {
    if( myFile->eof() ) return false;
    
    myFile->clear();
    myFile->seekg( byteOffset, std::ios_base::cur );
    myFile->read( buffer, byteSize );
    success = !myFile->fail();

    if( success ) {
      // Go back to where we were...
      myFile->seekg( -byteOffset - byteSize, std::ios_base::cur );
      if( myFile->fail() ) {
        throw( cseis_geolib::csException("csSegyReader::peekNextHeaderValue: Unknown problem occurred when trying to reset file pointer.") );
      }
    }
  }
  // Pick up bytes from stored data array
  else {
    memcpy( buffer, &myBigBuffer[myBufferCurrentTrace*myTraceByteSize + byteOffset], byteSize );
  }

  myCurrentPeekByteOffset = byteOffset;
  myCurrentPeekByteSize   = byteSize;

  return success;
}

//*******************************************************************
//
// Setting and changing SEGY headers before reading them in
//
//*******************************************************************

//-----------------------------------------------------------------------------------------
int csSegyReader::numTraceHeaders() const {
  return myTrcHdrMap->numHeaders();
}
csSegyHeaderInfo const* csSegyReader::header( int hdrIndex ) const {
  return myTrcHdrMap->header( hdrIndex );
}
bool csSegyReader::addHeader( int byteLoc, int byteSize, type_t inType, cseis_geolib::csHeaderInfo const& hdr ) {
  return myTrcHdrMap->addHeader( byteLoc, byteSize, inType, hdr );
}
bool csSegyReader::addHeader( cseis_geolib::csSegyHeaderInfo const& hdr ) {
  return myTrcHdrMap->addHeader( hdr );
}
bool csSegyReader::addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theType ){
  return myTrcHdrMap->addHeader( theName, theByteLoc, theByteSize, theType, theType, "" );
}
bool csSegyReader::addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theType, std::string const& theDesc ) {
  return myTrcHdrMap->addHeader( theName, theByteLoc, theByteSize, theType, theDesc );
}

bool csSegyReader::setSelection( std::string const& hdrValueSelectionText,
                                 std::string const& headerName,
                                 int sortOrder,
                                 int sortMethod )
{
  myIOSelection = new cseis_geolib::csIOSelection( headerName, sortOrder, sortMethod );
  return myIOSelection->initialize( this, hdrValueSelectionText );
}
void csSegyReader::setSelectionStep1( std::string const& hdrValueSelectionText,
                                      std::string const& headerName,
                                      int sortOrder,
                                      int sortMethod )
{
  myIOSelection = new cseis_geolib::csIOSelection( headerName, sortOrder, sortMethod );
  myIOSelection->step1( this, hdrValueSelectionText );
}
bool csSegyReader::setSelectionStep2( int& numTracesToRead )
{
  return myIOSelection->step2( this, numTracesToRead );
}
bool csSegyReader::setSelectionStep3()
{
  return myIOSelection->step3( this );
}


int csSegyReader::getNumSelectedTraces() const {
  return myIOSelection->getNumSelectedTraces();
}
cseis_geolib::csFlexNumber const* csSegyReader::getSelectedValue( int traceIndex ) const {
  return myIOSelection->getSelectedValue( traceIndex );
}
int csSegyReader::getSelectedIndex( int traceIndex ) const {
  return myIOSelection->getSelectedIndex( traceIndex );
}
