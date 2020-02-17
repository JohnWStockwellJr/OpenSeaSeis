/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegyWriter.h"
#include "csSegyDefines.h"
#include "csSegyHdrMap.h"
#include "csStandardHeaders.h"
#include "csSegyBinHeader.h"
#include "csSegyTraceHeader.h"
#include "csHeaderInfo.h"
#include "csSegyHeaderInfo.h"
#include "csException.h"
#include "csVector.h"
#include "csFlexNumber.h"
#include "csByteConversions.h"
#include "geolib_endian.h"
#include "geolib_math.h"
#include "methods_number_conversions.h"
#include "csFileUtils.h"
#include <string>
#include <cstring>

using namespace cseis_geolib;

csSegyWriter::csSegyWriter( std::string filename, int nTracesBuffer, bool reverseByteOrder, bool doAutoScaleHdrs, bool isSUFormat ) :
    NTRACES_BUFFER(nTracesBuffer)
{
  myHdrValues    = NULL;
  myCharHdrBlock = NULL;
  myBinHdrBlock  = NULL;
  myBinHdr       = NULL;
  myBigBuffer    = NULL;
  myHasBeenInitialized = false;
  myNumTraceHeaders = 0;
  myIsEBCDIC        = true;
  myIsAutoScaleHeaders = doAutoScaleHdrs;
  myForceToWrite = false;
  myIsSUFormat   = isSUFormat;

  myTrcHdr            = NULL;
  myTrcHdrMap         = NULL;

  myFile         = NULL;
  myFilename     = filename;
  myIsAtEOF      = false;
  myDoSwapEndian = isPlatformLittleEndian();
  if( reverseByteOrder ) myDoSwapEndian = !myDoSwapEndian;

  myCharHdrBlock = new char[cseis_segy::SIZE_CHARHDR];
  memset( myCharHdrBlock, 0, cseis_segy::SIZE_CHARHDR );
  myBinHdrBlock  = new byte_t[cseis_segy::SIZE_BINHDR];
  memset( myBinHdrBlock, 0, cseis_segy::SIZE_BINHDR );
  myBinHdr       = new csSegyBinHeader( myDoSwapEndian );

  // Create SEGY file if it does not exist yet.
  bool fileExists = csFileUtils::createDoNotOverwrite( myFilename );
  if( !fileExists ) throw( csException("Cannot open SEGY output file for writing: %s", myFilename.c_str() ) );
}
//-----------------------------------------------------------------------------------------
csSegyWriter::~csSegyWriter() {
  freeCharBinHdr();
  closeFile();
  if( myBigBuffer ) {
    delete [] myBigBuffer;
    myBigBuffer = NULL;
  }
  if( myFile ) {
    fclose( myFile );
    myFile = NULL;
  }
  if( myTrcHdrMap ) {
    delete myTrcHdrMap;
    myTrcHdrMap = NULL;
  }
  if( myBinHdr ) {
    delete myBinHdr;
    myBinHdr = NULL;
  }
  if( myHdrValues ) {
    delete [] myHdrValues;
    myHdrValues = NULL;
  }
  if( myTrcHdr ) {
    delete myTrcHdr;
    myTrcHdr = NULL;
  }
}
//-----------------------------------------------------------------------------------------
void csSegyWriter::initialize( csSegyHdrMap const* hdrMap, char const* newCharHdr ) {
  setCharHdr( newCharHdr );

  myBinHdr->encodeHeaders( myBinHdrBlock );
  myDataSampleFormat = myBinHdr->dataSampleFormat;
  myNumSamples = myBinHdr->numSamples;
  mySampleInt = myBinHdr->sampleIntUS * 0.001;
  mySampleByteSize = 4;   // assume 4 byte floating point
  myTotalTraceSize = myNumSamples*mySampleByteSize+cseis_segy::SIZE_TRCHDR;

  myBigBuffer = new char[ NTRACES_BUFFER * myTotalTraceSize ];
  memset( myBigBuffer, 0, NTRACES_BUFFER * myTotalTraceSize );

  if( !myBigBuffer ) {
    throw( csException("Not enough memory...") );
  }

  myTraceCounter   = 0;
  myCurrentTrace   = 0;
  myNumSavedTraces = 0;

  myTrcHdrMap = new csSegyHdrMap( hdrMap );
  myTrcHdrMap->initScalars();
  myTrcHdr = new csSegyTraceHeader( myTrcHdrMap );

  myHdrValues = new cseis_geolib::csFlexNumber[numTraceHeaders()];

  myHasBeenInitialized = true;
}
int csSegyWriter::numTraceHeaders() const {
  if( myTrcHdrMap == 0 ) return 0;
  else return myTrcHdrMap->numHeaders();
}

//-----------------------------------------------------------------------------------------
void csSegyWriter::openFile() {
  myFile = fopen( myFilename.c_str(), "wb" );
  if( myFile == NULL ) {
    throw csException("Could not open SEGY file %s", myFilename.c_str());
  }
  writeCharBinHdr();
}
void csSegyWriter::closeFile() {
  if( myFile != NULL ) {
    myForceToWrite = true;
    writeNextTrace( NULL, NULL, 0 );
    fclose( myFile );
    myFile = NULL;
  }
}

//*******************************************************************
//
// Char & bin headers
//
//*******************************************************************
void csSegyWriter::setCharHdr( char const* newCharHdr ) {
  int size = std::min( (int)strlen( newCharHdr ), cseis_segy::SIZE_CHARHDR );
  char space = ' ';
  if( myIsEBCDIC ) {
    for( int i = 0; i < size; i++ ) {
      myCharHdrBlock[i] = char2ebcdic( newCharHdr[i] );
    }
    space = char2ebcdic(' ');
  }
  else {
    memcpy( myCharHdrBlock, newCharHdr, size );
  }
  if( size < cseis_segy::SIZE_CHARHDR ) {
    memset( &myCharHdrBlock[size], space, cseis_segy::SIZE_CHARHDR-size );
  }
}

//-----------------------------------------------------------------------------------------
// Private method
void csSegyWriter::writeCharBinHdr() {
  int sizeWrite = 0;
  if( !myIsSUFormat ) {
    sizeWrite = fwrite( myCharHdrBlock, cseis_segy::SIZE_CHARHDR, 1, myFile );
    if( sizeWrite != 1 ) {
      throw csException("Unexpected error encountered while writing SEGY char header");
    }
  }

  if( !myIsSUFormat ) {
    sizeWrite = fwrite( (char*)myBinHdrBlock, cseis_segy::SIZE_BINHDR, 1, myFile );
    if( sizeWrite != 1 ) {
      throw csException("Unexpected end-of-file encountered while writing SEGY binary header");
    }
  }
}
//-----------------------------------------------------------------------------------------
void csSegyWriter::freeCharBinHdr() {
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
// The big method, writing one trace...
//
//*******************************************************************
//
void csSegyWriter::writeNextTrace( byte_t const* theBuffer, int nSamples ) {
  writeNextTrace( theBuffer, myTrcHdr, nSamples );
}
void csSegyWriter::writeNextTrace( byte_t const* theBuffer, csSegyTraceHeader const* trcHdr, int nSamples ) {
  if( !myHasBeenInitialized ) {
    throw( csException("Accessing method to read first trace before initializing SEGY Writer. This is a program bug in the calling method") );
  }
//  fprintf(stdout,"Has been initialized? %d %d %d %d\n", myHasBeenInitialized, myTotalTraceSize, myCurrentTrace, myNumSavedTraces );
  if( nSamples == 0 || nSamples > myNumSamples ) nSamples = myNumSamples;

  if( myCurrentTrace == NTRACES_BUFFER || myForceToWrite ) {
    if( myDataSampleFormat == cseis_segy::DATA_FORMAT_IBM ) {
      for( int itrc = 0; itrc < myNumSavedTraces; itrc++ ) {
        ieee2ibm( (unsigned char*)(myBigBuffer+myTotalTraceSize*itrc+cseis_segy::SIZE_TRCHDR), nSamples );
      }
    }
    // Convert trace header, store into trace header output array:
    if( myDoSwapEndian ) {
      for( int itrc = 0; itrc < myNumSavedTraces; itrc++ ) {
        swapEndian4( myBigBuffer+myTotalTraceSize*itrc+cseis_segy::SIZE_TRCHDR, nSamples*mySampleByteSize );
      }
    } // END doSwapEndian

    if( myNumSavedTraces > 0 ) {
      int sizeWrite;
      sizeWrite = fwrite( myBigBuffer, myTotalTraceSize, myNumSavedTraces, myFile );
      if( sizeWrite == 0 ) {
        fclose( myFile );
        myFile = NULL;
        throw( csException("Unexpected error occurred when writing to SEGY file") );
      }
    }
    myNumSavedTraces = 0;
    myCurrentTrace   = 0;
    //    fprintf(stdout,"Read... in traces: %d (saved: %d)\n", myCurrentTrace, myNumSavedTraces );
    // fprintf(stdout,"totalTraceSize: %d, %d (num traces buffer: %d)\n", myTotalTraceSize, 4*nSamples, NTRACES_BUFFER );
  }

  // buffer will be NULL when last traces shall be written out... Just return.
  if( theBuffer == NULL ) {
    return;
  }
  // Set input buffers
  //  fprintf(stdout,"Copy buffer: %d %d  (ntraces: %d, %d)\n", myCurrentTrace*myTotalTraceSize+cseis_segy::SIZE_TRCHDR, nSamples*mySampleByteSize, NTRACES_BUFFER, myTotalTraceSize);

  int indexCurrentTrace = myCurrentTrace*myTotalTraceSize;
  memcpy( &myBigBuffer[indexCurrentTrace+cseis_segy::SIZE_TRCHDR], theBuffer, nSamples*mySampleByteSize );
  
  byte_t* trcHdrPtr = reinterpret_cast<byte_t*>( &myBigBuffer[indexCurrentTrace] );
  trcHdr->writeHeaderValues( trcHdrPtr, myDoSwapEndian, myIsAutoScaleHeaders );

  myTraceCounter++;
  myCurrentTrace++;
  myNumSavedTraces++;
}
void csSegyWriter::setIntValue( int hdrIndex, int value ) {
  myTrcHdr->setIntValue( hdrIndex, value );
}
void csSegyWriter::setFloatValue( int hdrIndex, float value ) {
  myTrcHdr->setFloatValue( hdrIndex, value );
}
void csSegyWriter::setDoubleValue( int hdrIndex, double value ) {
  myTrcHdr->setDoubleValue( hdrIndex, value );
}
void csSegyWriter::setStringValue( int hdrIndex, std::string const& value ) {
  myTrcHdr->setStringValue( hdrIndex, value );
}

