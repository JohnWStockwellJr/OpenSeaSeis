/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSeismicReader_ver.h"
#include "csSeismicReader_ver00.h"
#include "csSeismicReader_ver01.h"
#include "csSeismicReader_ver02.h"
#include "csSeismicReader_ver03.h"
#include "csSeismicReader_ver04.h"
#include "csSeismicIOConfig.h"
#include "csException.h"
#include "csHeaderInfo.h"
#include "csIODefines.h"
#include "csFileUtils.h"
#include "csFlexHeader.h"
#include <cstring>
#include <limits>

extern "C" {
  #include <sys/stat.h>
  #include <unistd.h>
  #include <stdio.h>
}

using namespace cseis_io;

//
// numtracesbuffer should be specified as 1 in case specific trace selection is going to be set
//
csSeismicReader_ver::csSeismicReader_ver( std::string filename, bool enableRandomAccess, int numTracesBuffer ) {
  myFilename = filename;
  myEnableRandomAccess = enableRandomAccess;

  myTempBuffer   = NULL;
  myByteLoc = 0;
  myIsReadFileHeader = false;
  myHeaderByteSize = 0;
  myFileSize = cseis_geolib::csFileUtils::FILESIZE_UNKNOWN;

  myBufferCapacityNumTraces = numTracesBuffer;
  myByteSizeOneSample = 4;
  myByteSizeCompression = 0;

  myDataBuffer         = NULL;
  myTraceByteSize      = 0;
  myBufferNumTraces    = 0;
  myBufferCurrentTrace = 0;
  myCurrentTraceIndex  = 0;
  myLastTraceIndex     = 0;
  myBufferFirstTrace   = 0;

  myPeekIsInProgress = false;
  myCurrentPeekTraceIndex = 0;
  myCurrentPeekByteOffset = 0;
  myCurrentPeekByteSize   = 0;
  myNumSamples = 0;

  myFile = NULL;
  open();

  //  initialize();
}
csSeismicReader_ver* csSeismicReader_ver::createReaderObject( std::string filename, bool enableRandomAccess, int numTracesBuffer ) {
  std::string versionString;
  csSeismicReader_ver::extractVersionString( filename, versionString );
  if( !versionString.substr(5,3).compare("0.4") ) {
    return new csSeismicReader_ver04( filename, enableRandomAccess, numTracesBuffer );
  }
  else if( !versionString.substr(5,3).compare("0.3") ) {
    return new csSeismicReader_ver03( filename, enableRandomAccess, numTracesBuffer );
  }
  else if( !versionString.substr(5,3).compare("0.2") ) {
    return new csSeismicReader_ver02( filename, enableRandomAccess, numTracesBuffer );
  }
  else if( !versionString.substr(5,3).compare("0.1") ) {
    return new csSeismicReader_ver01( filename, enableRandomAccess, numTracesBuffer );
  }
  else if( !versionString.substr(5,3).compare("0.0") ) {
    return new csSeismicReader_ver00( filename, enableRandomAccess, numTracesBuffer );
  }
  throw( cseis_geolib::csException("Version not supported: %s", versionString.substr(5,3).c_str() ) );
}
void csSeismicReader_ver::extractVersionString( std::string filename, std::string& versionString ) {
  std::ifstream file;
  file.open( filename.c_str(), std::ios::in | std::ios::binary );
  if( file.fail() ) {
    throw cseis_geolib::csException("Could not open SeaSeis file " + filename);
  }

  int sizeText = 8;
  char text[9];

  file.read( text, sizeText );
  if( file.fail() ) {
    //    delete [] text;
    throw cseis_geolib::csException("Unexpected error occurred when reading SeaSeis header of file %s. Is file empty?", filename.c_str());
  }
  file.close();
  text[sizeText] ='\0';

  int lenCSEIS = strlen(ID_TEXT_CSEIS);
  int lenOSEIS = strlen(ID_TEXT_OSEIS);
  if( strncmp( text, ID_TEXT_CSEIS, lenCSEIS ) && strncmp( text, ID_TEXT_OSEIS, lenOSEIS ) ) {
    //    delete [] text;
    throw( cseis_geolib::csException("Format of input file '%s' unknown. This is not a SeaSeis file. File starts with '%s'",
                                     filename.c_str(), text ) );
  }

  versionString = text;
  //  delete [] text;
}
//----------------------------------------------------------------
csSeismicReader_ver::~csSeismicReader_ver() {
  if( myFile != NULL ) {
    closeFile();
    delete myFile;
    myFile = NULL;
  }
  if( myTempBuffer != NULL ) {
    delete [] myTempBuffer;
    myTempBuffer = NULL;
  }
  if( myDataBuffer ) {
    delete [] myDataBuffer;
    myDataBuffer = NULL;
  }
}
//----------------------------------------------------------------
void csSeismicReader_ver::open() {
  myFile = new std::ifstream();
  myFile->open( myFilename.c_str(), std::ios::in | std::ios::binary );
  if( myFile->fail() ) {
    throw cseis_geolib::csException("Could not open SeaSeis file " + myFilename);
  }
  try {
    myFileSize = cseis_geolib::csFileUtils::retrieveFileSize( myFilename );
  }
  catch( cseis_geolib::csException& e ) {
    throw cseis_geolib::csException("Error occurred while determining file size. System message: %s\n", e.getMessage() );
  }
}
//----------------------------------------------------------------
void csSeismicReader_ver::closeFile() {
  if( myFile != NULL ) {
    myFile->close();
    // Do not clear ios flags..
  }
}
//----------------------------------------------------------------
bool csSeismicReader_ver::initialize() {
  int sizeText = 5 + 3;
  char* text = new char[sizeText+1];

  myFile->read( text, sizeText );
  if( myFile->fail() ) {
    delete [] text;
    throw cseis_geolib::csException("Unexpected error occurred when reading SeaSeis header");
  }
  text[sizeText] ='\0';

  int lenCSEIS = strlen(ID_TEXT_CSEIS);
  int lenOSEIS = strlen(ID_TEXT_OSEIS);
  if( strncmp( text, ID_TEXT_CSEIS, lenCSEIS ) && strncmp( text, ID_TEXT_OSEIS, lenOSEIS ) ) {
    delete [] text;
    throw( cseis_geolib::csException("Format of input file '%s' unknown. This is not an SeaSeis file.", myFilename.c_str() ) );
  }
  int lenVersion = (int)strlen(myVersionText);
  if( strncmp( &(text[lenCSEIS]), myVersionText, lenVersion ) ) {
    delete [] text;
    throw( cseis_geolib::csException("Version '%s' of input file not supported. Supported version is '%s'.", &(text[lenCSEIS]), myVersionText ) );
  }
  delete [] text;

  myHeaderByteSize = sizeText;

  //  fprintf(stderr,"IN:  Byte size: %d %d %d\n", myByteSizeSamples, myByteSizeHdrValueBlock, myHeaderByteSize );
  
  return true;
}


// Reset file pointer to start of current trace, from whereever it is at the moment
bool csSeismicReader_ver::revertFromPeekPosition() {
  if (myPeekIsInProgress ) {
    csInt64_t bytePosRelative( (csInt64_t)(myCurrentTraceIndex - myCurrentPeekTraceIndex ) * (csInt64_t)myTraceByteSize - (csInt64_t)(myCurrentPeekByteOffset+myCurrentPeekByteSize) ); 
    if( !seekg_relative(bytePosRelative) ) return false;
  
    myPeekIsInProgress      = false;
    myCurrentPeekTraceIndex = 0;
  }    
  
  return true;
}

//--------------------------------------------------------------------
bool csSeismicReader_ver::peek( int byteOffset, int byteSize, char* buffer, int traceIndex ) {
  if( !myIsReadFileHeader ) {
    throw( cseis_geolib::csException("csSeismicReader_ver::peek: File header has not been read. This is a program bug in the calling function") );
  }
  if( myEnableRandomAccess == false ) {
    throw( cseis_geolib::csException("csSeismicReader_ver::peek: Random access not enabled. Set enableRandomAccess to true. This is a program bug in the calling function." ) );
  }
  if( myFileSize == cseis_geolib::csFileUtils::FILESIZE_UNKNOWN ) {
    throw( cseis_geolib::csException("csSeismicReader_ver::peek: File size unknown. This may be due to a compatibility problem of this compiled version of the program on the current platform." ) );
  }

  //  fprintf(stdout,"peekHeader: offset %d --> %d (inProgress: %s)   %d  or  %d\n", myCurrentPeekByteOffset, byteOffset, myPeekIsInProgress ? "yes" : "no",
  //        (traceIndex - myCurrentPeekTraceIndex) * myTraceByteSize  - myCurrentPeekByteSize + (byteOffset-myCurrentPeekByteOffset), (traceIndex - myCurrentTraceIndex) * myTraceByteSize + byteOffset);

  if( traceIndex < 0 && myCurrentTraceIndex >= myNumTraces && myBufferCurrentTrace == myBufferNumTraces ) {
    //    fprintf(stderr,"Seismic_ver: PEEK LAST TRACE  %d %d %d\n", myNumTraces, myBufferCurrentTrace, myBufferNumTraces );
    //    fflush(stderr);
    return false;  // Last trace has been reached. Cannot peek ahead.
  }

  bool success = true;
  // No traces are buffered: Read ahead in input file:
  if( traceIndex >= 0 ) {
    csInt64_t bytePosRelative = 0;
  
    if( myPeekIsInProgress ) {      // subsequent peeks, jump from header block offset     
      bytePosRelative = (csInt64_t)(traceIndex - myCurrentPeekTraceIndex) * (csInt64_t)myTraceByteSize - myCurrentPeekByteSize + (byteOffset-myCurrentPeekByteOffset);
    } 
    else {                          // 1st peek, sets file pointer to offset in header block    
      bytePosRelative = (csInt64_t)(traceIndex - myCurrentTraceIndex) * (csInt64_t)myTraceByteSize + (csInt64_t)byteOffset;
      myPeekIsInProgress = true;
    }  
    //    fprintf(stdout,"--> bytePosRelative:  %lld\n", bytePosRelative);
    if( bytePosRelative > (csInt64_t)std::numeric_limits<int>::max() ) {
      if( !seekg_relative(bytePosRelative) ) return false;  
    }
    else {
      if( !myFile->seekg( (int)bytePosRelative, std::ios_base::cur ) ) return false;
    }
   
    myCurrentPeekTraceIndex = traceIndex;    
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
        throw( cseis_geolib::csException("csSeismicReader_ver::peekNextHeaderValue: Unknown problem occurred when trying to reset file pointer.") );
      }
    }
  }
  // Pick up bytes from stored data array
  else {
    memcpy( buffer, &myDataBuffer[myBufferCurrentTrace*myTraceByteSize + byteOffset], byteSize );
  }
  myCurrentPeekByteOffset = byteOffset;
  myCurrentPeekByteSize   = byteSize;
  
//  fprintf(stderr, "buffer [%s] [%d]\n", buffer, *(int*)buffer);   

  return success;
}
//--------------------------------------------------------------------
bool csSeismicReader_ver::readFileHeader( csSeismicIOConfig* config ) {
  if( myFile == NULL ) return false;
  if( myIsReadFileHeader ) throw( cseis_geolib::csException("csSeismicReader_ver::readFileHeader: Attempt to re-read file header. This is probably a program bug in the calling function") );
  myIsReadFileHeader = true;
  
  char byteSizeChar[4];
  myFile->read( byteSizeChar, 4 );
  int byteSize = 0;
  memcpy( &byteSize, byteSizeChar, 4 );
  if( myFile->fail() ) {
    throw cseis_geolib::csException("Unexpected error occurred when reading SeaSeis header");
  }
  myHeaderByteSize += 4;

  myTempBuffer = new char[byteSize];
  myFile->read( myTempBuffer, byteSize );
  if( myFile->fail() ) {
    throw cseis_geolib::csException("Unexpected error occurred when reading SeaSeis header");
  }
  myHeaderByteSize += byteSize;

  //  fprintf(stdout,"Byte size read: %d\n", byteSize);
  // Set super header
  int byteLoc = 0;
  memcpy( &config->numSamples, &myTempBuffer[byteLoc], 4 );
  memcpy( &config->sampleInt, &myTempBuffer[byteLoc+4], 4 );
  memcpy( &config->domain, &myTempBuffer[byteLoc+8], 4 );
  myNumSamples = config->numSamples;
  int numEnsKeys = 0;
  memcpy( &numEnsKeys, &myTempBuffer[byteLoc+12], 4 );
  byteLoc += 16;

  //  fprintf(stdout,"Num ens keys: %d\n", numEnsKeys);

  for( int ikey = 0; ikey < numEnsKeys; ikey++ ) {
    int sizeText;
    memcpy( &sizeText, &myTempBuffer[byteLoc], 4 );
    char* text = new char[sizeText+1];
    text[sizeText] = '\0';
    memcpy( text, &myTempBuffer[byteLoc+4], sizeText );
    config->ensKeyNames.insertEnd( text );
    delete [] text;
    byteLoc += 4 + sizeText;
  }


  memcpy( &config->byteSizeSamples, &myTempBuffer[byteLoc], 4 );
  memcpy( &config->byteSizeHdrValueBlock, &myTempBuffer[byteLoc+4], 4 );
  //  fprintf(stdout,"Num bytes: %d %d\n", config->byteSizeSamples, config->byteSizeHdrValueBlock );

  int numTrcHdrs = 0;
  memcpy( &numTrcHdrs, &myTempBuffer[byteLoc+8], 4 );
  byteLoc += 12;
  //  fprintf(stdout,"Num trc hdrs: %d\n", numTrcHdrs);

  for( int ihdr = 0; ihdr < numTrcHdrs; ihdr++ ) {
    cseis_geolib::type_t type = (cseis_geolib::type_t)myTempBuffer[byteLoc];
    //fprintf(stdout,"Header %d: Type: %d, ", ihdr, type );
    byteLoc += 1;

    int nElements;
    memcpy( &nElements, &myTempBuffer[byteLoc], 4 );
  //  fprintf(stdout,"nElements: %d, ", nElements );
    byteLoc += 4;

    int sizeName;
    memcpy( &sizeName, &myTempBuffer[byteLoc], 4 );
    //    fprintf(stdout,"SizeName: %d, ", sizeName );
    byteLoc += 4;
    char* name = new char[sizeName+1];
    name[sizeName] = '\0';
    memcpy( name, &myTempBuffer[byteLoc], sizeName );
    byteLoc += sizeName;
    //fprintf(stdout,"Name: '%s', ", name);
    
    int sizeDesc;
    memcpy( &sizeDesc, &myTempBuffer[byteLoc], 4 );
    //fprintf(stdout,"SizeDesc: %d, ", sizeDesc );
    byteLoc += 4;
    char* desc = new char[sizeDesc+1];
    desc[sizeDesc] = '\0';
    memcpy( desc, &myTempBuffer[byteLoc], sizeDesc );
    byteLoc += sizeDesc;
    //fprintf(stdout,"Desc: '%s', ", desc);
    //fprintf(stdout,"'\n");

    config->addHeader( type, name, desc, nElements );

    delete [] name;
    delete [] desc;
  }

  myByteSizeSamples       = config->byteSizeSamples;
  myByteSizeHdrValueBlock = config->byteSizeHdrValueBlock;
  myTraceByteSize         = myByteSizeSamples + myByteSizeHdrValueBlock;

  if( myBufferCapacityNumTraces <= 0 ) {  // Set number of buffered traces if it wasn't set explicitely in constructor
    myBufferCapacityNumTraces = DEFAULT_BUFFERED_SAMPLES / config->numSamples;
    if( myBufferCapacityNumTraces <= 0 ) myBufferCapacityNumTraces = 1;
    else if( myBufferCapacityNumTraces > 20 ) myBufferCapacityNumTraces = 20;
  }
  if( myBufferCapacityNumTraces > 1 ) {
    resizeDataBuffer( myBufferCapacityNumTraces );
  }

  // Determine number of traces in file
  if( myFileSize != cseis_geolib::csFileUtils::FILESIZE_UNKNOWN ) {
    csInt64_t sizeTraces  = myFileSize - (csInt64_t)myHeaderByteSize;
    myNumTraces       = (int)(sizeTraces/(csInt64_t)(myByteSizeHdrValueBlock+myByteSizeSamples));  // Total number of full traces in input file
    config->numTraces = myNumTraces;
    myLastTraceIndex  = myNumTraces-1;
  }
  else {
    config->numTraces = 0;
    myLastTraceIndex = 0;
  }

  return true;
}
//----------------------------------------------------------------
//
bool csSeismicReader_ver::readSingleTrace( float* samples, char* hdrValueBlock, int numSamples ) {
  //  fprintf(stderr,"IN readSingleTrace, compression = %d, numSamples: %d (myNumSamples: %d)\n", myByteSizeOneSample, numSamples, myNumSamples);
  myFile->read( hdrValueBlock, myByteSizeHdrValueBlock );
  if( myFile->fail() ) {
    return false;
  }

  char* bufferPtr;
  if( myByteSizeOneSample == 4 ) {
    bufferPtr = (char*)samples;
  }
  else {
    bufferPtr = &myDataBuffer[myBufferCurrentTrace*myTraceByteSize+myByteSizeHdrValueBlock];
  }

  if( numSamples >= myNumSamples ) {
    //    myFile->read( (char*)samples, myByteSizeSamples );
    //    fprintf(stderr,"readSingleTrace, bytes to read: %d %d, current trace = %d, tracebytesize = %d, bytesizehdr = %d\n",
    //        myByteSizeSamples, myByteSizeCompression,
    //        myBufferCurrentTrace, myTraceByteSize, myByteSizeHdrValueBlock);
    // fprintf(stderr,"Data buffer size %d   %x\n", myBufferCapacityNumTraces * myTraceByteSize, myDataBuffer);
    myFile->read( bufferPtr, myByteSizeSamples+myByteSizeCompression );
    if( numSamples > myNumSamples ) {
      for( int i = numSamples; i < myNumSamples; i++ ) {
        samples[i] = 0.0;  // Zero out rest of buffer
      }
    }
  }
  else {
    int byteSize2Read = (myByteSizeSamples/myNumSamples) * numSamples + myByteSizeCompression;
    //    myFile->read( (char*)samples, byteSize2Read );
    myFile->read( bufferPtr, byteSize2Read );
    myFile->seekg( myByteSizeSamples-byteSize2Read );
  }
  if( myFile->fail() ) {
    return false;
  }

  if( myByteSizeOneSample != 4 ) {
    //    fprintf(stderr,"readSingleTrace, compression = %d, numSamples: %d (myNumSamples: %d)\n", myByteSizeOneSample, numSamples, myNumSamples);
    decompressBuffer( samples, numSamples );
  }
  myCurrentTraceIndex += 1;  // Advance current file pointer to next trace.

  return true;
}
//--------------------------------------------------------
bool csSeismicReader_ver::readDataBuffer() {
  myBufferCurrentTrace = 0;

  if( myFileSize != cseis_geolib::csFileUtils::FILESIZE_UNKNOWN ) {
    if( myCurrentTraceIndex == myNumTraces ) return false;

    // Set myBufferNumTraces: Number of traces to be read into buffer
    if( myCurrentTraceIndex < myLastTraceIndex ) {
      // Make sure only as many traces are read in as necessary:
      myBufferNumTraces = std::min( myBufferCapacityNumTraces, myLastTraceIndex-myCurrentTraceIndex+1 );
    }
    else {
      myBufferNumTraces = myBufferCapacityNumTraces;
    }
    if( myCurrentTraceIndex+myBufferNumTraces > myNumTraces ) myBufferNumTraces = myNumTraces - myCurrentTraceIndex;

    myFile->clear(); // Clear all flags
    myFile->read( myDataBuffer, myTraceByteSize*myBufferNumTraces );
    if( myFile->fail() ) {
      closeFile();
      throw( cseis_geolib::csException("csSeismicReader_ver::readDataBuffer: Unexpected error (1) occurred when reading in data from input file '%s'", myFilename.c_str()) );
    }
    else if( myFile->eof() ) {
      return false;
    }
    myCurrentTraceIndex += myBufferNumTraces;
  }
  else {
    myBufferCurrentTrace = 0;
    myBufferNumTraces = myBufferCapacityNumTraces;
    myFile->read( myDataBuffer, myTraceByteSize*myBufferNumTraces );
    if( myFile->fail() ) {
      if( myFile->eof() ) {
        int numBytesRead = (int)myFile->gcount();
        myBufferNumTraces = numBytesRead / myTraceByteSize;
        // BUGFIX 080709 Added following line. Without it, reader continued to read indefinitely when total number of traces in file was N x numTracesToBuffer
        if( myBufferNumTraces== 0 ) return false;
      }
      else {
        closeFile();
        throw( cseis_geolib::csException("csSeismicReader_ver::readDataBuffer: Unexpected error (2) occurred when reading in data from input file '%s'", myFilename.c_str()) );
      }
    }
  }
  
  return true;
}
//------------------------------------------------------------------------------------
void csSeismicReader_ver::resizeDataBuffer( int newNumTraces ) {
  myBufferCapacityNumTraces = newNumTraces;
  char* newBuffer = new char[myBufferCapacityNumTraces * myTraceByteSize];
  if( !newBuffer ) throw( cseis_geolib::csException("Not enough memory...") );
  if( myDataBuffer != NULL ) {
    delete [] myDataBuffer;
    myDataBuffer = NULL;
  }
  myDataBuffer     = newBuffer;
  //  fprintf(stderr,"Data buffer resized to %d   %x\n", myBufferCapacityNumTraces * myTraceByteSize, myDataBuffer);
}
//----------------------------------------------------------------
//
bool csSeismicReader_ver::moveToNextTrace() {
  if( myCurrentTraceIndex < myNumTraces-1 ) {
    return moveToTrace( myCurrentTraceIndex+1 );
  }
  else {
    return false;
  }
}
bool csSeismicReader_ver::moveToTrace( int traceIndex ) {
  return moveToTrace( traceIndex, myNumTraces-traceIndex );
}
bool csSeismicReader_ver::moveToTrace( int traceIndex, int numTracesToRead ) {
  if( !myIsReadFileHeader ) throw( cseis_geolib::csException("csSeismicReader_ver::moveToTrace: File header has not been read. This is a program bug in the calling function") );
  if( myEnableRandomAccess == false ) {
    throw( cseis_geolib::csException("csSeismicReader_ver::moveToTrace: Random access not enabled. Set enableRandomAccess to true. This is a program bug in the calling function." ) );
  }
  else if( myFileSize == cseis_geolib::csFileUtils::FILESIZE_UNKNOWN ) {
    throw( cseis_geolib::csException("csSeismicReader_ver::moveToTrace: File size unknown. This may be due to a compatibility problem of this compiled version of the program on the current platform." ) );
  }
  else if( traceIndex < 0 || traceIndex >= myNumTraces ) {
    throw( cseis_geolib::csException("csSeismicReader_ver::moveToTrace: Incorrect trace index: %d (number of traces in input file: %d). This is a program bug in the calling method",
                       traceIndex, myNumTraces) );
  }
  if (myPeekIsInProgress ) revertFromPeekPosition();

  if( myCurrentTraceIndex != traceIndex ) {
    csInt64_t bytePosRelative = (csInt64_t)(traceIndex-myCurrentTraceIndex) * (csInt64_t)myTraceByteSize;
    if( !seekg_relative( bytePosRelative ) ) return false;
    myBufferNumTraces = 0;
    myBufferCurrentTrace = 0;
    myCurrentTraceIndex  = traceIndex;
  }
  
  // Index of last trace that will be read in one go by consecutive calls to getNextTrace()
  myLastTraceIndex = std::min( traceIndex + numTracesToRead - 1, myNumTraces-1 );

  return true;
}

bool csSeismicReader_ver::seekg_relative( csInt64_t bytePosRelative ) {

  // Complex algorithm to be able to make step that is larger than 2Gb : Make several smaller steps instead
  
  int maxInt = std::numeric_limits<int>::max() - 1;   // -1 to be on the safe side, also for negative byte positions
  if( bytePosRelative < 0 ) maxInt *= -1;
  
  csInt64_t numSteps  = bytePosRelative / (csInt64_t)maxInt + 1LL;
  int bytePosResidual = (int)(bytePosRelative % (csInt64_t)maxInt);
  
  for( csInt64_t istep = 0; istep < numSteps - 1; ++istep ) { 
    myFile->clear(); // Clear all flags
    myFile->seekg( maxInt, std::ios_base::cur );
    if( myFile->fail() ) return false;
  }
  
  myFile->seekg( bytePosResidual, std::ios_base::cur );
  
  return true;
}

//----------------------------------------------------------------
bool csSeismicReader_ver::readTrace( float* samples, char* hdrValueBlock ) {
  return readTrace( samples, hdrValueBlock, myNumSamples );
}
bool csSeismicReader_ver::readTrace( float* samples, char* hdrValueBlock, int numSamples ) {
  //  fprintf(stderr,"IN readTrace, compression = %d, numSamples: %d (myNumSamples: %d)\n", myByteSizeOneSample, numSamples, myNumSamples);
  if( !myIsReadFileHeader ) throw( cseis_geolib::csException("csSeismicReader_ver::readTrace(): File header has not been read. This is a program bug in the calling function") );
  if (myPeekIsInProgress ) revertFromPeekPosition();
  if( myBufferCapacityNumTraces == 1 ) {
    return readSingleTrace( samples, hdrValueBlock, numSamples );
  }

  if( myBufferCurrentTrace == myBufferNumTraces ) {
    if( myFile->eof() ) fprintf(stderr,"EOF reached\n");
    if( myFile->eof() ) return false;
    bool success = readDataBuffer();
    if( !success ) return false;
  }

  memcpy( hdrValueBlock, &myDataBuffer[myBufferCurrentTrace*myTraceByteSize], myByteSizeHdrValueBlock );

  // UPDATE THE FOLLOWING:
  if( myByteSizeOneSample == 4 ) {
    if( numSamples >= myNumSamples ) {
      memcpy( samples, &myDataBuffer[myBufferCurrentTrace*myTraceByteSize+myByteSizeHdrValueBlock], myByteSizeSamples );
      if( numSamples > myNumSamples ) {
        for( int i = myNumSamples; i < numSamples; i++ ) {
          samples[i] = 0.0;  // Zero out rest of buffer
        }
      }
    }
    else {
      int byteSize2Read = (myByteSizeSamples/myNumSamples) * numSamples;
      memcpy( samples, &myDataBuffer[myBufferCurrentTrace*myTraceByteSize+myByteSizeHdrValueBlock], byteSize2Read );
    }
  }
  // Decompress data samples from N to 4 bytes
  else {
    //    fprintf(stderr,"readTrace, compression = %d, numSamples: %d (myNumSamples: %d)\n", myByteSizeOneSample, numSamples, myNumSamples);
    decompressBuffer( samples, numSamples );
  }

  myBufferCurrentTrace += 1;
  return true;
}
//----------------------------------------------------------------------
void csSeismicReader_ver::decompressBuffer( float* samples, int numSamples ) {
  float minValue;
  float rangeValue;
  int byteLoc = myBufferCurrentTrace*myTraceByteSize+myByteSizeHdrValueBlock;
  memcpy( &minValue,   &myDataBuffer[byteLoc], sizeof(float) );
  byteLoc += (int)sizeof(float);
  memcpy( &rangeValue, &myDataBuffer[byteLoc], sizeof(float) );
  byteLoc += (int)sizeof(float);
  int numSamplesToRead = std::min( numSamples, myNumSamples );
  if( myByteSizeOneSample == 2 ) {
    unsigned short maxShort = std::numeric_limits<unsigned short>::max();
    float stepValue  = rangeValue / (float)maxShort;
    unsigned short* shortPtr = (unsigned short*)&myDataBuffer[byteLoc];
    //    fprintf(stderr,"decompress, compression = %d, minValue = %.6e, rangeValue = %.6e, stepValue = %.6e, byteLoc = %d, tracebytesize: %d\n",
    //        myByteSizeOneSample, minValue, rangeValue, stepValue, byteLoc, myTraceByteSize);
    // fprintf(stderr,"tracebytesize: %d, all bytes to read: %d\n", myTraceByteSize, byteLoc+numSamplesToRead*2 );
    for( int isamp = 0; isamp < numSamplesToRead; isamp++ ) {
      samples[isamp] = (float)shortPtr[isamp] * stepValue + minValue;
    }
  }
  else {
    unsigned char maxChar = std::numeric_limits<unsigned char>::max();
    float stepValue  = rangeValue / (float)maxChar;
    unsigned char* charPtr = (unsigned char*)&myDataBuffer[byteLoc];
    for( int isamp = 0; isamp < numSamplesToRead; isamp++ ) {
      samples[isamp] = (float)charPtr[isamp] * stepValue + minValue;
    }
  }
  for( int isamp = myNumSamples; isamp < numSamples; isamp++ ) {
    samples[isamp] = 0.0f;
  }
}

