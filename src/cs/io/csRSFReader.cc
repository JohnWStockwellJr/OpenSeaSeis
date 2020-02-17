/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csRSFReader.h"
#include "csRSFHeader.h"
#include "csException.h"
#include "csVector.h"
#include "csFileUtils.h"
#include "csFlexHeader.h"
#include "geolib_endian.h"
#include "geolib_string_utils.h"
#include "csStandardHeaders.h"
#include "csIOSelection.h"
#include <string>
#include <cstring>
#include <cmath>
#include <limits>

using namespace cseis_io;
using namespace cseis_geolib;


csRSFReader::csRSFReader( std::string filenameRSF, int numTracesBuffer, bool reverseByteOrder ) {
  myDataBuffer         = NULL;
  myHasBeenInitialized = false;
  myIsAtEOF            = false;
  myHdr                = NULL;
  myPeekHdr            = csRSFHeader::HDR_NONE;
  myCopyBuffer         = NULL;
  myIOSelection        = NULL;

  myNumSamples = 0;
  mySampleInt  = 0;
  if( numTracesBuffer <= 0 ) {
    myBufferCapacityNumTraces = 20;
  }
  else {
    myBufferCapacityNumTraces = numTracesBuffer;
  }

  myFilenameRSF  = filenameRSF;
  myIsAtEOF      = false;
  myDoSwapEndian = reverseByteOrder;
  myNumDimensions = 3;

  mySampleByteSize = 0;
  myTraceByteSize = 0;
  myFileSize = cseis_geolib::csFileUtils::FILESIZE_UNKNOWN;

  myBufferNumTraces = 0;
  // Index pointer to current trace in buffer
  myBufferCurrentTrace = 0;
  // Total trace counter of all traces accessed via the getNextTrace() method
  myTotalTraceCounter = 0;
  // Total number of traces in input file
  myTotalNumTraces    = 0;
  // File pointer to current trace in input file
  myCurrentFileTraceIndex = 0;
  myCurrentTraceIndex = 0;

  myFileBin = NULL;
}
//-----------------------------------------------------------------------------------------
csRSFReader::~csRSFReader() {
  closeFile();
  if( myIOSelection != NULL ) {
    delete myIOSelection;
    myIOSelection = NULL;
  }
  if( myCopyBuffer != NULL ) {
    delete [] myCopyBuffer;
    myCopyBuffer = NULL;
  }
  if( myDataBuffer ) {
    delete [] myDataBuffer;
    myDataBuffer = NULL;
  }
  if( myFileBin ) {
    myFileBin->close();
    delete myFileBin;
    myFileBin = NULL;
  }
  if( myHdr ) {
    delete myHdr;
    myHdr = NULL;
  }
}

void csRSFReader::convert2standardUnit() {
  // Convert unit of sample interval to [ms] or [m] if necessary
  if( myHdr->unit1 == csRSFHeader::SAMPLE_UNIT_S ) {
    myHdr->unit1 = csRSFHeader::SAMPLE_UNIT_MS;
    myHdr->d1   *= 1000.0f;
  }
  else if( myHdr->unit1 == csRSFHeader::SAMPLE_UNIT_KM ) {
    myHdr->unit1 = csRSFHeader::SAMPLE_UNIT_M;
    myHdr->d1   *= 1000.0f;
  }

  if( myHdr->unit2 == csRSFHeader::SAMPLE_UNIT_S ) {
    myHdr->unit2 = csRSFHeader::SAMPLE_UNIT_MS;
    myHdr->d1   *= 1000.0f;
  }
  else if( myHdr->unit2 == csRSFHeader::SAMPLE_UNIT_KM ) {
    myHdr->unit2 = csRSFHeader::SAMPLE_UNIT_M;
    myHdr->d1   *= 1000.0f;
  }

  if( myHdr->unit3 == csRSFHeader::SAMPLE_UNIT_S ) {
    myHdr->unit3 = csRSFHeader::SAMPLE_UNIT_MS;
    myHdr->d1   *= 1000.0f;
  }
  else if( myHdr->unit3 == csRSFHeader::SAMPLE_UNIT_KM ) {
    myHdr->unit3 = csRSFHeader::SAMPLE_UNIT_M;
    myHdr->d1   *= 1000.0f;
  }
}

//-----------------------------------------------------------------------------------------
void csRSFReader::initialize( csRSFHeader* hdr ) {
  readRSFHdr();
  convert2standardUnit();
  hdr->set( *myHdr );

  myNumSamples = myHdr->n1;
  mySampleInt  = myHdr->d1;
  if( myHdr->n2 == 0 ) throw( csException("No traces in second dimension. n2=%d", myHdr->n2) );
  if( myHdr->n5 != 0 ) {
    if( myHdr->n4 == 0 ) throw( csException("No traces in forth dimension (..but fifth dimension exists). n3=%d, n4=%d, n5=%d", myHdr->n3, myHdr->n4, myHdr->n5) );
    myTotalNumTraces = myHdr->n2 * myHdr->n3 * myHdr->n4 * myHdr->n5;
    myNumDimensions  = 5;
  }
  else if( myHdr->n4 != 0 ) {
    if( myHdr->n3 == 0 ) throw( csException("No traces in third dimension (..but forth dimension exists). n3=%d, n4=%d", myHdr->n3, myHdr->n4) );
    myTotalNumTraces = myHdr->n2 * myHdr->n3 * myHdr->n4;
    myNumDimensions  = 4;
  }
  else if( myHdr->n3 != 0 ) {
    myTotalNumTraces = myHdr->n2 * myHdr->n3;
    myNumDimensions  = 3;
  }
  else {
    myTotalNumTraces = myHdr->n2; 
    myNumDimensions  = 2;
  }

  if( myHdr->esize == 4 ) {
    mySampleByteSize = 4;   // assume 4 byte floating point
  }
  else if( myHdr->esize == 8 ) {
    mySampleByteSize = 8;   // assume 8 byte floating point
  }
  else {
    throw( csException("Data format not supported. Byte size (esize): %d", myHdr->esize) );
  }
  myTraceByteSize = myNumSamples*mySampleByteSize;

  myBufferCapacityNumTraces = std::min( myTotalNumTraces, myBufferCapacityNumTraces );
  myDataBuffer = new char[ myBufferCapacityNumTraces * myTraceByteSize ];
  if( !myDataBuffer ) {
    throw( csException("Not enough memory. Too many buffered traces requested: %d.", myBufferCapacityNumTraces) );
  }
  memset( myDataBuffer, 0, myBufferCapacityNumTraces * myTraceByteSize );

  myTotalTraceCounter    = 0;
  myBufferCurrentTrace   = 0;

  openBinFile();

  myHasBeenInitialized = true;
}

//-----------------------------------------------------------------------------------------
void csRSFReader::openBinFile() {
  myFileBin = new std::ifstream();
  myFileBin->open( myHdr->filename_bin_full_path.c_str(), std::ios::in | std::ios::binary );
  if( myFileBin->fail() ) {
    throw csException("Cannot open RSF binary file for reading: %s", myHdr->filename_bin_full_path.c_str());
  }
  try {
    myFileSize = cseis_geolib::csFileUtils::retrieveFileSize( myHdr->filename_bin_full_path );
  }
  catch( cseis_geolib::csException& e ) {
    throw cseis_geolib::csException("Error occurred while determining file size. System message: %s\n", e.getMessage() );
  }
}
void csRSFReader::closeFile() {
  if( myFileBin != NULL ) {
    myFileBin->close();
    delete myFileBin;
    myFileBin = NULL;
  }
}

//*******************************************************************
//
// Read RSF header file
//
//*******************************************************************
void csRSFReader::readRSFHdr() {
  if( myHdr != NULL ) delete myHdr;
  myHdr = new csRSFHeader();

  FILE* fin_rsf = fopen(myFilenameRSF.c_str(),"r");
  if( fin_rsf == NULL ) {
    throw csException("Cannot open rsf header file for reading: %s", myFilenameRSF.c_str());
  }

  char buffer[256];
  char buffer2[256];
  while( fgets(buffer,256,fin_rsf) != NULL ) {
    int length = strlen(buffer);
    if( strlen(buffer) < 2 ) continue;
    for( int i = 0; i < length; i++ ) {
      if( buffer[i] == '\t' ) buffer[i] = ' ';
    }
    std::string bufferStr = trim(buffer);
    if( bufferStr.length() == 0 ) continue;
    // 1) Split into several pieces if more than one variable is defined in one line
    //    Example:  "o1=0    d1=6.25    n1=1597 label1=DEPTH"  -->  "o1=0"  "d1=6.25"  "n1=1597"  ...
    char* cPtrOuter = strtok(buffer," ");
    csVector<std::string> stringList;
    while( cPtrOuter != NULL ) {
      //      fprintf(stderr,"Text: '%s'\n", cPtrOuter);
      std::string newString( cPtrOuter );
      stringList.insertEnd( newString );
      cPtrOuter = strtok(NULL," ");
    }
    // 2) Split each item into name and value
    //    Example:  "n1=1597"  -->  name=n1  value=1597
    for( int i = 0; i < stringList.size(); i++ ) {
      char const* stringPtr = stringList.at(i).c_str();
      memcpy( buffer2, stringPtr, stringList.at(i).size() );
      buffer2[stringList.at(i).size()] = '\0';
      char* cPtr = strtok(buffer2,"=");
      if( cPtr != NULL ) {
        string name = trim(cPtr);
        cPtr = strtok(NULL,"=");
        if( name.length() > 0 && cPtr != NULL ) {
          string value = trim(cPtr);
          //          fprintf(stderr,"Name/value: '%s = %s'\n", name.c_str(), value.c_str() ); 
          if( value.length() > 0 ) {
            myHdr->setField( name.c_str(), value.c_str() );
          }
        }
      }
    } // END while cPtrOuter
  }
  if( fabs( myHdr->e1 - (myHdr->d1 * myHdr->n1 + myHdr->o1 ) ) > 0.01 ) {
    myHdr->e1 = (myHdr->d1 * myHdr->n1 + myHdr->o1 );
  }
  if( fabs( myHdr->e2 - (myHdr->d2 * myHdr->n2 + myHdr->o2 ) ) > 0.01 ) {
    myHdr->e2 = (myHdr->d2 * myHdr->n2 + myHdr->o2 );
  }
  if( fabs( myHdr->e3 - (myHdr->d3 * myHdr->n3 + myHdr->o3 ) ) > 0.01 ) {
    myHdr->e3 = (myHdr->d3 * myHdr->n3 + myHdr->o3 );
  }
  if( fabs( myHdr->e4 - (myHdr->d4 * myHdr->n4 + myHdr->o4 ) ) > 0.01 ) {
    myHdr->e4 = (myHdr->d4 * myHdr->n4 + myHdr->o4 );
  }
  if( myHdr->filename_bin.length() == 0 ) {
    throw( csException("Cannot extract filename for binary input file (in=... in rsf header file)") );
  }
  if( myHdr->filename_bin.at(0) == '/' ) {
    myHdr->filename_bin_full_path = myHdr->filename_bin;
  }
  else {
    // If binary file name in rsf file does not contain full path,
    // set binary file name with full path
    int counter = myFilenameRSF.length()-1;
    while( counter > 0 ) {
      if( myFilenameRSF.at(counter) == '/' ) {
        break;
      }
      counter -= 1;
    }
    if( counter == 0 ) {
      counter = -1;
    }
    int sizeFile = myHdr->filename_bin.length();
    int sizePath = counter + 1;
    char* tmpStr = new char[sizePath + sizeFile + 1];
    memcpy(tmpStr,myFilenameRSF.c_str(),sizePath);
    memcpy(&tmpStr[sizePath],myHdr->filename_bin.c_str(),sizeFile);
    tmpStr[sizePath+sizeFile] = '\0';
    myHdr->filename_bin_full_path = tmpStr;
    delete [] tmpStr;
  }

  fclose(fin_rsf);
}

//*******************************************************************
//
// Reading one trace...
//
//*******************************************************************
//
bool csRSFReader::getNextTrace( byte_t* theBuffer, int numSamplesToRead ) {
  if( !myHasBeenInitialized ) {
    throw( csException("Accessing method to read first trace before initializing RSF Reader. This is a program bug in the calling method") );
  }
  //  fprintf(stdout,"getNextTrace. bufferNumTrace/total: %d/%d, traceIndex: %d\n",myBufferCurrentTrace, myBufferNumTraces, myCurrentFileTraceIndex);
  if( myIOSelection ) {
    int traceIndex = myIOSelection->getNextTraceIndex();
    if( traceIndex < 0 ) return false;
    bool success = moveToTrace( traceIndex );
    if( !success ) return false;
  }

  if( myBufferCurrentTrace == myBufferNumTraces ) {
    if( myFileBin->eof() ) return false;
    bool success = readDataBuffer();
    if( !success ) return false;
  }

  int minNumSamples = std::min( myNumSamples, numSamplesToRead );
  if( mySampleByteSize == 4 ) {
    memcpy( theBuffer, &myDataBuffer[myBufferCurrentTrace*myTraceByteSize], minNumSamples*mySampleByteSize );
  }
  else if( mySampleByteSize == 8 ) {
    double val_dbl;
    for( int isamp = 0; isamp < minNumSamples; isamp++ ) {
      memcpy( &val_dbl, &myDataBuffer[myBufferCurrentTrace*myTraceByteSize + isamp*mySampleByteSize], mySampleByteSize );
      theBuffer[isamp] =  (float)val_dbl;
    }
  }

  if( numSamplesToRead > myNumSamples ) {
    // Zero out rest of output buffer if necessary:
    memset( &theBuffer[minNumSamples*mySampleByteSize], 0, (numSamplesToRead-minNumSamples)*mySampleByteSize );
  }

  myCurrentTraceIndex = myCurrentFileTraceIndex - (myBufferNumTraces-myBufferCurrentTrace);
  myBufferCurrentTrace += 1;

  return true;
}

//--------------------------------------------------------
bool csRSFReader::readDataBuffer() {
  if( myCurrentFileTraceIndex == myTotalNumTraces ) return false;

  // Set myBufferNumTraces: Number of traces to be read into buffer
  // Make sure only as many traces are read in as necessary:
  myBufferNumTraces = std::min( myBufferCapacityNumTraces, myTotalNumTraces - myCurrentFileTraceIndex );
  //  fprintf(stdout,"readDataBuffer, capacity: %d, totNumTraces: %d, bufferNumTraces: %d, total num bytes: %d\n", myBufferCapacityNumTraces, myTotalNumTraces, myBufferNumTraces, myTraceByteSize*myBufferNumTraces );

  myFileBin->clear(); // Clear all flags
  myFileBin->read( myDataBuffer, myTraceByteSize*myBufferNumTraces );
  if( myFileBin->fail() ) {
    closeFile();
    throw( cseis_geolib::csException("csRSFReader::readDataBuffer: Unexpected error (1) occurred when reading in data from input file '%s'", myFilenameRSF.c_str()) );
  }
  else if( myFileBin->eof() ) {
    return false;
  }
  if( myDoSwapEndian ) {
    for( int itrc = 0; itrc < myBufferNumTraces; itrc++ ) {
      swapEndian4( myDataBuffer+myTraceByteSize*itrc, myTraceByteSize );
    }
  }

  myBufferCurrentTrace = 0;
  myCurrentFileTraceIndex += myBufferNumTraces;
  
  return true;
}
void csRSFReader::dump( FILE* stream ) {
  myHdr->dump( stream );
}
int csRSFReader::numDim() const {
  return myNumDimensions;
}
int csRSFReader::computeTrace( double val_dim2, double val_dim3 ) const {
  int traceIndex = (int)round( (val_dim3-myHdr->o3) / myHdr->d3 ) * myHdr->n2 + (int)round( (val_dim2-myHdr->o2) / myHdr->d2 );
  return traceIndex;
}
int csRSFReader::computeTrace( double val_dim2, double val_dim3, double val_dim4 ) const {
  int traceIndex = ( (int)round( (val_dim4-myHdr->o4) / myHdr->d4 ) * myHdr->n3 + (int)round( (val_dim3-myHdr->o3) / myHdr->d3 )  ) * myHdr->n2 + (int)round( (val_dim2-myHdr->o2) / myHdr->d2 );
  return traceIndex;
}
int csRSFReader::computeTrace( double val_dim2, double val_dim3, double val_dim4, double val_dim5 ) const {
  int traceIndex = ( ( ( (int)round( (val_dim5-myHdr->o5) / myHdr->d5 ) * myHdr->n4 + (int)round( (val_dim4-myHdr->o4) / myHdr->d4 ) ) * myHdr->n3 + (int)round( (val_dim3-myHdr->o3) / myHdr->d3 ) ) ) * myHdr->n2 + (int)round( (val_dim2-myHdr->o2) / myHdr->d2 );
  return traceIndex;
}
double csRSFReader::computeDim2( int traceIndex ) const {
  if( myHdr == NULL ) throw( cseis_geolib::csException("csRSFReader::computeDim2:") );
  double value = 0;
  if( myNumDimensions == 2 ) {
    value = myHdr->o2 + myHdr->d2 * (double)traceIndex;
  }
  else if( myNumDimensions == 3 ) {
    int steps = (int)( traceIndex / myHdr->n2 );
    int remainder = traceIndex - ( steps * myHdr->n2 );
    value = myHdr->o2 + myHdr->d2 * (double)remainder;
  }
  else if( myNumDimensions == 4 ) {
    int steps4 = (int)( traceIndex / (myHdr->n2 * myHdr->n3) );
    int traceIndex_red4 = traceIndex - (steps4 * myHdr->n2 * myHdr->n3);
    int steps3 = (int)( traceIndex_red4 / myHdr->n2 );
    int remainder2 = traceIndex_red4 - ( steps3 * myHdr->n2 );
    value = myHdr->o2 + myHdr->d2 * (double)remainder2;
  }
  else { // if( myNumDimensions == 5 ) {
    int steps5 = (int)( traceIndex / (myHdr->n2 * myHdr->n3 * myHdr->n4) );
    int traceIndex_red5 = traceIndex - (steps5 * myHdr->n2 * myHdr->n3 * myHdr->n4);
    int steps4 = (int)( traceIndex_red5 / myHdr->n2 );
    int remainder2 = traceIndex_red5 - ( steps4 * myHdr->n2 );
    value = myHdr->o2 + myHdr->d2 * (double)remainder2;
  }
  return value;
}
double csRSFReader::computeDim3( int traceIndex ) const {
  if( myHdr == NULL ) throw( cseis_geolib::csException("csRSFReader::computeDim3:") );
  double value = 0;
  if( myNumDimensions == 3 ) {
    int steps = (int)( traceIndex / myHdr->n2 );
    value = myHdr->o3 + myHdr->d3 * (double)steps;
  }
  else if( myNumDimensions == 4 ) {
    int steps4 = (int)( traceIndex / (myHdr->n2 * myHdr->n3) );
    int traceIndex_red4 = traceIndex - (steps4 * myHdr->n2 * myHdr->n3);
    int steps3 = (int)( traceIndex_red4 / myHdr->n2 );
    value = myHdr->o3 + myHdr->d3 * (double)steps3;
  }
  else { // if( myNumDimensions == 5 ) {
    int steps5 = (int)( traceIndex / (myHdr->n2 * myHdr->n3 * myHdr->n4) );
    int traceIndex_red5 = traceIndex - (steps5 * myHdr->n2 * myHdr->n3 * myHdr->n4);
    int steps4 = (int)( traceIndex_red5 / myHdr->n2 );
    value = myHdr->o4 + myHdr->d4 * (double)steps4;
  }
  return value;
}
double csRSFReader::computeDim4( int traceIndex ) const {
  if( myHdr == NULL ) throw( cseis_geolib::csException("csRSFReader::computeDim4:") );
  int steps4 = (int)( traceIndex / (myHdr->n2 * myHdr->n3) );
  double value = myHdr->o4 + myHdr->d4 * (double)steps4;
  return value;
}
double csRSFReader::computeDim5( int traceIndex ) const {
  if( myHdr == NULL ) throw( cseis_geolib::csException("csRSFReader::computeDim5:") );
  int steps5 = (int)( traceIndex / (myHdr->n2 * myHdr->n3 * myHdr->n4) );
  double value = myHdr->o5 + myHdr->d5 * (double)steps5;
  return value;
}
bool csRSFReader::moveToComputedTrace( double val_dim2, double val_dim3 ) {
  return moveToTrace( computeTrace( val_dim2, val_dim3 ) );
}
bool csRSFReader::moveToComputedTrace( double val_dim2, double val_dim3, double val_dim4 ) {
  return moveToTrace( computeTrace( val_dim2, val_dim3, val_dim4 ) );
}
bool csRSFReader::moveToComputedTrace( double val_dim2, double val_dim3, double val_dim4, double val_dim5 ) {
  return moveToTrace( computeTrace( val_dim2, val_dim3, val_dim4, val_dim5 ) );
}
bool csRSFReader::moveToTrace( int traceIndex ) {
  return moveToTrace( traceIndex, myTotalNumTraces-traceIndex );
}
bool csRSFReader::moveToTrace( int traceIndex, int numTracesToRead ) {
  if( !myHasBeenInitialized ) {
    throw( cseis_geolib::csException("csRSFReader::moveToTrace: Input file has not been initialized yet. This is a program bug in the calling function") );
  }
  else if( myFileSize == cseis_geolib::csFileUtils::FILESIZE_UNKNOWN ) {
    throw( cseis_geolib::csException("csRSFReader::moveToTrace: File size unknown. This may be due to a compatibility problem of this compiled version of the program on the current platform." ) );
  }
  else if( traceIndex < 0 || traceIndex >= myTotalNumTraces ) {
    throw( cseis_geolib::csException("csRSFReader::moveToTrace: Incorrect trace index: %d (number of traces in input file: %d). This is a program bug in the calling method", traceIndex, myTotalNumTraces) );
  }

  if( myCurrentFileTraceIndex != traceIndex ) {
    csInt64_t bytePosRelative = (csInt64_t)(traceIndex-myCurrentFileTraceIndex) * (csInt64_t)myTraceByteSize;
    if( !seekg_relative( bytePosRelative ) ) return false;
    myBufferNumTraces    = 0;
    myBufferCurrentTrace = 0;
  }
  myCurrentFileTraceIndex  = traceIndex;
  myCurrentTraceIndex  = traceIndex;

  // Index of last trace that will be read in one go by consecutive calls to getNextTrace()
  // myLastTraceIndex = std::min( traceIndex + numTracesToRead - 1, myTotalNumTraces-1 );

  return true;
}

bool csRSFReader::setHeaderToPeek( std::string const& headerName ) {
  myPeekHdr = headerIndex( headerName );
  return( myPeekHdr != csRSFHeader::HDR_NONE ); 
}
bool csRSFReader::setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType ) {
  bool success = setHeaderToPeek( headerName );
  if( myPeekHdr == csRSFHeader::HDR_DIM2 || myPeekHdr == csRSFHeader::HDR_DIM3 || myPeekHdr == csRSFHeader::HDR_DIM4 || myPeekHdr == csRSFHeader::HDR_DIM5 ) {
    headerType = cseis_geolib::TYPE_DOUBLE;
  }
  else {
    headerType = cseis_geolib::TYPE_INT;
  }
  return success;
}

//--------------------------------------------------------------------

bool csRSFReader::peekHeaderValue( cseis_geolib::csFlexHeader* hdrValue, int traceIndex ) {
  if( !myHasBeenInitialized ) {
    throw( cseis_geolib::csException("csRSFReader::peek: Reader object has not been initialized. This is a program bug in the calling function") );
  }
  else if( myFileSize == cseis_geolib::csFileUtils::FILESIZE_UNKNOWN ) {
    throw( cseis_geolib::csException("csRSFReader::peek: File size unknown. This may be due to a compatibility problem of this compiled version of the program on the current platform." ) );
  }
  else if( myPeekHdr == csRSFHeader::HDR_NONE ) {
    throw(cseis_geolib::csException("csRSFReader::peekHeaderValue: No header has been set for checking. This is a program bug in the calling function"));
  }

  if( traceIndex < 0 || traceIndex >= myTotalNumTraces ) {
    return false;  // Trace beyond end of file. Cannot peek.
  }

  if( myPeekHdr == csRSFHeader::HDR_DIM2 ) {
    hdrValue->setDoubleValue( computeDim2(traceIndex) );
  }
  else if( myPeekHdr == csRSFHeader::HDR_DIM3 ) {
    hdrValue->setDoubleValue( computeDim3(traceIndex) );
  }
  else if( myPeekHdr == csRSFHeader::HDR_DIM4 ) {
    hdrValue->setDoubleValue( computeDim4(traceIndex) );
  }
  else if( myPeekHdr == csRSFHeader::HDR_DIM5 ) {
    hdrValue->setDoubleValue( computeDim5(traceIndex) );
  }
  else if( myPeekHdr == csRSFHeader::HDR_TRCNO ) {
    hdrValue->setIntValue( traceIndex+1 );
  }
  return true;;
}

bool csRSFReader::seekg_relative( csInt64_t bytePosRelative ) {
  // Complex algorithm to be able to make step that is larger than 2Gb : Make several smaller steps instead
  
  int maxInt = std::numeric_limits<int>::max() - 1;   // -1 to be on the safe side, also for negative byte positions
  if( bytePosRelative < 0 ) maxInt *= -1;
  
  csInt64_t numSteps  = bytePosRelative / (csInt64_t)maxInt + 1LL;
  int bytePosResidual = (int)(bytePosRelative % (csInt64_t)maxInt);
  
  for( csInt64_t istep = 0; istep < numSteps - 1; ++istep ) { 
    myFileBin->clear(); // Clear all flags
    myFileBin->seekg( maxInt, std::ios_base::cur );
    if( myFileBin->fail() ) return false;
  }
  
  myFileBin->seekg( bytePosResidual, std::ios_base::cur );
  
  return true;
}
//********************************************************************************

int csRSFReader::hdrIntValue( int hdrIndex ) const {
  switch( hdrIndex ) {
  case csRSFHeader::HDR_DIM2:
    return (int)computeDim2(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM3:
    return (int)computeDim3(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM4:
    return (int)computeDim4(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM5:
    return (int)computeDim5(myCurrentTraceIndex);
  case csRSFHeader::HDR_TRCNO:
    return( myCurrentTraceIndex+1 );
  default:
    return 0;
  }
}
float csRSFReader::hdrFloatValue( int hdrIndex ) const {
  switch( hdrIndex ) {
  case csRSFHeader::HDR_DIM2:
    return (float)computeDim2(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM3:
    return (float)computeDim3(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM4:
    return (float)computeDim4(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM5:
    return (float)computeDim5(myCurrentTraceIndex);
  case csRSFHeader::HDR_TRCNO:
    return (float)( myCurrentTraceIndex+1 );
  default:
    return 0;
  }
}
double csRSFReader::hdrDoubleValue( int hdrIndex ) const {
  switch( hdrIndex ) {
  case csRSFHeader::HDR_DIM2:
    return computeDim2(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM3:
    return computeDim3(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM4:
    return computeDim4(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM5:
    return computeDim5(myCurrentTraceIndex);
  case csRSFHeader::HDR_TRCNO:
    return (double)( myCurrentTraceIndex+1 );
  default:
    return 0;
  }
}
csInt64_t csRSFReader::hdrInt64Value( int hdrIndex ) const {
  switch( hdrIndex ) {
  case csRSFHeader::HDR_DIM2:
    return (csInt64_t)computeDim2(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM3:
    return (csInt64_t)computeDim3(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM4:
    return (csInt64_t)computeDim4(myCurrentTraceIndex);
  case csRSFHeader::HDR_DIM5:
    return (csInt64_t)computeDim5(myCurrentTraceIndex);
  case csRSFHeader::HDR_TRCNO:
    return (csInt64_t)( myCurrentTraceIndex+1 );
  default:
    return 0;
  }
}
std::string csRSFReader::hdrStringValue( int hdrIndex ) const {
  return "";
}

int csRSFReader::numTraceHeaders() const {
  return( myNumDimensions ); // Number of spatial dimensions plus trace counter
}
int csRSFReader::headerIndex( std::string const& headerName ) const {
  if( !headerName.compare( "dim2" ) ) {
    return csRSFHeader::HDR_DIM2;
  }
  else if( !headerName.compare( "dim3" ) ) {
    return csRSFHeader::HDR_DIM3;
  }
  else if( !headerName.compare( "dim4" ) ) {
    return csRSFHeader::HDR_DIM4;
  }
  else if( !headerName.compare( "dim5" ) ) {
    return csRSFHeader::HDR_DIM5;
  }
  else if( !headerName.compare(cseis_geolib::HDR_TRCNO.name) ) {
    return csRSFHeader::HDR_TRCNO;
  }
  else {
    return csRSFHeader::HDR_NONE;
  }
}
std::string csRSFReader::headerName( int hdrIndex ) const {
  if( hdrIndex == csRSFHeader::HDR_DIM2 ) {
    return "dim2";
  }
  else if( hdrIndex == csRSFHeader::HDR_DIM3 ) {
    return "dim3";
  }
  else if( hdrIndex == csRSFHeader::HDR_DIM4 ) {
    return "dim4";
  }
  else if( hdrIndex == csRSFHeader::HDR_DIM5 ) {
    return "dim5";
  }
  else if( hdrIndex == csRSFHeader::HDR_TRCNO ) {
    return cseis_geolib::HDR_TRCNO.name;
  }
  else {
    return "NONE";
  }
}
std::string csRSFReader::headerDesc( int hdrIndex ) const {
  if( hdrIndex == csRSFHeader::HDR_DIM2 ) {
    return "RSF data dimension 2";
  }
  else if( hdrIndex == csRSFHeader::HDR_DIM3 ) {
    return "RSF data dimension 3";
  }
  else if( hdrIndex == csRSFHeader::HDR_DIM4 ) {
    return "RSF data dimension 4";
  }
  else if( hdrIndex == csRSFHeader::HDR_DIM5 ) {
    return "RSF data dimension 5";
  }
  else if( hdrIndex == csRSFHeader::HDR_TRCNO ) {
    return cseis_geolib::HDR_TRCNO.description;
  }
  else {
    return "NONE";
  }
}
cseis_geolib::type_t csRSFReader::headerType( int hdrIndex ) const {
  if( hdrIndex == csRSFHeader::HDR_DIM2 ) {
    return cseis_geolib::TYPE_DOUBLE;
  }
  else if( hdrIndex == csRSFHeader::HDR_DIM3 ) {
    return cseis_geolib::TYPE_DOUBLE;
  }
  else if( hdrIndex == csRSFHeader::HDR_DIM4 ) {
    return cseis_geolib::TYPE_DOUBLE;
  }
  else if( hdrIndex == csRSFHeader::HDR_DIM5 ) {
    return cseis_geolib::TYPE_DOUBLE;
  }
  else if( hdrIndex == csRSFHeader::HDR_TRCNO ) {
    return cseis_geolib::HDR_TRCNO.type;
  }
  else {
    return cseis_geolib::TYPE_UNKNOWN;
  }
}
float const* csRSFReader::getNextTracePointer() {
  if( myCopyBuffer == NULL ) {
    myCopyBuffer = new float[numSamples()];
  }
  if( getNextTrace( (byte_t*)myCopyBuffer, numSamples() ) ) return myCopyBuffer;
  else return NULL;
}


int csRSFReader::getCSEISDomain() const {
  if( myHdr->unit1 == csRSFHeader::SAMPLE_UNIT_MS || myHdr->unit1 == csRSFHeader::SAMPLE_UNIT_S ) {
    return cseis_geolib::DOMAIN_XT;
  }
  else if( myHdr->unit1 == csRSFHeader::SAMPLE_UNIT_M || myHdr->unit1 == csRSFHeader::SAMPLE_UNIT_KM || myHdr->domain1 == cseis_geolib::DOMAIN_XD ) {
    return cseis_geolib::DOMAIN_XD;
  }
  return cseis_geolib::DOMAIN_XT;
}

bool csRSFReader::setSelection( std::string const& hdrValueSelectionText,
                                std::string const& headerName,
                                int sortOrder,
                                int sortMethod )
{
  myIOSelection = new cseis_geolib::csIOSelection( headerName, sortOrder, sortMethod );
  return myIOSelection->initialize( this, hdrValueSelectionText );
}
void csRSFReader::setSelectionStep1( std::string const& hdrValueSelectionText,
                                     std::string const& headerName,
                                     int sortOrder,
                                     int sortMethod )
{
  myIOSelection = new cseis_geolib::csIOSelection( headerName, sortOrder, sortMethod );
  myIOSelection->step1( this, hdrValueSelectionText );
}
bool csRSFReader::setSelectionStep2( int& numTracesToRead )
{
  return myIOSelection->step2( this, numTracesToRead );
}
bool csRSFReader::setSelectionStep3()
{
  return myIOSelection->step3( this );
}

int csRSFReader::getNumSelectedTraces() const {
  return myIOSelection->getNumSelectedTraces();
}
cseis_geolib::csFlexNumber const* csRSFReader::getSelectedValue( int traceIndex ) const {
  return myIOSelection->getSelectedValue( traceIndex );
}
int csRSFReader::getSelectedIndex( int traceIndex ) const {
  return myIOSelection->getSelectedIndex( traceIndex );
}
