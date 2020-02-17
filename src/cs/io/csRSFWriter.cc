/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csRSFWriter.h"
#include "csRSFHeader.h"
#include "csException.h"
#include "csVector.h"
#include "csFileUtils.h"
#include "geolib_endian.h"
#include "csByteConversions.h"
#include <string>
#include <cstring>

using namespace cseis_io;
using namespace cseis_geolib;

csRSFWriter::csRSFWriter( std::string filename, int nTracesBuffer, bool reverseByteOrder, int numTraces_swapDim3, bool outputGrid, double tolerance ) :
  myNumTracesBuffer(nTracesBuffer), TOLERANCE(tolerance)
{
  myBigBuffer    = NULL;
  myNumTraces_swapDim3 = numTraces_swapDim3;
  if( myNumTraces_swapDim3 > 0 ) myNumTracesBuffer = myNumTraces_swapDim3;
  myOutputGrid = outputGrid;
  myHdr       = NULL;
  myStage = csRSFWriter::STAGE_1_INIT;

  myFilename     = filename;
  myIsAtEOF      = false;
  myDoSwapEndian = reverseByteOrder;

  // Create output files if they do not exist yet. Do not overwrite.
  //  bool fileExists = csFileUtils::createDoNotOverwrite( myHdr->filename_bin_full_path );
  // if( !fileExists ) throw( csException("Cannot open RSF binary file for writing: %s", myHdr->filename_bin_full_path.c_str() ) );
  // fileExists = csFileUtils::createDoNotOverwrite( myFilename );
  // if( !fileExists ) throw csException("Cannot open RSF header file for writing: %s", myFilename.c_str());
}
//-----------------------------------------------------------------------------------------
csRSFWriter::~csRSFWriter() {
  closeFile();
  try {
    if( myStage == csRSFWriter::STAGE_6_COMPLETE ) {
      finalize();
    }
  }
  catch(...) {
    // nothing
  }
  if( myBigBuffer ) {
    delete [] myBigBuffer;
    myBigBuffer = NULL;
  }
  if( myFile ) {
    myFile->close();
    delete myFile;
    myFile = NULL;
  }
  if( myHdr ) {
    delete myHdr;
    myHdr = NULL;
  }
}
//-----------------------------------------------------------------------------------------
void csRSFWriter::initialize( csRSFHeader const* hdr ) {
  myHdr = new csRSFHeader();
  myHdr->set( *hdr );

  myNumSamples = myHdr->n1;
  mySampleInt  = myHdr->d1;

  /*
    if( myNumTraces_SwapDim3 > 0 ) {
    myHdr->n1 = hdr->n3;
    myHdr->o1 = hdr->o3;
    myHdr->d1 = hdr->d3;
    myHdr->e1 = hdr->e3;
    myHdr->n3 = hdr->n1;
    myHdr->o3 = hdr->o1;
    myHdr->d3 = hdr->d1;
    myHdr->e3 = hdr->e1;
    }
  */

  mySampleByteSize = 4;   // assume 4 byte floating point
  myTotalTraceSize = myNumSamples*mySampleByteSize;

  long long int totalNumBytes = myNumTracesBuffer * myTotalTraceSize;
  myBigBuffer = new char[ totalNumBytes ];
  if( !myBigBuffer ) {
    throw( csException("Not enough memory...") );
  }
  memset( myBigBuffer, 0, totalNumBytes );

  myTraceCounter   = 0;
  myCurrentTrace   = 0;
  myNumSavedTraces = 0;

  myStage = csRSFWriter::STAGE_2_SET_ORIG;
  //  myHasBeenInitialized = true;

  // Create output files if they do not exist yet. Do not overwrite.
  bool fileExists = csFileUtils::createDoNotOverwrite( myHdr->filename_bin_full_path );
  if( !fileExists ) throw( csException("Cannot open RSF binary file for writing: %s", myHdr->filename_bin_full_path.c_str() ) );
  fileExists = csFileUtils::createDoNotOverwrite( myFilename );
  if( !fileExists ) throw csException("Cannot open RSF header file for writing: %s", myFilename.c_str());
}
//-----------------------------------------------------------------------------------------
void csRSFWriter::finalize() {
  //  if( myStage == csRSFWriter::STAGE_6_COMPLETE ) return;
  myHdr->e2 = myHdr->o2 + myHdr->d2 * (myHdr->n2-1);
  myHdr->e3 = myHdr->o3 + myHdr->d3 * (myHdr->n3-1);

  if( myNumTraces_swapDim3 ) {
    finalizeSwapDim3();
  }
  writeRSFHdr();
  myStage = csRSFWriter::STAGE_6_COMPLETE;
}
void csRSFWriter::finalizeSwapDim3() {
  // Swap RSF headers:
  int ntmp    = myHdr->n1;
  double otmp = myHdr->o1;
  double etmp = myHdr->e1;
  double dtmp = myHdr->d1;

  myHdr->n1 = myHdr->n2;
  myHdr->o1 = myHdr->o2;
  myHdr->d1 = myHdr->d2;
  myHdr->e1 = myHdr->e2;
  myHdr->n2 = myHdr->n3;
  myHdr->o2 = myHdr->o3;
  myHdr->d2 = myHdr->d3;
  myHdr->e2 = myHdr->e3;
  myHdr->n3 = ntmp;
  myHdr->o3 = otmp;
  myHdr->d3 = dtmp;
  myHdr->e3 = etmp;

  if( myHdr->n1 * myHdr->n2 != myNumSavedTraces ) throw( csException("Inconsistent number of traces: %d(n1) x %d(n2) != %d(total ntrc)", myHdr->n1, myHdr->n2, myNumSavedTraces) );

  // Write swapped binary file:
  int totalTraceSize = myHdr->n1 * mySampleByteSize;;
  char* bufferOneTrace = new char[ totalTraceSize ];
  for( int isamp = 0; isamp < ntmp; isamp++ ) {
    for( int itrc2 = 0; itrc2 < myHdr->n2; itrc2++ ) {
      for( int itrc1 = 0; itrc1 < myHdr->n1; itrc1++ ) {
        memcpy( &bufferOneTrace[itrc1*mySampleByteSize], &myBigBuffer[ (itrc1+(itrc2*myHdr->n1)) * myTotalTraceSize + isamp*mySampleByteSize ], mySampleByteSize);
      }
      myFile->write( bufferOneTrace, totalTraceSize );
      if( myFile->fail() ) {
        myFile->close();
        delete myFile;
        myFile = NULL;
        throw( csException("Unexpected error occurred when writing to RSF binary file") );
      }
    }
  }
  delete bufferOneTrace;
}

//-----------------------------------------------------------------------------------------
void csRSFWriter::openFile() {
  myFile = new std::ofstream();
  myFile->open( myHdr->filename_bin_full_path.c_str(), std::ios::out | std::ios::binary );
  if( myFile->fail() ) {
    throw csException("Cannot open RSF binary file for writing: %s", myHdr->filename_bin_full_path.c_str());
  }
}
void csRSFWriter::closeFile() {
  if( myFile != NULL ) {
    if( myNumSavedTraces > 0 ) {
      writeNextTrace( NULL, 0, 0, 0 );
    }
    myFile->close();
    delete myFile;
    myFile = NULL;
  }
}
//--------------------------------------------------------------------------------
//
bool csRSFWriter::check( double valDim2, double valDim3 ) {
  if( myStage <= csRSFWriter::STAGE_4_DIM3_STEP ) {
    if( myStage <= csRSFWriter::STAGE_3_DIM2_STEP ) {

      // Stage 2: Set dimension 2 and 3 origin
      if( myStage == csRSFWriter::STAGE_2_SET_ORIG ) {
        myHdr->o2 = valDim2;
        myHdr->o3 = valDim3;
        myHdr->n2 = 1;
        myHdr->n3 = 1;
        myStage = csRSFWriter::STAGE_3_DIM2_STEP;
      }
      // Stage 3: Set dimension 2 step
      else if( myStage == csRSFWriter::STAGE_3_DIM2_STEP ) {
        myHdr->d2 = valDim2 - myHdr->o2;
        myHdr->n2 = 2;
        if( fabs(valDim2 - myHdr->o2) < TOLERANCE ) {
          throw( csException("Unexpected value in sequential trace #%d for dimension2: %f. Value of trace #1 was: %f.\nExpected dimension2 values to change from trace #1 to trace #2, since dimension2 is the fast dimension.",
                             myTraceCounter, valDim2, myHdr->o2) );
        }
        myStage = csRSFWriter::STAGE_4_DIM3_STEP;
      }
      else { 
        throw( csException("csRSFWriter::check(): Wrong 'stage' number: %d. This is a program bug in this class", myStage) );
      }
    }
    // Stage 4: Set dimension 3 step
    else {
      myHdr->d3 = valDim3 - myHdr->o3;
      if( fabs(myHdr->d3) > TOLERANCE ) {
        myHdr->n3 = 2;
        myCurrentCounterDim2 = 1;
        myStage = csRSFWriter::STAGE_5_WRITE_RSF;
      }
      else {
        if( fabs( (valDim2-myHdr->o2) - (myHdr->n2*myHdr->d2) ) > TOLERANCE ) {
          throw( csException("Unexpected value in sequential trace #%d for dimension2: %f. Expected value: %f. Incorrect sorting..?",
                             myTraceCounter, valDim2, myHdr->o2 + myHdr->n2*myHdr->d2 ) );
        }
        myHdr->n2 += 1;
      }
    }
  } // END if myStage <= STAGE_4
  // Beyond stage 4: Check that all traces are sorted correctly
  else {
    myCurrentCounterDim2 += 1;
    // Within a 'row' --> Check dim2 step & dim3 value
    if( myCurrentCounterDim2 <= myHdr->n2 ) {
      if( fabs( (valDim2-myPreviousValueDim2) - myHdr->d2 ) > TOLERANCE ) {
        throw( csException("Unexpected value in sequential trace #%d for dimension2: %f. Expected value: %f. Incorrect sorting..?",
                           myTraceCounter, valDim2, myPreviousValueDim2+myHdr->d2) );
      }
      else if( fabs( valDim3-myPreviousValueDim3 ) > TOLERANCE ) {
        throw( csException("Unexpected value in sequential trace #%d for dimension3: %f. Expected value: %f. Incorrect sorting..?",
                           myTraceCounter, valDim3, myPreviousValueDim3) );
      }
    }
    // New 'row' just starting --> Check dim2 origin, dim2 amount & dim3 step
    else {
      if( fabs( valDim2-myHdr->o2 ) > TOLERANCE ) {
        throw( csException("Unexpected value in sequential trace #%d for dimension2: %f. Expected value (=origin): %f. Incorrect sorting..?",
                           myTraceCounter, valDim2, myHdr->o2) );
      }
      else if( fabs( (valDim3-myPreviousValueDim3) - myHdr->d3 ) > TOLERANCE ) {
        throw( csException("Unexpected value in sequential trace #%d for dimension3: %f. Expected value: %f. Incorrect sorting..?",
                           myTraceCounter, valDim3, myPreviousValueDim3+myHdr->d3) );
      }
      myHdr->n3 += 1;
      myCurrentCounterDim2 = 1;
    }
  }

  myPreviousValueDim2 = valDim2;
  myPreviousValueDim3 = valDim3;
  return true;
}
//*******************************************************************
//
// Write RSF header file
//
//*******************************************************************
void csRSFWriter::writeRSFHdr() {
  FILE* fhdr = fopen(myFilename.c_str(),"w");
  if( fhdr == NULL ) {
    throw csException("Cannot open rsf header file for writing: %s", myFilename.c_str());
  }
  myHdr->dump(fhdr,myOutputGrid);
  fclose(fhdr);
}

//*******************************************************************
//
// The big method, writing one trace...
//
//*******************************************************************
//
void csRSFWriter::writeNextTrace( byte_t const* theBuffer, int nSamples, double valDim2, double valDim3 ) {
  //  fprintf(stdout,"Next trace: %d %d   %f %f\n", myNumSavedTraces, myCurrentTrace, valDim2, valDim3);
  // fflush(stdout);
  if( myStage == csRSFWriter::STAGE_1_INIT ) {
    throw( csException("Accessing method to read first trace before initializing RSF Writer. This is a program bug in the calling method") );
  }
  if( nSamples == 0 || nSamples > myNumSamples ) nSamples = myNumSamples;

  if( myCurrentTrace == myNumTracesBuffer || theBuffer == NULL ) {
    if( myDoSwapEndian ) {
      for( int itrc = 0; itrc < myNumSavedTraces; itrc++ ) {
        swapEndian4( myBigBuffer+myTotalTraceSize*itrc, nSamples*mySampleByteSize );
      }
    } // END doSwapEndian

    if( myNumTraces_swapDim3 == 0 ) {
      if( myNumSavedTraces > 0 ) {
        myFile->write( myBigBuffer, myTotalTraceSize*myNumSavedTraces );
        if( myFile->fail() ) {
          myFile->close();
          delete myFile;
          myFile = NULL;
          throw( csException("Unexpected error occurred when writing to RSF binary file") );
        }
      }
      myNumSavedTraces = 0;
      myCurrentTrace   = 0;
    }
  }

  // buffer will be NULL when last traces shall be written out... Just return.
  if( theBuffer == NULL ) {
    return;
  }
  // Set input buffers

  int indexCurrentTrace = myCurrentTrace*myTotalTraceSize;
  memcpy( &myBigBuffer[indexCurrentTrace], theBuffer, nSamples*mySampleByteSize );

  myTraceCounter++;
  myCurrentTrace++;
  myNumSavedTraces++;

  check( valDim2, valDim3 );
}

