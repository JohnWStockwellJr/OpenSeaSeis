/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csStackUtil.h"
#include "csVector.h"
#include "cseis_includes.h"
#include <cmath>

using namespace cseis_system;

csStackUtil::csStackUtil( int nSamples, float normFactor, int outputOption ) {
  init( nSamples, normFactor, outputOption );
}
//csStackUtil::csStackUtil( int nSamples, float normFactor, int outputOption, bool normTimeVariant ) {
//  init( nSamples, normFactor, outputOption, normTimeVariant );
//}
void csStackUtil::init( int nSamples, float normFactor, int outputOption ) {
  myNumSamples   = nSamples;
  myNormFactor   = normFactor;
  myOutputOption = outputOption;

  myNormTimeVariant = false;
  myOutputNormTrace = false;
  myNormTraceList = NULL;
  myNormTraceIndexMap = NULL;
  myHdrId_keyValue = -1;
  myZeroThreshold = 0;
}
csStackUtil::~csStackUtil() {
  if( myNormTraceList != NULL ) {
    for( int i = 0; i < myNormTraceList->size(); i++ ) {
      delete [] myNormTraceList->at(i);
    }
    delete myNormTraceList;
    myNormTraceList = NULL;
  }
  if( myNormTraceIndexMap != NULL ) {
    delete myNormTraceIndexMap;
    myNormTraceIndexMap = NULL;
  }
}
void csStackUtil::setOutputNormTrace( bool doOutputNormTrace ) {
  myOutputNormTrace = doOutputNormTrace;
}
void csStackUtil::setTimeVariantNorm( bool doTimeVariantNorm, int hdrId_keyValue, float zeroThreshold ) {
  myNormTimeVariant = doTimeVariantNorm;
  myHdrId_keyValue  = hdrId_keyValue;
  myZeroThreshold = zeroThreshold;
  if( myNormTimeVariant ) {
    myNormTraceList = new cseis_geolib::csVector<int*>();
    myNormTraceIndexMap = new std::map<int,int>();
  }
}
void csStackUtil::clearNormTrace_internal() {
  for( int i = 0; i < myNormTraceList->size(); i++ ) {
    int* normTrace = myNormTraceList->at(i);
    if( normTrace != NULL ) delete [] normTrace;
    //for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
    //  normTrace[isamp] = 0;
    // }
  }
  myNormTraceList->clear();
  myNormTraceIndexMap->clear();
}
int* csStackUtil::getNormTrace_internal( int keyValue, bool createNew ) {
  std::map<int,int>::iterator iter = myNormTraceIndexMap->find( keyValue );
  int* normTrace = NULL;
  if( iter == myNormTraceIndexMap->end() ) {
    if( createNew ) {
      // Create new coverage trace:
      normTrace = new int[myNumSamples];
      for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
        normTrace[isamp] = 0;
      }
      myNormTraceIndexMap->insert( std::pair<int,int>( keyValue, myNormTraceList->size() ) );
      myNormTraceList->insertEnd(normTrace);
    }
  }
  else {
    normTrace = myNormTraceList->at( iter->second );
  }
  return normTrace;
}
//-----------------------------------------------------------
// Stack all input traces into first input trace (=output trace)
//
void csStackUtil::stackAndNormalizeTraces( csTraceGather* traceGather ) {
  int nTraces = traceGather->numTraces();
  if( nTraces == 0 ) return;

  csTrace* traceOut = traceGather->trace(0);
  float* samplesOut = traceOut->getTraceSamples();

  if( !myNormTimeVariant ) {
    for( int itrc = 1; itrc < nTraces; itrc++ ) {
      float const* samplesIn = traceGather->trace(itrc)->getTraceSamples();
      for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
        samplesOut[isamp] += samplesIn[isamp];
      }
    }
  }
  else {
    int keyValue = 0;
    int* normTrace = getNormTrace_internal( keyValue, true );
    // 1) Loop through samples in output trace (= first input trace)
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      if( fabs(samplesOut[isamp]) > myZeroThreshold ) {
        normTrace[isamp] += 1;
      }
      else {
        samplesOut[isamp] = 0;
      }
    }
    for( int itrc = 1; itrc < nTraces; itrc++ ) {
      float const* samplesIn = traceGather->trace(itrc)->getTraceSamples();
      for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
        if( fabs(samplesIn[isamp]) > myZeroThreshold ) {
          samplesOut[isamp] += samplesIn[isamp];
          normTrace[isamp]  += 1;
        }
      }
    }
  }
  if( myOutputOption == csStackUtil::OUTPUT_FIRST ) {
    // Nothing to be done
  }
  else if( myOutputOption == csStackUtil::OUTPUT_LAST ) {
    traceOut->getTraceHeader()->copyFrom( traceGather->trace(nTraces-1)->getTraceHeader() );
  }
  else if( myOutputOption == csStackUtil::OUTPUT_AVERAGE ) {
    for( int itrc = 1; itrc < nTraces; itrc++ ) {
      stackHeaders( traceOut->getTraceHeader(), traceGather->trace(itrc)->getTraceHeader() ) ;
    }
  }
  normStackedTrace( traceOut, nTraces );
  traceGather->freeTraces( 1, nTraces-1 );
  if( myNormTimeVariant ) clearNormTrace_internal();
}
//-----------------------------------------------------------
//
void csStackUtil::stackTraceOld( csTrace* stackedTrace, csTrace const* traceIn ) {
  float* samplesOut = stackedTrace->getTraceSamples();
  float const* samplesIn = traceIn->getTraceSamples();
  for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
    samplesOut[isamp] += samplesIn[isamp];
  }

  if( myOutputOption == OUTPUT_FIRST ) {
    // Nothing to be done
  }
  else if( myOutputOption == csStackUtil::OUTPUT_LAST ) {
    stackedTrace->getTraceHeader()->copyFrom( traceIn->getTraceHeader() );
  }
  else if( myOutputOption == csStackUtil::OUTPUT_AVERAGE ) {
    stackHeaders( stackedTrace->getTraceHeader(), traceIn->getTraceHeader() ) ;
  }
}
void csStackUtil::stackTrace( csTrace* stackedTrace ) {
  stackTrace( stackedTrace, stackedTrace, true );
}
void csStackUtil::stackTrace( csTrace* stackedTrace, csTrace const* traceIn ) {
  stackTrace( stackedTrace, traceIn, false );
}
void csStackUtil::stackTrace( csTrace* stackedTrace, csTrace const* traceIn, bool newTrace ) {
  if( myNormTimeVariant ) {
    int keyValue = 0;
    if( myHdrId_keyValue >= 0 ) keyValue = traceIn->getTraceHeader()->intValue(myHdrId_keyValue);
    int* normTrace = getNormTrace_internal( keyValue, true );
    float const* samplesIn = traceIn->getTraceSamples();
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      if( fabs(samplesIn[isamp]) > myZeroThreshold ) normTrace[isamp] += 1;
    }
  }
  // If this is a new trace, no need to stack in data samples and set headers. Otherwise, yes:
  if( !newTrace ) stackTraceOld( stackedTrace, traceIn );
}
//-----------------------------------------------------------
//
void csStackUtil::normStackedTraceOld( csTrace* trace, int nTraces ) {
  if( !myNormTimeVariant ) {
    float* samplesOut = trace->getTraceSamples();
  
    float norm = 1.0;
    if( myNormFactor == 0.5 ) {
      norm = sqrt((float)nTraces);
    }
    else if( myNormFactor == 1.0 ) {
      norm = (float)nTraces;
    }
    else if( myNormFactor == 0.0 ) {
      norm = 1.0;
    }
    else {
      norm = pow( nTraces, myNormFactor );
    }
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      samplesOut[isamp] /= norm;
    }
  }
  
  if( myOutputOption == OUTPUT_AVERAGE ) {
    normHeaders( trace->getTraceHeader(), nTraces );
  }
}
void csStackUtil::normStackedTrace( csTrace* trace, int nTraces ) {
  normStackedTraceOld( trace, nTraces );
  if( myNormTimeVariant ) {
    float* samplesOut = trace->getTraceSamples();
    int keyValue = 0;
    if( myHdrId_keyValue >= 0 ) keyValue = trace->getTraceHeader()->intValue(myHdrId_keyValue);
    int* normTrace = getNormTrace_internal( keyValue, false );
    if( normTrace == NULL ) {
      throw( cseis_geolib::csException("csStackUtil::normStackedTrace(): Key value not found: %d\n", keyValue) );
    }

    //    std::map<int,int>::iterator iter = myNormTraceIndexMap->find( keyValue );
    //  if( iter == myNormTraceIndexMap->end() ) {
    //   throw( cseis_geolib::csException("csStackUtil::normStackedTrace(): Key value not found: %d\n", keyValue) );
    //   return; // Oops!
    //  }
    //  int* normTrace = myNormTraceList->at( iter->second );
    if( !myOutputNormTrace ) {
      for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
        if( normTrace[isamp] != 0 ) {
          samplesOut[isamp] /= pow(normTrace[isamp],myNormFactor);
        }
      }
    }
    else {
      for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
        samplesOut[isamp] = normTrace[isamp];
      }
    }
  } // END norm time variant
}
void csStackUtil::stackHeaders( csTraceHeader* trcHdrOut, csTraceHeader const* trcHdrIn ) {
  int nHeaders = trcHdrIn->numHeaders();
  for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
    cseis_geolib::type_t type = trcHdrIn->type(ihdr);
    switch( type ) {
      case cseis_geolib::TYPE_INT:
        trcHdrOut->setIntValue( ihdr, trcHdrIn->intValue( ihdr ) + trcHdrOut->intValue( ihdr ) );
        break;
      case cseis_geolib::TYPE_FLOAT:
        trcHdrOut->setFloatValue( ihdr, trcHdrIn->floatValue( ihdr ) + trcHdrOut->floatValue( ihdr ) );
        break;
      case cseis_geolib::TYPE_DOUBLE:
        trcHdrOut->setDoubleValue( ihdr, trcHdrIn->doubleValue( ihdr ) + trcHdrOut->doubleValue( ihdr ) );
        break;
    }
  }
}
void csStackUtil::normHeaders( csTraceHeader* trcHdrOut, int nTraces ) {
  int nHeaders = trcHdrOut->numHeaders();
  for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
    cseis_geolib::type_t type = trcHdrOut->type(ihdr);
    switch( type ) {
      case cseis_geolib::TYPE_INT:
        trcHdrOut->setIntValue( ihdr, trcHdrOut->intValue( ihdr ) / nTraces );
        break;
      case cseis_geolib::TYPE_FLOAT:
        trcHdrOut->setFloatValue( ihdr, trcHdrOut->floatValue( ihdr ) / nTraces );
        break;
      case cseis_geolib::TYPE_DOUBLE:
        trcHdrOut->setDoubleValue( ihdr, trcHdrOut->doubleValue( ihdr ) / nTraces );
        break;
    }
  }
}
//--------------------------------------------------------------------------------
//
// BUGFIX 080702 Commented out one line in method stackTrace that led to trace headers to be summed ALWAYS, not only for OUTPUT_AVERAGE
//


