/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csDespike.h"
#include "csVector.h"
#include "geolib_math.h"

using namespace cseis_geolib;

csDespike::csDespike( int numSamples, double sampleInt, DespikeConfig& config ) {
  init( numSamples, sampleInt );
  set( config );
}
void csDespike::init( int numSamples, double sampleInt ) {
  mySampleInt     = sampleInt;
  myNumSamples    = numSamples;
  myRatios = NULL;
  myValueRefWin = NULL;
}
csDespike::~csDespike() {
  if( myRatios != NULL ) {
    delete [] myRatios;
    myRatios = NULL;
  }
  if( myValueRefWin != NULL ) {
    delete [] myValueRefWin;
    myValueRefWin = NULL;
  }
}
//***************************************************************************************
//
//
void csDespike::set( DespikeConfig const& config ) {
  if( myRatios != NULL ) {
    delete [] myRatios;
    myRatios = NULL;
  }
  if( myValueRefWin != NULL ) {
    delete [] myValueRefWin;
    myValueRefWin = NULL;
  }
  myPerformDebias = config.performDebias;
  myMaxRatio = config.maxRatio;
  myMethod   = config.method;
  myMethodRefWin   = config.methodRefWin;
  myMethodSpikeWin = config.methodSpikeWin;

  if( config.incWin <= mySampleInt ) {
    myIncWin = 1;
  }
  else {
    throw( csException("Increments other than one sample interval are currently not supported.") );
    myIncWin = (int)( config.incWin / mySampleInt );
  }

  myWidthRefWin    = 2*(int)((int)( config.widthRefWin / mySampleInt + 0.5 )/2) + 1;  // Make sure this becomes an odd number
  myWidthSpikeWin  = (int)( config.minWidthSpikeWin / mySampleInt + 0.5 );
  myStartSample    = (int)( config.start / mySampleInt + 0.5 );
  if( config.stop > 0 ) {
    myStopSample     = (int)( config.stop / mySampleInt + 0.5 );
  }
  else {
    myStopSample     = myNumSamples - 1;
  }
  myIncWin         = 1;
  myWidthMeanWin   = config.advanced.widthMeanWin_inSamples;
  myWindowAmplifier= config.advanced.windowAmplifier;
  myRatioAmplifier = config.advanced.ratioAmplifier;

  if( myStartSample < 0 ) myStartSample = 0;
  if( myStopSample >= myNumSamples ) myStopSample = myNumSamples-1;
  int numSamplesToProcess = myStopSample - myStartSample + 1;
  myNumWindows = ( numSamplesToProcess - myWidthRefWin ) / myIncWin + 1;

  myRatios = new float[myNumWindows];
  myValueRefWin = new float[myNumWindows];
  //  fprintf(stderr,"#win: %d, inc: %d, widthMean: %d, widthSpike: %d, widthRef: %d,  start/stop: %d %d\n",
  //         myNumWindows, myIncWin, myWidthMeanWin, myWidthSpikeWin, myWidthRefWin, myStartSample, myStopSample);
}

//***************************************************************************************
//
//
void csDespike::apply( float* samples, int numSamples, int& numSpikesFound, int& numSamplesAffected ) {

  // Compute amplitude ratios between defined spike window and background
  if( !computeRatios( samples, numSamples ) ) {
    return;
  }

  numSpikesFound      = 0;
  numSamplesAffected  = 0;
  int countSampleLast = 0;

  int widthRefWinHalf = myWidthRefWin / 2;
  int winIndex = 0;
  while( winIndex < myNumWindows ) {
    // If computed ratio exceeds the defined maximum, a spike has been detected
    if( myRatios[winIndex] > myMaxRatio ) {
      int win2 = winIndex + 1;
      double maxRatio = myRatios[winIndex];
      while( win2 < myNumWindows && myRatios[win2] > myMaxRatio ) {
        if( myRatios[win2] > maxRatio ) {
          maxRatio = myRatios[win2];
        }
        win2 += 1;
      }
      maxRatio *= myRatioAmplifier;
      int width = myWindowAmplifier*(win2 - winIndex);
      if( width < myWidthSpikeWin/2 ) width = myWidthSpikeWin/2;
      int winMid = ( win2 + winIndex ) / 2;
      int sampleMid   = winMid * myIncWin + myStartSample + widthRefWinHalf;
      int sampleFirst = std::max( 0, sampleMid - width );
      int sampleLast  = MIN( sampleMid + width, myNumSamples-1 );
      if( myMethod == COSINE_TAPER ) {
        for( int isamp = sampleFirst; isamp <= sampleLast; isamp++ ) {
          double phase = ( ( (double)(isamp-sampleFirst) / (double)width ) - 1.0 ) * M_PI;
          double scalar = 0.5 * ( cos(phase) * (-1.0 + 1.0/maxRatio) + (1.0 + 1.0/maxRatio) );
          samples[isamp] *= (float)scalar;
        }
      }
      else if( myMethod == LINEAR_INTERPOLATION ) {
        float val1 = samples[sampleFirst];
        float val2 = samples[sampleLast];
        float length = (float)(sampleLast - sampleFirst);
        for( int isamp = sampleFirst+1; isamp < sampleLast; isamp++ ) {
          float valCurrent = val1 + (float)(isamp-sampleFirst)/length * ( val2 - val1 );
          samples[isamp] = valCurrent;
        }
      }
      else {
        for( int isamp = sampleFirst; isamp <= sampleLast; isamp++ ) {
          samples[isamp] = 0;
        }
      }
      // Count identified spikes, and number of affected samples
      numSpikesFound     += 1;
      numSamplesAffected += sampleLast - std::max( sampleFirst, countSampleLast ) + 1;
      countSampleLast    = sampleLast;

      winIndex = win2;
    }
    winIndex += 1;
  }

}
//***************************************************************************************
//
//
bool csDespike::computeRatios( float* samples, int numSamples ) {
  if( myWidthRefWin >= myNumSamples ) return false;
  if( numSamples != myNumSamples ) return false;

  // Sorted sample lists: These are only needed for median option
  csVector<DespikePoint*> sampleListRefWin( myWidthRefWin );     // List of sample values in reference window
  csVector<DespikePoint*> sampleListSpikeWin( myWidthSpikeWin ); // List of sample values in spike window
  int indexMidSpikeWin = myWidthSpikeWin / 2;
  int indexMidRefWin   = myWidthRefWin   / 2;

  int sampleShiftSpikeWin = indexMidRefWin - indexMidSpikeWin;

  //---------------------------------------------------------------------
  // For first window, add sample values to sampleListRefWin and compute mean amplitude over reference window
  // Omit last sample so that the following loop also works for first window
  double currentSumRefWin   = 0.0;
  double currentSumSpikeWin = 0.0;
  double currentDebiasSumRefWin   = 0.0;
  double currentDebiasSumSpikeWin = 0.0;

  if( myPerformDebias ) {
    for( int i = 0; i < myWidthRefWin; i++ ) {
      currentDebiasSumRefWin += samples[i + myStartSample];
    }
    for( int i = sampleShiftSpikeWin; i < myWidthSpikeWin; i++ ) {
      currentDebiasSumSpikeWin += samples[i + myStartSample];
    }
  }

  for( int i = 0; i < myWidthRefWin; i++ ) {
    int sampleIndex = i + myStartSample;
    float value = fabs(samples[sampleIndex]);
    if( myPerformDebias ) {
      value = fabs( samples[sampleIndex] - (currentDebiasSumRefWin/(double)myWidthRefWin) );
    }
    if( myMethodRefWin == csDespike::METHOD_WIN_MEDIAN ) {
      insertNewValue( new DespikePoint( sampleIndex, value ), &sampleListRefWin );
      /*
      DespikePoint* p = new DespikePoint( sampleIndex, value );
      int is = 0;
      for( ; is < sampleListRefWin.size(); is++ ) {
        DespikePoint* ptmp = sampleListRefWin.at(is);
        if( p->value < ptmp->value ) {
          sampleListRefWin.insert( p, is );
          break;
        }
      }
      if( is == sampleListRefWin.size() ) {
        sampleListRefWin.insert( p, sampleListRefWin.size() );
      }
      */
    }
    currentSumRefWin += value;
  }

  for( int i = 0; i < myWidthSpikeWin; i++ ) {
    int sampleIndex = i + myStartSample + sampleShiftSpikeWin;
    float value = fabs(samples[sampleIndex]);
    if( myPerformDebias ) {
      value = fabs( samples[sampleIndex] - (currentDebiasSumSpikeWin/(double)myWidthSpikeWin) );
    }
    if( myMethodSpikeWin == csDespike::METHOD_WIN_MEDIAN ) {
      insertNewValue( new DespikePoint( sampleIndex, value ), &sampleListSpikeWin );
      /*
      DespikePoint* p = new DespikePoint( sampleIndex, value );
      int is = 0;
      for( ; is < sampleListSpikeWin.size(); is++ ){
        DespikePoint* ptmp = sampleListSpikeWin.at(is);
        if( p->value < ptmp->value ) {
          sampleListSpikeWin.insert( p, is );
          break;
        }
      }
      if( is == sampleListSpikeWin.size() ) {
        sampleListSpikeWin.insert( p, sampleListSpikeWin.size() );
      }
      */
    }
    currentSumSpikeWin += value;
  }

//----------------------------------------------------------------------
// Main loop
//


  for( int iwin = 0; iwin < myNumWindows; iwin++ ) {
    int sampleFirst = iwin * myIncWin + myStartSample;
    int sampleLast  = sampleFirst + myWidthRefWin - 1;
    float valueRef   = fabs(samples[sampleLast]);
    float valueSpike = fabs(samples[sampleLast-sampleShiftSpikeWin]);

    if( iwin > 0 ) {
      int sampleDelRef   = sampleFirst - 1;
      int sampleDelSpike = sampleFirst - 1 + sampleShiftSpikeWin;
      if( myPerformDebias ) {
        currentDebiasSumRefWin = currentDebiasSumRefWin - samples[sampleDelRef] + samples[sampleLast];
        valueRef =  fabs(samples[sampleLast] - (currentDebiasSumRefWin / (double)myWidthRefWin) );
        currentDebiasSumSpikeWin = currentDebiasSumSpikeWin - samples[sampleDelSpike] + samples[sampleLast-sampleShiftSpikeWin];
        valueSpike = (float)( samples[sampleLast-sampleShiftSpikeWin] - (currentDebiasSumSpikeWin / (double)myWidthSpikeWin) );
      }
      currentSumRefWin -= fabs(samples[sampleDelRef]);
      currentSumRefWin += valueRef;

      currentSumSpikeWin -= fabs(samples[sampleDelSpike]);
      currentSumSpikeWin += valueSpike;
    }
    if( iwin > 0 && myMethodRefWin == csDespike::METHOD_WIN_MEDIAN ) {
// (1) Determine median value in reference window (sort all values, pick middle one)
//      fprintf(stderr,"Median start Index: %d %d %d\n", sampleFirst, sampleLast, sampleShiftSpikeWin );

      //-------------------- Reference WINDOW -------------------------
      // Remove previous last sample from sample List
      DespikePoint* p = NULL;
      int sampleDelRef = sampleFirst - 1;
      int length = sampleListRefWin.size();
      int is = 0; 
      for( ; is < length; is++ ) {   // Currently assumes that window increment is 1
        p = sampleListRefWin.at(is);
        if( p->index == sampleDelRef ) {
          sampleListRefWin.remove(is);
          p->index = sampleLast;
          p->value = valueRef;
          break;
        }
      }
      if( is == length ) { // ??? Previous sample not found ???
        p = new DespikePoint( sampleLast, valueRef );
      }
      // Add new point by sliding the window one sample forwards
      insertNewValue( p, &sampleListRefWin );
    }
    if( iwin > 0 && myMethodSpikeWin == csDespike::METHOD_WIN_MEDIAN ) {
      //-------------------- Spike WINDOW -------------------------
      // Remove previous last sample from sample List
      DespikePoint* p = NULL;
      int sampleDelSpike = sampleFirst - 1 + sampleShiftSpikeWin;
      int length = sampleListSpikeWin.size();
      int is = 0; 
      for( ; is < length; is++ ) {   // Currently assumes that window increment is 1
        p = sampleListSpikeWin.at(is);
        if( p->index == sampleDelSpike ) {
          sampleListSpikeWin.remove(is);
          p->index = sampleLast - sampleShiftSpikeWin;
          p->value = valueSpike;
          break;
        }
      }
      if( is == length ) { // ??? Previous sample not found ???
        p = new DespikePoint( sampleLast - sampleShiftSpikeWin, valueSpike );
      }
      // Add new point by sliding the window one sample forwards
      insertNewValue( p, &sampleListSpikeWin );
    }

    double resultRef;
    double resultSpike;
    if( myMethodRefWin == csDespike::METHOD_WIN_MEAN ) {
      resultRef = currentSumRefWin / (double)myWidthRefWin;
    }
    else { // MEDIAN
      resultRef = sampleListRefWin.at( indexMidRefWin )->value;
    }
    if( resultRef == 0.0 ) resultRef = 1.0;
    if( myMethodSpikeWin == csDespike::METHOD_WIN_MEAN ) {
      resultSpike = currentSumSpikeWin / (double)myWidthSpikeWin;
    }
    else { // MEDIAN
      resultSpike = sampleListRefWin.at( indexMidSpikeWin )->value;
    }
    myRatios[iwin] = resultSpike / resultRef;
    //    myValueRefWin[iwin] = currentDebiasSumRefWin;
    //    fprintf(stdout,"%d %e %e  %e  %d\n", iwin, resultRef, resultSpike, myRatios[iwin], sampleMid);
  }
  // Free memory
  for( int is = 0; is < sampleListRefWin.size(); is++ ) {
    delete sampleListRefWin.at( is );
  }
  for( int is = 0; is < sampleListSpikeWin.size(); is++ ) {
    delete sampleListSpikeWin.at( is );
  }
  return true;
}
//***************************************************************************************
//
//
void csDespike::insertNewValue( DespikePoint* p, csVector<DespikePoint*>* sortedList ) {
  int id1 = 0;
  int id2 = sortedList->size() - 1;
  float value = p->value;

  if( id2 < 0 ) {
    sortedList->insertEnd(p);
    return;
  }
  //  fprintf(stdout,"ID %d %d - num: %d\n", id1, id2, sortedList->size());
  //  fflush(stdout);

  while( true ) {
    int id0 = (id1+id2)/2;
    float val0 = sortedList->at(id0)->value;
    if( value > val0 ) {
      if( id2-id1 <= 1 ) {
        if( value > sortedList->at(id2)->value ) {
          sortedList->insert( p, id2+1 );
        }
        else {
          sortedList->insert( p, id1+1 );
        }
        break;
      }
      id1 = id0;
    }
    else if( value < val0 ) {
      if( id2-id1 <= 1 ) {
        sortedList->insert( p, id1 );
        break;
      }
      id2 = id0;
    }
    else {
      sortedList->insert( p, id0 );
      break;
    }
  }
}
void csDespike::getDefaultFrequencySpikeConfig( DespikeConfig& config ) {
  config.widthRefWin = 1.0;
  config.incWin      = 0;
  config.start       = 0;
  config.stop        = 0;
  config.minWidthSpikeWin = 0.1;
  config.advanced.widthMeanWin_inSamples = 1;
  config.advanced.ratioAmplifier  = 1;
  config.advanced.windowAmplifier = 1;
  config.performDebias = true;
  config.maxRatio      = 3.0;
  config.method        = csDespike::LINEAR_INTERPOLATION;
  config.methodRefWin     = csDespike::METHOD_WIN_MEDIAN;
  config.methodSpikeWin     = csDespike::METHOD_WIN_MEAN;
}
void csDespike::getDefaultTimeNoiseBurstConfig( DespikeConfig& config ) {
  config.widthRefWin = 3000;
  config.incWin      = 0;
  config.start       = 0;
  config.stop        = 0;
  config.minWidthSpikeWin = 300;
  config.advanced.widthMeanWin_inSamples = 11;
  config.advanced.ratioAmplifier  = 2;
  config.advanced.windowAmplifier = 4;
  config.performDebias = false;
  config.maxRatio   = 5.0;
  config.method     = csDespike::COSINE_TAPER;
  config.methodRefWin     = csDespike::METHOD_WIN_MEDIAN;
  config.methodSpikeWin     = csDespike::METHOD_WIN_MEAN;
}

