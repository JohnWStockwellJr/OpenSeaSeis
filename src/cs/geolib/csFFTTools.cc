/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csFFTTools.h"
#include "csFFT.h"
#include "csException.h"
#include "csTimeFunction.h"
#include "csInterpolation.h"
#include "geolib_math.h"
#include "geolib_defines.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

using namespace cseis_geolib;

csFFTTools::csFFTTools( int numSamples, float sampleInt ) {
  init( numSamples, sampleInt );

  myFFT    = new csFFT( myNumSamplesIn );
  myBuffer = new float[2*myFFT->numFreqValues()];
}

csFFTTools::csFFTTools( int numSamplesIn, int numSamplesOut, float sampleIntIn, float sampleIntOut ) {
  init( numSamplesIn, sampleIntIn );

  myNumSamplesOut = numSamplesOut;
  mySampleIntOut  = sampleIntOut;
  myFFTOut        = new csFFT( myNumSamplesOut );
  myBufferOut     = new float[2*myFFTOut->numFreqValues()];

  int numFreqOut  = myFFTOut->numFreqValues();
  float ratio     = mySampleIntOut / mySampleIntIn;
  int ratio_int   = (int)round(ratio);
  int numFFTValuesIn_forced = 2 * ratio_int * ( numFreqOut - 1 );

  myFFT    = new csFFT( myNumSamplesIn, numFFTValuesIn_forced );
  myBuffer = new float[2*myFFT->numFreqValues()];
}
int csFFTTools::numFFTSamples() const {
  return myFFT->numFFTValues();
}
int csFFTTools::numFFTSamplesOut() const {
  return myFFTOut->numFFTValues();
}
//--------------------------------------------------------------------------------
//
//
void csFFTTools::init( int numSamples_time, float sampleInt_time ) {
  mySampleIntIn  = sampleInt_time;
  myNumSamplesIn = numSamples_time;
  mySampleIntOut  = mySampleIntIn;
  myNumSamplesOut = myNumSamplesIn;

  myFFT    = NULL;
  myFFTOut = NULL;

  myBuffer = NULL;
  myBufferOut = NULL;

  myFilterWavelet = NULL;
  myFilterScalars = NULL;
  myLengthFilterWavelet = 0;
  myIsFilterWavelet     = false;

  myNotchFilterScalars = NULL;
}

//--------------------------------------------------------------------------------
//
//
csFFTTools::~csFFTTools() {
  if( myBuffer != NULL ) {
    delete [] myBuffer;
    myBuffer = NULL;
  }
  if( myBufferOut != NULL ) {
    delete [] myBufferOut;
    myBufferOut = NULL;
  }
  if( myNotchFilterScalars == NULL ) {
    delete []   myNotchFilterScalars;
    myNotchFilterScalars = NULL;
  }

  if( myFilterScalars != NULL ) {
    delete [] myFilterScalars;
    myFilterScalars = NULL;
  }
  if( myFilterWavelet != NULL ) {
    delete [] myFilterWavelet;
    myFilterWavelet = NULL;
  }
  if( myFFT != NULL ) {
    delete myFFT;
    myFFT = NULL;
  }
}
//--------------------------------------------------------------------------------
//
//
void csFFTTools::setFilter( float cutOffFreqHz, float slope ) {
  prepareBandpassFilter( 1, &cutOffFreqHz, &slope, 0, NULL, NULL );
}
void csFFTTools::setFilterWavelet( int length ) {
  myIsFilterWavelet = true;
  if( myFilterWavelet != NULL ) {
    delete [] myFilterWavelet;
  }
  myLengthFilterWavelet = length;
  myFilterWavelet = new double[myLengthFilterWavelet];
}
//--------------------------------------------------------------------------------
//
//
double csFFTTools::sampleIntFreq() const {
  double freqNyquist = 1.0 / ( 2.0 * ( (double)mySampleIntIn / 1000.0 ) );
  return( ( freqNyquist / (double)(myFFT->numFreqValues()-1) ) );
}
//--------------------------------------------------------------------------------
//
//
void csFFTTools::prepareBandpassFilter( int numLowPassFilters, float const* lowPassFreq, float const* lowPassSlope, int numHighPassFilters, float const* highPassFreq, float const* highPassSlope ) {
  double G0 = 1.0;
  double df = sampleIntFreq();

  if( myFilterScalars != NULL ) delete [] myFilterScalars;

  int numFreq = myFFT->numFreqValues();
  myFilterScalars = new float[numFreq];

  for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
    myFilterScalars[ifreq] = 1.0f;
  }
  for( int i = 0; i < numLowPassFilters; i++ ) {
    float slope   = lowPassSlope[i];
    double power  = fabs(slope/3.0f);
    float cutOffFreq = lowPassFreq[i];
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      double freq = (double)ifreq * df;
      float dampG = (float)sqrt(G0 / (1.0 + pow(freq/cutOffFreq,power) ));
      if( slope <= 0 && dampG == 0.0 ) dampG = 1.0;
      myFilterScalars[ifreq] *= dampG;
    }
  }
  for( int i = 0; i < numHighPassFilters; i++ ) {
    float slope   = highPassSlope[i];
    double power  = fabs(slope/3.0f);
    float cutOffFreq = highPassFreq[i];
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      double freq = (double)ifreq * df;
      float dampG = (float)sqrt(G0 / (1.0 + pow(cutOffFreq/freq,power) ));
      if( freq == 0 ) {
        if( slope > 0 ) {
          dampG = 0.0f;
        }
        else {
          dampG = 1.0f;
        }
      }
      else if( slope <= 0 ) {
        if( dampG == 0.0f ) dampG = 1.0f;
      }
      myFilterScalars[ifreq] *= dampG;
    }
  }
}

void csFFTTools::applyBandpassFilter( float* samples, bool outputImpulseResponse ) {
  if( myFilterScalars == NULL ) throw( csException("csFFTTools::applyBandpassFilter: Filters have not been initialized. Call prepareBandpassFilter() first") );
  if( outputImpulseResponse ) {
    memset( samples, 0, myNumSamplesIn * sizeof(float) );
    samples[myNumSamplesIn/2] = 1.0f/mySampleIntIn;
  }
  myFFT->applyAmpFilter( samples, myFilterScalars );
}

float csFFTTools::resample( float* samples, bool applyFilter, bool applyRMSNorm ) {
  float rmsIn  = 1.0;
  if( applyRMSNorm ) rmsIn = compute_rms( samples, myNumSamplesIn );

  myFFT->forwardTransform( samples, myBuffer, cseis_geolib::FX_AMP_PHASE );

  int numFreqIn  = myFFT->numFreqValues();
  int numFreqOut = myFFTOut->numFreqValues();

  myBuffer[0] = 0.0;         // Set DC amp to 0
  myBuffer[numFreqIn] = 0.0; // Set DC phase to 0

  if( applyFilter ) {
    if( myFilterScalars == NULL ) throw( csException("csFFTTools::resample: Anti-alias filter has not been initialized.") );
    for( int ifreq = 0; ifreq < numFreqIn; ifreq++ ) {
      myBuffer[ifreq] *= myFilterScalars[ifreq];
    }
  }

  for( int ifreqOut = 0; ifreqOut < numFreqOut; ifreqOut++ ) {
    myBufferOut[ifreqOut]            = myBuffer[ifreqOut];    // Amplitude
    myBufferOut[ifreqOut+numFreqOut] = myBuffer[ifreqOut+numFreqIn]; // Phase
  }

  myFFTOut->inverseTransform( myBufferOut, samples, cseis_geolib::FX_AMP_PHASE );

  if( applyRMSNorm ) {
    float rmsOut = compute_rms( samples, myNumSamplesOut );
    if( rmsOut != 0.0 ) {
      float ratio = rmsIn / rmsOut;
      for( int i = 0; i < myNumSamplesOut; i++ ) {
        samples[i] *= ratio;
      }
      return ratio;
    }
  }
  return 1.0f;
}

//--------------------------------------------------------------------------------
// Hilbert transform == Shift phase by 90deg
//
void csFFTTools::hilbertTransform( float* samplesInOut ) {
  hilbertTransform( samplesInOut, samplesInOut );
}
void csFFTTools::hilbertTransform( float const* samplesIn, float* samplesOut) {

  myFFT->forwardTransform( samplesIn, myBuffer, cseis_geolib::FX_AMP_PHASE );

  int numFreq    = myFFT->numFreqValues();
  int numFreq_x2 = 2*numFreq;
  float pi_half  = (float)( M_PI / 2.0 );
  for( int ifreq = numFreq; ifreq < numFreq_x2; ifreq++ ) {
    myBuffer[ifreq] -= pi_half; // Hilbert transform == quadrature phase == 90deg phase shift
  }

  myFFT->inverseTransform( myBuffer, samplesOut, cseis_geolib::FX_AMP_PHASE );

  //  for( int isamp = 0; isamp < myNumSamplesOut; isamp++ ) {
  //   samplesOut[isamp] *= 0.5f;
  //  }
}

//--------------------------------------------------------------------------------
void csFFTTools::envelope( float* samples ) {
  if( myBufferOut == NULL ) {
    myBufferOut = new float[myNumSamplesIn];
  }

  // 1) Hilbert transform == Shift phase by 90deg, stored in myBuffer
  hilbertTransform( samples, myBufferOut );

  // 2) Combine input trace and 90deg shifted trace to form envelope
  for( int isamp = 0; isamp < myNumSamplesIn; isamp++ ) {
    float sampleInput = samples[isamp];
    float sampleHilb  = myBufferOut[isamp];
    samples[isamp] = (float)sqrt( sampleInput*sampleInput + sampleHilb*sampleHilb );
  } 
}
//--------------------------------------------------------------------------------
void csFFTTools::inst_phase( float* samples ) {
  if( myBufferOut == NULL ) {
    myBufferOut = new float[myNumSamplesIn];
  }

  hilbertTransform( samples, myBufferOut );

  // Compute instantaneous phase
  // analytic_signal(x) = signal(x) + i * hilbert(x)
  for( int isamp = 0; isamp < myNumSamplesIn; isamp++ ) {
    float phase_inst = (float) atan2( -myBufferOut[isamp], samples[isamp] );
    samples[isamp] = phase_inst;
  }
}
//--------------------------------------------------------------------------------
void csFFTTools::inst_freq( float* samples ) {
  inst_phase( samples );

  // Compute instantaneous frequency
  double phaseThreshold = M_PI - 1e-5;
  double phasePrev = 0;
  for( int isamp = 0; isamp < myNumSamplesIn; isamp++ ) {
    double phase     = samples[isamp];
    double phaseDiff = phase - phasePrev;
    if( phaseDiff > phaseThreshold ) phaseDiff -= 2.0 * M_PI;
    else if( phaseDiff < -phaseThreshold ) phaseDiff += 2.0 * M_PI;
    samples[isamp] = (float)( fabs(phaseDiff) / ( 2.0 * M_PI ) * 1000.0 / mySampleIntIn );
    phasePrev = phase;
  }
}

void csFFTTools::unwrap_phase() {
  int numFreqIn   = myFFT->numFreqValues();
  float* phasePtr = &myBuffer[numFreqIn];
  float* phaseTmp = new float[numFreqIn];
  double ph0  = phasePtr[0];
  phaseTmp[0] = (float)ph0;
  double threshold = M_PI - 1e-5;
  double po = 0.0;
  for( int ifreq = 1; ifreq < numFreqIn; ifreq++ ) {
    double cp = phasePtr[ifreq] + po;
    double dp = cp - ph0;
    ph0 = cp;
    if( dp > threshold ) {
      while( dp > threshold ) {
        po -= 2.0*M_PI;
        dp -= 2.0*M_PI;
      }
    }
    else if( fabs(dp) > threshold ) {
      while( fabs(dp) > threshold ) {
        po += 2.0*M_PI;
        dp += 2.0*M_PI;
      }
    }
    phaseTmp[ifreq] = phasePtr[ifreq] + (float)po;
    ph0 = phaseTmp[ifreq];
  }
  memcpy( phasePtr, phaseTmp, sizeof(float) * numFreqIn );
  delete [] phaseTmp;
}

float const* csFFTTools::prepareNotchFilter( float notchFreqHz, float notchWidthHz, float slope, bool isCosineTaper ) {
  int numFreq = myFFT->numFreqValues();
  double df = sampleIntFreq();
  double freq1 = notchFreqHz - 0.5*notchWidthHz;
  double freq2 = notchFreqHz + 0.5*notchWidthHz;
  int indexFirst = (int)round( freq1 / df );
  int indexLast  = (int)round( freq2 / df );
  int width = indexLast - indexFirst;
  if( width == 0 ) {
    indexFirst = std::max(0,indexFirst-1);
    indexLast  = std::min(indexFirst + 2,numFreq-1);
    width = indexLast - indexFirst;
  }
  int indexFirstRed = std::max(0,indexFirst);
  int indexLastRed  = std::min(numFreq-1,indexLast);

  myNotchFilterScalars = new float[numFreq];
  //  myNotchFilterScalars = new double[myNumFFTSamplesIn];
  //  for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
  //  myNotchFilterScalars[ifreq] = 0.0;
  // }
  memset( myNotchFilterScalars, 0, numFreq * sizeof(float) );

  if( isCosineTaper ) {
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      myNotchFilterScalars[ifreq] = 1.0;
    }
    for( int ifreq = indexFirstRed; ifreq <= indexLastRed; ifreq++ ) {
      double phase  = 2.0*( ( (double)(ifreq-indexFirst) / (double)width ) - 1.0 ) * M_PI;
      double scalar = 0.5 * (cos(phase) + 1.0);
      myNotchFilterScalars[ifreq] = (float)scalar;
    }
  }
  else { // Apply Butterworth notch filter
    double G0 = 1.0;
    double power = fabs(slope) / 3.0f;

    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      double freq = (double)ifreq * df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/freq1,power) ));
      if( ifreq > 0 ) {
        dampG += sqrt(G0 / (1.0 + pow(freq2/freq,power) ));    // HIGHPASS
      }
      myNotchFilterScalars[ifreq] += (float)dampG;
    }
  }
  
  // Dampen potential side lobes of filter close to notch (must be done for positive and negative frequency portion of complex spectrum)
  double maxValue1 = 0;
  int maxIndex1 = 0;
  int indexNotch  = (int)round( notchFreqHz / df );
  for( int ifreq = 0; ifreq < indexNotch; ifreq++ ) {
    if( myNotchFilterScalars[ifreq] > maxValue1 ) {
      maxValue1 = myNotchFilterScalars[ifreq];
      maxIndex1 = ifreq;
    }
  }
  if( maxIndex1 != 0 ) {
    int maxIndex2 = 0;
    double maxValue2 = 0;
    for( int ifreq = indexNotch; ifreq < 2*indexNotch; ifreq++ ) {
      if( myNotchFilterScalars[ifreq] > maxValue2 ) {
        maxValue2 = myNotchFilterScalars[ifreq];
        maxIndex2 = ifreq;
      }
    }
    float correctionTerm = (float)( 1.0 - 0.5*( maxValue1 + maxValue2 ) );
    for( int ifreq = 0; ifreq < maxIndex1; ifreq++ ) {
      if( myNotchFilterScalars[ifreq] > 1.0f ) myNotchFilterScalars[ifreq] = 1.0f;
    }
    for( int ifreq = maxIndex1; ifreq <= maxIndex2; ifreq++ ) {
      myNotchFilterScalars[ifreq] += correctionTerm;
    }
    int maxIndex3 = numFreq - maxIndex2;
    int maxIndex4 = numFreq - maxIndex1;
    for( int ifreq = maxIndex2+1; ifreq < maxIndex3; ifreq++ ) {
      if( myNotchFilterScalars[ifreq] > 1.0f ) myNotchFilterScalars[ifreq] = 1.0f;
    }
    for( int ifreq = maxIndex3; ifreq <= maxIndex4; ifreq++ ) {
      myNotchFilterScalars[ifreq] += correctionTerm;
    }
    for( int ifreq = maxIndex4+1; ifreq < numFreq; ifreq++ ) {
      if( myNotchFilterScalars[ifreq] > 1.0f ) myNotchFilterScalars[ifreq] = 1.0f;
    }
  }
  return myNotchFilterScalars;
}
//--------------------------------------------------------------------------------
// Apply notch filter
//
void csFFTTools::applyNotchFilter( float* samples, bool outputImpulseResponse ) {
  if( myNotchFilterScalars == NULL ) throw( csException("csFFTTools::applyNotchFilter: Filters have not been initialized. Call prepareNotchFilter() first") );
  if( outputImpulseResponse ) {
    memset( samples, 0, myNumSamplesIn * sizeof(float) );
    samples[myNumSamplesIn/2] = 1.0f/mySampleIntIn;  // 
  }
  myFFT->applyAmpFilter( samples, myNotchFilterScalars );
}

//--------------------------------------------------------------------------------
// Apply filter defined by function
//
void csFFTTools::applyFilter( float* samples, csTimeFunction<double> const* freqFunc ) {
  int numFreq = myFFT->numFreqValues();
  if( myFilterScalars == NULL ) {
    myFilterScalars = new float[numFreq];
  }
  double df = sampleIntFreq();
  for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
    double freq = df * (double)ifreq;
    myFilterScalars[ifreq] = (float)freqFunc->valueAt( freq );
  }
  myFFT->applyAmpFilter( samples, myFilterScalars );
}
//--------------------------------------------------------------------------------
// Apply filter defined by coefficients
//
void csFFTTools::applyAmpFilter( float* samplesInOut, float const* filterCoefficients ) {
  myFFT->applyAmpFilter( samplesInOut, filterCoefficients );
}

