/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csFFTDesignature.h"
#include "csFFT.h"
#include "csException.h"
#include "geolib_math.h"
#include "geolib_defines.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace cseis_geolib;

/*
 * numSamples (input wavelet)
 */

csFFTDesignature::csFFTDesignature( int numSamples, float sampleInt ) :
  csFFTTools( numSamples, sampleInt )
{
  myDesigAmpFilter  = NULL;
  myDesigPhaseShift = NULL;
  myAmpPhaseSpecIn   = NULL;
  myFilterType = AMP_PHASE;
  // Nothing more to do. User needs to call initialize next
}

csFFTDesignature::csFFTDesignature( int numSamples, float sampleInt, float const* input_wavelet, float timeZero_s, float percWhiteNoise, float const* output_wavelet ) :
  csFFTTools( numSamples, sampleInt )
{
  myDesigAmpFilter  = NULL;
  myDesigPhaseShift = NULL;
  myAmpPhaseSpecIn  = NULL;
  myFilterType = AMP_PHASE;
  initDesignature( input_wavelet, timeZero_s, percWhiteNoise, output_wavelet );
}
void csFFTDesignature::initialize( float const* input_wavelet, float timeZero_s, float percWhiteNoise, float const* output_wavelet ) {
  initDesignature( input_wavelet, timeZero_s, percWhiteNoise, output_wavelet );
}

csFFTDesignature::~csFFTDesignature() {
  if( myDesigAmpFilter != NULL ) {
    delete [] myDesigAmpFilter;
    myDesigAmpFilter  = NULL;
  }
  if( myDesigPhaseShift != NULL ) {
    delete [] myDesigPhaseShift;
    myDesigPhaseShift  = NULL;
  }
  if( myAmpPhaseSpecIn != NULL ) {
    delete [] myAmpPhaseSpecIn;
    myAmpPhaseSpecIn = NULL;
  }
}

//--------------------------------------------------------------------------------
//
//
//
void csFFTDesignature::initDesignature( float const* input_wavelet, float timeZero_s, float percWhiteNoise, float const* output_wavelet ) {
  int numFreq = myFFT->numFreqValues();

  myAmpPhaseSpecIn  = new float[2*numFreq];
  myDesigAmpFilter  = new float[numFreq];
  myDesigPhaseShift = new float[numFreq];

  // Forward transform wavelet
  myFFT->forwardTransform( input_wavelet, myAmpPhaseSpecIn, cseis_geolib::FX_AMP_PHASE );

  // Compute maximum amplitude, for application of percentage white noise
  float maxAmp = fabs( myAmpPhaseSpecIn[0] );
  for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
    float valueAbs = fabs( myAmpPhaseSpecIn[ifreq] );
    if( valueAbs > maxAmp ) maxAmp = valueAbs;
  }
  double whiteNoise = (percWhiteNoise/100.0) * maxAmp;
  if( whiteNoise <= 0.0 ) whiteNoise = 1e-50;

  double df = sampleIntFreq();
  // Spiking filter:
  if( output_wavelet == NULL ) {
    // Compute inverse (designature) filter = Inverse amplitude spectrum after adding white noise
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      myDesigAmpFilter[ifreq] = (float)( maxAmp / (myAmpPhaseSpecIn[ifreq] + whiteNoise) );
    }
    // Compute phase shift for zero-phasing
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      double freq          = df * (double)(ifreq);
      double phaseShift    = 2.0 * M_PI * freq * timeZero_s;
      myDesigPhaseShift[ifreq] = (float)fmod( phaseShift - myAmpPhaseSpecIn[ifreq+numFreq] , 2.0*M_PI );
      if( myDesigPhaseShift[ifreq] > M_PI ) {
        myDesigPhaseShift[ifreq] -= (float)( 2.0 * M_PI );
      }
    }
  }
  // 'Transfer' filter from input to output wavelet
  else {
    float* ampPhaseSpecOut = new float[2*numFreq];
    myFFT->forwardTransform( output_wavelet, ampPhaseSpecOut, cseis_geolib::FX_AMP_PHASE );

    /* 
       float maxAmpOut = fabs( ampPhaseSpecOut[0] );
       for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
       float valueAbs = fabs( ampPhaseSpecOut[ifreq] );
       if( valueAbs > maxAmpOut ) maxAmpOut = valueAbs;
       }
    */
    // Compute inverse (designature) filter = Inverse amplitude spectrum after adding white noise
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      myDesigAmpFilter[ifreq] =  (float)( ampPhaseSpecOut[ifreq] / (myAmpPhaseSpecIn[ifreq] + whiteNoise) );
    }

    // Compute phase shift
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      int indexPhase = ifreq+numFreq;
      myDesigPhaseShift[ifreq] = (float)fmod( ampPhaseSpecOut[indexPhase] - myAmpPhaseSpecIn[indexPhase] , (float)(2.0*M_PI) );
      if( myDesigPhaseShift[ifreq] > M_PI ) {
        myDesigPhaseShift[ifreq] -= (float)( 2.0*M_PI );
      }
    }
    delete [] ampPhaseSpecOut;
  }
}
//--------------------------------------------------------------------------------
//
// TODO: Apply optional taper to data before FFT transform
//
void csFFTDesignature::applyFilter( float* samples, int numSamples ) {
  if( numSamples != myNumSamplesIn ) {
    throw(csException("csFFTDesignature::applyFilter(): Inconsistent number of samples: %d != %d", numSamples, myNumSamplesOut));
  }

  // Forward transform input data and compute amplitude & phase spectrum
  myFFT->forwardTransform( samples, myAmpPhaseSpecIn, cseis_geolib::FX_AMP_PHASE );

  int numFreq = myFFT->numFreqValues();
  if( myFilterType == csFFTDesignature::AMP_PHASE ) {
    // Apply inverse filter & zero-phasing to amplitude & phase spectrum
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      myAmpPhaseSpecIn[ifreq]         *= myDesigAmpFilter[ifreq];
      myAmpPhaseSpecIn[ifreq+numFreq] += myDesigPhaseShift[ifreq];
    }
  }
  else if( myFilterType == csFFTDesignature::AMP_ONLY ) {
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      myAmpPhaseSpecIn[ifreq]   *= myDesigAmpFilter[ifreq];
    }
  }
  else if( myFilterType == csFFTDesignature::PHASE_ONLY ) {
    for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
      myAmpPhaseSpecIn[ifreq+numFreq] += myDesigPhaseShift[ifreq];
    }
  }

  // Transform back to X-T. Do not perform normalisation in fft ???
  myFFT->inverseTransform( myAmpPhaseSpecIn, samples, cseis_geolib::FX_AMP_PHASE );
}

void csFFTDesignature::setDesigFilterType( int filterType ) {
  myFilterType = filterType;
}

void csFFTDesignature::setDesigLowPass( float cutOffHz, float slope ) {
  double G0 = 1.0;
  double power = slope / 3.0f;
  double df = sampleIntFreq();

  int numFreq = myFFT->numFreqValues();
  for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
    double freq = (double)ifreq * df;
    float dampG = (float)sqrt( G0 / (1.0 + pow(freq/cutOffHz,power) ) );
    myDesigAmpFilter[ifreq] *= (float)dampG;
  }
}

void csFFTDesignature::setDesigHighPass( float cutOffHz, float slope ) {
  double G0 = 1.0;
  double power = slope / 3.0f;
  double df = sampleIntFreq();

  myDesigAmpFilter[0] = 0.0;

  int numFreq = myFFT->numFreqValues();
  for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
    double freq = (double)ifreq * df;
    float dampG = (float)sqrt( G0 / (1.0 + pow(cutOffHz/freq,power) ) );
    myDesigAmpFilter[ifreq] *= dampG;
  }
}

void csFFTDesignature::setDesigHighEnd( float freq ) {
  double df = sampleIntFreq();
  int numFreq = myFFT->numFreqValues();
  int index = (int)round( freq / df );
  if( index <= 0 || index >= numFreq ) return;
  float amplitude = myDesigAmpFilter[index];
  float phase     = myDesigPhaseShift[index];

  for( int ifreq = index+1; ifreq < numFreq; ifreq++ ) {
    myDesigAmpFilter[ifreq]  = amplitude;
    myDesigPhaseShift[ifreq] = phase;
  }
}

// Apply cosine taper around notch frequency
void csFFTDesignature::setNotchSuppression( float notchFreq, float notchWidth ) {
  double df = sampleIntFreq();
  int numFreq = myFFT->numFreqValues();
  int indexFirst = (int)round( ( notchFreq - 0.5*notchWidth ) / df );
  int indexLast  = (int)round( ( notchFreq + 0.5*notchWidth ) / df );
  int width = indexLast - indexFirst;

  int indexFirstRed = std::max(0,indexFirst);
  int indexLastRed  = std::min(numFreq-1,indexLast);
  /*
  myIsNotchSuppression = true;
  myNotchIndexFirst    = indexFirst;
  myNotchWidth = width;
  myNotchIndexFirstRed = indexFirstRed;
  myNotchIndexLastRed  = indexLastRed;
  */
  for( int ifreq = indexFirstRed; ifreq <= indexLastRed; ifreq++ ) {
    double phase  = 2.0*( ( (double)(ifreq-indexFirst) / (double)width ) - 1.0 ) * M_PI;
    float scalar  = (float)( 0.5 * (cos(phase) + 1.0) );
    myDesigAmpFilter[ifreq] = myDesigAmpFilter[ifreq] * scalar;
  }
}

//--------------------------------------------------------------------------------

void csFFTDesignature::dump_spectrum( FILE* stream ) const {
  double df = sampleIntFreq();
  int numFreq = myFFT->numFreqValues();
  for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
    double freq = df * (double)(ifreq);
    fprintf(stream,"%.6f  %.10e %.10e\n", freq, myDesigAmpFilter[ifreq], myDesigPhaseShift[ifreq] );
  }
}

void csFFTDesignature::dump_wavelet( FILE* stream, bool doNormalize, float timeShift_s ) {
  double df = sampleIntFreq();
  int numFreq = myFFT->numFreqValues();
  float* desigAmpPhaseFilter  = new float[2*numFreq];
  float* samples              = new float[myNumSamplesIn];

  for( int ifreq = 0; ifreq < numFreq; ifreq++ ) {
    double freq = df * (float)(ifreq);
    float phaseShift    = (float)( 2.0 * freq * M_PI * timeShift_s );
    desigAmpPhaseFilter[ifreq]         = myDesigAmpFilter[ifreq];
    desigAmpPhaseFilter[ifreq+numFreq] = myDesigPhaseShift[ifreq] + phaseShift;
  }
  // Transform back to X-T
  myFFT->inverseTransform( desigAmpPhaseFilter, samples, cseis_geolib::FX_AMP_PHASE, false );

  float normScalar = 1.0f;
  if( doNormalize ) normScalar = 1.0f / (float)numFreq;
  for( int isamp = 0; isamp < myNumSamplesIn; isamp++ ) {
    float time = mySampleIntIn * (float)isamp;
    fprintf(stream,"%.6f  %.10e\n", time, normScalar*samples[isamp]  );
  }
  delete [] desigAmpPhaseFilter;
  delete [] samples;
}
