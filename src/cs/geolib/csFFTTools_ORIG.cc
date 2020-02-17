/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csFFTTools_ORIG.h"
#include "csException.h"
#include "csTimeFunction.h"
#include "geolib_math.h"
#include "geolib_defines.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

using namespace cseis_geolib;

csFFTTools_ORIG::csFFTTools_ORIG( int numSamples ) {
  myNumSamplesIn  = numSamples;
  myNumSamplesOut = numSamples;
  myOrder = 0;
  myCutOffFreqHz = 0;
  myAmpSpecConsistencyScalar = 2.0f;
  init();
}

csFFTTools_ORIG::csFFTTools_ORIG( int numSamples, float sampleInt ) {
  myNumSamplesIn  = numSamples;
  mySampleIntIn   = sampleInt;
  myNumSamplesOut = myNumSamplesIn;
  mySampleIntOut  = mySampleIntIn;
  myAmpSpecConsistencyScalar = sampleInt;

  myOrder = 4;
  float freqNyquist = (float)(500.0/mySampleIntIn);
  myCutOffFreqHz = freqNyquist;

  init();
}

csFFTTools_ORIG::csFFTTools_ORIG( int theNumSamplesIn, int theNumSamplesOut, float theSampleIntIn, float theSampleIntOut ) {
  myNumSamplesIn  = theNumSamplesIn;
  mySampleIntIn   = theSampleIntIn;
  myNumSamplesOut = theNumSamplesOut;
  mySampleIntOut  = theSampleIntOut;
  myAmpSpecConsistencyScalar = 2.0f;

  myOrder  = 20;
  float freqNyquist = (float)(500.0/mySampleIntIn);
  double ratio = (double)mySampleIntOut / (double)mySampleIntIn;
  myCutOffFreqHz = (float)( freqNyquist / ratio ) * 0.8f;

  init();
}

void csFFTTools_ORIG::setVintageAmpSpecNormalization() {
  myAmpSpecConsistencyScalar = 2.0f;  
}
//--------------------------------------------------------------------------------
//
//
void csFFTTools_ORIG::init() {
  myNumFFTSamplesIn  = myNumSamplesIn;
  myNumFFTSamplesOut = myNumSamplesOut;

  myFilterWavelet = NULL;
  myLengthFilterWavelet = 0;
  myIsFilterWavelet     = false;
  myFilterScalars = NULL;

  int two_power_m;
  csFFTTools_ORIG::Powerof2( myNumFFTSamplesIn, &myTwoPowerIn, &two_power_m );
  if( two_power_m != myNumFFTSamplesIn ) {
    myNumFFTSamplesIn = two_power_m * 2;
    myTwoPowerIn += 1;
  }

  myNumFFTSamplesOut = myNumSamplesOut;
  csFFTTools_ORIG::Powerof2( myNumFFTSamplesOut, &myTwoPowerOut, &two_power_m );
  if( two_power_m != myNumFFTSamplesOut ) {
    myNumFFTSamplesOut = two_power_m * 2;
    myTwoPowerOut += 1;
  }
  myBufferReal = new double[myNumFFTSamplesIn];
  myBufferImag = new double[myNumFFTSamplesIn];
  myNotchFilter = NULL;

  myOutputImpulseResponse = false;
}

//--------------------------------------------------------------------------------
//
//
csFFTTools_ORIG::~csFFTTools_ORIG() {
  if( myBufferReal != NULL ) {
    delete [] myBufferReal;
    myBufferReal = NULL;
  }
  if( myBufferImag != NULL ) {
    delete [] myBufferImag;
    myBufferImag = NULL;
  }
  if( myNotchFilter == NULL ) {
    delete []   myNotchFilter;
    myNotchFilter = NULL;
  }

  if( myFilterScalars != NULL ) {
    delete [] myFilterScalars;
    myFilterScalars = NULL;
  }
  if( myFilterWavelet != NULL ) {
    delete [] myFilterWavelet;
    myFilterWavelet = NULL;
  }
}
//--------------------------------------------------------------------------------
//
//
void csFFTTools_ORIG::setFilter( float cutOffFreqHz, float order, bool outputImpulseResponse ) {
  myCutOffFreqHz = cutOffFreqHz;
  myOrder = order;
  myOutputImpulseResponse = outputImpulseResponse;
}
void csFFTTools_ORIG::setFilterWavelet( int length ) {
  myIsFilterWavelet = true;
  if( myFilterWavelet != NULL ) {
    delete [] myFilterWavelet;
  }
  myLengthFilterWavelet = length;
  myFilterWavelet = new double[myLengthFilterWavelet];
}
//--------------------------------------------------------------------------------
bool csFFTTools_ORIG::fft_forward( float const* samples ) {
  setBuffer( samples );
  return csFFTTools_ORIG::fft( csFFTTools_ORIG::FORWARD, myTwoPowerIn, myBufferReal, myBufferImag );
}
//--------------------------------------------------------------------------------
bool csFFTTools_ORIG::fft_forward( float const* samples, float* ampSpec ) {
  return fft_forward( samples, ampSpec, NULL );
}
bool csFFTTools_ORIG::fft_forward( float const* samples, float* ampSpec, float* phaseSpec ) {
  setBuffer( samples );
  bool success = csFFTTools_ORIG::fft( csFFTTools_ORIG::FORWARD, myTwoPowerIn, myBufferReal, myBufferImag );
  if( !success ) return false;
  
  convertToAmpPhase( ampSpec, phaseSpec );
  return true;
}
//--------------------------------------------------------------------------------
// Convert real/imaginary spectrum to to amp/phase
//
void csFFTTools_ORIG::convertToAmpPhase( float* ampSpec, float* phaseSpec ) {
  if( ampSpec != NULL ) {
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      ampSpec[i] = (float)( myAmpSpecConsistencyScalar * sqrt(myBufferReal[i]*myBufferReal[i] + myBufferImag[i]*myBufferImag[i]) );
    }
  }
  if( phaseSpec != NULL ) {
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      phaseSpec[i] = (float)atan2( -myBufferImag[i], myBufferReal[i] );
    }
  }
}
//--------------------------------------------------------------------------------
// Convert amp/phase spectrum to to real/imaginary
//
void csFFTTools_ORIG::convertFromAmpPhase( float const* ampSpec, float const* phaseSpec ) {
  float scalar = 1.0f / myAmpSpecConsistencyScalar;
  for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
    myBufferReal[i] = scalar * ampSpec[i] * cos(phaseSpec[i]);
    myBufferImag[i] = -scalar * ampSpec[i] * sin(phaseSpec[i]);
  }
  for( int i = 1; i < myNumFFTSamplesIn/2; i++ ) {
    myBufferReal[myNumFFTSamplesIn-i] = myBufferReal[i];
    myBufferImag[myNumFFTSamplesIn-i] = -myBufferImag[i];
  }

  // Phase at 0Hz and numFFTSamplesNyquist is always set to zero. If not, time series will have imaginary part.
  myBufferImag[0] = 0.0;
  myBufferImag[myNumFFTSamplesIn/2] = 0.0;
}
//--------------------------------------------------------------------------------
//
bool csFFTTools_ORIG::fft_inverse() {
  return csFFTTools_ORIG::fft( csFFTTools_ORIG::INVERSE, myTwoPowerIn, myBufferReal, myBufferImag );
}
bool csFFTTools_ORIG::fft_inverse( float const* samples, int fftDataType ) {
  if( fftDataType == FX_REAL_IMAG ) {
    for( int i = 0; i < myNumFFTSamplesIn; i++ ) {
      myBufferReal[i] = samples[i];
      myBufferImag[i] = samples[i+myNumFFTSamplesIn];
    }
  }
  else if( fftDataType == FX_AMP_PHASE ) {
    // Input float array contains both amplitude and phase spectrum
    convertFromAmpPhase( &samples[0], &samples[myNumFFTSamplesIn/2+1] );
  }
  else {
    // ...just use currently stored real/imag buffers for inverse FFT
//    return false;
  }
  return csFFTTools_ORIG::fft( csFFTTools_ORIG::INVERSE, myTwoPowerIn, myBufferReal, myBufferImag );
}

bool csFFTTools_ORIG::fft_inverse( float const* ampSpec, float const* phaseSpec, bool doNormalize ) {
  convertFromAmpPhase( ampSpec, phaseSpec );
  return csFFTTools_ORIG::fft( csFFTTools_ORIG::INVERSE, myTwoPowerIn, myBufferReal, myBufferImag, doNormalize );
}
//--------------------------------------------------------------------------------
//
//
void csFFTTools_ORIG::lowPass( float* samples, float cutOffFreqHz, float order, bool outputImpulseResponse ) {
  myCutOffFreqHz = cutOffFreqHz;
  myOrder = order;
  myOutputImpulseResponse = outputImpulseResponse;
  filter( samples, LOWPASS );
}
//--------------------------------------------------------------------------------
//
//
void csFFTTools_ORIG::highPass( float* samples, float cutOffFreqHz, float order, bool outputImpulseResponse ) {
  myCutOffFreqHz = cutOffFreqHz;
  myOrder = order;
  myOutputImpulseResponse = outputImpulseResponse;
  filter( samples, HIGHPASS );
}
//--------------------------------------------------------------------------------
//
//
void csFFTTools_ORIG::applyFilter( float* samples, csTimeFunction<double> const* timeFunc ) {
  setBuffer( samples );

  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::FORWARD, myTwoPowerIn, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    throw( csException("csFFTTools_ORIG::filter(): Unknown error occurred during forward FFT transform.") );
  }
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );
  for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
    double freq = (double)is * df;
    double scalar = timeFunc->valueAt( freq );
    myBufferReal[is] *= scalar;
    myBufferImag[is] *= scalar;
  }
  for( int is = myNumFFTSamplesIn/2+1; is < myNumFFTSamplesIn; is++ ) {
    double freq = -(double)(is-myNumFFTSamplesIn)*df;
    double scalar = timeFunc->valueAt( freq );
    myBufferReal[is] *= scalar;
    myBufferImag[is] *= scalar;
  }

  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::INVERSE, myTwoPowerOut, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    myBufferReal = NULL;
    myBufferImag = NULL;
    throw( csException("filter: Unknown error occurred during inverse FFT transform.") );
  }
  for( int i = 0; i < myNumSamplesOut; i++ ) {
    samples[i] = (float)myBufferReal[i];
  }
}
//--------------------------------------------------------------------------------
// Apply cosine taper around notch frequency
//
void csFFTTools_ORIG::notchFilter( float* samples, bool addNoise ) {

  setBuffer( samples );
  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::FORWARD, myTwoPowerIn, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    throw( csException("csFFTTools_ORIG::notchFilter(): Unknown error occurred during forward FFT transform.") );
  }

  for( int is = 0; is < myNumFFTSamplesIn; is++ ) {
    myBufferReal[is] *= myNotchFilter[is];
    myBufferImag[is] *= myNotchFilter[is];
  }

  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::INVERSE, myTwoPowerOut, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    myBufferReal = NULL;
    myBufferImag = NULL;
    throw( csException("csFFTTools_ORIG::notchFilter(): Unknown error occurred during inverse FFT transform.") );
  }
  for( int i = 0; i < myNumSamplesOut; i++ ) {
    samples[i] = (float)myBufferReal[i];
  }
}

float csFFTTools_ORIG::sampleIntFreqHz() const {
  return 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );
}

//    vars->fftTool->prepareMultipleFilters( vars->numLowPassFilters, vars->freqLowPass, vars->orderLowPass, vars->numHighPassFilters, vars->freqHighPass, vars->orderHighPass );

double const* csFFTTools_ORIG::setupNotchFilter( float notchFreqHz, float notchWidthHz, float slope, bool isCosineTaper ) {
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );
  double freq1 = notchFreqHz - 0.5*notchWidthHz;
  double freq2 = notchFreqHz + 0.5*notchWidthHz;
  int indexFirst = (int)round( freq1 / df );
  int indexLast  = (int)round( freq2 / df );
  int width = indexLast - indexFirst;
  if( width == 0 ) {
    indexFirst = std::max(0,indexFirst-1);
    indexLast  = std::min(indexFirst + 2,myNumFFTSamplesIn);
    width = indexLast - indexFirst;
  }
  int indexFirstRed = std::max(0,indexFirst);
  int indexLastRed  = std::min(myNumFFTSamplesIn,indexLast);

  myNotchFilter = new double[myNumFFTSamplesIn];
  for( int is = 0; is < myNumFFTSamplesIn; is++ ) {
    myNotchFilter[is] = 0.0;
  }

  if( isCosineTaper ) {
    for( int is = 0; is < myNumFFTSamplesIn; is++ ) {
      myNotchFilter[is] = 1.0;
    }
    for( int isamp = indexFirstRed; isamp <= indexLastRed; isamp++ ) {
      double phase  = 2.0*( ( (double)(isamp-indexFirst) / (double)width ) - 1.0 ) * M_PI;
      double scalar = 0.5 * (cos(phase) + 1.0);
      myNotchFilter[isamp] = scalar;
      myNotchFilter[myNumFFTSamplesIn-isamp] = scalar;
    }
  }
  else { // Apply Butterworth notch filter
    double G0 = 1.0;
    double power = fabs(slope) / 3.0f;
    double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );

    for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
      double freq = (double)is * df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/freq1,power) ));
      myNotchFilter[is] += dampG;
    }
    for( int is = myNumFFTSamplesIn/2+1; is < myNumFFTSamplesIn; is++ ) {
      double freq = -(double)(is-myNumFFTSamplesIn)*df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/freq1,power) ));
      myNotchFilter[is] += dampG;
    }
    // HIGHPASS
    for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
      double freq = (double)is * df;
      double dampG;
      if( freq == 0 ) {
        dampG = 0;
      }
      else {
        dampG = sqrt(G0 / (1.0 + pow(freq2/freq,power) ));
      }
      myNotchFilter[is] += dampG;
    }
    for( int is = myNumFFTSamplesIn/2+1; is < myNumFFTSamplesIn; is++ ) {
      double freq = -(double)(is-myNumFFTSamplesIn)*df;
      double dampG;
      if( freq == 0 ) {
        dampG = 0.0;
      }
      else {
        dampG = sqrt(G0 / (1.0 + pow(freq2/freq,power) ));
      }
      myNotchFilter[is] += dampG;
    }
  }
  
  // Dampen potential side lobes of filter close to notch (must be done for positive and negative frequency portion of complex spectrum)
  double maxValue1 = 0;
  int maxIndex1 = 0;
  int indexNotch  = (int)round( notchFreqHz / df );
  for( int is = 0; is < indexNotch; is++ ) {
    if( myNotchFilter[is] > maxValue1 ) {
      maxValue1 = myNotchFilter[is];
      maxIndex1 = is;
    }
  }
  if( maxIndex1 != 0 ) {
    int maxIndex2 = 0;
    double maxValue2 = 0;
    for( int is = indexNotch; is < 2*indexNotch; is++ ) {
      if( myNotchFilter[is] > maxValue2 ) {
        maxValue2 = myNotchFilter[is];
        maxIndex2 = is;
      }
    }
    double correctionTerm = 1.0 - 0.5*( maxValue1 + maxValue2 );
    for( int is = 0; is < maxIndex1; is++ ) {
      if( myNotchFilter[is] > 1.0 ) myNotchFilter[is] = 1.0;
    }
    for( int is = maxIndex1; is <= maxIndex2; is++ ) {
      myNotchFilter[is] += correctionTerm;
    }
    double maxIndex3 = myNumFFTSamplesIn - maxIndex2;
    double maxIndex4 = myNumFFTSamplesIn - maxIndex1;
    for( int is = maxIndex2+1; is < maxIndex3; is++ ) {
      if( myNotchFilter[is] > 1.0 ) myNotchFilter[is] = 1.0;
    }
    for( int is = maxIndex3; is <= maxIndex4; is++ ) {
      myNotchFilter[is] += correctionTerm;
    }
    for( int is = maxIndex4+1; is < myNumFFTSamplesIn; is++ ) {
      if( myNotchFilter[is] > 1.0 ) myNotchFilter[is] = 1.0;
    }
  }
  
  return myNotchFilter;
}
//--------------------------------------------------------------------------------
void csFFTTools_ORIG::envelope( float* samples ) {
  // 1) Hilbert transform == Shift phase by 90deg
  hilbertTransform( samples );

  // 2) Combine input trace and 90deg shifted trace to form envelope
  for( int isamp = 0; isamp < myNumSamplesIn; isamp++ ) {
    float sampleInput = samples[isamp];
    float sampleHilb  = (float)myBufferReal[isamp];
    samples[isamp] = sqrt( sampleInput*sampleInput + sampleHilb*sampleHilb );
  } 
}
//--------------------------------------------------------------------------------
// Hilbert transform == Shift phase by 90deg
//
void csFFTTools_ORIG::hilbertTransform( float* samples ) {
  fft_forward( samples );
  for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
    float amp   = (float)( sqrt(myBufferReal[i]*myBufferReal[i] + myBufferImag[i]*myBufferImag[i]) );
    float phase = (float)atan2( -myBufferImag[i], myBufferReal[i] );
    phase -= M_PI/2; // Hilbert transform == quadrature phase == 90deg phase shift

    myBufferReal[i] = amp * cos(phase);
    myBufferImag[i] = -amp * sin(phase);
  }
  for( int i = 1; i < myNumFFTSamplesIn/2; i++ ) {
    myBufferReal[myNumFFTSamplesIn-i] = myBufferReal[i];
    myBufferImag[myNumFFTSamplesIn-i] = -myBufferImag[i];
  }
  myBufferImag[0] = 0.0;
  myBufferImag[myNumFFTSamplesIn/2] = 0.0;

  fft_inverse();
}
//--------------------------------------------------------------------------------
void csFFTTools_ORIG::inst_phase( float* samples ) {
  hilbertTransform( samples );

  // Compute instantaneous phase
  // analytic_signal(x) = signal(x) + i * hilbert(x)
  for( int isamp = 0; isamp < myNumSamplesIn; isamp++ ) {
    float phase_inst = (float) atan2( -myBufferReal[isamp], (double)samples[isamp] );
    samples[isamp] = phase_inst;
  }
}
//--------------------------------------------------------------------------------
void csFFTTools_ORIG::inst_freq( float* samples ) {
  inst_phase( samples );

  // Compute instantaneous frequency
  float phaseThreshold = M_PI - 1e-5;
  float phasePrev = 0;
  for( int isamp = 0; isamp < myNumSamplesIn; isamp++ ) {
    float phase = samples[isamp];
    float phaseDiff = phase - phasePrev;
    if( phaseDiff > phaseThreshold ) phaseDiff -= 2.0 * M_PI;
    else if( phaseDiff < -phaseThreshold ) phaseDiff += 2.0 * M_PI;
    samples[isamp] = fabs(phaseDiff) / ( 2 * M_PI ) * 1000 / mySampleIntIn;
    phasePrev = phase;
  }
}
float csFFTTools_ORIG::resample( float* samples ) {
  return resample( samples, false, false );
}
float csFFTTools_ORIG::resample( float* samples, bool applyFilter, bool applyNorm ) {
  return resample( samples, myCutOffFreqHz, myOrder, applyFilter, applyNorm );
}
float csFFTTools_ORIG::resample( float* samples, float cutOffFreqHz, float order, bool applyFilter, bool applyNorm ) {
  myCutOffFreqHz = cutOffFreqHz;
  myOrder = order;
  myOutputImpulseResponse = false;

  //  filter( samples, RESAMPLE );

  setBuffer( samples );

  float rmsIn  = 1.0;
  float rmsOut = 1.0;
  if( applyNorm ) rmsIn = compute_rms( samples, myNumSamplesIn );

  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::FORWARD, myTwoPowerIn, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    throw( csException("csFFTTools_ORIG::filter(): Unknown error occurred during forward FFT transform.") );
  }

  double G0 = 1.0;
  double power = fabs(myOrder) * 2.0;
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );

  myBufferReal[myNumFFTSamplesOut/2] = 0.0; //myBufferReal[numFFTSamples/2];
  myBufferImag[myNumFFTSamplesOut/2] = 0.0; //myBufferImag[numFFTSamples/2];

  if( applyFilter ) {
    for( int is = 0; is < myNumFFTSamplesOut/2; is++ ) {
      double freq = (double)is * df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/myCutOffFreqHz,power) ));
      myBufferReal[is] *= dampG;
      myBufferImag[is] *= dampG;
    }
    for( int is = 1; is < myNumFFTSamplesOut/2; is++ ) {
      int indexOut  = myNumFFTSamplesOut/2 + is;
      int indexIn = 3*myNumFFTSamplesIn/4 + is;
      double freq = (double)fabs((double)(indexIn-myNumFFTSamplesIn))*df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/myCutOffFreqHz,power) ));
      myBufferReal[indexOut] = myBufferReal[indexIn] * dampG;
      myBufferImag[indexOut] = myBufferImag[indexIn] * dampG;
    }
  }
  else {
    for( int is = 1; is < myNumFFTSamplesOut/2; is++ ) {
      int indexOut  = myNumFFTSamplesOut/2 + is;
      int indexIn = 3*myNumFFTSamplesIn/4 + is;
      myBufferReal[indexOut] = myBufferReal[indexIn];
      myBufferImag[indexOut] = myBufferImag[indexIn];
    }
  }

  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::INVERSE, myTwoPowerOut, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    throw( csException("filter: Unknown error occurred during inverse FFT transform.") );
  }
  for( int i = 0; i < myNumSamplesOut; i++ ) {
    samples[i] = (float)myBufferReal[i];
  }

  if( applyNorm ) {
    compute_rms( samples, myNumSamplesOut );
    if( rmsOut != 0.0 ) {
      float ratio = rmsIn / rmsOut;
      for( int i = 0; i < myNumSamplesOut; i++ ) {
        samples[i] *= ratio;
      }
      return ratio;
    }
  }
  //  fprintf(stderr,"Num samples in/out: %d %d, order: %d %f\n", myNumSamplesIn, myNumSamplesOut, myOrder, myCutOffFreqHz );
  return 1.0;
}
//--------------------------------------------------------------------------------
//
//
void csFFTTools_ORIG::prepareMultipleFilters( int numLowPassFilters, float const* lowPassFreq, float const* lowPassOrder, int numHighPassFilters, float const* highPassFreq, float const* highPassOrder ) {
  double G0 = 1.0;
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );

  if( myFilterScalars != NULL ) delete [] myFilterScalars;
  myFilterScalars = new float[myNumFFTSamplesIn];

  for( int is = 0; is < myNumFFTSamplesIn; is++ ) {
    myFilterScalars[is] = 1.0f;
  }
  for( int i = 0; i < numLowPassFilters; i++ ) {
    float order      = lowPassOrder[i];
    double power     = fabs(order) * 2.0;
    float cutOffFreq = lowPassFreq[i];
    for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
      double freq = (double)is * df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/cutOffFreq,power) ));
      if( order <= 0 && dampG == 0.0 ) dampG = 1.0;
      myFilterScalars[is] = dampG;
    }
    for( int is = myNumFFTSamplesIn/2+1; is < myNumFFTSamplesIn; is++ ) {
      double freq = -(double)(is-myNumFFTSamplesIn)*df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/cutOffFreq,power) ));
      if( order <= 0 && dampG == 0.0 ) dampG = 1.0;
      myFilterScalars[is] = dampG;
    }
  }
  for( int i = 0; i < numHighPassFilters; i++ ) {
    float order      = highPassOrder[i];
    double power     = fabs(order) * 2.0;
    float cutOffFreq = highPassFreq[i];
    for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
      double freq = (double)is * df;
      double dampG = sqrt(G0 / (1.0 + pow(cutOffFreq/freq,power) ));
      if( freq == 0 ) {
        if( order > 0 ) {
          dampG = 0;
        }
        else {
          dampG = 1.0;
        }
      }
      else if( order <= 0 ) {
        if( dampG == 0.0 ) dampG = 1.0;
      }
      myFilterScalars[is] *= dampG;
    }
    for( int is = myNumFFTSamplesIn/2+1; is < myNumFFTSamplesIn; is++ ) {
      double freq = -(double)(is-myNumFFTSamplesIn)*df;
      double dampG = sqrt(G0 / (1.0 + pow(cutOffFreq/freq,power) ));
      if( freq == 0 ) {
        if( order > 0 ) {
          dampG = 0;
        }
        else {
          dampG = 1.0;
        }
      }
      else if( order <= 0 ) {
        if( dampG == 0.0 ) dampG = 1.0;
      }
      myFilterScalars[is] *= dampG;
    }
  }
}

void csFFTTools_ORIG::applyLowHighPassFilters( float* samples, bool outputImpulseResponse ) {
  if( myFilterScalars == NULL ) throw( csException("csFFTTools_ORIG::applyLowHighPassFilters: Filters have not been initialized. Call prepareMultipleFilters() first") );

  myOutputImpulseResponse = outputImpulseResponse;
  setBuffer( samples );
  
  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::FORWARD, myTwoPowerIn, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    throw( csException("csFFTTools_ORIG::applyLowHighPassFilter(): Unknown error occurred during forward FFT transform.") );
  }

  for( int is = 0; is < myNumFFTSamplesIn; is++ ) {
    myBufferReal[is] *= myFilterScalars[is];
    myBufferImag[is] *= myFilterScalars[is];
  }

  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::INVERSE, myTwoPowerOut, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    myBufferReal = NULL;
    myBufferImag = NULL;
    throw( csException("applyLowHighPassFilters: Unknown error occurred during inverse FFT transform.") );
  }
  for( int i = 0; i < myNumSamplesOut; i++ ) {
    samples[i] = (float)myBufferReal[i];
  }
}
//--------------------------------------------------------------------------------
//
void csFFTTools_ORIG::applyLowPassFilter( float cutOffFreq, float order ) {
  double G0 = 1.0;
  double power = fabs(order) * 2.0;
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );

  if( order > 0 ) {
    for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
      double freq = (double)is * df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/cutOffFreq,power) ));
      myBufferReal[is] *= dampG;
      myBufferImag[is] *= dampG;
    }
    for( int is = myNumFFTSamplesIn/2+1; is < myNumFFTSamplesIn; is++ ) {
      double freq = -(double)(is-myNumFFTSamplesIn)*df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/cutOffFreq,power) ));
      myBufferReal[is] *= dampG;
      myBufferImag[is] *= dampG;
    }
  }
  else { //if( order < 0 ) {
    for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
      double freq = (double)is * df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/cutOffFreq,power) ));
      if( dampG == 0.0 ) dampG = 1.0;
      myBufferReal[is] /= dampG;
      myBufferImag[is] /= dampG;
    }
    for( int is = myNumFFTSamplesIn/2+1; is < myNumFFTSamplesIn; is++ ) {
      double freq = -(double)(is-myNumFFTSamplesIn)*df;
      double dampG = sqrt(G0 / (1.0 + pow(freq/cutOffFreq,power) ));
      if( dampG == 0.0 ) dampG = 1.0;
      myBufferReal[is] /= dampG;
      myBufferImag[is] /= dampG;
    }
  }
}

void csFFTTools_ORIG::applyHighPassFilter( float cutOffFreq, float order ) {
  double G0 = 1.0;
  double power = fabs(order) * 2.0;
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );

  if( order > 0 ) {
    for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
      double freq = (double)is * df;
      double dampG;
      if( freq == 0 ) {
        dampG = 0;
      }
      else {
        dampG = sqrt(G0 / (1.0 + pow(cutOffFreq/freq,power) ));
      }
      myBufferReal[is] *= dampG;
      myBufferImag[is] *= dampG;
    }
    for( int is = myNumFFTSamplesIn/2+1; is < myNumFFTSamplesIn; is++ ) {
      double freq = -(double)(is-myNumFFTSamplesIn)*df;
      double dampG;
      if( freq == 0 ) {
        dampG = 0.0;
      }
      else {
        dampG = sqrt(G0 / (1.0 + pow(cutOffFreq/freq,power) ));
      }
      myBufferReal[is] *= dampG;
      myBufferImag[is] *= dampG;
    }
  }
  else { //if( order < 0 ) {
    for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
      double freq = (double)is * df;
      double dampG;
      if( freq == 0 ) {
        dampG = 1.0;
      }
      else {
        dampG = sqrt(G0 / (1.0 + pow(cutOffFreq/freq,power) ));
        if( dampG == 0.0 ) dampG = 1.0;
      }
      myBufferReal[is] /= dampG;
      myBufferImag[is] /= dampG;
    }
    for( int is = myNumFFTSamplesIn/2+1; is < myNumFFTSamplesIn; is++ ) {
      double freq = -(double)(is-myNumFFTSamplesIn)*df;
      double dampG;
      if( freq == 0 ) {
        dampG = 1.0;
      }
      else {
        dampG = sqrt(G0 / (1.0 + pow(cutOffFreq/freq,power) ));
        if( dampG == 0.0 ) dampG = 1.0;
      }
      myBufferReal[is] /= dampG;
      myBufferImag[is] /= dampG;
    }
  }
}

void csFFTTools_ORIG::filter( float* samples, int filterType ) {
  setBuffer( samples );

  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::FORWARD, myTwoPowerIn, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    throw( csException("csFFTTools_ORIG::filter(): Unknown error occurred during forward FFT transform.") );
  }

  if( filterType == LOWPASS ) {
    applyLowPassFilter( myCutOffFreqHz, myOrder );
  }
  else if( filterType == HIGHPASS ) {
    applyHighPassFilter( myCutOffFreqHz, myOrder );
  }
  else {
    // ?
  }

  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::INVERSE, myTwoPowerOut, myBufferReal, myBufferImag ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    myBufferReal = NULL;
    myBufferImag = NULL;
    throw( csException("filter: Unknown error occurred during inverse FFT transform.") );
  }
  for( int i = 0; i < myNumSamplesOut; i++ ) {
    samples[i] = (float)myBufferReal[i];
  }
}
//--------------------------------------------------------------------------------
//
//
/*void csFFTTools_ORIG::applyQCompensation( float* samples, float qvalue, float freqRef, bool applyAmp, bool applyPhase ) {
  setBuffer( samples );

  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::FORWARD, myTwoPowerIn, myBufferReal, myBufferImag, false ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    throw( csException("csFFTTools_ORIG::applyQCompensation(): Unknown error occurred during forward FFT transform.") );
  }


  if( !csFFTTools_ORIG::fft( csFFTTools_ORIG::INVERSE, myTwoPowerOut, myBufferReal, myBufferImag, true ) ) {
    delete [] myBufferReal;
    delete [] myBufferImag;
    throw( csException("filter: Unknown error occurred during inverse FFT transform.") );
  }
  for( int i = 0; i < myNumSamplesOut; i++ ) {
    samples[i] = (float)myBufferReal[i];
  }

}
*/
//--------------------------------------------------------------------------------
// Set FFT coefficients
//
void csFFTTools_ORIG::setBuffer( float const* samples ) {
  if( !myOutputImpulseResponse ) {
    for( int i = 0; i < myNumSamplesIn; i++ ) {
      myBufferReal[i] = samples[i];
      myBufferImag[i] = 0;
    }
    for( int i = myNumSamplesIn; i < myNumFFTSamplesIn; i++ ) {
      myBufferReal[i] = 0;
      myBufferImag[i] = 0;
    }
  }
  //----------------------------------------------------------
  // Create filter impulse response
  else {
    double amplitude = 1.0;
    for( int i = 0; i < myNumFFTSamplesIn; i++ ) {
      myBufferReal[i] = 0;
      myBufferImag[i] = 0;
    }
    myBufferReal[myNumFFTSamplesIn/2] = amplitude;
  }
}


bool csFFTTools_ORIG::Powerof2( int numFFTSamplesX, int* m, int* twopm ) {
  int value = numFFTSamplesX;
  *m = 0;
  while( (value = (int)(value / 2)) > 0 ) {
    *m += 1;
  }
  *twopm = (int)pow( 2.0, *m );
  return( *twopm == numFFTSamplesX );
}


/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
 * This fft implementation is pulled somewhere from the Internet, I forgot where...
 */
bool csFFTTools_ORIG::fft( int dir, int power_of_two, double *realValues, double *imagValues ) {
  bool doNormalize = ( dir == csFFTTools_ORIG::INVERSE );
  return csFFTTools_ORIG::fft( dir, power_of_two, realValues, imagValues, doNormalize );
}
bool csFFTTools_ORIG::fft( int dir, int power_of_two, double *realValues, double *imagValues, bool doNormalize )
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<power_of_two;i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = realValues[i];
         ty = imagValues[i];
         realValues[i] = realValues[j];
         imagValues[i] = imagValues[j];
         realValues[j] = tx;
         imagValues[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<power_of_two;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * realValues[i1] - u2 * imagValues[i1];
            t2 = u1 * imagValues[i1] + u2 * realValues[i1];
            realValues[i1] = realValues[i] - t1;
            imagValues[i1] = imagValues[i] - t2;
            realValues[i] += t1;
            imagValues[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   // Normalisation should be done only for inverse transform
   if( doNormalize ) {
     for( i = 0; i < nn; i++ ) {
       realValues[i] /= (double)nn;
       imagValues[i] /= (double)nn;
     }
   }

   return(true);
}

// Not tested yet...
bool csFFTTools_ORIG::fft_2d( int dir, double** realValues, double** imagValues, int numFFTSamplesX, int numFFTSamplesY )
{
  int mx,my,twopm;

  // Transform the rows
  int max_numFFTSamples = std::max(numFFTSamplesX,numFFTSamplesY);
  double* real = new double[max_numFFTSamples];
  double* imag = new double[max_numFFTSamples];

  if( real == NULL || imag == NULL ) return false;
  if( !csFFTTools_ORIG::Powerof2( numFFTSamplesX, &mx, &twopm ) || twopm != numFFTSamplesX ) {
    return false;
  }
  if( !csFFTTools_ORIG::Powerof2( numFFTSamplesY, &my, &twopm ) || twopm != numFFTSamplesY ) {
    return false;
  }
  for( int j = 0; j < numFFTSamplesY; j++ ) {
    for( int i = 0; i < numFFTSamplesX; i++ ) {
      real[i] = realValues[i][j];
      imag[i] = imagValues[i][j];
    }
    csFFTTools_ORIG::fft(dir,mx,real,imag);
    for( int i = 0; i < numFFTSamplesX; i++ ) {
      realValues[i][j] = real[i];
      imagValues[i][j] = imag[i];
    }
  }

   // Transform the columns
  for( int i = 0; i < numFFTSamplesX; i++) {
//    memcpy( real, &realValues[i][0], numFFTSamplesY*sizeof(double) );
//    memcpy( imag, &imagValues[i][0], numFFTSamplesY*sizeof(double) );
    for( int j = 0; j < numFFTSamplesY; j++ ) {
      real[j] = realValues[i][j];
      imag[j] = imagValues[i][j];
    }
    csFFTTools_ORIG::fft(dir,my,real,imag);
//    memcpy( &realValues[i][0], real, numFFTSamplesY*sizeof(double) );
//    memcpy( &imagValues[i][0], imag, numFFTSamplesY*sizeof(double) );
    for( int j = 0; j < numFFTSamplesY; j++ ) {
      realValues[i][j] = real[j];
      imagValues[i][j] = imag[j];
    }
  }
  delete [] real;
  delete [] imag;
  
  return true;
}

void csFFTTools_ORIG::applyTaper( int taperType, int taperLengthInSamples, int numSamplesIn, float* samples ) {
  if( taperType == csFFTTools_ORIG::TAPER_COSINE || taperType == csFFTTools_ORIG::TAPER_HANNING ) {
    for( int i = 0; i < taperLengthInSamples; i++ ) {
      float scalar = cos( M_PI_2 * (float)(taperLengthInSamples-i)/(float)taperLengthInSamples );
      samples[i] *= scalar;
    }
    for( int i = numSamplesIn-taperLengthInSamples; i < numSamplesIn; i++ ) {
      float scalar = cos( M_PI_2 * (float)(taperLengthInSamples-numSamplesIn+i+1)/(float)taperLengthInSamples );
      samples[i] *= scalar;
    }
  }
  else if( taperType == csFFTTools_ORIG::TAPER_BLACKMAN ) {
    float alpha = 0.16;
    float a0 = 0.5 * (1.0 - alpha);
    float a1 = 0.5;
    float a2 = 0.5 * alpha;
    for( int i = 0; i < numSamplesIn; i++ ) {
      float piFactor = (2.0 * M_PI) * (float)i / (float)(numSamplesIn - 1);
      float weight = a0 - a1*cos( piFactor ) + a2*cos( 2 * piFactor );
      samples[i] *= weight;
    }
  }
}    
