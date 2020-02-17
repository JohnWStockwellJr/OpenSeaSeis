/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csFFT.h"
#include <cmath>
#include <limits>
#include <cstring>
#include <cstdio>
#include "geolib_defines.h"
#include "csException.h"

using namespace cseis_geolib;
using namespace std;

csFFT::csFFT( int numSamples_time, int numSamples_fft ) {
#ifdef USE_FFTW
  // Use FFTW library:
  myFFTW_floatBuffer   = NULL;
  myFFTW_complexBuffer = NULL;
  myFFTW_plan_forward  = NULL;
  myFFTW_plan_inverse  = NULL;
#else
  // Non-FFTW implementation:
  myBufferReal = NULL;
  myBufferImag = NULL;
#endif

  init( numSamples_time, numSamples_fft );
}
csFFT::~csFFT() {
#ifdef USE_FFTW
  // Use FFTW library:
  if( myFFTW_floatBuffer != NULL ) {
    fftwf_free( myFFTW_floatBuffer );
    myFFTW_floatBuffer = NULL;
  }
  if( myFFTW_complexBuffer != NULL ) {
    fftwf_free( myFFTW_complexBuffer );
    myFFTW_complexBuffer = NULL;
  }
  if( myFFTW_plan_forward != NULL ) {
    fftwf_destroy_plan(myFFTW_plan_forward);
    myFFTW_plan_forward = NULL;
  }
  if( myFFTW_plan_inverse != NULL ) {
    fftwf_destroy_plan(myFFTW_plan_inverse);
    myFFTW_plan_inverse = NULL;
  }
#else
  // Non-FFTW implementation:
  if( myBufferReal != NULL ) {
    delete [] myBufferReal;
    myBufferReal = NULL;
  }
  if( myBufferImag != NULL ) {
    delete [] myBufferImag;
    myBufferImag = NULL;
  }
#endif
}
void csFFT::init( int numSamples_time, int numSamples_fft ) {
  myNumTimeSamples    = numSamples_time;

#ifdef USE_FFTW
  // Use FFTW library:
  if( numSamples_fft <= 0 ) {
    myNumFFTSamples = csFFT::factor_2357( myNumTimeSamples ); 
    myNumFFTSamples = (int)( ( myNumFFTSamples + 1 ) / 2 ) * 2;
    if ( myNumFFTSamples <= 0  ) {
      //    throw( cseis_geolib::csException("Failed optimizing FFTW length %d",myNumTimeSamples) );
      myNumFFTSamples = myNumTimeSamples;
    }
  }
  else {
    myNumFFTSamples = numSamples_fft; // Override number of FFT samples
  }
  myNumFFTFreq    = myNumFFTSamples / 2 + 1;
  myFFTW_floatBuffer   = (float*)fftwf_malloc( sizeof(float) * myNumFFTSamples );
  myFFTW_complexBuffer = (fftwf_complex*)fftwf_malloc( sizeof(fftwf_complex) * myNumFFTFreq );
  memset( myFFTW_floatBuffer, 0, myNumFFTSamples * sizeof(float) );
  memset( myFFTW_complexBuffer, 0, sizeof(fftwf_complex) * myNumFFTFreq );
  myFFTW_plan_forward = fftwf_plan_dft_r2c_1d( myNumFFTSamples, myFFTW_floatBuffer, myFFTW_complexBuffer, FFTW_MEASURE );
  myFFTW_plan_inverse = fftwf_plan_dft_c2r_1d( myNumFFTSamples, myFFTW_complexBuffer, myFFTW_floatBuffer, FFTW_MEASURE );
#else
  // Non-FFTW implementation:
  myNumFFTSamples = myNumTimeSamples;
  int two_power_m;
  csFFT::powerof2( myNumFFTSamples, &myTwoPower, &two_power_m );
  if( two_power_m != myNumFFTSamples ) {
    myNumFFTSamples = two_power_m * 2;
    myTwoPower += 1;
  }
  myBufferReal = new double[myNumFFTSamples];
  myBufferImag = new double[myNumFFTSamples];
  myNumFFTFreq = myNumFFTSamples / 2 + 1;
#endif
}

int csFFT::numFreqValues() const {
  return myNumFFTFreq;
}
int csFFT::numFFTValues() const {
  return myNumFFTSamples;
}
double csFFT::sampleIntFreq( double sampleInt_ms ) const {
  double freqNyquist = 1.0 / ( 2.0 * ( sampleInt_ms / 1000.0 ) );
  return( ( freqNyquist / (double)(numFreqValues()-1) ) );
}
#ifdef USE_FFTW
void csFFT::buffer2RealImag( float* realValues, float* imagValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    realValues[ifreq] = (float)( myFFTW_complexBuffer[ifreq][0] );
    imagValues[ifreq] = (float)( myFFTW_complexBuffer[ifreq][1] );
  }
}
void csFFT::buffer2AmpPhase( float* ampValues, float* phaseValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float val1 = myFFTW_complexBuffer[ifreq][0];
    float val2 = myFFTW_complexBuffer[ifreq][1];
    ampValues[ifreq]   = sqrt( val1*val1 + val2*val2 );
    phaseValues[ifreq] = atan2( -val2, val1 );
  }
}
void csFFT::buffer2Amp( float* ampValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float val1 = myFFTW_complexBuffer[ifreq][0];
    float val2 = myFFTW_complexBuffer[ifreq][1];
    ampValues[ifreq]   = sqrt( val1*val1 + val2*val2 );
  }
}
void csFFT::buffer2PSD( float* psdValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float val1 = myFFTW_complexBuffer[ifreq][0];
    float val2 = myFFTW_complexBuffer[ifreq][1];
    psdValues[ifreq] = val1*val1 + val2*val2;
  }
}
void csFFT::buffer2Phase( float* phaseValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float val1 = myFFTW_complexBuffer[ifreq][0];
    float val2 = myFFTW_complexBuffer[ifreq][1];
    phaseValues[ifreq] = atan2( -val2, val1 );
  }
}
void csFFT::ampPhase2Buffer( float const* ampValues, float const* phaseValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float amplitude = ampValues[ifreq];
    float phase     = phaseValues[ifreq];
    myFFTW_complexBuffer[ifreq][0] =  amplitude * cos(phase);
    myFFTW_complexBuffer[ifreq][1] = -amplitude * sin(phase);
  }
}
void csFFT::realImag2Buffer( float const* realValues, float const* imagValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    myFFTW_complexBuffer[ifreq][0] = realValues[ifreq];
    myFFTW_complexBuffer[ifreq][1] = imagValues[ifreq];
  }
}

void csFFT::applyAmpFilter( float* samplesInOut, float const* filterCoefficients ) {
  forwardTransform( samplesInOut );
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    myFFTW_complexBuffer[ifreq][0] *= filterCoefficients[ifreq];
    myFFTW_complexBuffer[ifreq][1] *= filterCoefficients[ifreq];
  }
  inverseTransform( samplesInOut );
}
void csFFT::forwardTransform( float const* samplesIn ) {
  if( myNumFFTSamples > myNumTimeSamples ) {
    memset( &myFFTW_floatBuffer[myNumTimeSamples], 0, sizeof(float) * (myNumFFTSamples-myNumTimeSamples) );
  }
  memcpy( myFFTW_floatBuffer, samplesIn, myNumTimeSamples * sizeof(float) );
  //  memset( myFFTW_complexBuffer, 0, sizeof(fftwf_complex) * myNumFFTFreq );
  fftwf_execute( myFFTW_plan_forward );
}
void csFFT::inverseTransform( float* samplesOut, bool doNormalize ) {
  //  memset( myFFTW_floatBuffer, 0, sizeof(float) * myNumFFTSamples );
  fftwf_execute( myFFTW_plan_inverse );

  float scalar = !doNormalize ? 1.0f : ( 1.0f / (float)numFFTValues() );
  for( int isamp = 0; isamp < myNumTimeSamples; isamp++ ) {
    samplesOut[isamp] = scalar * myFFTW_floatBuffer[isamp];
  }
}

#else
//************************************************************************************************************************************
// Non-FFTW implementation:
//
//--------------------------------------------------------------------------------
// Set FFT coefficients
//
void csFFT::setBuffer( float const* samples ) {
  for( int i = 0; i < myNumTimeSamples; i++ ) {
    myBufferReal[i] = samples[i];
    myBufferImag[i] = 0;
  }
  for( int i = myNumTimeSamples; i < myNumFFTSamples; i++ ) {
    myBufferReal[i] = 0;
    myBufferImag[i] = 0;
  }
  //----------------------------------------------------------
  // Create filter impulse response
  //  else {
  //   for( int i = 0; i < myNumFFTSamples; i++ ) {
  //    myBufferReal[i] = 0;
  //     myBufferImag[i] = 0;
  //  }
  //  myBufferReal[myNumTimeSamples/2] = 1.0;
  // }
}
void csFFT::buffer2RealImag( float* realValues, float* imagValues ) {
  // Only first numFreq values in Real & Imag buffers are unique. Negative frequencies are redundant
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    realValues[ifreq] = (float)myBufferReal[ifreq];
    imagValues[ifreq] = (float)myBufferImag[ifreq];
  }
}
void csFFT::buffer2AmpPhase( float* ampValues, float* phaseValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float val1 = (float)myBufferReal[ifreq];
    float val2 = (float)myBufferImag[ifreq];
    ampValues[ifreq]   = sqrt( val1*val1 + val2*val2 );
    phaseValues[ifreq] = atan2( -val2, val1 );
  }
}
void csFFT::buffer2Amp( float* ampValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float val1 = (float)myBufferReal[ifreq];
    float val2 = (float)myBufferImag[ifreq];
    ampValues[ifreq]   = sqrt( val1*val1 + val2*val2 );
  }
}
void csFFT::buffer2PSD( float* psdValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float val1 = (float)myBufferReal[ifreq];
    float val2 = (float)myBufferImag[ifreq];
    psdValues[ifreq] = val1*val1 + val2*val2;
  }
}
void csFFT::buffer2Phase( float* phaseValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float val1 = (float)myBufferReal[ifreq];
    float val2 = (float)myBufferImag[ifreq];
    phaseValues[ifreq] = atan2( -val2, val1 );
  }
}
void csFFT::ampPhase2Buffer( float const* ampValues, float const* phaseValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    float amplitude = ampValues[ifreq];
    float phase     = phaseValues[ifreq];
    myBufferReal[ifreq] =  amplitude * cos(phase);
    myBufferImag[ifreq] = -amplitude * sin(phase);
  }
  for( int i = 1; i < myNumFFTFreq-1; i++ ) {
    myBufferReal[myNumFFTSamples-i] = myBufferReal[i];
    myBufferImag[myNumFFTSamples-i] = -myBufferImag[i];
  }
  // Phase at 0Hz and numFFTSamplesNyquist is always set to zero. If not, time series will have imaginary part.
  myBufferImag[0] = 0.0;
  myBufferImag[myNumFFTFreq-1] = 0.0;
}
void csFFT::realImag2Buffer( float const* realValues, float const* imagValues ) {
  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    myBufferReal[ifreq] = realValues[ifreq];
    myBufferImag[ifreq] = imagValues[ifreq];
  }
  // Negative frequencies: Redundant information
  for( int is = myNumFFTFreq; is < myNumFFTSamples; is++ ) {
    int index = myNumFFTSamples - is;
    myBufferReal[is] = realValues[index];
    myBufferImag[is] = -imagValues[index];  // Negative
  }
//  for( int is = 0; is < myNumFFTSamples; is++ ) {
  //  fprintf(stdout,"%d %.6e %.6e   CHECK\n", is, myBufferReal[is], myBufferImag[is]);
  // }
}

//--------------------------------------------------------------------------------
void csFFT::applyAmpFilter( float* samplesInOut, float const* filterCoefficients ) {
  forwardTransform( samplesInOut );

  for( int ifreq = 0; ifreq < myNumFFTFreq; ifreq++ ) {
    myBufferReal[ifreq] *= filterCoefficients[ifreq];
    myBufferImag[ifreq] *= filterCoefficients[ifreq];
  }
  for( int ifreq = myNumFFTFreq; ifreq < myNumFFTSamples; ifreq++ ) {
    myBufferReal[ifreq] *= filterCoefficients[myNumFFTSamples-ifreq];
    myBufferImag[ifreq] *= filterCoefficients[myNumFFTSamples-ifreq];
  }
  inverseTransform( samplesInOut );
}

void csFFT::forwardTransform( float const* samplesIn ) {
  setBuffer( samplesIn );
  csFFT::fft( csFFT::FORWARD, myTwoPower, myBufferReal, myBufferImag, false );
}

void csFFT::inverseTransform( float* samplesOut, bool doNormalize ) {
  csFFT::fft( csFFT::INVERSE, myTwoPower, myBufferReal, myBufferImag, true );

  for( int isamp = 0; isamp < myNumTimeSamples; isamp++ ) {
    samplesOut[isamp] = myBufferReal[isamp];
  }
}

#endif
//************************************************************************************************************************************
// Methods that work for both FFTW and non-FFTW implementations:
//

void csFFT::forwardTransform( float* samplesInOut, int fftDataType ) {
  forwardTransform( samplesInOut, samplesInOut, fftDataType );
}
void csFFT::forwardTransform( float const* samplesIn, float* samplesOut, int fftDataType ) {
  forwardTransform( samplesIn );

  if( samplesOut == NULL) return;

  // Convert to desired output
  if ( fftDataType == cseis_geolib::FX_AMP_PHASE ) {
    buffer2AmpPhase( &samplesOut[0], &samplesOut[myNumFFTFreq] );
  }
  else if ( fftDataType == cseis_geolib::FX_AMP ) {
    buffer2Amp( samplesOut );
  }
  else if ( fftDataType == cseis_geolib::FX_PHASE ) {
    buffer2Phase( samplesOut );
  }
  else if ( fftDataType == cseis_geolib::FX_PSD ) {
    buffer2PSD( samplesOut );
  }
  else { // if ( fftDataType == cseis_geolib::FX_REAL_IMAG ) {
    buffer2RealImag( &samplesOut[0], &samplesOut[myNumFFTFreq] );
  }
}

void csFFT::inverseTransform( float* samplesInOut, int fftDataType, bool doNormalize ) {
  inverseTransform( samplesInOut, samplesInOut, fftDataType, doNormalize );
}
void csFFT::inverseTransform( float const* samplesIn, float* samplesOut, int fftDataType, bool doNormalize ) {
  if ( fftDataType == cseis_geolib::FX_AMP_PHASE ) {
    ampPhase2Buffer( &samplesIn[0], &samplesIn[myNumFFTFreq] );
  }
  else if ( fftDataType == FX_REAL_IMAG ) {
    realImag2Buffer( &samplesIn[0], &samplesIn[myNumFFTFreq] );
  }
  else {
    throw( csException("csFFT::inverseTransform(): Input data has wrong FX data type. Must be FX_REAL_IMAG of FX_AMP_PHASE") );
  }

  inverseTransform( samplesOut, doNormalize );
}


////////////////////////////////////////////////////////////////////////////////
// Return the first multiple of small primes (2,3,5 and 7) greater than or 
// equal to the input value. 
//
// The return value is -1 if the value can't be reduced to a factor of small
// primes. That is, if it is less than 0 or too large (*).
//
// (*)NOTE: The largest 32-bit integer that is still a factor of small primes 
// appears to be 
//
//          INT_MAX - 3330622  = 2147483647 - 3330622
//                             = 2144153025
//
// Numbers larger than this are still valid integers up to INT_MAX but can't 
// be factored.
// 
// It is up to the calling program to decide what to do in that case. The return
// value will still be -1.
//
////////////////////////////////////////////////////////////////////////////////
int csFFT::factor_2357( int in ){
  int  test,  n2, n3, n5, n7, count; 
  int  out; 

  out = -1;

  if ( in < 0 || in > (std::numeric_limits<int>::max()-1) ) return out; // Fail if input bad  
  if ( in <=1 ){ out = 1; return out; } // Skip the loop for the easy one. 

  // Factor-out small primes from the input number. It should reduce to 1 once 
  // it is a perfect factor of 2,3,5 and 7. If not, increment the input and repeat.
  // Fail if the number of iterations or the test value gets to large. 
  test  = 0;
  count = 0;
  while ( test != 1 && count < std::numeric_limits<int>::max() ){
     test = in+count;
     if ( test == std::numeric_limits<int>::max() ) break;
     count++;     

     n2=n3=n5=n7=0;
     while ( test%7 == 0 ){ test = test/7; n7++; }
     while ( test%5 == 0 ){ test = test/5; n5++; }
     while ( test%3 == 0 ){ test = test/3; n3++; }
     while ( test%2 == 0 ){ test = test/2; n2++; }

  }
  if ( test != 1 ) return out;
  out = (int)round( pow(2.0,n2) * pow(3.0,n3) * pow(5.0,n5) * pow(7.0,n7) );

  return out;
}

bool csFFT::powerof2( int numFFTSamples, int* m, int* twopm ) {
  int value = numFFTSamples;
  *m = 0;
  while( (value = (int)(value / 2)) > 0 ) {
    *m += 1;
  }
  *twopm = (int)pow( 2.0, *m );
  return( *twopm == numFFTSamples );
}

bool csFFT::fft( int dir, int power_of_two, double *realValues, double *imagValues, bool doNormalize )
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
      if (dir == csFFT::FORWARD)
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
