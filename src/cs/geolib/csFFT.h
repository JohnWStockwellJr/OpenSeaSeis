/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_FFT_H
#define CS_FFT_H

#ifdef USE_FFTW
 extern "C" {
   #include <fftw3.h>
 }
#endif

namespace cseis_geolib {

/**
   Forward transform:
   - Compute single-sided amplitude spectra (and phase...)
   Inverse transform:
   - Restore original time trace
 */
class csFFT {
public:
  static int const FORWARD = 1;
  static int const INVERSE = 2;

public:
  csFFT( int numSamples_time, int numSamples_fft = -1 );
  ~csFFT();

  void init( int numSamples_time, int numSamples_fft );
  void forwardTransform( float* samplesInOut, int fftDataType );
  void forwardTransform( float const* samplesIn, float* samplesOut, int fftDataType );
  void inverseTransform( float* samplesInOut, int fftDataType, bool doNormalize = true );
  void inverseTransform( float const* samplesIn, float* samplesOut, int fftDataType, bool doNormalize = true );
  void applyAmpFilter( float* samplesInOut, float const* filterCoefficients );

  int numFreqValues() const;
  int numFFTValues() const;
  double sampleIntFreq( double sampleInt_ms ) const;

  static bool fft( int dir, int power_of_two, double *realValues, double *imagValues, bool doNormalize );
  static int factor_2357( int in );
  static bool powerof2( int nx, int* m, int* twopm );

 private:
  static int const INIT_FLAG_NONE    = -1;
  static int const INIT_FLAG_FORWARD = 1;
  static int const INIT_FLAG_INVERSE = 2;
  csFFT();

  void forwardTransform( float const* samplesIn );
  void inverseTransform( float* samplesOut, bool doNormalize = true );

  void buffer2RealImag( float* realValues, float* imagValues );
  void buffer2AmpPhase( float* ampValues, float* phaseValues );
  void buffer2Amp( float* ampValues );
  void buffer2PSD( float* psdValues );
  void buffer2Phase( float* phaseValues );

  void realImag2Buffer( float const* realValues, float const* imagValues );
  void ampPhase2Buffer( float const* ampValues, float const* phaseValues );

  int myNumTimeSamples;
  int myNumFFTSamples;
  int myNumFFTFreq;

#ifdef USE_FFTW
  // Use FFTW library:
  float*         myFFTW_floatBuffer;   // Float buffer 
  fftwf_complex* myFFTW_complexBuffer; // Complex buffer
  fftwf_plan     myFFTW_plan_forward;  // FFTW "plan" for forward transform
  fftwf_plan     myFFTW_plan_inverse;  // FFTW "plan" for inverse transform
#else
  // Non-FFTW implementation:
  void setBuffer( float const* samples );
  double* myBufferReal;
  double* myBufferImag;
  int myTwoPower;
#endif
};

} // namespace
#endif
