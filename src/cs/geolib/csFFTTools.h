/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_FFT_TOOLS_H
#define CS_FFT_TOOLS_H

namespace cseis_geolib {

template <typename T> class csTimeFunction;
class csFFT;

class csFFTTools {
public:
  static int const FORWARD = 1;
  static int const INVERSE = 2;

  static int const LOWPASS  = 50;
  static int const HIGHPASS = 51;
  static int const RESAMPLE = 52;
  static int const DESIGNATURE = 53;

  static int const TAPER_NONE     = -1;
  static int const TAPER_COSINE   = 1;
  static int const TAPER_HANNING  = 2;
  static int const TAPER_BLACKMAN = 3;

public:
  /**
   * Constructor - General use
   * @param numSamples Number of samples in input data.
   * @param sampleInt  Sample interval in input data in milliseconds [ms].
   */
  csFFTTools( int numSamplesIn, float sampleInt );
  /**
   * Constructor - Use for resampling operation
   * @param numSamplesIn  Number of samples in input data.
   * @param numSamplesOut Number of samples in output data.
   * @param sampleIntIn   Sample interval in input data in milliseconds [ms].
   * @param sampleIntOut  Sample interval in output data in milliseconds [ms].
   */
  csFFTTools( int numSamplesIn, int numSamplesOut, float sampleIntIn, float sampleIntOut );
  ~csFFTTools();

  void prepareBandpassFilter( int numLowPassFilters, float const* lowPassFreq, float const* lowPassSlope, int numHighPassFilters, float const* highPassFreq, float const* highPassSlope );
  float const* prepareNotchFilter( float notchFreqHz, float notchWidthHz, float slope, bool isCosineTaper );

  void applyBandpassFilter( float* samples, bool outputImpulseResponse );
  void applyNotchFilter( float* samples, bool outputImpulseResponse );
  void applyFilter( float* samples, csTimeFunction<double> const* freqFunc );
  void applyAmpFilter( float* samplesInOut, float const* filterCoefficients );

  void setFilter( float cutOffFreqHz, float slope );
  void setFilterWavelet( int length );

  float resample( float* samples, bool applyFilter, bool applyRMSNorm = false );
 
  /**
   * Apply Hilbert transform == Shift phase by 90deg
   * @param samples    Input/output samples, x-t
   */
  void hilbertTransform( float* samplesInOut );
  void hilbertTransform( float const* samplesIn, float* samplesOut );
  /**
   * Compute envelope
   * @param samples      Input/output samples, x-t
   */
  void envelope( float* samples );
  /**
   * Compute instantaneous phase
   * @param samples    Input/output samples, x-t
   */
  void inst_phase( float* samples );
  /**
   * Compute instantaneous frequency
   * @param samples    Input/output samples, x-t
   */
  void inst_freq( float* samples );

  /**
   * @return Number of samples in FFT transform of input data. Equals nearest power of 2 to numInputSamples().
   */
  int numFFTSamples() const;
  int numFFTSamplesOut() const;
  double sampleIntFreq() const;
  /**
   * @return Number of samples in input data.
   */
  int numInputSamples() const { return myNumSamplesIn; }
  
  static void applyTaper( int taperType, int taperLengthInSamples, int numSamplesIn, float* samples );

protected:
  void init( int numSamples_time, float sampleInt_time );
  void convertFromAmpPhase( float const* ampSpec, float const* phaseSpec );
  void unwrap_phase();
  cseis_geolib::csFFT* myFFT;
  cseis_geolib::csFFT* myFFTOut;

  /// Number of samples in input data
  int myNumSamplesIn;
  /// Number of samples in output data. Only different for resample operation
  int myNumSamplesOut;
  /// Sample interval of input data
  float mySampleIntIn;  // [ms]
  /// Sample interval of output data
  float mySampleIntOut; // [ms]

  float* myBuffer;
  float* myBufferOut;

  double* myFilterWavelet;
  int     myLengthFilterWavelet;

  float* myFilterScalars;
  float* myNotchFilterScalars;

  bool  myIsFilterWavelet;
};

} // namespace
#endif
