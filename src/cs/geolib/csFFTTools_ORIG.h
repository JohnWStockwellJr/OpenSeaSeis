/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_FFT_TOOLS_ORIG_H
#define CS_FFT_TOOLS_ORIG_H

namespace cseis_geolib {

template <typename T> class csTimeFunction;

class csFFTTools_ORIG {
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

  static bool Powerof2( int nx, int* m, int* twopm );
  static bool fft( int dir, int power_of_two, double *realValues, double *imagValues );
  static bool fft( int dir, int power_of_two, double *realValues, double *imagValues, bool doNormalize );
  static bool fft_2d( int dir, double** realValues, double** imagValues, int numFFTSamplesX, int numFFTSamplesY );

public:
  /**
   * @param numSamples Number of samples in input data.
   */
  csFFTTools_ORIG( int numSamples );
  /**
   * @param numSamples Number of samples in input data.
   * @param sampleInt  Sample interval in input data in milliseconds [ms].
   */
  csFFTTools_ORIG( int numSamples, float sampleInt );
  /**
   * @param numSamplesIn  Number of samples in input data.
   * @param numSamplesOut Number of samples in output data.
   * @param sampleIntIn   Sample interval in input data in milliseconds [ms].
   * @param sampleIntOut  Sample interval in output data in milliseconds [ms].
   */
  csFFTTools_ORIG( int numSamplesIn, int numSamplesOut, float sampleIntIn, float sampleIntOut );
  ~csFFTTools_ORIG();

  void setFilter( float cutOffFreqHz, float order, bool outputImpulseResponse );
  void setFilterWavelet( int length );

  void highPass( float* samples, float cutOffFreqHz, float order, bool outputImpulseResponse );
  void lowPass( float* samples, float cutOffFreqHz, float order, bool outputImpulseResponse );
  void applyFilter( float* samples, csTimeFunction<double> const* freqFunc );
  //  void applyMultipleFilters( float* samples, int numLowPassFilters, float const* lowPassFreq, float const* lowPassOrder, int numHighPassFilters, float const* highPassFreq, float const* highPassOrder, bool outputImpulseResponse );
  void prepareMultipleFilters( int numLowPassFilters, float const* lowPassFreq, float const* lowPassOrder, int numHighPassFilters, float const* highPassFreq, float const* highPassOrder );
  void applyLowHighPassFilters( float* samples, bool outputImpulseResponse );

// Notch filter:
  void notchFilter( float* samples, bool addNoise );
  double const* setupNotchFilter( float notchFreqHz, float notchWidthHz, float slope, bool isCosineTaper );
  float resample( float* samples );
  float resample( float* samples, bool applyFilter, bool applyNorm );
  float resample( float* samples, float cutOffFreqHz, float order, bool applyFilter, bool applyNorm );

  /**
   * Apply Hilbert transform == Shift phase by 90deg
   * @param samples    Input/output samples, x-t
   */
  void hilbertTransform( float* samples );
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
  //  void applyQCompensation( float* samples, float qvalue, float freqRef, bool applyAmp, bool applyPhase );

  /**
   * Apply forward FFT transform
   * @param samples          Input data. Must have numInputSamples() number of samples.
   */
  bool fft_forward( float const* samples );
  /**
   * Apply forward FFT transform. Compute amplitude spectrum.
   * @param samples          Input data. Must have numInputSamples() number of samples.
   * @param ampSpec          Output amplitude spectrum. Array must have numFFTSamples()/2+1 number of samples.
   */
  bool fft_forward( float const* samples, float* ampSpec );
  /**
   * Apply forward FFT transform. Compute amplitude & phase spectra.
   * @param samples          Input data. Must have myNumSamplesIn number of samples.
   * @param ampSpec          Output amplitude spectrum. Array must have numFFTSamples()/2+1 number of samples.
   * @param phaseSpec        Output phase spectrum. Array must have numFFTSamples()/2+1 number of samples.
   */
  bool fft_forward( float const* samples, float* ampSpec, float* phaseSpec );

  /**
   * Apply inverse FFT transform
   */
  bool fft_inverse( );
  /**
   * Apply inverse FFT transform from given amplitude & phase spectrum.
   * @param samples      Input data. Consists of either
   *                      a) Concatenated amplitude and phase spectra, numFFTSamples() number of samples, or
   *                      b) Concatenated real and imaginary spectra, 2*numFFTSamples() number of samples.
   * @param fftDataType  Data 'type': Either FX_REAL_IMAG or FX_AMP_PHASE.
   */
  bool fft_inverse( float const* samples, int fftDataType );
  /**
   * Apply inverse FFT transform from given amplitude & phase spectrum.
   * @param amplitude  Input amplitude spectrum. Array must consist of numFFTSamples()/2+1 number of samples.
   * @param phase      Input phase spectrum. Array must consist of numFFTSamples()/2+1 number of samples.
   */
  bool fft_inverse( float const* amplitude, float const* phase, bool doNormalize = true );

  double const* realData() const { return myBufferReal; }
  double const* imagData() const { return myBufferImag; }

  double* getRealDataPointer() { return myBufferReal; }
  double* getImagDataPointer() { return myBufferImag; }
  void convertToAmpPhase( float* ampSpec, float* phaseSpec );
  
  /**
   * @return Number of samples in FFT transform of input data. Equals nearest power of 2 to numInputSamples().
   */
  int numFFTSamples() const { return myNumFFTSamplesIn; }
  int numFFTSamplesOut() const { return myNumFFTSamplesOut; }
  float sampleIntFreqHz() const;
  /**
   * @return Number of samples in input data.
   */
  int numInputSamples() const { return myNumSamplesIn; }
  
  /**
   * ...for backward compatability: Set amplitude spectrum normalization as in vintage implementation
   */
  void setVintageAmpSpecNormalization();
  static void applyTaper( int taperType, int taperLengthInSamples, int numSamplesIn, float* samples );

protected:
  void filter( float* samples, int filterType );
  void applyLowPassFilter( float cutOffFreq, float order );
  void applyHighPassFilter( float cutOffFreq, float order );
  void init();
  void setBuffer( float const* samples );
  void convertFromAmpPhase( float const* ampSpec, float const* phaseSpec );

  /// Number of samples in input data
  int myNumSamplesIn;
  /// Number of samples in output data. Only different for resample operation
  int myNumSamplesOut;
  /// Sample interval of input data
  float mySampleIntIn;  // [ms]
  /// Sample interval of output data
  float mySampleIntOut; // [ms]
  /**
   * Number of samples used in FFT transform for input data.
   * This may be the same as myNumSamplesIn except it is the nearest power of 2 (y^2).
   */
  int myNumFFTSamplesIn;
  /// Number of samples used in FFT transform for output data. Required for resample operation
  int myNumFFTSamplesOut;
  int myTwoPowerIn;
  int myTwoPowerOut;
  double* myBufferReal;
  double* myBufferImag;
  double* myNotchFilter;

  double* myFilterWavelet;
  int     myLengthFilterWavelet;

  float* myFilterScalars;

  float myOrder;
  float myCutOffFreqHz;
  bool  myOutputImpulseResponse;
  bool  myIsFilterWavelet;
  float myAmpSpecConsistencyScalar;
};

} // namespace
#endif
