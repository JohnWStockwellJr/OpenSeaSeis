/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_FFT_DESIGNATURE_H
#define CS_FFT_DESIGNATURE_H

#include "csFFTTools.h"
#include <cstdio>

namespace cseis_geolib {

/**
 * Provide simple functionality for source designature and similar processes
 */
class csFFTDesignature : public csFFTTools {
public:
  static int const AMP_ONLY   = 101;
  static int const PHASE_ONLY = 102;
  static int const AMP_PHASE  = 103;

  csFFTDesignature( int numSamples, float sampleInt );
  /**
   * @param numSamples     Number of samples in input wavelet.
   * @param sampleInt      Sample interval of input wavelet, in milliseconds [ms].
   * @param inputWavelet   Input wavelet. For example a source signature, sensor response etc..
   * @param timeZero_s     Time in input wavelet which corresponds to time zero, in seconds [s].
   * @param percWhiteNoise Percent white noise to apply to wavelet spectrum before computing inverse filter.
   */
  csFFTDesignature( int numSamples, float sampleInt, float const* inputWavelet, float timeZero_s, float percWhiteNoise,
                    float const* output_wavelet = NULL );
  ~csFFTDesignature();

  /**
   * @param inputWavelet   Input wavelet. For example a source signature, sensor response etc..
   * @param timeZero_s     Time in input wavelet which corresponds to time zero, in seconds [s].
   * @param percWhiteNoise Percent white noise to apply to wavelet spectrum before computing inverse filter.
   */
  void initialize( float const* inputWavelet, float timeZero_s, float percWhiteNoise, float const* output_wavelet = NULL );

  /**
   * Set desired output wavelet.
   * If no output wavelet is specified, the desired output wavelet is assumed to be a spike/white spectrum.
   */
  // void setOutputWavelet( int numSamples, float sampleInt, float const* outputWavelet );
  
  /**
   * Apply low pass filter to designature filter.
   * @param cutOffHz Cut-off frequency [Hz].
   * @param slope    Filter slope [dB/oct].
   */
  void setDesigLowPass( float cutOffHz, float slope );
  /**
   * Apply high pass filter to designature filter.
   * @param cutOffHz Cut-off freuquency [Hz].
   * @param slope    Filter slope [dB/oct].
   */
  void setDesigHighPass( float cutOffHz, float slope );

  void setDesigHighEnd( float freq );
  /**
   * Apply cosine taper around notch frequency
   */
  void setNotchSuppression( float notchFreq, float notchWidth );
  /**
   * Set designature filter type.
   * @param filterType  AMP_ONLY for amplitude only, PHASE_ONLY for phase only, or AMP_PHASE for both amplitude and phase.
   */
  void setDesigFilterType( int filterType );

  /**
   * Apply designature filter to input data.
   * @param samples    Data to be filtered --> filtered output data.
   * @param numSamples Number of samples in input data. Must be the same as for the signature wavelet.
   */
  void applyFilter( float* samples, int numSamples );

  /**
   * Write designature filter as spectrum to output stream
   */
  void dump_spectrum( std::FILE* stream ) const;
  /**
   * Write designature filter as wavelet to output stream
   */
  void dump_wavelet( std::FILE* stream, bool doNormalize, float timeShift_s );

private:
  void initDesignature( float const* input_wavelet, float percWhiteNoise, float timeZero_s, float const* output_wavelet );

  float* myDesigAmpFilter;
  float* myDesigPhaseShift;

  float* myAmpPhaseSpecIn;

  int myFilterType;
};

} // namespace
#endif

