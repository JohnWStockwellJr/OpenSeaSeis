/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_INTERPOLATION_H
#define CS_INTERPOLATION_H

namespace cseis_geolib {

/**
 * Provides interpolation methods
 */
class csInterpolation {
 public:
  static const int METHOD_LIN  = 30;
  static const int METHOD_QUAD = 31;
  static const int METHOD_SINC = 32;

 public:
  csInterpolation( int numSamples, float sampleInt );
  csInterpolation( int numSamples, float sampleInt, int numCoefficients );
  csInterpolation( int numSamples, int interpolationMethod, float sampleInt );
  ~csInterpolation();

  void setExtrapolation( float valLeft, float valRight );
  void static_shift( float shift_ms, float const* samplesIn, float* samplesOut );
  /**
   * @param shift_ms   Time shift to be applied [ms]
   * @param samplesIn  Input samples
   * @param samplesOut Output samples, shifted
   */
  void static_shift( float shift_ms, float const* samplesIn, float* samplesOut, bool sameShift );
  void process( float sampleIntSkew, float xVal1, float const* samplesIn, float const* sIndexOut, float* samplesOut, int numSamplesOut );
  void process( float sampleIntSkew, float xVal1, float const* samplesIn, float const* sIndexOut, float* samplesOut ) {
    process( sampleIntSkew, xVal1, samplesIn, sIndexOut, samplesOut, myNumSamples );
  }
  float valueAt( float time_ms, float const* samplesIn );
  float timeAt_lin( float value, float const* samplesIn );
  void timeAt_lin( float const* samplesIn, float sampleInt_out, int numSamplesOut, float* timesOut );

  static double sincFunction( double value );
  static bool toeplitzSolver( int numDimensions, double const* topRow, double const* vecRight, double* vecLeft, double* vecSolve );
  static void createSincFunction( float ratio, int numCoef, float* sinc );

  /**
   * Linear interpolation between knee points
   * @param numSamplesIn   Number of samples in input array
   * @param xin            Input samples
   * @param yin            Input sample values
   * @param numSamplesOut  Number of samples in input array
   * @param sampleIntOut   Output sample interval (must be in same unit as input samples xin)
   * @param yout           Output sample values, to be computed
   */
  static void linearInterpolation( int numSamplesIn, float const* xin, float const* yin,
                                   int numSamplesOut, float sampleIntOut, float* yout );

  static void xy2yxInterpolation( float const* arrayIn, float* arrayOut, int numSamples );
  static void xy2yxInterpolation( float const* arrayIn, float* arrayOut, int numSamples, float sampleInt );
  static void xy2yxInterpolation( float const* arrayIn, float* arrayOut, int numSamples, float sampleInt, float valueInvalid );

 private:
  void init( int numSamples, int interpolationMethod, float sampleInt, int numCoefficients );
  float valueAt_sinc( float time_ms, float const* samplesIn );
  float valueAt_lin( float time_ms, float const* samplesIn );
  float valueAt_quad( float time_ms, float const* samplesIn );

  int myNumCoefficients;
  int myNumValues;
  float myExtrapolValLeft;
  float myExtrapolValRight;
  float** myCoefficients;
  float mySampleInt;
  int myNumSamples;
  int myMethod; /// Interpolation method
  float* myIndexBuffer;
  float myPreviousShift_ms; // Saved, previously applied time shift
};


}

#endif

