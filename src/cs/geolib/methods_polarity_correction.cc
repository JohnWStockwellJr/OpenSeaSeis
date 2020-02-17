/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

/*-------------------------------------------------------------------
 *
 * Estimate polarity on input data
 * Apply polarity correction to input vector
 *
 * Polarity of input data is estimated by picking first peak and first trough
 * inside data window, assuming minimum phase wavelet.
 *
 *-------------------------------------------------------------------*/

#include "geolib_math.h"
#include "geolib_methods.h"

namespace {
#ifndef MAXFLOAT
#ifndef FLT_MAX
  const float MAXFLOAT = 1.e+20;
#else
  const float MAXFLOAT = FLT_MAX;
#endif
#endif
}

namespace cseis_geolib {

void polarity_correction( float* xSamples,
        float* ySamples,
        float* zSamples,
        int firstSample,
        int lastSample,
        bool isMinPhase,
        float* vec_inout )
{
  const int NCOLS = 3;
  float pol_x, pol_y, pol_z, polarity;
  int nSamples;
  int icol;

  //------------------------------------------------------------
  nSamples = lastSample - firstSample + 1;

  if( isMinPhase ) {
    pol_x = get_polarity_minphase( xSamples+firstSample, nSamples );
    pol_y = get_polarity_minphase( ySamples+firstSample, nSamples );
    pol_z = get_polarity_minphase( zSamples+firstSample, nSamples );
  }
  else {
    pol_x = get_polarity_zerophase( xSamples+firstSample, nSamples );
    pol_y = get_polarity_zerophase( ySamples+firstSample, nSamples );
    pol_z = get_polarity_zerophase( zSamples+firstSample, nSamples );
  }

  // To estimate the polarity of the input vector, weight the polarity of individual components by their amplitudes
  polarity = SIGN( pol_x * vec_inout[0] + pol_y * vec_inout[1] + pol_z * vec_inout[2] );

  for( icol = 0; icol < NCOLS; icol++ ) {
    vec_inout[icol] *= polarity;
  }

  //  printf("POLXYZ: %10.5f %10.5f %10.5f ++ %10.5f   %10.5f\n", pol_x, pol_y, pol_z, polarity, vec_inout[0]);

  return;
}

void polarity_correction_p( float* xSamples,
        float* ySamples,
        float* zSamples,
        float* pSamples,
        int firstSample,
        int lastSample,
        float* vec_inout )
{
  const int NCOLS = 3;
  float pol_x, pol_y, pol_z, polarity;
  int nSamples;
  int icol;

  //------------------------------------------------------------
  nSamples = lastSample - firstSample + 1;

  pol_x = get_polarity_p( pSamples+firstSample, xSamples+firstSample, nSamples );
  pol_y = get_polarity_p( pSamples+firstSample, ySamples+firstSample, nSamples );
  pol_z = get_polarity_p( pSamples+firstSample, zSamples+firstSample, nSamples );

  // To estimate the polarity of the input vector, weight the polarity of individual components by their amplitudes
  polarity = SIGN( pol_x * vec_inout[0] + pol_y * vec_inout[1] + pol_z * vec_inout[2] );

  for( icol = 0; icol < NCOLS; icol++ ) {
    vec_inout[icol] *= polarity;
  }

  //  printf("POLXYZ: %10.5f %10.5f %10.5f ++ %10.5f   %10.5f\n", pol_x, pol_y, pol_z, polarity, vec_inout[0]);

  return;
}

/**
 * Returns polarity (+/-1) of wavelet inside specified time window
 *
 * Polarity is determined in the following way (optimal for minimum phase wavelets):
 * In the specified time window, the minimum and maximum amplitudes are determined.
 * If the maximum amplitude precedes the minimum amplitude, the polarity is positive.
 */
float get_polarity_minphase( float* samples, int nSamples ) {

  float min_value = MAXFLOAT;
  float max_value = -MAXFLOAT;
  int isamp, min_index, max_index;
  float polarity;

  min_index = 0;
  max_index = 0;
  for( isamp = nSamples-1; isamp >= 0; isamp-- ) {
    if( samples[isamp] <= min_value ) {
      min_value = samples[isamp];
      min_index = isamp;
    }
    if( samples[isamp] >= max_value ) {
      max_value = samples[isamp];
      max_index = isamp;
    }
  }
  polarity = ( max_index < min_index ) ? 1.0 : -1.0;

  return polarity;
}

/**
 * Returns polarity (+/-1) of wavelet inside specified time window
 *
 * Polarity is determined in the following way (optimal for zero phase wavelets):
 * In the specified time window, the minimum and maximum amplitudes are determined.
 * If the maximum amplitude is posituve and larger than the absolute of the minimum amplitude, the polarity is positive.
 */
float get_polarity_zerophase( float* samples, int nSamples ) {

  float min_value = MAXFLOAT;
  float max_value = -MAXFLOAT;
  int isamp;
  float polarity;

  for( isamp = nSamples-1; isamp >= 0; isamp-- ) {
    if( samples[isamp] <= min_value ) {
      min_value = samples[isamp];
    }
    if( samples[isamp] >= max_value ) {
      max_value = samples[isamp];
    }
  }
  polarity = ( max_value > -min_value && max_value > 0.0 ) ? 1.0 : -1.0;

  return polarity;
}
//-----------------------------------------------------------------------
float get_polarity_p( float* pSamples, float* samples, int nSamples ) {
  int isamp;
  float polarity;
  float sum;

  sum = 0.0;
  for( isamp = nSamples-1; isamp >= 0; isamp-- ) {
    sum += pSamples[isamp]*samples[isamp];
  }
  // Expected polarity of hydrophone for direct arrival is opposite to geophones
  polarity = ( sum >= 0.0 ) ? -1.0 : 1.0;

  return polarity;
}

} // namespace

