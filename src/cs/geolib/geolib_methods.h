/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef GEOLIB_METHODS_H
#define GEOLIB_METHODS_H

namespace cseis_geolib {

/**
* SVD decomposition
* A = U V W
*
* @param mat_a   (io) Input matrix, NROWS x NCOLS. Output matrix 'U', containing eigen vectors
* @param nrows   (i) Number of rows in input/output matrix
* @param ncols   (i) Number of columns in input/output matrix
* @param vec_w   (o) Output vector 'W', NCOLS
* @param mat_v   (o) Output matrix 'V', NCOLS x NCOLS
* @param vec_work(io) Help vector, needed during SVD decomposition, NROWS
*/
bool svd_decomposition( float** mat_a, int nrows, int ncols, float vec_w[], float** mat_v, float vec_work[] );
bool svd_decomposition( double** mat_a, int nrows, int ncols, double vec_w[], double** mat_v, double vec_work[] );

/**
* Solve linear equation system via SVD decomposition
*
* @param mat_u   (i) NROWS x NCOLS. Output matrix 'U' from SVD decomposition
* @param vec_w   (i) NCOLS, Output vector 'W' from SVD decomposition
* @param mat_v   (i) NCOLS x NCOLS, Output matrix 'V' from SVD decomposition
* @param nrows   (i) Number of rows in input/output matrix
* @param ncols   (i) Number of columns in input/output matrix
* @param vec_b   (i) Known solution vector, NROWS
* @param vec_x   (o) Unknown solution vector, NCOLS
* @param vec_work(io) Help vector, needed during SVD decomposition, NROWS
*/
void svd_linsolve(float **mat_u, float vec_w[], float **mat_v, int nrows, int ncols, float vec_b[], float x[], float vec_work[]);
void svd_linsolve(double **mat_u, double vec_w[], double **mat_v, int nrows, int ncols, double vec_b[], double x[], double vec_work[]);

/**
* Fit 3D line to data points
* @param xSamples    (i) X data values
* @param ySamples    (i) Y data values
* @param zSamples    (i) Z data values
* @param firstSample (i) Sample index of first sample to use (starting at 0)
* @param lastSample  (i) Sample index of last sample to use (starting at 0)
* @param forceOrigin (i) 'True' if 3D line shall pass through origin point (0,0,0)
* @param vec_out     (o) Output (XYZ) vector giving 3D line direction
*/
bool linefit_3d( float* xSamples,
                 float* ySamples,
                 float* zSamples,
                 int firstSample,
                 int lastSample,
                 int force_origin,
                 float* vec_out );

/**
* Fit 3D line to data points
* @param xSamples    (i) X data values
* @param ySamples    (i) Y data values
* @param zSamples    (i) Z data values
* @param firstSample (i) Sample index of first sample to use (starting at 0)
* @param lastSample  (i) Sample index of last sample to use (starting at 0)
* @param forceOrigin (i) 'True' if 3D line shall pass through origin point (0,0,0)
* @param vec_out     (o) Output (XYZ) vector giving 3D line direction
* @param vec_axes    (o) Size of major, minor, and second minor axes of best fit spheroid
* @param vec_minor   (o) 3D vector of minor spheroid axis
*/
bool linefit_3d_all( float* xSamples,
                     float* ySamples,
                     float* zSamples,
                     int firstSample,
                     int lastSample,
                     int force_origin,
                     float* vec_out,
                     float* vec_axes,
                     float* vec_minor );  // Size of major, minor, and second minor axes of best fit spheroid
/**
* ...not working
*/
bool linefit_3d_2step( float* xSamples,
                       float* ySamples,
                       float* zSamples,
                       int nSamples,
                       int method,
                       int force_origin,
                       float* vec_out );

/**
 * Returns polarity (+/-1) of wavelet inside specified time window
 *
 * Polarity is determined in the following way (optimal for minimum phase wavelets):
 * In the specified time window, the minimum and maximum amplitudes are determined.
 * If the maximum amplitude precedes the minimum amplitude, the polarity is positive.
 *
 * @param samples  (i) Input data values
 * @param nSamples (i) Number samples in input trace
 * @return -1.0 or +1.0
 */
float get_polarity_minphase( float* samples, int nSamples );

/**
 * Returns polarity (+/-1) of wavelet inside specified time window
 *
 * Polarity is determined in the following way (optimal for zero phase wavelets):
 * In the specified time window, the minimum and maximum amplitudes are determined.
 * If the maximum amplitude is posituve and larger than the absolute of the minimum amplitude, the polarity is positive.
 *
 * @param samples  (i) Input data values
 * @param nSamples (i) Number samples in input trace
 * @return -1.0 or +1.0
 */
float get_polarity_zerophase( float* samples, int nSamples );

/**
* Estimate polarity of input data by cross-correlation with second trace
* Determine cross-correlation value. Output polarity of cross-correlation value
* @param psamples (i) Input data values, additional trace
* @param samples  (i) Input data values
* @param nSamples (i) Number samples in input trace
* @return -1.0 or +1.0
*/
float get_polarity_p( float *pSamples, float* samples, int nSamples );

/**
* Correct polarity of input vector
* Checks polarity and multiplies input vector by -1 if necessary.
* @param xSamples    (i)
* @param ySamples    (i)
* @param zSamples    (i)
* @param firstSample (i) Sample index of first sample to use (starting at 0)
* @param lastSample  (i) Sample index of last sample to use (starting at 0)
* @param isMinPhase  (i) true if wavelet is assumed to be minimum phase, false if zero phase
* @param vec_inout   (io) Input vector. Will be polarity flipped if necessary.
*/
void polarity_correction( float* xSamples,
                          float* ySamples,
                          float* zSamples,
                          int firstSample,
                          int lastSample,
                          bool isMinPhase,
                          float* vec_inout );

/**
* Correct polarity of input vector
* Checks polarity and multiplies input vector by -1 if necessary.
* @param xSamples    (i)
* @param ySamples    (i)
* @param ZSamples    (i)
* @param pSamples    (i)
* @param firstSample (i) Sample index of first sample to use (starting at 0)
* @param lastSample  (i) Sample index of last sample to use (starting at 0)
* @param vec_inout   (io) Input vector. Will be polarity flipped if necessary.
*/
void polarity_correction_p( float* xSamples,
                            float* ySamples,
                            float* zSamples,
                            float* pSamples,
                            int firstSample,
                            int lastSample,
                            float* vec_inout );


/*
* Fix: Quad interpolation, computation of isamp, 18.04.06
*
*/

/**
* Estimate amplitude at the given sample index, using linear interpolation
* @param traceData   (i) Input data values/trace
* @param sample      (i) Sample index at which amplitude shall be computed (starting at 0)
* @param numSamples  (i) Number of samples in input trace
* @return estimated amplitude
*/
float getLinAmplitudeAtSample( float const* traceData, double sample, int numSamples );

/**
* Estimate amplitude at the given sample index, using quadratic interpolation
* @param traceData   (i) Input data values/trace
* @param sample      (i) Sample index at which amplitude shall be computed (starting at 0)
* @param numSamples  (i) Number of samples in input trace
* @return estimated amplitude
*/
float getQuadAmplitudeAtSample( float const* traceData, double sample, int numSamples );


/**
* Find maximum amplitude around given sample index, using quadratic interpolation
* @param traceData    (i) Input data values/trace
* @param midSamp      (i) Sample index around which maximum amplitude shall be computed
* @param numSamples   (i) Number of samples in input trace
* @param maxAmplitude (o) Maximum amplitude
* @return sample index giving maximum amplitude
*/
float getQuadMaxSample( float const* traceData, int midSamp, int numSamples, float* maxAmplitude );
float getQuadZeroSample( float const* traceData, int midSamp, int numSamples );

void orientation_compute_roll_tilt (
  float const* vec_polarisation, // Unit polarisation vector in as-laid coordinate system
  double rcvx,     // Receiver X position [m]
  double rcvy,     // Receiver Y position [m]
  double rcvz,     // Receiver Z position [m]
  double srcx,     // Source X position [m]
  double srcy,     // Source Y position [m]
  double srcz,     // Source Z position [m]
  float v1v2_ratio,         // Velocity ratio at seabed
  int   isOBCMode,          // true(1) if source is above receiver, false(0) if source is beneath receiver
  float angle_azim,         // Receiver azimuth [rad]
  float* angle_tilt,        // Tilt angle [rad]
  float* angle_roll );       // Roll angle [rad]

 void orientation_compute_azimuth(
  float const* vec_polarisation,  // Unit polarisation vector in as-laid coordinate system
  double rcvx,     // Receiver X position [m]
  double rcvy,     // Receiver Y position [m]
  double srcx,     // Source X position [m]
  double srcy,     // Source Y position [m]
  float* angle_azim // Computed azimuth [rad]
  );


/**
 * Solves for isotropic CCP position using Tessmer & Behle's (1988) formula
 * (i) vpvs : vp/vs ratio
 * (i) offset: src-rcv offset [m]
 * (i) depth : Reflector depth [m]
 * (o) offset_pside: src-ccp offset [m]
*/
void ccp_offset_iso( double vpvs, double offset, double depth, double* offset_pside );


void pz_analysis( float const* const* samples_p,  // (i) Hydrophone samples, for all input traces
                  float const* const* samples_z,  // (i) Geophone samples, for all input traces
                  int    nTraces,           // (i) Number of input traces
                  int    nSamples,          // (i) Number of input samples
                  float* refl_coef,         // (o) Seabed reflection coefficient
                  float* pz_scalar,         // (o) PZ scalar
                  float* wdep,              // (o) Water depth [m]
                  float* z_phase );         // (o) Relative phase of geophone [rad]

void pz_analysis2( float const* const* p_samples,  // (i) Hydrophone samples, for all input traces
                   float const* const* z_samples,  // (i) Geophone samples, for all input traces
                   int    nTraces,           // (i) Number of input traces
                   int    nSamples,          // (i) Number of input samples
                   float  waterVel,         // (i) Water velocity [m/s]
                   float* refl_coef,         // (o) Seabed reflection coefficient
                   float* pz_scalar,         // (o) PZ scalar
                   float* wdep,              // (i/o) Water depth [m]
                   float* z_phase );         // (o) Relative phase of geophone [rad]


} // end namespace

#endif


