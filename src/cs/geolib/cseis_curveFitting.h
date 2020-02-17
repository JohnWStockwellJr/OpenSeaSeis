/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CSEIS_CURVE_FITTING_H
#define CSEIS_CURVE_FITTING_H

namespace cseis_geolib {

/**
* Fit cosine curve of given periodicity by cross-correlation of data point with cosine curve
*
* @param angles        (i) Angles in radians
* @param values        (i) Values (amplitudes)
* @param nAngles       (i) Number of angles/values in array
* @param periodicity   (i) Cosine periodicity in 360deg (for example 2 to look for azimuthal anisotropy)
* @param angleInc_deg  (i) Angle increment to search for [deg]
* @param computeStddev (i) If true, compute standard deviation of fitted cosine curve
* @param result        (o) Resultant angle of maximum value
* @param stddev        (o) Standard deviation of result
* @param amplitude     (o) Resultant amplitude of best-fit cosine function
*/
void computeXcorCos(
  double const* angles,
  double const* values,
  int nAngles,
  int periodicity,
  double  angleInc_deg,
  bool computeStddev,
  double& angleOut_deg,
  double& stddev,
  double& amplitude );

void computeLSCos(
                  double const* angles,
                  double const* values,
                  int nAngles,
                  int periodicity,
                  double  angleInc_deg,
                  bool computeStddev,
                  double  assumedAmplitude,
                  double& angleOut_deg,
                  double& stddev,
                  double& amplitude );

 bool polynom_fit( double* xValues, double* yValues, int nValues, int order, double* coefficients );

} // end namespace

#endif

