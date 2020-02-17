/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_curveFitting.h"
#include <cmath>
#include <cstdio>
#include "geolib_math.h"
#include "geolib_methods.h"

namespace cseis_geolib {

  void computeXcorCos(
                      double const* angles,
                      double const* values,
                      int     nAngles,
                      int     periodicity,
                      double  angleInc_deg,
                      bool    computeStddev,
                      double& angleOut_deg,
                      double& stddev,
                      double& amplitude )
  {
    if( periodicity < 1 ) periodicity = 1;
    if( angleInc_deg < 1.0e-5 ) angleInc_deg = 1.0e-5;
    double periodicity_dbl = (double)periodicity;
    double maxTestAngle = 360.0 / periodicity_dbl;
    int numTestAngles   = (int)round( maxTestAngle / angleInc_deg ) + 1;
    double angleInc_rad = angleInc_deg * M_PI / 180.0;

    // 1) Find best fit angle (angleOut_rad)
    double angleOut_rad = 0;
    double maxValue = -9999999.0;
    for( int itest = 0; itest < numTestAngles; itest++ ) {
      double angle_rad = double(itest) * angleInc_rad;
      double sum = 0.0;
      for( int i = 0; i < nAngles; i++ ) {
        sum += values[i] * cos( periodicity_dbl*( angles[i] - angle_rad) );
      }
      if( sum > maxValue ) {
        maxValue    = sum;
        angleOut_rad = angle_rad;
      }
      //      fprintf(stdout,"%d %.2f %e  XCOR\n", itest, angle_rad*180/M_PI, sum);
    }
    angleOut_deg = angleOut_rad * 180.0 / M_PI;

    /// Do quadratic interpolation to get closer to real maximum
    //  float getQuadAmplitudeAtSample( float const* values, angleOut_deg, int numSamples )

    // 2) Find amplitude of cosine function, by least-square's solution for f(A,a) = A*cos(a - phi)
    //    a(i): Specified angles
    //    y(i): Specified amplitude
    //    A   : Amplitude of best-fit cosine function
    //    phi : Phase shift of best-fit cosine function
    //
    //    g(a)  = SUM{ y - A * cos(a-phi) }^2 = MINIMUM
    // -> dg/dA = SUM{ 2*( y - A * cos(a-phi) ) * -cos(a-phi) } = 0
    // ->         2 * SUM{ -y * cos(a-phi) } + 2*A * SUM{ cos^2(a-phi) } = 0
    // ->     A = SUM{ y * cos(a-phi) } / SUM{ cos^2(a-phi) }
    //
    double sum = 0.0;
    for( int i = 0; i < nAngles; i++ ) {
      double tmp_cos = cos( periodicity_dbl*( angles[i] - angleOut_rad) );
      sum += tmp_cos*tmp_cos;
    }
    if( sum != 0.0 ) {
      amplitude = maxValue/sum;
    }
    else {
      amplitude = 0.0;
    }

    if( computeStddev ) {
      sum = 0.0;
      for( int i = 0; i < nAngles; i++) {
        sum += CS_SQR( values[i] - amplitude*cos( periodicity_dbl*( angles[i] - angleOut_rad) ) );
        //printf("SUM %d %f %f   %f  (%f)\n", i, RAD2DEG(angles[i]), values[i], sum, amplitude*cos( periodicity_dbl*( angles[i] - angleOut_rad) ) );
      }
      stddev = sqrt( sum/(double)nAngles );
    }
  }


  void computeLSCos(
                    double const* angles,
                    double const* values,
                    int     nAngles,
                    int     periodicity,
                    double  angleInc_deg,
                    bool    computeStddev,
                    double  assumedAmplitude,
                    double& angleOut_deg,
                    double& stddev,
                    double& amplitude )
  {
    if( periodicity < 1 ) periodicity = 1;
    if( angleInc_deg < 1.0e-5 ) angleInc_deg = 1.0e-5;
    double periodicity_dbl = (double)periodicity;
    double maxTestAngle = 360.0 / periodicity_dbl;
    int numTestAngles   = (int)round( maxTestAngle / angleInc_deg ) + 1;
    double angleInc_rad = angleInc_deg * M_PI / 180.0;

    // 1) Find best fit angle (angleOut_rad)
    double angleOut_rad = 0;
    double minValue = 1.0e30;
    for( int itest = 0; itest < numTestAngles; itest++ ) {
      double angle_rad = double(itest) * angleInc_rad;
      double sum = 0.0;
      for( int i = 0; i < nAngles; i++ ) {
        sum += pow( values[i] - assumedAmplitude * cos( periodicity_dbl*( angles[i] - angle_rad) ), 2 );
      }
      sum = sqrt( sum / (double)nAngles );
      if( sum < minValue ) {
        minValue    = sum;
        angleOut_rad = angle_rad;
      }
      //      fprintf(stdout,"%d %.2f %e  LS\n", itest, angle_rad*180/M_PI, sum);
    }
    angleOut_deg = angleOut_rad * 180.0 / M_PI;

    /*
    // Old way of estimating amplitude
    double sum1 = 0.0;
    double sum2 = 0.0;
    for( int i = 0; i < nAngles; i++ ) {
      double tmp_cos = cos( periodicity_dbl*( angles[i] - angleOut_rad) );
      sum1 += values[i]*tmp_cos;
      sum2 += tmp_cos*tmp_cos;
    }
    if( sum2 != 0.0 ) {
      amplitude = sum1/sum2;
    }
    else {
      amplitude = 0.0;
    }
    */

    // Step 2: LS fit to find best-fit amplitude
    int numDim = 2;
    double* AT_A_inv = new double[numDim*numDim];
    double sum = 0.0;
    double sum2 = 0.0;
    double ysum1 = 0.0;
    double ysum2 = 0.0;
    for( int i = 0; i < nAngles; i++ ) {
      double value = cos( periodicity_dbl*( angles[i] - angleOut_rad ) );
      sum += value;
      sum2 += value*value;
      ysum1 += 1 * values[i];
      ysum2 += value * values[i];
    }
    //  AT_A[0] = nAngles;
    //  AT_A[1] = sum;
    //  AT_A[2] = sum;
    //  AT_A[3] = sum2;

    double det = nAngles * sum2 - sum*sum;
    AT_A_inv[0] = sum2 / det;
    AT_A_inv[1] = -sum / det;
    AT_A_inv[2] = -sum / det;
    AT_A_inv[3] = nAngles / det;

    //    double c1 = AT_A_inv[0] * ysum1 + AT_A_inv[1] * ysum2; // Bias
    double c2 = AT_A_inv[2] * ysum1 + AT_A_inv[3] * ysum2;

    delete [] AT_A_inv;

    //    fprintf(stderr,"Amplitude %f --> %f   Bias: %f\n", amplitude, c2, c1);
    amplitude = c2;

    if( computeStddev ) {
      double sum = 0.0;
      for( int i = 0; i < nAngles; i++) {
        sum += CS_SQR( values[i] - amplitude*cos( periodicity_dbl*( angles[i] - angleOut_rad) ) );
      }
      stddev = sqrt( sum/(double)nAngles );
    }
  }

  //================================================================================
  //

bool polynom_fit( double* xValues, double* yValues, int nValues, int order, double* coefficients ) {
  int nUnknowns = order + 1;

  if( nUnknowns > nValues ) {
    //    fprintf(stderr,"Fewer input values (%d) than unknowns (%d)\n", nValues, nUnknowns );
    return false;
  }
  
  double** mat_v  = NULL;
  double** mat_u  = NULL;
  double*  vec_w  = new double[nValues];
  double*  vec_work = new double[nValues];

  mat_u = new double*[nValues];
  mat_v = new double*[nUnknowns];
  for( int irow = 0; irow < nValues; irow++ ) {
    mat_u[irow] = new double[nUnknowns];
  }
  for( int icol = 0; icol < nUnknowns; icol++ ) {
    mat_v[icol] = new double[nUnknowns];
  }
  
  
  for( int irow = 0; irow < nValues; irow++ ) {
    double xval = xValues[irow];
    mat_u[irow][0] = 1.0;
    for( int icol = 1; icol < nUnknowns; icol++ ) {
      mat_u[irow][icol] = pow( xval, icol );
    }
  }
  bool ret = svd_decomposition( mat_u, nValues, nUnknowns, vec_w, mat_v, vec_work );
  if( !ret ) {
    fprintf(stderr,"Error occurred in SVD\n");
    return false;
  }

  svd_linsolve( mat_u, vec_w, mat_v, nValues, nUnknowns, yValues, coefficients, vec_work );

  for( int irow = 0; irow < nValues; irow++ ) {
    delete [] mat_u[irow];
  }
  for( int icol = 0; icol < nUnknowns; icol++ ) {
    delete [] mat_v[icol];
  }
  delete [] mat_u;
  delete [] mat_v;
  delete [] vec_w;
  delete [] vec_work;

  return true;
}

} // end namespace

