/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <cstdio>

namespace cseis_geolib {

//------------------------------------------------------
// Returns an estimated amplitude at the given sample index, using linear interpolation
//
float getLinAmplitudeAtSample( float const* traceData, double sample, int numSamples ) {
  float amplitude = 0.0;
  int isamp = ( int )sample;
  double y1 = 0.0;
  double y2 = 0.0;
  double weight = 0.0;
  if( isamp < 0 || isamp >= numSamples-1 ) {
    amplitude = 0.0;
  }
  else {
    y1 = traceData[isamp];
    y2 = traceData[isamp + 1];
    weight = (sample - (double)isamp)/1.0;
    amplitude = (float)(y1 + weight * (y2 - y1));
  }
  return amplitude;
}

//------------------------------------------------------
// Returns an estimated amplitude at the given sample index, using quadratic interpolation
//
float getQuadAmplitudeAtSample( float const* traceData, double sample, int numSamples )
{
  float amplitude = 0.0;
  int isamp = ( int )( sample );
  if( isamp < 1 ) {
    if( isamp < 0 ) {
      amplitude = 0.0;
    }
    else {
      amplitude = traceData[0];
    }
  }
  else if( isamp >= (numSamples - 1) ) {
    if( isamp == (numSamples - 1) ) {
      amplitude = traceData[numSamples-1];
    }
    else {
      amplitude = 0.0;
    }
  }
  else {
    double y1 = traceData[isamp - 1];
    double y2 = traceData[isamp];
    double y3 = traceData[isamp + 1];

    // Set x values to similar range as y, to avoid rounding errors
    double min = (y1 < y2 ? y1 : y2);  // Math.min( Math.min( y1, y2 ), y3 );
    min = (min < y3 ? min : y3);

    double max = (y1 > y2 ? y1 : y2);  // Math.max( Math.max( y1, y2 ), y3 );
    max = (max > y3 ? max : y3);

    double x1 = (y1 >= 0.0 ? y1 : -y1);  // = Math.abs( y1 );
    double x2 = x1 + ( max - min );
    double x3 = x1 + 2.0 * ( max - min );

    double xAmp = x2;
    if( (double)isamp < sample ) {
      xAmp = x2 + ( (sample-(double)isamp) * (x3 - x2) );
    }
    else {
      xAmp = x1 + ( (sample-(double)(isamp-1)) * (x2 - x1) );
    }
    double absolute = x2-x1;
    if( absolute < 0.0 ) absolute = -absolute;  // absolute = abs(x2-x1)
    if( absolute != 0.0 ) {
      double tmp1 = ( y2 - y1 ) / ( x2 - x1 );
      double tmp2 = tmp1 * ( x3 - x1 ) + y1 - y3;
      double tmp3 = ( x2 + x1 ) * ( x3 - x1 ) + x1 * x1 - x3 * x3;

      double a = 0.0;
      absolute = (tmp3 >= 0.0 ? tmp3 : -tmp3);  // absolute = abs(tmp3)
      if( absolute != 0.0 ) {
        a = tmp2 / tmp3;
      }
      double b = tmp1 - a * ( x2 + x1 );
      double c = y1 - ( a * x1 + b ) * x1;
      amplitude = (float)(( a * xAmp + b ) * xAmp + c);
    }
  }
  return amplitude;
}
//================================================================================
//
//
float getQuadMaxSample( float const* traceData, int midSamp, int numSamples, float* maxAmplitude )
{
  float maxSample = (float)midSamp;
  *maxAmplitude = traceData[midSamp];
  //  fprintf(stdout,"midSamp: %d %d\n",midSamp,numSamples);
  if( midSamp < 1 ) {
    if( midSamp < 0 ) {
      maxSample = 0.0;
      *maxAmplitude = traceData[0];
    }
    else {
      maxSample = 0.0;
      *maxAmplitude = traceData[0];
    }
  }
  else if( midSamp >= numSamples - 1 ) {
    maxSample = (float)(numSamples-1);
    *maxAmplitude = traceData[numSamples-1];
  }
  else {
    double y1 = traceData[midSamp - 1];
    double y2 = traceData[midSamp];
    double y3 = traceData[midSamp + 1];

    //    fprintf(stdout,"Amplitudes: %f %f %f\n",y1,y2,y3);

    // Set x values to similar range as y, to avoid rounding errors
    double min = (y1 < y2 ? y1 : y2);  // Math.min( Math.min( y1, y2 ), y3 );
    min = (min < y3 ? min : y3);

    double max = (y1 > y2 ? y1 : y2);  // Math.max( Math.max( y1, y2 ), y3 );
    max = (max > y3 ? max : y3);

    double x1 = (y1 >= 0.0 ? y1 : -y1);  // = Math.abs( y1 );
    double x2 = x1 + ( max - min );
    double x3 = x1 + 2.0 * ( max - min );

    double absolute = x2-x1;
    if( absolute < 0.0 ) absolute = -absolute;  // absolute = abs(x2-x1)
    if( absolute != 0.0 ) {
      double tmp1 = ( y2 - y1 ) / ( x2 - x1 );
      double tmp2 = tmp1 * ( x3 - x1 ) + y1 - y3;
      double tmp3 = ( x2 + x1 ) * ( x3 - x1 ) + x1 * x1 - x3 * x3;

      double a = 0.0;
      absolute = (tmp3 >= 0.0 ? tmp3 : -tmp3);  // absolute = abs(tmp3)
      if( absolute != 0.0 ) {
        a = tmp2 / tmp3;
      }
      double b = tmp1 - a * ( x2 + x1 );
      double c = y1 - ( a * x1 + b ) * x1;
      if( a != 0.0 ) {
        maxSample = -b/(2.0*a);
        *maxAmplitude = (float)(( a * maxSample + b ) * maxSample + c);
        // Convert max sample back to original range
        maxSample = (maxSample-x2)/(x3-x2)+midSamp;
      }
    }
  }
  return maxSample;
}

//================================================================================
//
//
float getQuadZeroSample( float const* traceData, int midSamp, int numSamples )
{
  float zeroSample = (float)midSamp;
  //  fprintf(stdout,"midSamp: %d %d\n",midSamp,numSamples);
  if( midSamp < 1 ) {
    if( midSamp < 0 ) {
      zeroSample = 0.0;
    }
    else {
      zeroSample = 0.0;
    }
  }
  else if( midSamp >= numSamples - 1 ) {
    zeroSample = (float)(numSamples-1);
  }
  else {
    double y1 = traceData[midSamp - 1];
    double y2 = traceData[midSamp];
    double y3 = traceData[midSamp + 1];

    // Set x values to similar range as y, to avoid rounding errors
    double min = (y1 < y2 ? y1 : y2);  // Math.min( Math.min( y1, y2 ), y3 );
    min = (min < y3 ? min : y3);

    double max = (y1 > y2 ? y1 : y2);  // Math.max( Math.max( y1, y2 ), y3 );
    max = (max > y3 ? max : y3);

    double x1 = (y1 >= 0.0 ? y1 : -y1);  // = Math.abs( y1 );
    double x2 = x1 + ( max - min );
    double x3 = x1 + 2.0 * ( max - min );

    double absolute = x2-x1;
    if( absolute < 0.0 ) absolute = -absolute;  // absolute = abs(x2-x1)
    if( absolute != 0.0 ) {
      double tmp1 = ( y2 - y1 ) / ( x2 - x1 );
      double tmp2 = tmp1 * ( x3 - x1 ) + y1 - y3;
      double tmp3 = ( x2 + x1 ) * ( x3 - x1 ) + x1 * x1 - x3 * x3;

      double a = 0.0;
      absolute = (tmp3 >= 0.0 ? tmp3 : -tmp3);  // absolute = abs(tmp3)
      if( absolute != 0.0 ) {
        a = tmp2 / tmp3;
      }
      double b = tmp1 - a * ( x2 + x1 );
      double c = y1 - ( a * x1 + b ) * x1;
      tmp1 = b*b - 4.0*c*a;

      if( tmp1 >= 0 && a != 0 ) {
        tmp2 = sqrt(tmp1);
        double result1 = 0.5/a * ( tmp2 - b );
        double result2 = 0.5/a * ( -tmp2 - b );
        result1 = (result1-x2)/(x3-x2)+midSamp;
        result2 = (result2-x2)/(x3-x2)+midSamp;

        if( result1 >= midSamp-1 && result1 <= midSamp ) {
          zeroSample = result1;
        }
        else if( result2 >= midSamp-1 && result2 <= midSamp ) {
          zeroSample = result2;
        }
      }
    }
  }
  return zeroSample;
}

/**
 * Sinc function
 * sinc(x) = sin(pi*x)/(pi*x)
 */
template <typename T> T sincFunction( T value ) {
  if( value != 0.0 ) {
    T tmp = M_PI * value;
    return sin(tmp)/tmp;
  }
  else {
    return 1.0;
  }
}

/**
 * Toeplitz equation solver
 *   M x vecLeft  = vecRight
 *   M x vecSolve = val0
 */

template <typename T> bool toeplitz_solver( int numDimensions, T const* topRow, T const* vecRight, T* vecLeft, T* vecSolve ) {
  T val0 = topRow[0];
  if( val0 == 0.0 ) return false;

  int numDimHalf = numDimensions/2;
  vecLeft[0]  = vecRight[0] / topRow[0];
  vecSolve[0] = 1.0;

  for( int idim = 1; idim < numDimensions; idim++ ) {
    vecSolve[idim] = 0.0;
    vecLeft[idim]  = 0.0;

    T solve1 = 0.0;
    for( int i = 0; i < idim; i++ ) {
      solve1 += vecSolve[i]*topRow[idim-i];
    }

    T solve2 = solve1/val0;
    val0  -= solve2*solve1;
    for( int i = 0; i <= numDimHalf; i++ ) {
      T tmp = vecSolve[idim-i] - solve2 * vecSolve[i];
      vecSolve[i] -= solve2 * vecSolve[idim-i];
      vecSolve[idim-i] = tmp;
    }

    T solve3 = 0.0;
    for( int i = 0; i < idim; i++ ) {
      solve3 += vecLeft[i] * topRow[idim-i];
    }

    solve2 = ( solve3 - vecRight[idim] ) / val0;
    for( int i = 0; i <= idim; i++ ) {
      vecLeft[i] -= solve2 * vecSolve[idim-i];
    }
  }
}

} // end namespace

