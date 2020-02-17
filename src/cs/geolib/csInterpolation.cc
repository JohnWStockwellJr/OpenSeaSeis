/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csInterpolation.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "geolib_math.h"
#include "csException.h"

using namespace cseis_geolib;

csInterpolation::csInterpolation( int numSamples, float sampleInt ) {
  init( numSamples, csInterpolation::METHOD_SINC, sampleInt, 16 );
}
csInterpolation::csInterpolation( int numSamples, float sampleInt, int numCoefficients ) {
  init( numSamples, csInterpolation::METHOD_SINC, sampleInt, numCoefficients );
}
csInterpolation::csInterpolation( int numSamples, int interpolationMethod, float sampleInt ) {
  init( numSamples, interpolationMethod, sampleInt, 16 );
}

void csInterpolation::init( int numSamples, int interpolationMethod, float sampleInt, int numCoefficients ) {
  if( interpolationMethod != csInterpolation::METHOD_LIN && interpolationMethod != csInterpolation::METHOD_QUAD && interpolationMethod != csInterpolation::METHOD_SINC ) {
    throw( csException("csInterpolation::setInterpolationMethod: unknown method index provided (=%d). This is a bug in the calling function", interpolationMethod) );
  }
  myMethod = interpolationMethod;
  mySampleInt        = sampleInt;
  myNumSamples       = numSamples;
  myNumCoefficients  = numCoefficients;
  myNumValues        = 0;
  myPreviousShift_ms = 0;
  myIndexBuffer      = NULL;
  myCoefficients     = NULL;
  myExtrapolValLeft  = 0.0;
  myExtrapolValRight = 0.0;

  if( myMethod == csInterpolation::METHOD_SINC ) {
    myNumValues = 513;
    myCoefficients = new float*[myNumValues];
    for( int ival = 0; ival < myNumValues; ival++ ) {
      myCoefficients[ival] = new float[myNumCoefficients];
    }
    for( int ival = 1; ival < myNumValues-1; ival++ ) {
      float ratio = (float)ival/(float)(myNumValues-1);
      csInterpolation::createSincFunction( ratio, myNumCoefficients, &myCoefficients[ival][0] );
    }
    for( int icoef = 0; icoef < myNumCoefficients; icoef++ ) {
      myCoefficients[0][icoef] = 0.0;
      myCoefficients[myNumValues-1][icoef] = 0.0;
    }
    myCoefficients[0][myNumCoefficients/2-1] = 1.0;
    myCoefficients[myNumValues-1][myNumCoefficients/2] = 1.0;
  }
}

csInterpolation::~csInterpolation() {
  if( myCoefficients != NULL ) {
    for( int ival = 0; ival < myNumValues; ival++ ) {
      delete [] myCoefficients[ival];
    }
    delete [] myCoefficients;
  }
  if( myIndexBuffer != NULL ) {
    delete [] myIndexBuffer;
    myIndexBuffer = NULL;
  }
}

//void csInterpolation::setInterpolationMethod( int method ) {
//  if( method != csInterpolation::METHOD_LIN && method != csInterpolation::METHOD_QUAD && method != csInterpolation::METHOD_SINC ) {
//    throw( csException("csInterpolation::setInterpolationMethod: unknown method index provided (=%d). This is a bug in the calling function", method) );
//  }
//  myMethod = method;
//}
void csInterpolation::setExtrapolation( float valLeft, float valRight ) {
  myExtrapolValLeft  = valLeft;
  myExtrapolValRight = valRight;
}

//--------------------------------------------------------------------------
//
//
void csInterpolation::static_shift( float shift_ms, float const* samplesIn, float* samplesOut ) {
  if( myMethod == csInterpolation::METHOD_SINC ) {
    if( (shift_ms != myPreviousShift_ms) || myIndexBuffer == NULL ) {
      if( myIndexBuffer == NULL ) {
        myIndexBuffer = new float[myNumSamples];
      }
      for( int i = 0; i < myNumSamples; i++ ) {
        myIndexBuffer[i] = (float)(i) - shift_ms/mySampleInt;
      }
      myPreviousShift_ms = shift_ms;
    }
    process( 1.0, 0.0, samplesIn, myIndexBuffer, samplesOut );
  }
  else {
    for( int i = 0; i < myNumSamples; i++ ) {
      float time = (float)i * mySampleInt - shift_ms;
      samplesOut[i] = valueAt( time, samplesIn );
    }

  }

}
//--------------------------------------------------------------------------
//
//
float csInterpolation::valueAt( float time_ms, float const* samplesIn ) {
  switch( myMethod ) {
  case csInterpolation::METHOD_SINC:
    return valueAt_sinc( time_ms, samplesIn );
  case csInterpolation::METHOD_LIN:
    return valueAt_lin( time_ms, samplesIn );
  case csInterpolation::METHOD_QUAD:
    return valueAt_quad( time_ms, samplesIn );
  default:
    return valueAt_sinc( time_ms, samplesIn );
  }
}

float csInterpolation::valueAt_sinc( float time_ms, float const* samplesIn ) {
  float indexValue[1];
  indexValue[0] = time_ms/mySampleInt;

  float sampleOut[1];
  process( 1.0, 0.0, samplesIn, indexValue, sampleOut, 1 );
  return sampleOut[0];
}

float csInterpolation::valueAt_lin( float time_ms, float const* samplesIn ) {
  float amplitude = 0.0;
  double sampleIndex = time_ms / mySampleInt;
  int isamp = (int)sampleIndex;
  int isamp1 = isamp;
  if( isamp < 0 ) isamp1 = 0;
  if( isamp >= myNumSamples-1 ) isamp1 = myNumSamples-2;
  double y1 = samplesIn[isamp1];
  double y2 = samplesIn[isamp1 + 1];
  double weight = ( sampleIndex - (double)isamp ) / 1.0;
  amplitude = (float)( y1 + weight * (y2 - y1) );
  return amplitude;
}
// This works only for monotonically increasing functions using linear interpolation
float csInterpolation::timeAt_lin( float value, float const* samplesIn ) {
  // 1) Search two samples encompassing search value
  int sampLeft  = 0;
  int sampRight = myNumSamples-1;
  do {
    int sampMid = ( sampLeft + sampRight ) / 2;
    float valMid = samplesIn[sampMid];
    if( value >= valMid ) {
      sampLeft = sampMid;
    }
    else if( value <= valMid ) {
      sampRight = sampMid;
    }
    else if( sampRight-sampLeft > 1 ) {
      return -1; // Input sample values are not monotonically increasing
    }
  } while( sampRight-sampLeft > 1 );

  double y1 = samplesIn[sampLeft];
  double y2 = samplesIn[sampRight];
  double y2y1 = y2 - y1;
  float sampleIndex = (float)sampLeft;
  if( abs(y2y1) > 1e-6 ) {
    sampleIndex = ( ( value - y1 ) / y2y1 ) + (float)sampLeft;
  }
  return( sampleIndex * mySampleInt );
}

// This works only for monotonically increasing functions using linear interpolation
void csInterpolation::timeAt_lin( float const* samplesIn, float sampleIntOut, int numSamplesOut, float* timesOut ) {
  // 1) Search two samples encompassing search value
  int sampLeft  = 0;
  int sampRight = myNumSamples-1;
  float valLeft  = samplesIn[sampLeft];
  float valRight = samplesIn[sampRight];
  int samp1 = 0;
  while( samp1 < numSamplesOut && valLeft > (float)samp1 * sampleIntOut ) {
    timesOut[samp1] = sampLeft * mySampleInt;
    samp1 += 1;
  }
  int samp2 = numSamplesOut-1;
  while( samp2 > sampLeft && valRight < (float)samp2 * sampleIntOut ) {
    timesOut[samp2] = sampRight * mySampleInt;
    samp2 -= 1;
  }
  for( int isampOut = samp1; isampOut <= samp2; isampOut++ ) {
    float valueOut = (float)isampOut * sampleIntOut;
    do {
      int sampMid = ( sampLeft + sampRight ) / 2;
      float valMid = samplesIn[sampMid];
      if( valueOut >= valMid ) {
        sampLeft = sampMid;
      }
      else if( valueOut <= valMid ) {
        sampRight = sampMid;
      }
      else if( sampRight-sampLeft > 1 ) {
        return; // Input sample values are not monotonically increasing
      }
    } while( sampRight-sampLeft > 1 );
    double y1 = samplesIn[sampLeft];
    double y2 = samplesIn[sampRight];
    double y2y1 = y2 - y1;
    if( abs(y2y1) > 1e-6 ) {
      timesOut[isampOut] = mySampleInt * ( ( valueOut - y1 ) / y2y1 + (float)sampLeft );
    }
    else {
      timesOut[isampOut] = mySampleInt * (float)sampLeft;
    }
    //    fprintf(stderr,"%d %.1f  %d %d   %.2f\n", isampOut, valueOut, sampLeft, sampRight, timesOut[isampOut] );
    sampRight = myNumSamples-1;
  }
}


float csInterpolation::valueAt_quad( float time_ms, float const* samplesIn ) {
  float amplitude = 0.0;
  double sampleIndex = time_ms / mySampleInt;
  int isamp = (int)sampleIndex;
  int isampMid = isamp;
  if( isamp < 1 ) {
    isampMid = 1;
  }
  else if( isamp >= myNumSamples-2 ) {
    isampMid = myNumSamples-2;
  }

  double y1 = samplesIn[isampMid - 1];
  double y2 = samplesIn[isampMid];
  double y3 = samplesIn[isampMid + 1];

  // Set x values to similar range as y, to avoid rounding errors
  double min = (y1 < y2 ? y1 : y2);  // Math.min( Math.min( y1, y2 ), y3 );
  min = (min < y3 ? min : y3);

  double max = (y1 > y2 ? y1 : y2);  // Math.max( Math.max( y1, y2 ), y3 );
  max = (max > y3 ? max : y3);

  double x1 = (y1 >= 0.0 ? y1 : -y1);  // = Math.abs( y1 );
  double x2 = x1 + ( max - min );
  double x3 = x1 + 2.0 * ( max - min );

  double xAmp = x2;
  if( (double)isamp < sampleIndex ) {
    xAmp = x2 + ( (sampleIndex-(double)isamp) * (x3 - x2) );
  }
  else {
    xAmp = x1 + ( (sampleIndex-(double)(isamp-1)) * (x2 - x1) );
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
  return amplitude;
}

/**
* Estimate amplitude at the given sample index, using quadratic interpolation
* @param traceData   (i) Input data values/trace
* @param sample      (i) Sample index at which amplitude shall be computed (starting at 0)
* @param numSamples  (i) Number of samples in input trace
* @return estimated amplitude
*/
float getQuadAmplitudeAtSample( float const* traceData, double sample, int numSamples );


//--------------------------------------------------------------------------
//
//
void csInterpolation::process( float sampleIntSkew, float xVal1, float const* samplesIn, float const* sIndexOut, float* samplesOut, int numSamplesOut ) {
  int sampleOffset  = 1 - 3*myNumCoefficients/2;
  float sampleRate  = 1.0f / sampleIntSkew;
  float helpIndex   = (float)myNumCoefficients - xVal1 * sampleRate;
  float numValMin   = (float)(myNumValues-1);
  int numSamplesMin = myNumSamples - myNumCoefficients;

  for( int isamp = 0; isamp < numSamplesOut; isamp++ ) {
    float sampleOut   = helpIndex + sIndexOut[isamp] * sampleRate;
    int isampleOut    = (int)sampleOut;
    int currentSample = sampleOffset+isampleOut;
    float const* ptrSample = &samplesIn[currentSample];
    float remainder  = sampleOut-(float)isampleOut;
    int valIndex     = (int)( remainder >= 0.0 ? remainder*numValMin+0.5 : (remainder+1.0)*numValMin-0.5 );
    if( valIndex < 0 || valIndex >= myNumValues ) throw(csException("csInterpolation::process(): Inconsistent value index (%d != [0,%d]). This is probably a program bug in the calling function", valIndex, myNumValues) );
    float* coefRed   = myCoefficients[ valIndex ];

    if( currentSample >= 0 && currentSample <= numSamplesMin ) {
      float sum = 0.0;
      for( int icoef = 0; icoef < myNumCoefficients; icoef++ ) {
        sum += ptrSample[icoef] * coefRed[icoef];
      }
      samplesOut[isamp] = sum;
    }
    else {
      float sum = 0.0;
      float valOut = 0.0;
      for( int icoef = 0; icoef < myNumCoefficients; icoef++,currentSample++ ) {
        if( currentSample < 0 ) {
          valOut = myExtrapolValLeft;
        }
        else if( currentSample >= myNumSamples ) {
          valOut = myExtrapolValRight;
        }
        else {
          valOut = samplesIn[currentSample];
        }
        sum += valOut * coefRed[icoef];
      }
      samplesOut[isamp] = sum;
    }
  }

}

/**
 * Sinc function
 * sinc(x) = sin(pi*x)/(pi*x)
 */
double csInterpolation::sincFunction( double value ) {
  if( value != 0.0 ) {
    double tmp = M_PI * value;
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

bool csInterpolation::toeplitzSolver( int numDimensions, double const* topRow, double const* vecRight, double* vecLeft, double* vecSolve ) {
  double val0 = topRow[0];
  if( val0 == 0.0 ) return false;

  vecLeft[0]  = vecRight[0] / topRow[0];
  vecSolve[0] = 1.0;

  for( int idim = 1; idim < numDimensions; idim++ ) {
    vecSolve[idim] = 0.0;
    vecLeft[idim]  = 0.0;

    double solve1 = 0.0;
    for( int i = 0; i < idim; i++ ) {
      solve1 += vecSolve[i]*topRow[idim-i];
    }

    double solve2 = solve1/val0;
    val0  -= solve2*solve1;
    int idimHalf = idim/2;
    for( int i = 0; i <= idimHalf; i++ ) {
      double tmp = vecSolve[idim-i] - solve2 * vecSolve[i];
      vecSolve[i] -= solve2 * vecSolve[idim-i];
      vecSolve[idim-i] = tmp;
    }

    double solve3 = 0.0;
    for( int i = 0; i < idim; i++ ) {
      solve3 += vecLeft[i] * topRow[idim-i];
    }

    solve2 = ( solve3 - vecRight[idim] ) / val0;
    for( int i = 0; i <= idim; i++ ) {
      vecLeft[i] -= solve2 * vecSolve[idim-i];
    }
  }
  return true;
}

void csInterpolation::createSincFunction( float ratio, int numCoef, float* sincPulse ) {
  double* vecLeft = new double[numCoef];
  double* topRow = new double[numCoef];
  double* vecRight = new double[numCoef];
  double* vecSolve = new double[numCoef];

  double fmax = 0.066 + 0.265 * log((double)numCoef);
  fmax = ( fmax < 1.0 ) ? fmax : 1.0;
  for( int i = 0; i < numCoef; i++ ) {
    topRow[i]   = sincFunction( fmax * i );
    vecRight[i] = sincFunction( fmax * (numCoef/2-i-1+ratio));
  }

  csInterpolation::toeplitzSolver( numCoef, topRow, vecRight, vecLeft, vecSolve );

  for( int i = 0; i < numCoef; i++ ) {
    sincPulse[i] = (float)vecLeft[i];
  }
  delete [] topRow;
  delete [] vecRight;
  delete [] vecLeft;
  delete [] vecSolve;
}

void csInterpolation::xy2yxInterpolation( float const* arrayIn, float* arrayOut, int numSamples ) {
  csInterpolation::xy2yxInterpolation( arrayIn, arrayOut, numSamples, 1.0, -1.0 );
}
void csInterpolation::xy2yxInterpolation( float const* arrayIn, float* arrayOut, int numSamples, float sampleInt ) {
  csInterpolation::xy2yxInterpolation( arrayIn, arrayOut, numSamples, sampleInt, -1.0 );
}
void csInterpolation::xy2yxInterpolation( float const* arrayIn, float* arrayOut, int numSamples, float sampleInt, float valueInvalid ) {
  int isamp1 = (int)ceil(arrayIn[0]/sampleInt);
  int isamp2 = (int)(arrayIn[numSamples-1]/sampleInt);
  if( isamp1 < 0 ) isamp1 = 0;
  if( isamp2 > numSamples-1 ) isamp2 = numSamples-1;
  if( isamp1 >= numSamples ) {
    isamp1 = numSamples-1;
  }
  for( int isamp = 0; isamp <= isamp1; isamp++ ) {
    arrayOut[isamp] = valueInvalid;
  }
  for( int isamp = isamp2; isamp < numSamples; isamp++ ) {
    arrayOut[isamp] = valueInvalid;
  }
  int sIndex = isamp1;
  int counter = 1;
  while( sIndex <= isamp2 ) {
    while( arrayIn[counter] < (float)sIndex*sampleInt && counter < numSamples ) counter += 1;
    if( counter == numSamples ) break;
    float diff = arrayIn[counter] - arrayIn[counter-1];
    if( fabs(diff) > 0.0 ) {
      arrayOut[sIndex] = ( (float)counter - 1.0f + ((float)sIndex*sampleInt - arrayIn[counter-1]) / diff ) * sampleInt;
    }
    else {
      arrayOut[sIndex] = ((float)counter - 1.0f) * sampleInt;
    }
    sIndex += 1;
  }
}

void csInterpolation::linearInterpolation( int numSamplesIn, float const* xin, float const* yin, int numSamplesOut, float sampleIntOut, float* yout )
{
  float x0 = xin[0];
  float y0 = yin[0];
  float xN = xin[numSamplesIn-1];
  float yN = yin[numSamplesIn-1];

  int sample0 = MIN( (int)(x0 / sampleIntOut), numSamplesOut-1 );
  int sampleN = (int)(xN / sampleIntOut) + 1;
  if( sampleN > numSamplesOut ) sampleN = numSamplesOut;

  for( int isamp = 0; isamp <= sample0; isamp++ ) {
    yout[isamp] = y0;
  }
  for( int isamp = sampleN; isamp < numSamplesOut; isamp++ ) {
    yout[isamp] = yN;
  }
  
  float x1  = x0;
  float x2  = numSamplesIn > 1 ? xin[1] : x0;
  float y1  = y0;
  float y2  = numSamplesIn > 1 ? yin[1] : y0;
  float dx = x2-x1;
  float dy = y2-y1;
  float dydx = dx > 0.0f ? dy/dx : 0.0f;
  int counter = 0;
  for( int isamp = sample0+1; isamp < sampleN; isamp++ ) {
    float x = (float)isamp*sampleIntOut;
    while( x > x2 ) {
      x1 = x2;
      y1 = y2;
      counter += 1;
      if( counter == numSamplesIn ) {
        x2 = x; // Avoids entering this loop again
        break;
      }
      x2 = xin[counter];
      y2 = yin[counter];
      dx = x2-x1;
      dy = y2-y1;
      dydx = dx > 0.0f ? dy/dx : 0.0f;
    }
    yout[isamp] = y1 + (x - x1) * dydx;
  }
}
