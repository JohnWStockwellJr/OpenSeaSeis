/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <cstdio>
#include <limits>
#include <cstring>
#include "geolib_math.h"
#include "geolib_methods.h"
#include "csSort.h"

#include <algorithm>

namespace cseis_geolib {

//---------------------------------------------
  int compute_correlation_length( int maxlag ) {
    return( 2*maxlag+1 );
  }

  /*
   *
   * auto[f(i)] = a[tau] = a[dt * i] = 1 / (m - i + 1) * SUM_j=0-m[ f(j) * f(i+j) ]
   */
  void compute_onesided_auto_correlation( float const* samples, int nSampIn, float* autocorr ) {
    compute_onesided_auto_correlation( samples, nSampIn, autocorr, nSampIn-1, false );
  }
  void compute_onesided_auto_correlation( float const* samples, int nSampIn, float* autocorr, int maxlag_in_num_samples ) {
    compute_onesided_auto_correlation( samples, nSampIn, autocorr, maxlag_in_num_samples, false );
  }
  void compute_onesided_auto_correlation( float const* samples, int nSampIn, float* autocorr, int maxlag_in_num_samples, bool dampen ) {
    compute_onesided_correlation( samples, samples, nSampIn, autocorr, maxlag_in_num_samples, dampen );
  }
  void compute_onesided_correlation( float const* samplesLeft, float const* samplesRight, int nSampIn, float* autocorr, int maxlag_in_num_samples, bool dampen ) {
    for( int ilag = 0; ilag <= maxlag_in_num_samples; ilag++ ) {
      float sum = 0;
      int sampEnd = nSampIn-ilag;
      for( int isamp = 0; isamp < sampEnd; isamp++ ) {
        sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
      }
      if( !dampen ) {
        autocorr[ilag] = sum;
      }
      else {
        autocorr[ilag] = (float)sampEnd * sum / (float)nSampIn;
      }
    }
  }
  
  void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight, int nSampIn, float* corr ) {
    compute_twosided_correlation( samplesLeft, samplesRight, nSampIn, corr, nSampIn-1, false );
  }
  void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight, int nSampIn, float* corr, int maxlag_in_num_samples ) {
    compute_twosided_correlation( samplesLeft, samplesRight, nSampIn, corr, maxlag_in_num_samples, false );
  }
  void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight, int nSampIn, float* corr, int maxlag_in_num_samples, bool dampen ) {

    if( !dampen ) {
      //---------------------------------------
      // Compute negative lags
      //
      int sampEnd    = nSampIn;
      for( int ilag = -maxlag_in_num_samples; ilag < 0; ilag++ ) {
        int sampStart  = -ilag;
        float sum = 0;
        for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
          sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
        }
        corr[ilag+maxlag_in_num_samples] = sum;
      }
      int sampStart = 0;
      //---------------------------------------
      // Compute positive lags
      //
      for( int ilag = 0; ilag <= maxlag_in_num_samples; ilag++ ) {
        int sampEnd    = nSampIn-ilag;
        float sum = 0;
        for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
          sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
        }
        corr[ilag+maxlag_in_num_samples] = sum;
      }
    }
    else { // dampen
      //---------------------------------------
      // Compute negative lags
      //
      int sampEnd    = nSampIn;
      for( int ilag = -maxlag_in_num_samples; ilag < 0; ilag++ ) {
        int sampStart  = -ilag;
        float sum = 0;
        for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
          sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
        }
        int nSamp = sampEnd-sampStart;
        corr[ilag+maxlag_in_num_samples] = (float)nSamp * sum / (float)nSampIn;
      }
      int sampStart  = 0;
      //---------------------------------------
      // Compute positive lags
      //
      for( int ilag = 0; ilag <= maxlag_in_num_samples; ilag++ ) {
        int sampEnd    = nSampIn-ilag;
        float sum = 0;
        for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
          sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
        }
        int nSamp = sampEnd-sampStart;
        corr[ilag+maxlag_in_num_samples] = (float)nSamp * sum / (float)nSampIn;
      }
    }
  }

  void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight,
                                     int nSampIn, float* corr, int maxlag_in_num_samples, float* sampleIndex_maxAmp, float* maxAmp ) {

    compute_twosided_correlation( samplesLeft, samplesRight, nSampIn, corr, maxlag_in_num_samples );

    // Determine maximum cross-correlation lag time & amplitude
    int nSampCorr = maxlag_in_num_samples*2 + 1;
    int sampleIndex_int = 0;
    *maxAmp = corr[sampleIndex_int];
    for( int isamp = 0; isamp < nSampCorr; isamp++ ) {
      if( corr[isamp] > *maxAmp ) {
        *maxAmp = corr[isamp];
        sampleIndex_int = isamp;
      }
    }
    *sampleIndex_maxAmp = getQuadMaxSample( corr, sampleIndex_int, nSampCorr, maxAmp );
  }

  float compute_correlation_coefficient( float const* samplesLeft, float const* samplesRight, int nSampIn ) {
    float x_mean = 0;
    float y_mean = 0;
    float x2_sum = 0;
    float y2_sum = 0;
    float xy_sum = 0;
    for( int isamp = 0; isamp < nSampIn; isamp++ ) {
      x_mean += samplesLeft[isamp];
      y_mean += samplesRight[isamp];
      x2_sum += samplesLeft[isamp]  * samplesLeft[isamp];
      y2_sum += samplesRight[isamp] * samplesRight[isamp];
      xy_sum += samplesLeft[isamp]  * samplesRight[isamp];
    }
    x_mean /= (float)nSampIn;
    y_mean /= (float)nSampIn;

    float sum1 = xy_sum - (float)nSampIn * x_mean * y_mean;;
    float sum2 = x2_sum - (float)nSampIn * x_mean * x_mean;;
    float sum3 = y2_sum - (float)nSampIn * y_mean * y_mean;;

    float tmp = sum2 * sum3;
    if( tmp <= 0 ) return 0;
    return(  sum1 / sqrt( tmp ) );
  }

  float compute_rms( float const* samples, int nSamples ) {
    if( nSamples > 0 ) {
      double sum_sqr = 0.0;
      for( int isamp = 0; isamp < nSamples; isamp++ ) {
        sum_sqr += (double)samples[isamp]*(double)samples[isamp];
      }
      return( sqrt( (float)sum_sqr/(float)nSamples ) );
    }
    else {
      return 0.0;
    }
  }

  bool factor_2357( int valIn, int& valOut ) {
    valOut = cseis_geolib::factor_2357( valIn );
    return( valOut != -1 );
  }
  int factor_2357( int valIn ) {
    if( valIn < 0 || valIn > (std::numeric_limits<int>::max()-1) ) {
      return -1; // Fail if input bad
    }
    if( valIn <= 1 ) {
      return 1;
    } // Skip the loop for the easy one. 
    
    // Factor-out small primes from the input number. It should reduce to 1 once 
    // it is a perfect factor of 2,3,5 and 7. If not, increment the input and repeat.
    // Fail if the number of iterations or the test value gets to large. 
    int test    = 0;
    int counter = 0;
    int n2,n3,n5,n7;
    while( test != 1 && counter < std::numeric_limits<int>::max() ){
      test = valIn + counter;
      if ( test == std::numeric_limits<int>::max() ) break;
      counter++;
      
      n2 = n3 = n5 = n7 = 0;
      while( test%7 == 0 ){ test = test/7; n7++; }
      while( test%5 == 0 ){ test = test/5; n5++; }
      while( test%3 == 0 ){ test = test/3; n3++; }
      while( test%2 == 0 ){ test = test/2; n2++; }
    }
    if( test != 1 ) return -1;
    int valOut = (int)round( pow(2.0,n2) * pow(3.0,n3) * pow(5.0,n5) * pow(7.0,n7) );
    
    return valOut;
  }

//---------------------------------------

  bool is_nan_inf( float value ) {
    if( value == +1.0f/0.0f || value == -1.0f/0.0f ) return true; // +/-inf
    if( value != value ) return true;      // nan
    if( value == +0.0f/0.0f ) return true; // nan (in case previous line did not work)
    return false;
  }

  
  void meanFilter( double const* valuesIn, double* valuesOut, int numValues, int numValuesFilter ) {
    int numValuesFilterHalf = numValuesFilter/2;
    numValuesFilter = 2*numValuesFilterHalf + 1;

    for( int i1 = 0; i1 < numValues; i1++ ) {
      int start = std::max( 0, i1-numValuesFilterHalf );
      int end   = std::min( numValues-1, i1+numValuesFilterHalf );
      double sum = 0.0;
      int numSummed = end-start+1;
      for( int i2 = start; i2 <= end; i2++ ) {
        sum += valuesIn[i2];
      }
      // Mitigate edge effects - left-hand side
      int numValuesEdge = numValuesFilterHalf-i1;
      if( numValuesEdge > 0 ) {
        sum += numValuesEdge*valuesIn[0];
        numSummed += numValuesEdge;
      }
      // Mitigate edge effects - right-hand side
      numValuesEdge = (i1+numValuesFilterHalf)-(numValues-1);
      if( numValuesEdge > 0 ) {
        sum += numValuesEdge*valuesIn[numValues-1];
        numSummed += numValuesEdge;
      }
      double weight = 2.0*(double)numSummed/(double)numValuesFilter - 1.0;
      sum = (sum / (double)numSummed) * weight;
      valuesOut[i1] = sum;
    }
    
  }

  void medianFilter( double const* valuesIn, double* valuesOut, int numValues, int numValuesFilterALL, bool applyAtEdges ) {
    if( numValues <= 0 || numValuesFilterALL < 2 ) return;
    int numValuesFilter = numValuesFilterALL / 2;

    double* buffer_filter = new double[numValues];
    cseis_geolib::csSort<double> sortObj;

    /*
    int val1_out = 0;
    int val2_out = numValues-1;
    if( !applyAtEdges ) {
      val1_out = numValuesFilter;
      val2_out = numValues - numValuesFilter - 1;
    }
    */
    //--------------------------------------------------------------------------------
    for( int ival = numValuesFilter; ival < numValues-numValuesFilter; ival++ ) {
      int index1 = ival - numValuesFilter;
      memcpy( buffer_filter, &valuesIn[index1], numValuesFilterALL*sizeof(double) );
      sortObj.treeSort( buffer_filter, numValuesFilterALL );
      valuesOut[ival] = buffer_filter[numValuesFilter];
    }
    // Apply filter at edges:
    if( applyAtEdges ) {
      for( int ival = 0; ival < numValuesFilter; ival++ ) {
        int nval2 = std::max( ival + numValuesFilter, 3 );
        memcpy( buffer_filter, valuesIn, nval2*sizeof(double) );
        sortObj.treeSort( buffer_filter, nval2 );
        valuesOut[ival] = buffer_filter[nval2/2];
      }
      for( int ival = std::max(0,numValues-numValuesFilter); ival < numValues; ival++ ) {
        int index1 = std::max( 0, ival - numValuesFilter );
        int nval2  = numValues - index1 - 1;
        memcpy( buffer_filter, &valuesIn[index1], nval2*sizeof(double) );
        sortObj.treeSort( buffer_filter, nval2 );
        valuesOut[ival] = buffer_filter[nval2/2];
      }
    }
    delete [] buffer_filter;
  }

  //--------------------------------------------------------------------------------
  void padTraceCosineTaper( float const* samplesIn, float* samplesOut, int numSamplesIn, int numSampPad ) {
    int numSamplesInclPad = numSamplesIn + 2*numSampPad;
    // Padding: Pad samples and extrapolate data samples. Apply cosine taper to padded data.
    for( int isamp = 0; isamp < numSamplesInclPad; isamp++ ) {
      samplesOut[isamp] = 0;
    }
    memcpy( &samplesOut[numSampPad], samplesIn, numSamplesIn * sizeof(float) );
  
    // Extrapolate start of trace
    double value2x = 2 * samplesIn[0];
    for( int isamp = numSampPad-1; isamp >= 0; isamp-- ) {
      double ratio = (double)isamp/(double)(numSampPad-1);
      double taper = 0.5 * ( 1 + cos( M_PI*(ratio-1.0) ) );
      samplesOut[isamp] = taper * ( value2x - samplesIn[numSampPad-isamp] );
    }
    // Extrapolate end of trace
    value2x = 2 * samplesIn[numSamplesIn - 1];
    int num = 2 * numSamplesIn + numSampPad-2;
    for( int isamp = numSamplesIn + numSampPad; isamp < numSamplesInclPad; isamp++ ) {
      double ratio = (double)(numSamplesInclPad-isamp-1)/(double)(numSampPad-1);
      double taper = 0.5 * ( 1 + cos( M_PI*(ratio-1.0) ) );
      samplesOut[isamp] = taper * ( value2x - samplesIn[num-isamp] );
    }
  }

  bool is_inf( float value ) {
    if( value == +1.0f/0.0f || value == -1.0f/0.0f ) return true; // +/-inf
    return false;
  }
  bool is_nan( float value ) {
    if( value != value ) return true;      // nan
    if( value == +0.0f/0.0f ) return true; // nan (in case previous line did not work)
    return false;
  }

}

