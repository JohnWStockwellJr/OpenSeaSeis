/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <algorithm>
#include "geolib_methods.h"
#include "cseis_jni_csNativeFFTTransform.h"
//#include "csFFT.h"

bool Powerof2( int nx, int* m, int* twopm);
bool fft(int dir,int m,double *x,double *y);

using namespace std;

/**
 * JNI Native FFT transform
 *
 * @author Bjorn Olofsson
 * @date 2007
 */

/*
 * Class:     cseis_jni_csNativeFFTTransform
 * Method:    native_performFFT
 * Signature: (I[DI)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeFFTTransform_native_1performFFT
(JNIEnv *env, jobject obj, jint direction, jdoubleArray samples, jint numSamples )
{
  /*
  int fftDataType = cseis_geolib::FX_AMP;
  cseis_geolib::csFFT fft( numSamples );

  float* bufferInOut = new float[numSamples];
  env->GetDoubleArrayRegion( samples, 0, numSamples, bufferInOut );
  fft.forwardTransform( bufferInOut, bufferInOut, fftDataType );
  */
    
  int twoPower;
  int two_power_m;
  Powerof2( numSamples, &twoPower, &two_power_m );
  if( two_power_m != numSamples ) {
    fprintf(stderr,"Wrong number of samples input. Must be power of 2\n");
    return;
  }

  double* bufferReal = new double[numSamples];
  double* bufferImag = new double[numSamples];

  env->GetDoubleArrayRegion( samples, 0, numSamples, bufferReal );
  for( int i = 0; i < numSamples; i++ ) {
    bufferImag[i] = 0;
  }

  if( !fft( direction, twoPower, bufferReal, bufferImag ) ) {
    fprintf(stderr,"FFT transform failed for unknown reasons...\n");
    delete [] bufferReal;
    delete [] bufferImag;
    return;
  }

  // Convert to amplitude and phase
  for( int i = 1; i < numSamples/2; i++ ) {
    double amp  = sqrt( bufferReal[i]*bufferReal[i] + bufferImag[i]*bufferImag[i] );
    double ph   = atan2( -bufferImag[i], bufferReal[i] );
    bufferReal[i] = amp;
    bufferImag[i] = ph;
  }
  bufferReal[0] = fabs( bufferReal[0] );
  bufferImag[0] = 0.0;

  env->SetDoubleArrayRegion( samples, 0, numSamples/2, bufferReal );
  env->SetDoubleArrayRegion( samples, numSamples/2, numSamples/2, bufferImag );

  delete [] bufferReal;
  delete [] bufferImag;
}

