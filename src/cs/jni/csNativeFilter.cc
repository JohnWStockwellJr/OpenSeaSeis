/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <algorithm>
#include "geolib_methods.h"
#include "cseis_jni_csNativeFilter.h"
#include "csFilterTool.h"

using namespace std;

/**
 * JNI Native Filter
 *
 * @author Felipe Punto
 * @date 2015
 */

/*
 * Class:     cseis_jni_csNativeFilter
 * Method:    native_createInstance
 * Signature: (IFFFFF)J
 */
JNIEXPORT jlong JNICALL Java_cseis_jni_csNativeFilter_native_1createInstance
(JNIEnv *env, jobject obj, jint numSamples_in, jfloat sampleInt_in, jfloat freqLowPass_in, jfloat slopeLowPass_in, jfloat freqHighPass_in, jfloat slopeHighPass_in ) {
  cseis_jni::csFilterTool* filterTool = new cseis_jni::csFilterTool( numSamples_in, sampleInt_in );
  filterTool->setParam( freqLowPass_in, slopeLowPass_in, freqHighPass_in, slopeHighPass_in );
  return( reinterpret_cast<jlong>( filterTool ) );
}

/*
 * Class:     cseis_jni_csNativeFilter
 * Method:    native_performFilter
 * Signature: (J[F[F)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeFilter_native_1performFilter
( JNIEnv *env, jobject obj, jlong ptr_in, jfloatArray samples_in, jfloatArray samples_out ) {
  cseis_jni::csFilterTool* filterTool = reinterpret_cast<cseis_jni::csFilterTool*>(ptr_in);

  if( filterTool->isFilter() ) {
    float* samples_c = filterTool->retrieveSamplesPointer();
    env->GetFloatArrayRegion( samples_in, 0, filterTool->numInputSamples(), samples_c );

    filterTool->applyFilter();

    env->SetFloatArrayRegion( samples_out, 0, filterTool->numInputSamples(), samples_c );
  }
}

/*
 * Class:     cseis_jni_csNativeFilter
 * Method:    native_freeInstance
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeFilter_native_1freeInstance
(JNIEnv *env, jobject obj, jlong ptr_in) {
  cseis_jni::csFilterTool* ptr = reinterpret_cast<cseis_jni::csFilterTool*>(ptr_in);
  delete ptr;
}
