/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.jni;

/**
 * JNI interface to Filter (C++).
 * @author Felipe Punto
 */
public class csNativeFilter {
  private native long native_createInstance( int numSamples, float sampleInt, float freqCutOffLowPass, float slopeLowPass_db, float freqCutOffHighPass, float slopeHighPass_db );
  private native void native_performFilter( long nativePtr, float[] samplesIn, float[] samplesOut );
  private native void native_freeInstance( long nativePtr );
  protected long myNativePtr = 0;
//  private static native void native_createFilterWavelet(
//    float cutOffFreqHz,
//    float slope_dbPerOct,
//    float sampleInt_ms,
//    int filterLength_samples,
//    double[] filterSamples );

  public csNativeFilter( int numSamples, float sampleInt_ms, float freqLowPass, float slopeLowPass, float freqHighPass, float slopeHighPass ) {
    myNativePtr = native_createInstance( numSamples, sampleInt_ms, freqLowPass, slopeLowPass, freqHighPass, slopeHighPass );
    if( myNativePtr == 0 ) {
      // ..should never happen
    }
  }
  public void performFilter( float[] samplesIn, float[] samplesOut ) {
    if( myNativePtr != 0 ) native_performFilter( myNativePtr, samplesIn, samplesOut );
  }
  @Override
  public void finalize() {
    if( myNativePtr != 0 ) {
      native_freeInstance( myNativePtr );
      myNativePtr = 0;
    }
  }
//  public static csFilterObject createFilterWavelet( float cutOffFreqHz, float slope_dbPerOct,
//          float sampleInt_ms, int filterLength_samples ) {
//    double[] filterSamples = new double[filterLength_samples];
//    native_createFilterWavelet( cutOffFreqHz, slope_dbPerOct, sampleInt_ms, filterLength_samples,
//            filterSamples );
//    csFilterObject filterObject = new csFilterObject();
//    filterObject.samples = filterSamples;
//    return filterObject;
//  }
}

