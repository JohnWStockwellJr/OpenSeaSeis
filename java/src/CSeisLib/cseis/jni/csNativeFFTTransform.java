/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

/**
 * JNI interface to FFT transform (C++).
 * @author Bjorn Olofsson
 */
public class csNativeFFTTransform {
  public static final int FORWARD = 1;
  public static final int INVERSE = 2;

  private native void native_performFFT(
      int direction,
      double[] samples,
      int numSamples );

  private int myNumSamplesOrig;
  private int myNumFFTSamples;

  public csNativeFFTTransform( int numSamples ) {
    myNumSamplesOrig = numSamples;
    myNumFFTSamples  = myNumSamplesOrig;
    int twoPower = 0;
    int ns = myNumSamplesOrig;
    while( (ns = (int)(ns / 2)) > 0 ) {
      twoPower += 1;
    }
    int two_power_m = (int)Math.pow( 2.0, twoPower );
    
    Powerof2( myNumSamplesOrig, twoPower, two_power_m );
    if( two_power_m != ns ) {
      myNumFFTSamples = two_power_m * 2;
      twoPower = twoPower + 1;
    }

  }
  /**
   * Make sure number of input samples is power of 2!
   * @param samples
   * @param direction
   * @param amplitude
   * @param phase
   */
  public csFFTObject performFFT( float[] samples, int direction, float sampleInt_ms ) {
    if( samples.length != myNumSamplesOrig ) {
      return null;
    }
//    vars->fftTool->lowPass( samplesInOut, vars->order, vars->cutoffLow, (vars->output == OUTPUT_IMPULSE) );
    int numSamplesHalf = myNumFFTSamples/2;
    csFFTObject fftObject = new csFFTObject( numSamplesHalf );
    double[] samplesDouble = new double[myNumFFTSamples];
    for( int i = 0; i < myNumSamplesOrig; i++ ) {
      samplesDouble[i] = samples[i];
    }
    for( int i = myNumSamplesOrig; i < myNumFFTSamples; i++ ) {
      samplesDouble[i] = 0.0;
    }
    native_performFFT( direction, samplesDouble, myNumFFTSamples );
    System.arraycopy( samplesDouble, 0, fftObject.amplitude, 0, numSamplesHalf );
    System.arraycopy( samplesDouble, numSamplesHalf, fftObject.phase, 0, numSamplesHalf );

    fftObject.freqInc = 1000.0/(myNumFFTSamples*sampleInt_ms);

    return fftObject;
  }
  public static boolean Powerof2( int nx, Integer m, Integer twopm ) {
    int value = nx;
    m = 0;
    while( (value = (int)(value / 2)) > 0 ) {
      m += 1;
    }
    twopm = (int)Math.pow( 2.0, m );
    return( twopm == nx );
  }
}


