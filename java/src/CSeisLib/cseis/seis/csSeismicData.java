/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seis;

/**
 * Seismic data.
 * @author Bjorn Olofsson
 */
public class csSeismicData {
  private float[]  mySampleValues;
  private int myNumSamples;
  protected float myMeanAmplitude;     // Mean absolute amplitude = sum( abs(sample[i]) )/nsamp
  protected float myDCAmplitude;       // Average amplitude = DC = sum( sample[i] )/nsamp
  protected float myMaxAmplitude;      // Maximum amplitude
  protected float myMinAmplitude;      // Minimum amplitude
  protected boolean myAreAmplitudesComputed;

  public csSeismicData() {
    mySampleValues = null;
    myNumSamples   = 0;
    myMeanAmplitude = 0;
    myDCAmplitude   = 0;
    myMaxAmplitude  = 0;
    myMinAmplitude  = 0;
    myAreAmplitudesComputed = false;
  }
  public csSeismicData( int nSamples ) {
    this();
    mySampleValues = new float[nSamples];
    myNumSamples   = nSamples;
  }
  public csSeismicData( float[] samples ) {
    this();
    mySampleValues = samples;
    myNumSamples   = samples.length;
  }
  public csSeismicData( csSeismicData data ) {
    mySampleValues  = data.mySampleValues;
    myNumSamples    = data.myNumSamples;
    myMeanAmplitude = data.myMeanAmplitude;
    myDCAmplitude   = data.myDCAmplitude;
    myMaxAmplitude  = data.myMaxAmplitude;
    myMinAmplitude  = data.myMinAmplitude;
    myAreAmplitudesComputed = data.myAreAmplitudesComputed;
  }
  public int numSamples() {
    return myNumSamples;
  }
  public float[] samples() {
    return mySampleValues;
  }
  public void setSamples( float[] samples ) {
    mySampleValues = samples;
    myNumSamples   = mySampleValues.length;
    myAreAmplitudesComputed = false;
  }
  /**
   * 
   * @param traceIndex
   * @return Maximum amplitude in specified trace
   */
  public float maxAmplitude() {
    if( myAreAmplitudesComputed ) return myMaxAmplitude;
    computeTraceAmplitudes();
    return myMaxAmplitude;
  }
  /**
   * 
   * @param traceIndex
   * @return Mean (absolute) amplitude of specified trace
   */
  public float meanAmplitude() {
    if( myAreAmplitudesComputed ) return myMeanAmplitude;
    computeTraceAmplitudes();
    return myMeanAmplitude;
  }
  /**
   * 
   * @param traceIndex
   * @return DC amplitude (= mean of signed amplitude) of specified trace
   */
  public float dcAmplitude() {
    if( myAreAmplitudesComputed ) return myDCAmplitude;
    computeTraceAmplitudes();
    return myDCAmplitude;
  }
  /**
   * 
   * @param traceIndex
   * @return Minimum amplitude in specified trace
   */
  public float minAmplitude() {
    if( myAreAmplitudesComputed ) return myMinAmplitude;
    computeTraceAmplitudes();
    return myMinAmplitude;
  }
  private void computeTraceAmplitudes() {
    myMeanAmplitude = 0;
    myDCAmplitude   = 0;
    myMaxAmplitude  = 0;
    myMinAmplitude  = 0;

    float sum = 0.0f;
    float sum_abs = 0.0f;
    float min = Float.MAX_VALUE;
    float max = -Float.MAX_VALUE;
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      float amp = mySampleValues[isamp];
      if( amp > max ) {
        max = amp;
      }
      if( amp < min ) {
        min = amp;
      }
      if( amp >= 0 ) sum_abs += amp;
      else sum_abs -= amp;
      sum += amp;
    }
    myMeanAmplitude = sum_abs / myNumSamples;
    myDCAmplitude   = sum / myNumSamples;
    myMaxAmplitude  = max;
    myMinAmplitude  = min;
    myAreAmplitudesComputed = true;
  }
}


