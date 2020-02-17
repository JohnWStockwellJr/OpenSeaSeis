/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

/**
 * FFT object.
 * @author 2009 Bjorn Olofsson
 */
public class csFFTObject {
  public csFFTObject( int numSamples_in ) {
    numSamples = numSamples_in;
    amplitude = new double [numSamples];
    phase     = new double [numSamples];
  }
  public int numSamples;
  public double[] amplitude;
  public double[] phase;
  public double freqInc;
}


