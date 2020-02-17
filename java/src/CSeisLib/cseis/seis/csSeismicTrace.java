/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seis;

/**
 * Seismic trace.
 * @author Bjorn Olofsson
 */
public class csSeismicTrace {
  private csHeader[] myHeaderValues;
  private csSeismicData myData;
  private int myNumHeaders;
  private int myOriginalTraceNumber;

  public csSeismicTrace( int nSamples, int nHeaders ) {
    myNumHeaders = nHeaders;
    myHeaderValues  = new csHeader[myNumHeaders];
    for( int i = 0; i < myNumHeaders; i++ ) {
      myHeaderValues[i] = new csHeader();
    }
    myData      = new csSeismicData( nSamples );
    myOriginalTraceNumber = 1;
  }
  public csSeismicTrace( float[] samples, csHeader[] headers ) {
    this( samples, headers, 1 );
  }
  public csSeismicTrace( float[] samples, csHeader[] headers, int originalTraceNumber ) {
    myNumHeaders   = headers.length;
    myHeaderValues = headers;
    myData = new csSeismicData( samples );
    myOriginalTraceNumber = originalTraceNumber;
  }
  public csSeismicTrace( csSeismicData data, csHeader[] headers, int originalTraceNumber ) {
    myNumHeaders   = headers.length;
    myHeaderValues = headers;
    myData = new csSeismicData( data );
    myOriginalTraceNumber = originalTraceNumber;
  }
  public csSeismicData data() {
    return myData;
  }
  public float[] samples() {
    return myData.samples();
  }
  public void setSamples( float[] samples ) {
    myData.setSamples(samples);
  }
  public csHeader[] headerValues() {
    return myHeaderValues;
  }
  public int numSamples() {
    return myData.numSamples();
  }
  public int numHeaders() {
    return myNumHeaders;
  }
  /**
   * Get original trace number (starting at 1 for first trace)
   * @return
   */
  public int originalTraceNumber() {
    return myOriginalTraceNumber;
  }
  /**
   * Set original trace number (starting at 1 for first trace)
   * @param traceNumber
   */
  public void setOriginalTraceNumber( int traceNumber ) {
    myOriginalTraceNumber = traceNumber;
  }
/*
  public void setHeader( int index, double value ) {
    if( index >= 0 && index < myNumHeaders ) {
      myHeaderValues[index].setValue( value );
    }
  }
  */
  
    public void setDoubleHeader( int index, double value ) {
    if( index >= 0 && index < myNumHeaders ) {
      myHeaderValues[index].setValue( value );
    }
  }
  public void setLongHeader( int index, long value ) {
    if( index >= 0 && index < myNumHeaders ) {
      myHeaderValues[index].setValue( value );
    }
  }
  public void setFloatHeader( int index, float value ) {
    if( index >= 0 && index < myNumHeaders ) {
      myHeaderValues[index].setValue( value );
    }
  }
  public void setIntHeader( int index, int value ) {
    if( index >= 0 && index < myNumHeaders ) {
      myHeaderValues[index].setValue( value );
    }
  }
  public void setStringHeader( int index, String value ) {
    if( index >= 0 && index < myNumHeaders ) {
      myHeaderValues[index].setValue( value );
    }
  }
 
  
  /**
   * 
   * @param traceIndex
   * @return Maximum amplitude in specified trace
   */
  public float maxAmplitude() {
    return myData.maxAmplitude();
  }
  /**
   * 
   * @param traceIndex
   * @return Mean (absolute) amplitude of specified trace
   */
  public float meanAmplitude() {
    return myData.meanAmplitude();
  }
  /**
   * 
   * @param traceIndex
   * @return DC amplitude (= mean of signed amplitude) of specified trace
   */
  public float dcAmplitude() {
    return myData.dcAmplitude();
  }
  /**
   * 
   * @param traceIndex
   * @return Minimum amplitude in specified trace
   */
  public float minAmplitude() {
    return myData.minAmplitude();
  }
  
  public static void main( String[] args ) {
    csSeismicTrace t = new csSeismicTrace(100,20);
    t.setIntHeader( 0, 123 );
  }
}


