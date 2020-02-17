/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

import cseis.seis.csHeader;

/**
 * Bundles methods and fields needed across JNI interface in order to perform trace header 'selection'.<br>
 * Trace header 'selection' refers to a way to scan and select traces with certain trace header values
 * in a seismic data set.
 * @author 2013 Felipe Punto
 */
public class csSelectedHeaderBundle {
  /// Value of each header
  public csHeader[] hdrValues;
  public int[] traceIndexList;
  public csSelectedHeaderBundle( int numTraces ) {
    hdrValues = new csHeader[numTraces];
    traceIndexList = new int[numTraces];
    for( int i = 0; i < numTraces; i++ ) {
      hdrValues[i] = new csHeader();
      traceIndexList[i] = i;
    }
  }
  public void setSortedTraceIndex( int traceIndex, int sortedTraceIndex ) {
    traceIndexList[traceIndex] = sortedTraceIndex;
  }
  public void setIntValue( int traceIndex, int value ) {
    hdrValues[traceIndex].setValue(value);
  }
  public void setDoubleValue( int traceIndex, double value ) {
    hdrValues[traceIndex].setValue(value);
  }
  public void setFloatValue( int traceIndex, float value ) {
    hdrValues[traceIndex].setValue(value);
  }
  public void setLongValue( int traceIndex, long value ) {
    hdrValues[traceIndex].setValue(value);
  }
}

