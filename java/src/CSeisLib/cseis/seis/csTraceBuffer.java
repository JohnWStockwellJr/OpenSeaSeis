/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seis;

import cseis.processing.csIProcessing;
import java.util.*;

/**
 * Trace buffer.<br>
 * Stores any number of seismic traces.
 * @author Bjorn Olofsson
 */
public class csTraceBuffer implements csISeismicTraceBuffer {
  protected int myNumSamples;
  protected int myNumHeaders;
  protected ArrayList<csSeismicTrace> myTraces;
  protected csDataBuffer myProcessedDataBuffer;
  protected boolean myIsProcessed;
  
  protected boolean myAreAmplitudesComputed;
  protected float myTotalMinAmplitude;
  protected float myTotalMaxAmplitude;
  protected float myTotalMeanAmplitude;

  /**
   * Constructor.
   * 
   * @param numSamples Number of samples
   * @param numHeaders Number of headers
   */
  public csTraceBuffer( int numSamples, int numHeaders ) {
    myTraces = new ArrayList<csSeismicTrace>();
    myNumSamples = numSamples;
    myNumHeaders = numHeaders;
    myAreAmplitudesComputed = false;
    myProcessedDataBuffer = new csDataBuffer();
    myIsProcessed = false;
  }
  /**
   * Add trace to trace buffer
   * 
   * @param samples Sample values (size of array must be the same as number of samples specified in constructor)
   * @param headers Header values (size of array must be the same as number of headers specified in constructor)
   */
  public void addTrace( float[] samples, csHeader[] headers ) {
    addTrace( new csSeismicTrace( samples, headers ) );
  }
  /**
   * Add trace to trace buffer
   *
   * @param trace Seismic trace to add
   */
  public void addTrace( csSeismicTrace trace ) {
    myAreAmplitudesComputed = false;
    myIsProcessed = false;
    myTraces.add( trace );
  }
  /**
   * Add trace to trace buffer
   *
   * @param trace Seismic trace to add
   */
  public void addTraceAtFront( csSeismicTrace trace ) {
    myAreAmplitudesComputed = false;
    myIsProcessed = false;
    myTraces.add( 0, trace );
  }
  /**
   * Add trace to beginning of trace buffer
   *
   * @param trace Seismic trace to add
   */
  public void addTraceAtStart( csSeismicTrace trace ) {
    myAreAmplitudesComputed = false;
    myIsProcessed = false;
    myTraces.add( 0, trace );
  }
  /**
   * Clear all processing steps
   */
  public void clearProcessing() {
    myIsProcessed = false;
    myProcessedDataBuffer.clear();
  }
  /**
   * @return true if data traces have some type of processing applied
   */
  public boolean isProcessed() {
    return myIsProcessed;
  }
  /**
   * Apply processing steps to data
   * @param processingSteps List of processing steps to go through
   */
  public void applyProcessing( java.util.List<csIProcessing> processingSteps ) {
    clearProcessing();
    if( processingSteps.size() > 0 ) {
      for( int i = 0; i < processingSteps.size(); i++ ) {
        csIProcessing proc = processingSteps.get(i);
        if( proc.isActive() ) proc.apply( this, myProcessedDataBuffer );
      }
      myIsProcessed = true;
    }
  }
  /**
   * Apply processing steps to data
   * @param processingSteps List of processing steps to go through
   */
  public void applyProcessing( csIProcessing processingStep ) {
    clearProcessing();
    processingStep.apply( this, myProcessedDataBuffer );
    myIsProcessed = true;
  }
  /**
   * 
   * @return Number of traces in trace buffer
   */
  @Override
  public int numTraces() {
    return myTraces.size();
  }
  /**
   * 
   * @return Number of samples in each trace
   */
  @Override
  public int numSamples() {
    return myNumSamples;
  }
  /**
   * 
   * @return Number of headers in each trace
   */
  @Override
  public int numHeaders() {
    return myNumHeaders;
  }
  /**
   * 
   * @param traceIndex
   * @return Handle to sample values of specified trace
   */
  @Override
  public float[] samples( int traceIndex ) {
    return originalTrace( traceIndex ).samples();
  }
  /**
   * 
   * @param traceIndex
   * @return Handle to header values of specified trace
   */
  @Override
  public csHeader[] headerValues( int traceIndex ) {
    return myTraces.get( traceIndex ).headerValues();
  }
  /**
   * 
   * @param traceIndex
   * @return Specified seismic trace
   */
  public csSeismicTrace originalTrace( int traceIndex ) {
    csSeismicTrace trace = myTraces.get( traceIndex );
    if( !myIsProcessed ) {
      return trace;
    }
    else {
      return( new csSeismicTrace( myProcessedDataBuffer.getDataTrace(traceIndex), trace.headerValues(), trace.originalTraceNumber() ) );
    }
  }
  @Override
  public int originalTraceNumber( int traceIndex ) {
    return myTraces.get( traceIndex ).originalTraceNumber();
  }
  /**
   * Clear trace buffer. Remove all traces.
   */
  @Override
  public void clear() {
    myTraces.clear();
    myIsProcessed = false;
    myProcessedDataBuffer.clear();
  }
  /**
   * 
   * @param traceIndex
   * @return Maximum amplitude in specified trace
   */
  @Override
  public float maxAmplitude( int traceIndex ) {
    return originalTrace(traceIndex).maxAmplitude();
  }
  /**
   * 
   * @param traceIndex
   * @return Mean (absolute) amplitude of specified trace
   */
  @Override
  public float meanAmplitude( int traceIndex ) {
    return originalTrace(traceIndex).meanAmplitude();
  }
  /**
   * 
   * @param traceIndex
   * @return DC amplitude (= mean of signed amplitude) of specified trace
   */
  @Override
  public float dcAmplitude( int traceIndex ) {
    return originalTrace(traceIndex).dcAmplitude();
  }
  /**
   * 
   * @param traceIndex
   * @return Minimum amplitude in specified trace
   */
  @Override
  public float minAmplitude( int traceIndex ) {
    return originalTrace(traceIndex).minAmplitude();
  }
  /**
   * 
   * @return Maximum amplitude in all traces
   */
  @Override
  public float maxTotalAmplitude() {
    if( myAreAmplitudesComputed ) return myTotalMaxAmplitude;
    computeTotalAmplitudes();
    return myTotalMaxAmplitude;
  }
  /**
   * 
   * @return Minimum amplitude in all traces
   */
  @Override
  public float minTotalAmplitude() {
    if( myAreAmplitudesComputed ) return myTotalMinAmplitude;
    computeTotalAmplitudes();
    return myTotalMinAmplitude;
  }
  /**
   * 
   * @return Mean amplitude in all traces
   */
  @Override
  public float meanTotalAmplitude() {
    if( myAreAmplitudesComputed ) return myTotalMeanAmplitude;
    computeTotalAmplitudes();
    return myTotalMeanAmplitude;
  }
  public void computeTotalAmplitudes() {
    int nTraces = numTraces();
    myTotalMaxAmplitude = -Float.MAX_VALUE;
    myTotalMinAmplitude = Float.MAX_VALUE;
    myTotalMeanAmplitude = 0.0f;
    for( int itrc = 0; itrc < nTraces; itrc++ ) {
      csSeismicData data = originalTrace(itrc).data();
      myTotalMeanAmplitude += data.meanAmplitude();
      float max = data.maxAmplitude();
      float min = data.minAmplitude();
      if( max > myTotalMaxAmplitude ) myTotalMaxAmplitude = max;
      if( min < myTotalMinAmplitude ) myTotalMinAmplitude = min;
    }
    myTotalMeanAmplitude /= myNumSamples;

    myAreAmplitudesComputed = true;
  }
}


