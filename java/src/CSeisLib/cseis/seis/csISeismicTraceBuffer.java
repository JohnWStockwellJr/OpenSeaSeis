/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seis;

/**
 * Interface to be implemented by seismic trace buffers
 * @author Bjorn Olofsson
 */
public interface csISeismicTraceBuffer {
  /**
   * @return Number of traces in buffer
   */
  public int numTraces();
  /**
   * @return Number of samples in each trace
   */
  public int numSamples();
  /**
   * @return Number of trace headers
   */
  public int numHeaders();
  /**
   * Gain access to trace samples of one trace.<br>
   * Any changes made to the retrieved array are permanent in the trace buffer.
   * @param traceIndex Index of trace
   * @return Array of trace samples
   */
  public float[] samples( int traceIndex );
  /**
   * Gain access to trace header values of one trace.
   * Any changes made to the retrieved array are permanent in the trace buffer.
   * @param traceIndex Index of trace
   * @return Array of trace header values
   */
  public csHeader[] headerValues( int traceIndex );
  /**
   * Clear trace buffer. Delete all traces.
   */
  public void clear();
  /**
   * @param traceIndex Index of trace
   * @return Maximum amplitude in trace
   */
  public float maxAmplitude( int traceIndex );
  /**
   * @param traceIndex Index of trace
   * @return Mean (absolute) amplitude in trace
   */
  public float meanAmplitude( int traceIndex );
  /**
   * @param traceIndex Index of trace
   * @return Minimum amplitude in trace
   */
  public float minAmplitude( int traceIndex );
  /**
   * @param traceIndex Index of trace
   * @return DC amplitude (= mean of signed amplitude) of specified trace
   */
  public float dcAmplitude( int traceIndex );
  /**
   * @return Maximum amplitude in whole trace buffer.
   */
  public float maxTotalAmplitude();
  /**
   * @return Minimum amplitude in whole trace buffer.
   */
  public float minTotalAmplitude();
  /**
   * @return Mean (absolute) amplitude in whole buffer.
   */
  public float meanTotalAmplitude();
  /**
   * Retrieve original trace number.<br>
   * The original trace number typically differs from the trace index.
   * @param traceIndex Index of trace
   * @return Original trace number
   */
  public int originalTraceNumber( int traceIndex );
}


