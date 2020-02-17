/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import cseis.general.csStandard;

/**
 * Trace/sample information for current trace/sample.
 * @author Bjorn Olofsson
 */
public class csSampleInfo {
  /// Seismic amplitude
  public double amplitude;
  /// Depth [m]
  public double depth;

  public int domainType;
  /// Offset
  public double offset;
  /// Time [s]
  public double time;
  /// Sample index, starting at 0
  public int sample;
  /// Sample index, including decimal places, starting at 0
  public double sampleDouble;
  /// Trace number (interval counting, starting with 0 from the left)
  public int trace;
  /// Trace number, including decimal places
  public double traceDouble;
  
  public csSampleInfo() {
    time        = csStandard.ABSENT_VALUE;
    depth       = csStandard.ABSENT_VALUE;
    amplitude   = 0.0;
    offset      = csStandard.ABSENT_VALUE;
    trace       = csStandard.ABSENT_VALUE_INT;
    traceDouble = csStandard.ABSENT_VALUE;
    domainType  = 0;
  }
  public csSampleInfo( csSampleInfo sInfo ) {
    time        = sInfo.time;
    depth       = sInfo.depth;
    amplitude   = sInfo.amplitude;
    offset      = sInfo.offset;
    trace       = sInfo.trace;
    traceDouble = sInfo.traceDouble;
    domainType  = sInfo.domainType;
  }
}


