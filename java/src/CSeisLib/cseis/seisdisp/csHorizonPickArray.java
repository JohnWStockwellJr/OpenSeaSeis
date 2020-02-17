/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.seisdisp;
/**
 * 
 */
public class csHorizonPickArray {
  public int trace1;
  public int trace2;
  public int numTraces;
  public float[] picks;
  
  public csHorizonPickArray( int traceFrom, int traceTo, float[] pickedTimesIn ) {
    trace1 = traceFrom;
    trace2 = traceTo;
    numTraces = trace2 - trace1 + 1;
    picks = pickedTimesIn;
  }
  public csHorizonPickArray( int traceFrom, int traceTo ) {
    trace1 = traceFrom;
    trace2 = traceTo;
    numTraces = trace2 - trace1 + 1;
    picks = new float[numTraces];
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      picks[itrc] = csPickOverlay.NO_VALUE;
    } 
  }
  public csHorizonPickArray createSubArray( int traceFrom, int traceTo ) {
    int numTracesNew = traceTo - traceFrom + 1;
    float[] pickedTimesNew = new float[numTracesNew];
    int startArrayIndex = traceFrom - trace1;
    System.arraycopy(picks, startArrayIndex, pickedTimesNew, 0, numTracesNew);
    csHorizonPickArray array = new csHorizonPickArray( traceFrom, traceTo, pickedTimesNew );
    return array;
  }
  void removeAboveAndIncluding( int traceFrom ) {
//      if( traceFrom <= traceIndex1 ) throw( new Exception("") );
    int numTracesNew = traceFrom - trace1;
    float[] pickedTimesNew = new float[numTracesNew];
    System.arraycopy(picks, 0, pickedTimesNew, 0, numTracesNew);
    trace2 = traceFrom-1;
  }
  void removeBelowAndIncluding( int traceTo ) {
//      if( traceFrom <= traceIndex1 ) throw( new Exception("") );
    int numTracesNew = trace2 - traceTo;
    float[] pickedTimesNew = new float[numTracesNew];
    int startArrayIndex = traceTo - trace1 + 1;
    System.arraycopy(picks, startArrayIndex, pickedTimesNew, 0, numTracesNew);
    trace1 = traceTo+1;
  }
}
