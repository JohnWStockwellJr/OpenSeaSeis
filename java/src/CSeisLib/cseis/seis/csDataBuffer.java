/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seis;

import java.util.ArrayList;

/**
 * Data buffer. Basically a trace buffer (csTraceBuffer) stripped of its headers.
 * @author 2013 Felipe Punto
 */
public class csDataBuffer {
  private ArrayList<csSeismicData> data;
  public csDataBuffer() {
    data = new ArrayList<csSeismicData>();
  }
  public void addDataTrace( csSeismicData dataTrace ) {
    data.add( dataTrace );
  }
  public csSeismicData getDataTrace( int traceIndex ) {
    return data.get( traceIndex );
  }
  public void clear() {
    data.clear();
  }
}

