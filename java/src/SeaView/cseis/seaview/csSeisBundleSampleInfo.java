/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.seisdisp.csSampleInfo;

/**
 * Class used to convey information about a seismic sample (typically at current mouse location).
 * @author 2013 Felipe Punto
 */
public class csSeisBundleSampleInfo {
  public csSampleInfo info;
  public double timeFullSample;
  public boolean isFrequencyDomain;
  public long time_samp1_ms;
  public float zoomVert;
  public int originalTrace;
  
  public csSeisBundleSampleInfo() {
    info = new csSampleInfo();
  }
  public csSeisBundleSampleInfo( csSampleInfo info_in ) {
    info = info_in;
  }
}

