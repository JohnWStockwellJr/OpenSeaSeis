/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.seis.csHeaderDef;

/**
 * Listener interface which is used to communicate updates related to trace header scan process.
 * @author 2013 Felipe Punto
 */
public interface csITraceHeaderScanListener {
  /**
   * Indicates that scan operation has completed.
   * @param headerDef Trace header definition object defining scanned trace header.
   */
  public void scanCompleted( csHeaderDef headerDef );
}

