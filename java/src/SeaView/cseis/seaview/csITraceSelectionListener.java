/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.seis.csHeaderDef;

/**
 * Listener interface which is used to communicate updates related to trace selection.
 * @author 2013 Felipe Punto
 */
public interface csITraceSelectionListener {
  /**
   * Update trace selection.
   * @param param Parameters defining current trace selection settings.
   */
  public void updateTraceSelection( csTraceSelectionParam param );
  /**
   * Start new trace header selection scan operation.
   * @param headerDef Trace header definition object defining the trace header which shall be scanned.
   */
  public void startScanTraceSelection( csHeaderDef headerDef );
  /**
   * Cancel current scan operation.
   */
  public void cancelScanTraceSelection();
}

