/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

/**
 * Listener interface which is used to communicate updates while file read or scan operation is in progress.
 * @author 2013 Felipe Punto
 */
public interface csIDataSetListener {
  /**
   * Update data read progress.
   * @param currentTraceIndex Index of currently read trace.
   */
  public void updateTrace( int currentTraceIndex );
  /**
   * Update data read progress.
   * @param numTracesRead       Number of traces which have been read in so far.
   * @param currentTraceIndex   Index of currently read trace.
   */
  public void updateTrace( int numTracesRead, int currentTraceIndex );
  /**
   * Indicates that read operation is complete.
   * @param success true if read operation was successful, false otherwise.
   */
  public void isComplete( boolean success );
  /**
   * Indicates that current scan operation shall be stopped.
   * @param bundle  Seismic bundle for which scan operation was run.
   * @param success true if scan operation was successfully completed, false otherwise.
   */
  public void stopScan( csSeisPaneBundle bundle, boolean success );
  /**
   * Update progress of scan operation.
   * @param bundle  Seismic bundle for which scan operation is being run.
   * @param currentTraceIndex   Index of currently scanned trace.
   */
  public void updateScan( csSeisPaneBundle bundle, int currentTraceIndex );
}

