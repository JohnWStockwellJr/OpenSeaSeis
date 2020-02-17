/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

/**
 * Listener interface which is used to communicate updates related to seismic bundle.
 * @author 2013 Felipe Punto
 */
public interface csISeisPaneBundleListener {
  /**
   * Update seismic display scalar.
   * @param bundle Seismic bundle which called this function.
   */
  public void updateBundleDisplayScalar( csSeisPaneBundle bundle );
  /**
   * Update seismic display settings.
   * @param bundle Seismic bundle which called this function.
   */
  public void updateBundleDisplaySettings( csSeisPaneBundle bundle );
  /**
   * Update sample information at current mouse location.
   * @param bundle Seismic bundle which called this function
   * @param sampleInfo Information about sample at current mouse location.
   */
  public void updateBundleSampleInfo( csSeisPaneBundle bundle, csSeisBundleSampleInfo sampleInfo );
}

