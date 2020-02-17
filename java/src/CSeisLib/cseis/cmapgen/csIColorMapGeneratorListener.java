/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.cmapgen;

import cseis.general.csCustomColorMap;

/**
 * Listener interface which is used to communicate events from color map generator dialog.
 * @author 2013 Felipe Punto
 */
public interface csIColorMapGeneratorListener {
  /**
   * Apply/test new color map.
   * @param cmap The color map.
   */
  public void applyColorMap( csCustomColorMap cmap );
  /**
   * Update list of custom color maps.
   * @param list List of custom color maps.
   */
  public void updateColorMaps( java.util.List<csCustomColorMap> list );
}

