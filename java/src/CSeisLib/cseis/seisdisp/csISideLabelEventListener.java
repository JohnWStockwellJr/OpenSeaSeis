/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.Point;

/**
 * Event listener for side label
 * @author Bjorn Olofsson
 */
public interface csISideLabelEventListener {
  public void zoomEvent( Point p1, Point p2, int zoomMode );
}


