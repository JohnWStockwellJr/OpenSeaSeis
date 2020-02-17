/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.Graphics2D;

/**
 * Interface for seismic display overlay
 * @author Bjorn Olofsson
 */
public interface csISeisOverlay {
  public void draw( csSeisView seisview, Graphics2D g );
}


