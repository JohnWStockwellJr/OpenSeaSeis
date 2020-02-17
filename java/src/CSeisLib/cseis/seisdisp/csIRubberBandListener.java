/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

/**
 * Listener interface for 'rubberband' events
 * @author Bjorn Olofsson
 */
public interface csIRubberBandListener {
  public void rubberBandCompleted( csSampleInfo posStart, csSampleInfo posEnd );
}


