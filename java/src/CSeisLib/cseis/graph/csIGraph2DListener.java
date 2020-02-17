/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

/**
 * Listener interface for graph events.
 * @author Bjorn Olofsson
 */
public interface csIGraph2DListener {
  public void graph2DValues( float xModel, float yModel );
  public void graphChanged();
}
