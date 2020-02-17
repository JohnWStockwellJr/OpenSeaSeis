/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

/**
 * Listener interface which is used to communicate updates relevant to SeaView's toolbars.
 * @author 2013 Felipe Punto
 */
public interface csISeaViewToolBarListener {
  /**
   * Update seismic display scalar.
   * @param scalar  New value for display scalar.
   */
  public void updateDisplayScalar( float scalar );
  /**
   * Increment/increase seismic display scalar.
   */
  public void incDisplayScalar();
  /**
   * Decrement/decrease seismic display scalar.
   */
  public void decDisplayScalar();
}

