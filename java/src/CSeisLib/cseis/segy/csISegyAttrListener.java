/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.segy;

/**
 * Listener interface communicating user defined SEGY attributes.
 * @author 2013 Felipe Punto
 */
public interface csISegyAttrListener {
  /**
   * Update attributes defining how SEGY files are read in.
   * @param attr SEGY attributes
   */
  public void updateSegyAttr( csSegyAttr attr );
}

