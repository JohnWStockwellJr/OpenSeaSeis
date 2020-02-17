/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.Color;

public class csHorizonAttr {
  private final int idNumber;
  public String name;
  public Color color;
  public int mode;
  public csHorizonAttr( int id ) {
    idNumber = id;
    name = "";
    color = Color.black;
    mode = csPickOverlay.MODE_PEAK;
  }
  public csHorizonAttr( int id, String name_in, Color color_in, int mode_in ) {
    idNumber = id;
    name  = name_in;
    color = color_in;
    mode  = mode_in;
  }
  public int getID() {
    return idNumber;
  }
}
