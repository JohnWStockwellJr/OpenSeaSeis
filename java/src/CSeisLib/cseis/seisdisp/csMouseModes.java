/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.Cursor;
import java.awt.Toolkit;
import java.awt.Point;

/**
 * Mouse mode definitions
 * @author Bjorn Olofsson
 */
public class csMouseModes {
  public static final int NO_MODE   = 1;
  public static final int ZOOM_MODE = 2;
  public static final int PAN_MODE  = 3;
  public static final int KILL_MODE = 4;
  public static final int SPECTRUM_MODE = 5;
  public static final int PICK_MODE = 6;
  public static final int PAINT_MODE  = 7;
  public static final int VELPICK_MODE = 8;

  public static final Cursor ZOOM_CURSOR;
  public static final Cursor PAN_CURSOR;
  public static final Cursor KILL_CURSOR;
  public static final Cursor SPECTRUM_CURSOR;
  public static final Cursor PICK_CURSOR;
  public static final Cursor PAINT_CURSOR;
  public static final Cursor VELPICK_CURSOR;

  static {
    ZOOM_CURSOR = Toolkit.getDefaultToolkit().createCustomCursor(
      cseis.resources.csResources.getIcon("csZoomCursor.png").getImage(), new Point( 15, 15 ), "ZOOM_CURSOR" );
    PICK_CURSOR = Toolkit.getDefaultToolkit().createCustomCursor(
      cseis.resources.csResources.getIcon("csPickCursor.png").getImage(), new Point( 15, 15 ), "PICK_CURSOR" );
    VELPICK_CURSOR = Toolkit.getDefaultToolkit().createCustomCursor(
      cseis.resources.csResources.getIcon("csVelPickCursor.png").getImage(), new Point( 15, 15 ), "VELPICK_CURSOR" );
    SPECTRUM_CURSOR = Toolkit.getDefaultToolkit().createCustomCursor(
      cseis.resources.csResources.getIcon("csSpectrumCursor.png").getImage(), new Point( 15, 15 ), "SPECTRUM_CURSOR" );
    PAN_CURSOR = Cursor.getPredefinedCursor( Cursor.HAND_CURSOR );
    KILL_CURSOR = Cursor.getPredefinedCursor( Cursor.CROSSHAIR_CURSOR );
    PAINT_CURSOR = Toolkit.getDefaultToolkit().createCustomCursor(
      cseis.resources.csResources.getIcon("csPaintCursor.png").getImage(), new Point( 15, 15 ), "PAINT_CURSOR" );
  }
}


