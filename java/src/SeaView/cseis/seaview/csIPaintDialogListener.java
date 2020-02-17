/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import java.awt.Color;

public interface csIPaintDialogListener {
  public void closePaintDialog();
  public void clearPaintOverlay();
  public void updatePaintOverlay( int lineSize, Color lineColor, int pointSize, Color pointColor );
}
