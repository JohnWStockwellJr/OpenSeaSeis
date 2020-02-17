/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.general.*;

/**
 * Color bar used for image plot.
 * @author 2009 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csPlotImageColorBar extends csColorBarPanel {

  private int myHeight;
  private int myWidth;

  public csPlotImageColorBar( csColorMap map, int width, int height ) {
    this( map, csColorBarPanel.ORIENT_VERTICAL, width, height );
  }
  public csPlotImageColorBar( csColorMap map, int orientation, int width, int height ) {
    this( map, orientation, width, height, csColorBarPanel.ANNOTATION_ADVANCED);
  }
  public csPlotImageColorBar( csColorMap map, int orientation, int width, int height, int annotation ) {
    super( map, orientation, annotation );
    myHeight = height;
    myWidth = width;
  }
  @Override
  public int getHeight() {
    return myHeight;
  }
  @Override
  public int getWidth() {
    return myWidth;
  }
}


