/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import cseis.general.csColorMap;
import cseis.general.csSquareIcon;
import cseis.general.csStandard;

/**
 * List item for color maps.
 * @author Bjorn Olofsson
 */
public class csColorMapListItem {
  public csSquareIcon icon;
  public csColorMap map;
  public csColorMapListItem( csSquareIcon iconIn, csColorMap mapIn ) {
    icon = iconIn;
    map  = mapIn;
  }
  public static csColorMapListItem createStandardItem( csColorMap cmap ) {
    csColorMapListItem item = new csColorMapListItem( new csSquareIcon( csStandard.ICON_SIZE, cmap ), cmap );
    return item;
  }
  @Override
  public String toString() {
    return map.toString();
  }
}

