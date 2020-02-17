/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.seis.csHeaderDef;
import cseis.seisdisp.csHorizon;

public interface csIPickDialogListener {
  public void closePickDialog();
  public void updatePicks();
  public void setPicksFromHeader( csHeaderDef hdrDef );  
  public void savePicks( csHorizon horizon, String name );  
}
