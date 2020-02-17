/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.velanal;

import cseis.velocity.csVelFunction;

public interface csIVelPickDialogListener {
  public void velPickDialogClosing();
  public void velPickDialogReset();
  public void savePicks( java.util.Collection<csVelFunction> velFuncList, String name );
  public void loadPicks();
}
