/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.velanal;

import cseis.seaview.csSeisPaneBundle;
import cseis.velocity.csVelEnsembleInfo;

public interface csIVelSetupDialogListener {
  public void setupVelocityAnalysis( csSeisPaneBundle bundle, csVelEnsembleInfo ensInfo );
}
