/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.io.File;
import java.util.EventObject;

/**
 * File menu event
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csFileMenuEvent extends EventObject {
  private File myFile;

  public csFileMenuEvent( Object source, File file ) {
    super( source );
    myFile = file;
  }
  public File file() {
    return myFile;
  }

}


