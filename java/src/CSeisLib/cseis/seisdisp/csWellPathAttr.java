/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.seisdisp;

import java.awt.Color;
import java.util.ArrayList;
import cseis.math.csPoint3D;

public class csWellPathAttr {
  private static int ID_COUNTER = 0;

  public final int id;
  public String filename;
  public String path;
  public ArrayList<csPoint3D> posXYZList;
  public Color   color;
  public boolean showWell;

  public csWellPathAttr() {
    this( null, null );
  }
  public csWellPathAttr( String filename_in, String path_in ) {
    id = ++ID_COUNTER;
    posXYZList = new ArrayList<csPoint3D>();
    filename  = filename_in;
    path      = path_in;
    color     = Color.red;
    showWell  = false;
  }
  public csWellPathAttr( csWellPathAttr attr ) {
    id = attr.id;
    posXYZList = new ArrayList<csPoint3D>( attr.posXYZList );
    filename  = attr.filename;
    path      = attr.path;
    color     = attr.color;
    showWell  = attr.showWell;
  }
  public boolean isDifferent( csWellPathAttr attr ) {
    return( showWell != attr.showWell || color != attr.color || filename != attr.filename ); 
  }
  public void setAttributes( csWellPathAttr attr ) {
    color     = attr.color;
    showWell  = attr.showWell;
  }
}
