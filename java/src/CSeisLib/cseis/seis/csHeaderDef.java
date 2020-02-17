/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seis;

import cseis.jni.csJNIDef;

/**
 * Header definition
 * @author Bjorn Olofsson
 */
public class csHeaderDef implements Comparable<csHeaderDef> {
  public static final int TYPE_INT    = csJNIDef.TYPE_INT;
  public static final int TYPE_LONG   = csJNIDef.TYPE_LONG;
  public static final int TYPE_DOUBLE = csJNIDef.TYPE_DOUBLE;
  public static final int TYPE_FLOAT  = csJNIDef.TYPE_FLOAT;
  
  /// Header type. Currently supported are TYPE_INT, TYPE_FLOAT and TYPE_DOUBLE
  public int type;
  /// Header name
  public String name;
  /// Header description
  public String desc;
  /// Index, user specifiable. May for example be used if some kind of sorting is required.
  public int index;

  /**
   * Default constructor
   */
  public csHeaderDef() {
    name = "";
    desc = "";
    type = TYPE_FLOAT;
    index = 0;
  }
  /**
   * 'Copy' constructor
   * @param h
   */
  public csHeaderDef( csHeaderDef h ) {
    name = h.name;
    desc = h.desc;
    type = h.type;
    index = h.index;
  }
  /**
   * Constructor
   * @param theName  Header name
   * @param theDesc  Header description
   * @param theType  Header type
   */
  public csHeaderDef( String theName, String theDesc, int theType ) {
    this( theName, theDesc, theType, 0 );
  }
  /**
   * Constructor
   * @param theName  Header name
   * @param theDesc  Header description
   * @param theType  Header type
   * @param theIndex Index, user specifiable. May for example be used if some kind of sorting is required.

   */
  public csHeaderDef( String theName, String theDesc, int theType, int theIndex ) {
    name = theName;
    type = theType;
    desc = theDesc;
    index = theIndex;
  }
  @Override
  public String toString() {
    return name;
  }
  @Override
  public int compareTo( csHeaderDef otherHeader ) {
    return( name.compareTo(otherHeader.name) );
  }
  @Override
  public boolean equals( Object obj ) {
    if( obj instanceof csHeaderDef ) {
      return( ((csHeaderDef)obj).name.compareTo(name) == 0 );
    }
    else {
      return false;
    }
  }
}


