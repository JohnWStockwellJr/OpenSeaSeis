/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.cmapgen;

import cseis.general.csColorMap;

/**
 * Object used in event interface of csIColorMapChangeListener.
 * @author 2013 Felipe Punto
 */
public class csColorMapChangeEvent {
  public static final int TYPE_UPDATE   = 0;
  public static final int TYPE_ADD      = 1;
  public static final int TYPE_REMOVE   = 2;
  public static final int TYPE_NEW      = 3;

  public java.util.List<csColorMap> cmapList;
  public csColorMap cmap;
  public int type;
  public int itemIndex;
  
  public csColorMapChangeEvent( java.util.List<csColorMap> list ) {
    cmapList = list;
    type = TYPE_NEW;
  }
  public csColorMapChangeEvent( csColorMap cmap_in, int eventType_in ) {
    cmap = cmap_in;
    type = eventType_in;
  }
  public csColorMapChangeEvent( csColorMap cmap_in, int eventType_in, int itemIndex_in ) {
    cmap = cmap_in;
    itemIndex = itemIndex_in;
    type = TYPE_UPDATE;
  }
}

