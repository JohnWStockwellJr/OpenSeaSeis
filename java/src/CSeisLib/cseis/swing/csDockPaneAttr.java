/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

/**
 * Attributes defining state of csDockPane, as well as its title and description.<br>
 * A csDockPane may be (un-)docked, (un-)synched and visible/hidden from view.
 * @author 2013 Felipe Punto
 */
public class csDockPaneAttr {
  public boolean isSync;
  public boolean isScrollVisible;
  public String title;
  public String titleDescription;
  public boolean isDocked;
  public boolean isActive;

  public csDockPaneAttr() {
    this("Panel");
  }
  public csDockPaneAttr( String title_in ) {
    isSync = false;
    isScrollVisible = true;
    title = title_in;
    titleDescription = title_in;
    isDocked = true;
    isActive = false;
  }
  public csDockPaneAttr( csDockPaneAttr attr ) {
    isSync = attr.isSync;
    isScrollVisible = attr.isScrollVisible;
    title = attr.title;
    titleDescription = attr.titleDescription;
    isDocked = attr.isDocked;
    isActive = attr.isActive;
  }
}

