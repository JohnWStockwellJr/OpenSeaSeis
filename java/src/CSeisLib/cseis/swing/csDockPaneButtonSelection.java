/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

/**
 * Attributes defining buttons on top ribbon of csDockPane.<br>
 * Required when creating a new csDockPane object. Each button may be shown or hidden.
 * @author 2013 Felipe Punto
 */
public class csDockPaneButtonSelection {
  public boolean sync;
  public boolean dock;
  public boolean scrollbars;
  public boolean hide;
  public boolean close;
  public boolean maximize;
  public csDockPaneButtonSelection() {
    sync = true;
    dock = true;
    scrollbars = true;
    hide = true;
    close = true;
    maximize = true;
  }
}

