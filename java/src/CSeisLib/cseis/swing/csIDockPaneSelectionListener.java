/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

/**
 * Listener interface which is used to communicate events related to dock pane selection.
 * @author 2013 Felipe Punto
 */
public interface csIDockPaneSelectionListener {
  /**
   * Dock pane has been selected.
   * @param pane The pane.
   */
  public void dockPaneSelected( csDockPane pane );
  /**
   * Dock pane has been closed.
   * @param pane The pane.
   */
  public void dockPaneClosed( csDockPane pane );
  /**
   * Number of hidden dock panes changed.
   * @param hiddenPanesExist True if at least one pane is hidden from view
   */
  public void hideStateChanged( boolean hiddenPanesExist );
}

