/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

/**
 * Listener interface which is used to communicate events related to dock panes.
 * @author 2013 Felipe Punto
 */
public interface csIDockPaneListener {
  /**
   * Hide pane.
   * @param pane The pane
   */
  public void hidePane( csDockPane pane );
  /**
   * Close pane.
   * @param pane The pane
   */
  public void closePane( csDockPane pane );
  /**
   * Dock pane.
   * @param pane The pane
   */
  public void dockPane( csDockPane pane, boolean dock );
  /**
   * Maximize pane. (=hide all other panes)
   * @param pane The pane
   */
  public void maximizePane( csDockPane pane );
  /**
   * Show or hide scrollbars.
   * @param showScrollbars Pass 'true' if scrollbars shall be shown.
   */
  public void scrollbars( boolean showScrollbars );
}

