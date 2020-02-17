/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

/**
 * Listener interface to communicate events related to dock pane synchronization
 * @author 2013 Felip Punto
 */
public interface csIDockPaneSyncListener {
  /**
   * Sync state changed.
   * @param isSync true if pane shall now be synced
   */
  public void syncStateChanged( csDockPane pane, boolean isSync );
}

