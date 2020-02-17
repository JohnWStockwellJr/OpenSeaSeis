/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import javax.swing.ImageIcon;

/**
 * Provides static sub-classes defining different csDockPane actions.<br>
 * Action sub-classes really only define the action name, tooltip and icon.
 * @author 2013 Felipe Punto
 */
public class csDockPaneActions {
  //------------------------------------------------------------------------
  public static final int DockAction     = 0;
  public static final int UndockAction   = 1;
  public static final int SyncAction     = 2;
  public static final int UnsyncAction   = 3;
  public static final int HideAction     = 4;
  public static final int CloseAction    = 5;
  public static final int MaximizeAction = 6;
  public static final int ScrollbarAction= 7;
  public static final int CloseTabAction= 8;

  public static final String [] ACTION_DESC = {
    "Dock pane to main window" ,
    "Undock pane to separate window" ,
    "Synchronize this pane with other panes" ,
    "Unsync this pane" ,
    "Hide pane to background"   ,
    "Close this pane"  ,
    "Maximize pane"    ,
    "Hide scroll bars" ,
    "Close this tab"
  };
  public static final int [] ACTION = {
    DockAction     ,
    UndockAction   ,
    SyncAction     ,
    UnsyncAction   ,
    HideAction     ,
    CloseAction    ,
    MaximizeAction ,
    ScrollbarAction,
    CloseTabAction
  };
  public static final String [] ACTION_TITLE = {
    "Dock pane",
    "Undock pane",
    "Sync pane",
    "Unsync pane",
    "Hide pane",
    "Close pane",
    "Maximize pane",
    "Hide scroll bars",
    "Close tab"
  };
  public static final String [] ACTION_ICON = {
    "csDockPane_dock.png",
    "csDockPane_undock.png",
    "csDockPane_unsynced.png",
    "csDockPane_synced.png",
    "csDockPane_hide.png",
    "csDockPane_close.png",
    "csDockPane_maximize.png",
    "csDockPane_scrollbars.png",
    "csDockPane_close_tab.png"
  };
  
  public static ImageIcon getIcon( int actionIndex ) {
    String iconName = csDockPaneActions.ACTION_ICON[actionIndex];
    if( iconName == null ) return null;
    return cseis.resources.csResources.getIcon( iconName );
  }

}


