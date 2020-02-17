/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import javax.swing.ImageIcon;

/**
 * Actions that can be applied from SeaView application.<br>
 * This class defines standard names, descriptions and icons for standard actions performed by different
 * GUI entities within SeaView application.
 * 
 * @author 2007 Bjorn Olofsson
 * @author 2013 Felipe Punto
 */
public class csSeaViewActions {
  public static final int ForwardSeismicAction  = 0;
  public static final int BackwardSeismicAction = 1;
  public static final int BeginSeismicAction    = 2;
  public static final int EndSeismicAction      = 3;
  public static final int SetAnnotationAction   = 4;
  public static final int ZoomInAction          = 5;
  public static final int ZoomOutAction         = 6;
  public static final int ZoomInHorzAction      = 7;
  public static final int ZoomOutHorzAction     = 8;
  public static final int ZoomInVertAction      = 9;
  public static final int ZoomOutVertAction     = 10;
  public static final int SnapShotAction        = 11;
  public static final int SnapShotPaneAction    = 12;
  public static final int IncreaseScalingAction = 13;
  public static final int DecreaseScalingAction = 14;
  public static final int KillTraceAction       = 15;
  public static final int ShowGraphAction       = 16;
  public static final int ShowOverlayAction     = 17;
  public static final int SpectrumAction        = 18;
  public static final int RubberBandZoomAction  = 19;
  public static final int PanModeAction         = 20;
  public static final int PickModeAction        = 21;
  public static final int FitToScreenAction     = 22;
  public static final int SelectPanesAction     = 23;
  public static final int SelectPanesAlertAction = 24;
  public static final int PaintModeAction       = 25;

  public static final String [] ACTION_DESC = {
    "Move to next seismic section" ,
    "Move to previous seismic section" ,
    "Move to first seismic section" ,
    "Move to last seismic section" ,
    "Select trace headers for trace annotation" ,
    "Zoom in" ,
    "Zoom out" ,
    "Zoom in horizontally" ,
    "Zoom out horizontally" ,
    "Zoom in vertically" ,
    "Zoom out vertically" ,
    "Create snapshot" ,
    "Create snapshot including side labels" ,
    "Increase scaling" ,
    "Decrease scaling" ,
    "Trace-to-zero-setter. Click on trace to set trace samples to zero." ,
    "Show/hide trace header graph" ,
    "Plot trace header values on top of seismic data" ,
    "Compute spectrum in selected window" ,
    "Rubberband zoom" ,
    "Pan mode" ,
    "Horizon picking" ,
    "Fit pane content to screen size",
    "Hide/show panes",
    "Hide/show panes. Hidden panes exist.",
    "Paint mode"
  };
  public static final int [] ACTION = {
    ForwardSeismicAction  ,
    BackwardSeismicAction ,
    BeginSeismicAction    ,
    EndSeismicAction      ,
    SetAnnotationAction   ,
    ZoomInAction          ,
    ZoomOutAction         ,
    ZoomInHorzAction      ,
    ZoomOutHorzAction     ,
    ZoomInVertAction      ,
    ZoomOutVertAction     ,
    SnapShotAction        ,
    SnapShotPaneAction    ,
    IncreaseScalingAction ,
    DecreaseScalingAction ,
    KillTraceAction       ,
    ShowGraphAction       ,
    ShowOverlayAction     ,
    SpectrumAction        ,
    RubberBandZoomAction  ,
    PanModeAction         ,
    PickModeAction        ,
    FitToScreenAction     ,
    SelectPanesAction     ,
    SelectPanesAlertAction,
    PaintModeAction
  };
  public static final String [] ACTION_TITLE = {
    "Next seismic section",
    "Previous seismic section",
    "Go to beginning",
    "Go to end",
    "Set trace annotation...",
    "Zoom in",
    "Zoom out",
    "Zoom in horz",
    "Zoom out horz",
    "Zoom in vert",
    "Zoom out vert",
    "Snapshot",
    "Snapshot incl. side labels",
    "Increase",
    "Decrease",
    "Kill trace",
    "Set trace header graph...",
    "Set trace header overlay...",
    "Compute spectrum",
    "Rubberband zoom",
    "Pan mode",
    "Horizon picking...",
    "Fit to screen",
    "Hide/show panes",
    "Hide/show panes",
    "Paint mode..."
  };
  public static final String [] ACTION_ICON = {
    "csArrowRight.png",
    "csArrowLeft.png",
    "csArrowLeftLeft.png",
    "csArrowRightRight.png",
    null,  // Trace annotation: Does not have icon
    "csZoomIn.gif",
    "csZoomOut.gif",
    "csZoomInHorz.gif",
    "csZoomOutHorz.gif",
    "csZoomInVert.gif",
    "csZoomOutVert.gif",
    "seaview_snapshot.png",
    "seaview_snapshot_pane.png",
    "icon_plus.png",
    "icon_minus.png",
    "seaview_kill_trace.png",
    "seaview_graph.png",
    "seaview_overlay.png",
    "seaview_spectrum.png",
    "csZoomRubberBand.gif",
    "csPanMode.png",
    "csPickMode.png",
    "seaview_fitToScreen.png",
    "seaview_selectPanes.png",
    "seaview_selectPanesAlert.png",
    "csPaintMode.png"
  };
  
  public static ImageIcon getIcon( int actionIndex ) {
    String iconName = csSeaViewActions.ACTION_ICON[actionIndex];
    if( iconName == null ) return null;
    return cseis.resources.csResources.getIcon( iconName );
  }
  public static void main( String[] args ) {
    System.out.println(" Index list:       " + csSeaViewActions.ACTION.length);
    System.out.println(" Title list:       " + csSeaViewActions.ACTION_TITLE.length);
    System.out.println(" Description list: " + csSeaViewActions.ACTION_DESC.length);
    System.out.println(" Icon name list:   " + csSeaViewActions.ACTION_ICON.length);
  }

}


