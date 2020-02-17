/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.event.ActionEvent;

import javax.swing.AbstractAction;

/**
 * Seismic pane actions.<br>
 * Different actions that can be applied to seismic display pane.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csSeisPaneActions {
  private csSetAnnotationAction mySetAnnotationAction;
  private csZoomInAction myZoomInAction;
  private csZoomOutAction myZoomOutAction;
  private csZoomInHorzAction myZoomInHorzAction;
  private csZoomOutHorzAction myZoomOutHorzAction;
  private csZoomInVertAction myZoomInVertAction;
  private csZoomOutVertAction myZoomOutVertAction;
  private csForwardSeismicAction myForwardSeismicAction;
  private csBackwardSeismicAction myBackwardSeismicAction;
  private csBeginSeismicAction myBeginSeismicAction;
  private csEndSeismicAction myEndSeismicAction;
  private csSeisPane myPane;

  public csSeisPaneActions( csSeisPane pane ) {
    myPane = pane;
    mySetAnnotationAction = new csSetAnnotationAction();
    myZoomInAction        = new csZoomInAction();
    myZoomOutAction       = new csZoomOutAction();
    myZoomInHorzAction    = new csZoomInHorzAction();
    myZoomOutHorzAction   = new csZoomOutHorzAction();
    myZoomInVertAction    = new csZoomInVertAction();
    myZoomOutVertAction   = new csZoomOutVertAction();

    myForwardSeismicAction  = new csForwardSeismicAction();
    myBackwardSeismicAction = new csBackwardSeismicAction();
    myBeginSeismicAction    = new csBeginSeismicAction();
    myEndSeismicAction      = new csEndSeismicAction();
  }
  public csSetAnnotationAction getSetAnnotationAction() {
    return mySetAnnotationAction;
  }
  public csZoomInAction getZoomInAction() {
    return myZoomInAction;
  }
  public csZoomOutAction getZoomOutAction() {
    return myZoomOutAction;
  }
  public csZoomInHorzAction getZoomInHorzAction() {
    return myZoomInHorzAction;
  }
  public csZoomOutHorzAction getZoomOutHorzAction() {
    return myZoomOutHorzAction;
  }
  public csZoomInVertAction getZoomInVertAction() {
    return myZoomInVertAction;
  }
  public csZoomOutVertAction getZoomOutVertAction() {
    return myZoomOutVertAction;
  }
  public csForwardSeismicAction getForwardSeismicAction() {
    return myForwardSeismicAction;
  }
  public csBackwardSeismicAction getBackwardSeismicAction() {
    return myBackwardSeismicAction;
  }
  public csBeginSeismicAction getBeginSeismicAction() {
    return myBeginSeismicAction;
  }
  public csEndSeismicAction getEndSeismicAction() {
    return myEndSeismicAction;
  }
  //------------------------------------------------------------------------
  class csForwardSeismicAction extends AbstractAction {
    public csForwardSeismicAction() {
      super( "Next seismic section", cseis.resources.csResources.getIcon("csArrowRight.png") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Move to next seismic section" );
    }
    public void actionPerformed( ActionEvent e ) {
//      myTremor.seismicForward();
    }
  }
  //------------------------------------------------------------------------
  class csBackwardSeismicAction extends AbstractAction {
    public csBackwardSeismicAction() {
      super( "Previous seismic section", cseis.resources.csResources.getIcon("csArrowLeft.png") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Move to previous seismic section" );
    }
    public void actionPerformed( ActionEvent e ) {
//      myTremor.seismicBackward();
    }
  }
  //------------------------------------------------------------------------
  class csBeginSeismicAction extends AbstractAction {
    public csBeginSeismicAction() {
      super( "Go to beginning", cseis.resources.csResources.getIcon("csArrowLeftLeft.png") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Move to first seismic section" );
    }
    public void actionPerformed( ActionEvent e ) {
//      myTremor.seismicBegin();
    }
  }
  //------------------------------------------------------------------------
  class csEndSeismicAction extends AbstractAction {
    public csEndSeismicAction() {
      super( "Go to end", cseis.resources.csResources.getIcon("csArrowRightRight.png") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Move to last seismic section" );
    }
    public void actionPerformed( ActionEvent e ) {
//      myTremor.seismicEnd();
    }
  }
  //------------------------------------------------------------------------
  class csSetAnnotationAction extends AbstractAction {
    public csSetAnnotationAction() {
      super( "Set trace annotation..." );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Select trace headers for trace annotation" );
    }
    public void actionPerformed( ActionEvent e ) {
      myPane.showAnnotationDialog();
    }
  }
  //------------------------------------------------------------------------
  class csZoomInAction extends AbstractAction {
    public csZoomInAction() {
      super( "Zoom in", cseis.resources.csResources.getIcon("osZoomIn.gif") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Zoom in" );
    }
    public void actionPerformed( ActionEvent e ) {
      myPane.zoom( csSeisPane.ZOOM_IN, csSeisPane.ZOOM_BOTH );
    }
  }
  //------------------------------------------------------------------------
  class csZoomOutAction extends AbstractAction {
    public csZoomOutAction() {
      super( "Zoom out", cseis.resources.csResources.getIcon("osZoomOut.gif") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Zoom out" );
    }
    public void actionPerformed( ActionEvent e ) {
      myPane.zoom( csSeisPane.ZOOM_OUT, csSeisPane.ZOOM_BOTH );
    }
  }
  //------------------------------------------------------------------------
  class csZoomInHorzAction extends AbstractAction {
    public csZoomInHorzAction() {
      super( "Zoom in horz", cseis.resources.csResources.getIcon("osZoomInHorz.gif") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Zoom in horizontally" );
    }
    public void actionPerformed( ActionEvent e ) {
      myPane.zoom( csSeisPane.ZOOM_IN, csSeisPane.ZOOM_HORZ );
    }
  }
  //------------------------------------------------------------------------
  class csZoomOutHorzAction extends AbstractAction {
    public csZoomOutHorzAction() {
      super( "Zoom out horz", cseis.resources.csResources.getIcon("osZoomOutHorz.gif") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Zoom out horizontally" );
    }
    public void actionPerformed( ActionEvent e ) {
      myPane.zoom( csSeisPane.ZOOM_OUT, csSeisPane.ZOOM_HORZ );
    }
  }
  //------------------------------------------------------------------------
  class csZoomInVertAction extends AbstractAction {
    public csZoomInVertAction() {
      super( "Zoom in vert", cseis.resources.csResources.getIcon("osZoomInVert.gif") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Zoom in vertically" );
    }
    public void actionPerformed( ActionEvent e ) {
      myPane.zoom( csSeisPane.ZOOM_IN, csSeisPane.ZOOM_VERT );
    }
  }
  //------------------------------------------------------------------------
  class csZoomOutVertAction extends AbstractAction {
    public csZoomOutVertAction() {
      super( "Zoom out vert", cseis.resources.csResources.getIcon("osZoomOutVert.gif") );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Zoom out vertically" );
    }
    public void actionPerformed( ActionEvent e ) {
      myPane.zoom( csSeisPane.ZOOM_OUT, csSeisPane.ZOOM_VERT );
    }
  }
}


