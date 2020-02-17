/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.Point;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import javax.swing.JPanel;
import javax.swing.SwingUtilities;

/**
 * Mouse handler for side label
 * @author Bjorn Olofsson
 */
public class csSideLabelMouseHandler implements MouseListener, MouseMotionListener {
  public static final int MOUSE_NONE     = 0;
  public static final int MOUSE_MOVED    = 1;
  public static final int MOUSE_DRAGGED  = 2;
  public static final int MOUSE_RELEASED = 3;
  public int myMouseMode;
  public Point myMouseMovedPos;
  public Point myMouseDraggedPos;
  public int myZoomMode;
  private boolean myIsMouseInside;
  private boolean myIsMousePressed;
  private JPanel myParentPanel;
  private csISideLabelEventListener myListener;

  public csSideLabelMouseHandler( JPanel panel, csISideLabelEventListener listener ) {
    myParentPanel = panel;
    myListener    = listener;
    myMouseMode = MOUSE_NONE;
    myZoomMode  = csSeisPane.ZOOM_OUT;
    myMouseMovedPos   = new Point();
    myMouseDraggedPos = new Point();
    myIsMouseInside  = false;
    myIsMousePressed = false;
  }
  public void mouseClicked( MouseEvent e ) {
  }
  public void mouseEntered( MouseEvent e ) {
    myIsMouseInside = true;
    if( myMouseMode != MOUSE_DRAGGED ) {
      myMouseMode = MOUSE_MOVED;
      myMouseMovedPos = e.getPoint();
    }
  }
  public void mouseExited( MouseEvent e ) {
    myIsMouseInside = false;
    if( !myIsMousePressed ) myMouseMode = MOUSE_NONE;
    myParentPanel.repaint();
  }
  public void mousePressed( MouseEvent e ) {
    if( myMouseMode == MOUSE_DRAGGED ) {
      return;
    }
    myIsMousePressed = true;
    myMouseMovedPos   = e.getPoint();
    myMouseDraggedPos = myMouseMovedPos;
  }
  public void mouseReleased( MouseEvent e ) {
    if( myIsMousePressed && myMouseMode != MOUSE_NONE ) {
//      int button = e.getButton();
      if( SwingUtilities.isLeftMouseButton(e) ) {
        myZoomMode = csSeisPane.ZOOM_IN;
      }
      else if( SwingUtilities.isMiddleMouseButton(e) ) {
        myZoomMode = csSeisPane.ZOOM_OUT;
      }
      else if( SwingUtilities.isRightMouseButton(e) ) {
        myZoomMode = csSeisPane.ZOOM_OUT;
      }
      fireZoomEvent();
    }
    myMouseMode = MOUSE_NONE;
    myMouseMovedPos = e.getPoint();
    myIsMousePressed = false;
    myParentPanel.repaint();
    if( myIsMouseInside ) myMouseMode = MOUSE_MOVED;  // Allow subsequent clicks
  }
  public void mouseDragged( MouseEvent e ) {
    if( !myIsMousePressed ) return;
    myMouseDraggedPos = e.getPoint();
    myMouseMode = MOUSE_DRAGGED;
    myParentPanel.repaint();
  }
  public void mouseMoved( MouseEvent e ) {
    myMouseMovedPos = e.getPoint();
    myMouseMode = MOUSE_MOVED;
    myParentPanel.repaint();
  }
  private void fireZoomEvent() {
    myListener.zoomEvent( myMouseMovedPos, myMouseDraggedPos, myZoomMode );
  }
}


