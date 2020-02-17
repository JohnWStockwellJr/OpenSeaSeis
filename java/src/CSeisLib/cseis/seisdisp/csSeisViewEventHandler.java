/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.event.*;
import java.util.*;
import java.awt.Point;
import java.awt.Cursor;
import javax.swing.SwingUtilities;

/**
 * Event handler for csSeisView.<br>
 * Catches mouse and key events for csSeisView objects.
 * @author Bjorn Olofsson
 */
public class csSeisViewEventHandler implements MouseMotionListener, MouseListener, KeyListener {
  private csSeisView myView;
  private final ArrayList<csISampleInfoListener> mySampleInfoListeners;
  private final ArrayList<csIKeyListener> myKeyListeners;
  private final ArrayList<csIPanningListener> myPanningListeners;
  private final ArrayList<csIRubberBandListener> myRubberBandListeners;
//  private boolean myIsSHIFTPressed = false;
//  private boolean myIsCTRLPressed = false;
  private Point myPressedPoint;
  private boolean myIsMouseInWindow;
  private final csRubberBandOverlay myRubberBand;

  public csSeisViewEventHandler( csSeisView view ) {
    myView = view;
    myView.addMouseListener( this );
    myView.addMouseMotionListener( this );
    myView.addKeyListener(this);
    mySampleInfoListeners = new ArrayList<csISampleInfoListener>();
    myKeyListeners = new ArrayList<csIKeyListener>();
    myPanningListeners = new ArrayList<csIPanningListener>();
    myRubberBandListeners = new ArrayList<csIRubberBandListener>();
    myIsMouseInWindow = false;
    myRubberBand = new csRubberBandOverlay();
  }
//------------------------------------------------------------
  @Override
  public void mouseDragged( MouseEvent e ) {
    int xpos = e.getX();
    int ypos = e.getY();
    csSampleInfo sInfo = myView.getSampleInfo( xpos, ypos );
    fireSampleInfoEvent( sInfo );

    if( myView.getMouseMode() == csMouseModes.ZOOM_MODE || myView.getMouseMode() == csMouseModes.SPECTRUM_MODE ) {
      if( SwingUtilities.isLeftMouseButton(e) ) {
        myRubberBand.updateEndPoint(new Point(xpos,ypos));
        myView.repaint();
      }
    }
    else if( myView.getMouseMode() == csMouseModes.PAN_MODE ) {
      if( SwingUtilities.isLeftMouseButton(e) || SwingUtilities.isMiddleMouseButton(e) ) {
        int dx = xpos - myPressedPoint.x;
        int dy = ypos - myPressedPoint.y;
        if( dx != 0 || dy != 0 ) firePanningEvent( -dx, -dy );
      }
    }
    myPressedPoint = e.getPoint();

//    if( !myView.hasFocus() ) {
//      myIsSHIFTPressed = false;
//      myIsCTRLPressed  = false;
//    }
  }
//----------------------------------------------------------------
  @Override
  public void mouseMoved( MouseEvent e ) {
    int xpos = e.getX();
    int ypos = e.getY();

    csSampleInfo sInfo = myView.getSampleInfo( xpos, ypos );
    
    fireSampleInfoEvent( sInfo );
    if( sInfo != null ) {
      myView.mouseMoved( sInfo );
    }
  }
//----------------------------------------------------------------
  @Override
  public void mousePressed( MouseEvent e ) {
    int xpos = e.getX();
    int ypos = e.getY();
    myPressedPoint = e.getPoint();
//    csSampleInfo sInfo = myView.getSampleInfo( xpos, ypos );

    if(SwingUtilities.isRightMouseButton(e)) {
      myView.getPopupMenu().show( myView, xpos, ypos );
    }
    // Initiate dragging/panning of seismic view:
    else if( myView.getMouseMode() == csMouseModes.PAN_MODE ) {
      if( SwingUtilities.isLeftMouseButton(e) || SwingUtilities.isMiddleMouseButton(e) ) {
        myView.setCursor(new Cursor(Cursor.MOVE_CURSOR));
      }
    }
    else if( SwingUtilities.isLeftMouseButton(e) && myView.getMouseMode() == csMouseModes.ZOOM_MODE || myView.getMouseMode() == csMouseModes.SPECTRUM_MODE ) {
      myRubberBand.updateStartPoint(myPressedPoint);
      myView.addOverlay(myRubberBand);
    }
  }
//----------------------------------------------------------------
  @Override
  public void mouseReleased( MouseEvent e ) {
    if( myView.getMouseMode() == csMouseModes.PAN_MODE ) {
      myView.setCursor( csMouseModes.PAN_CURSOR );
    }
    else if( SwingUtilities.isLeftMouseButton(e) && ( myView.getMouseMode() == csMouseModes.ZOOM_MODE || myView.getMouseMode() == csMouseModes.SPECTRUM_MODE ) ) {
      csSampleInfo sInfoStart = myView.getSampleInfo( myRubberBand.startPos().x, myRubberBand.startPos().y );
      csSampleInfo sInfoEnd   = myView.getSampleInfo( myRubberBand.endPos().x, myRubberBand.endPos().y );
      fireRubberBandEvent(sInfoStart,sInfoEnd);
      myView.removeOverlay(myRubberBand);
      myView.repaint();
    }
//    else if( SwingUtilities.isMiddleMouseButton(e) && myView.getMouseMode() == csMouseModes.ZOOM_MODE ) {
//    }
//    if( !myView.hasFocus() ) {
//      myIsSHIFTPressed = false;
//      myIsCTRLPressed = false;
//    }
  }

//----------------------------------------------------------------
  @Override
  public void mouseEntered( MouseEvent e ) {
    myIsMouseInWindow = true;
    if( myView.getMouseMode() == csMouseModes.PAN_MODE ) {
      myView.setCursor( csMouseModes.PAN_CURSOR );
    }
    else if( myView.getMouseMode() == csMouseModes.ZOOM_MODE ) {
      myView.setCursor( csMouseModes.ZOOM_CURSOR );
    }
    else if( myView.getMouseMode() == csMouseModes.SPECTRUM_MODE ) {
      myView.setCursor( csMouseModes.SPECTRUM_CURSOR );
    }
    else if( myView.getMouseMode() == csMouseModes.KILL_MODE ) {
      myView.setCursor( csMouseModes.KILL_CURSOR );
    }
    else if( myView.getMouseMode() == csMouseModes.PICK_MODE ) {
      myView.setCursor( csMouseModes.PICK_CURSOR );
    }
    else if( myView.getMouseMode() == csMouseModes.PAINT_MODE ) {
      myView.setCursor( csMouseModes.PAINT_CURSOR );
    }
    else if( myView.getMouseMode() == csMouseModes.VELPICK_MODE ) {
      myView.setCursor( csMouseModes.VELPICK_CURSOR );
    }
    else { //( myView.getMouseMode() == csMouseModes.NO_MODE ) {
      myView.setCursor( Cursor.getDefaultCursor() );
    }
  }
//----------------------------------------------------------------
  @Override
  public void mouseExited( MouseEvent e ) {
    myIsMouseInWindow = false;
    if( myView.hasFocus() ) {
      myView.transferFocus(); // release focus
    }
//    myIsSHIFTPressed = false;
//    myIsCTRLPressed  = false;
    myView.mouseExited();
  }

//----------------------------------------------------------------
  @Override
  public void mouseClicked( MouseEvent e ) {
    int xpos = e.getX();
    int ypos = e.getY();
    csSampleInfo sInfo = myView.getSampleInfo( xpos, ypos );
    fireSampleInfoClickEvent( sInfo );
  }

//----------------------------------------------------------------
  @Override
  public void keyPressed( KeyEvent e ) {
    if( !myIsMouseInWindow ) return;
//    myIsSHIFTPressed = e.isShiftDown();
//    myIsCTRLPressed  = e.isControlDown();
    fireKeyPressedEvent( e );
  }

//----------------------------------------------------------------
  @Override
  public void keyReleased( KeyEvent e ) {
    if( !myIsMouseInWindow ) return;
//    myIsSHIFTPressed = e.isShiftDown();
//    myIsCTRLPressed  = e.isControlDown();
    fireKeyReleasedEvent( e );
  }

//----------------------------------------------------------------
  @Override
  public void keyTyped( KeyEvent e ) {
    if( !myIsMouseInWindow ) return;
  }
//-------------------------------------------------------------------
  public void addRubberBandListener( csIRubberBandListener listener ) {
    if( listener != null && !myRubberBandListeners.contains( listener ) ) {
      myRubberBandListeners.add( listener );
    }
  }
//-------------------------------------------------------------------
  public void addSampleInfoListener( csISampleInfoListener listener ) {
    if( listener != null && !mySampleInfoListeners.contains( listener ) ) {
      mySampleInfoListeners.add( listener );
    }
  }
//-------------------------------------------------------------------
  public void removeSampleInfoListener( csISampleInfoListener listener ) {
    if( mySampleInfoListeners.contains( listener ) ) {
      mySampleInfoListeners.remove( listener );
    }
  }
//-------------------------------------------------------------------
  public void addKeyListener( csIKeyListener listener ) {
    if( listener != null && !myKeyListeners.contains( listener ) ) {
      myKeyListeners.add( listener );
    }
  }
//-------------------------------------------------------------------
  public void addPanningListener( csIPanningListener listener ) {
    if( listener != null && !myPanningListeners.contains( listener ) ) {
      myPanningListeners.add( listener );
    }
  }
//-------------------------------------------------------------------
  public void removeKeyListener( csIKeyListener listener ) {
    if( myKeyListeners.contains( listener ) ) {
      myKeyListeners.remove( listener );
    }
  }

//-------------------------------------------------------------------
  private void fireRubberBandEvent( csSampleInfo infoStart, csSampleInfo infoEnd ) {
    for( int i = 0; i < myRubberBandListeners.size(); i++ ) {
      csIRubberBandListener target = myRubberBandListeners.get( i );
      target.rubberBandCompleted( infoStart, infoEnd );
    }
  }
//-------------------------------------------------------------------
  private void fireSampleInfoEvent( csSampleInfo info ) {
    for( int i = 0; i < mySampleInfoListeners.size(); i++ ) {
      csISampleInfoListener target = mySampleInfoListeners.get( i );
      target.mouseMoved( myView, info );
    }
  }
//-------------------------------------------------------------------
  private void fireSampleInfoClickEvent( csSampleInfo info ) {
    for( int i = 0; i < mySampleInfoListeners.size(); i++ ) {
      csISampleInfoListener target = mySampleInfoListeners.get( i );
      target.mouseClicked( myView, info );
    }
  }
//-------------------------------------------------------------------
  private void fireKeyReleasedEvent( KeyEvent event) {
    for( int i = 0; i < myKeyListeners.size(); i++ ) {
      csIKeyListener target = myKeyListeners.get( i );
      target.keyReleased( event );
    }
  }
//-------------------------------------------------------------------
  private void fireKeyPressedEvent( KeyEvent event) {
    for( int i = 0; i < myKeyListeners.size(); i++ ) {
      csIKeyListener target = myKeyListeners.get( i );
      target.keyPressed( event );
    }
  }
//-------------------------------------------------------------------
  private void firePanningEvent( int dx, int dy ) {
    for( int i = 0; i < myPanningListeners.size(); i++ ) {
      csIPanningListener target = myPanningListeners.get( i );
      target.hasPanned( dx, dy );
    }
  }
  //------------------------------------------------------------
  public void finalize() {
    myView      = null;
  }
}


