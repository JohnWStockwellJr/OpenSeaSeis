/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import javax.swing.SwingUtilities;

/**
 */
public class csPaintOverlay implements csISeisOverlay, MouseMotionListener, MouseListener {
  final private csSeisView mySeisView;
  final private ArrayList<Point2D.Float> myLinePoints;
  final private ArrayList<Point2D.Float> myPoints;
  
  private Color myColorLine;
  private Color myColorPoint;
  private int mySizeLine;
  private int mySizePoint;

  public csPaintOverlay( csSeisView seisView ) {
    mySeisView   = seisView;
    myLinePoints = new ArrayList<Point2D.Float>();
    myPoints     = new ArrayList<Point2D.Float>();
    mySizeLine   = 6;
    mySizePoint  = 11;
    myColorLine  = Color.red;
    myColorPoint = Color.blue;
  }
  public void activate() {
    mySeisView.addMouseListener( this );
    mySeisView.addMouseMotionListener( this );
  }
  public void deactivate() {
    mySeisView.removeMouseListener( this );
    mySeisView.removeMouseMotionListener( this );
  }
  public void update( int lineSize, Color lineColor, int pointSize, Color pointColor ) {
    mySizeLine   = lineSize;
    mySizePoint  = pointSize;
    myColorLine  = lineColor;
    myColorPoint = pointColor;
    mySeisView.repaint();
  }
  public void clearOverlay() {
    myPoints.clear();
    myLinePoints.clear();
    mySeisView.repaint();
  }
  @Override
  public void draw( csSeisView seisview, Graphics2D g ) {
    if( seisview.getTraceBuffer() == null ) return;
    RenderingHints rhints_save = g.getRenderingHints();
    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

    // 1) Draw circles
    int sizeHalf  = mySizePoint/2;
    int sizePoint = mySizeLine/2;
    if( sizePoint < 1 ) sizePoint = 1;
    g.setColor( myColorPoint );
    g.setStroke( new BasicStroke(sizePoint) );
    for( Point2D.Float pointModel : myPoints ) {
      float xView = seisview.xModel2View( pointModel.x );
      float yView = seisview.yModel2View( pointModel.y );
      g.draw( new Ellipse2D.Double(xView-sizeHalf, yView-sizeHalf, mySizePoint, mySizePoint) );
		}
    // 2) Draw lines
    GeneralPath path = new GeneralPath();
    boolean isGap = true;
    for( Point2D.Float pointModel : myLinePoints ) {
      if( pointModel.x == -1 && pointModel.y == -1 ) isGap = true;
      else {
        float xView = seisview.xModel2View( pointModel.x );
        float yView = seisview.yModel2View( pointModel.y );
        if( isGap ) {
          path.moveTo( xView, yView );
          isGap = false;
        }
        else {
          path.lineTo( xView, yView );
        }
      }
		}
    g.setColor( myColorLine );
    g.setStroke( new BasicStroke(mySizeLine) );
    g.draw( path );

    g.setRenderingHints(rhints_save);
  }
  //-----------------------------------------------------------
  // Mouse is dragged in seismic view
  //
  @Override
  public void mouseDragged( MouseEvent e ) {
    if( SwingUtilities.isLeftMouseButton(e) ) {
      Point pointView = e.getPoint();
      csSampleInfo sInfo = mySeisView.getSampleInfo( pointView.x, pointView.y );
      Point2D.Float pointModel = new Point2D.Float( (float)sInfo.traceDouble, (float)sInfo.sampleDouble );
      if( !myLinePoints.isEmpty() ) {
        Point2D.Float pointPrevious = myLinePoints.get( myLinePoints.size()-1 );
        if( pointPrevious.x == pointModel.x && pointPrevious.y == pointModel.y ) return;
      }
      myLinePoints.add( pointModel );
      mySeisView.repaint();
    }
  }
  //-----------------------------------------------------------
  // Mouse button is clicked in seismic view
  //
  @Override
  public void mouseClicked(MouseEvent e) {
    if( SwingUtilities.isMiddleMouseButton(e) ) {
      clearOverlay();
    }
    else if( SwingUtilities.isLeftMouseButton(e) ) {
      Point pointView = e.getPoint();
      csSampleInfo sInfo = mySeisView.getSampleInfo( pointView.x, pointView.y );
      myPoints.add( new Point2D.Float( (float)sInfo.traceDouble, (float)sInfo.sampleDouble ) );
      mySeisView.repaint();
    }
  }
  //-----------------------------------------------------------
  // Mouse button is released in seismic view
  //
  @Override
  public void mouseReleased(MouseEvent e) {
    if( myLinePoints.isEmpty() ) return;
    Point2D.Float pointPrevious = myLinePoints.get( myLinePoints.size()-1 );
    if( pointPrevious.x != -1 || pointPrevious.y != -1 ) {
      Point2D.Float pointModel = new Point2D.Float( -1, -1 );
      myLinePoints.add( pointModel );
    }
  }

  @Override
  public void mouseMoved(MouseEvent e) {
  }
  @Override
  public void mousePressed(MouseEvent e) {
  }
  @Override
  public void mouseEntered(MouseEvent e) {
  }
  @Override
  public void mouseExited(MouseEvent e) {
  }
  
//  public class PointAttr {
//    int lineSize;
//    int pointSize;
//    Color lineColor;
//    Color pointColor;
//    Point2D.Float point;
//    PointAttr( int theLineSize, Color theLineColor, int thePointSize, Color thePointColor, float x, float y ) {
//      lineSize   = theLineSize;
//      pointSize  = thePointSize;
//      lineColor  = theLineColor;
//      pointColor = thePointColor;
//      point = new Point2D.Float( x, y );
//    }
//  }
}

