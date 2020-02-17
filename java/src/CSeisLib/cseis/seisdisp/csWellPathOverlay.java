/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.seisdisp;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.util.ArrayList;

import cseis.math.csPoint3D;
import cseis.seis.csHeader;
import cseis.seis.csISeismicTraceBuffer;

/**
 */
public class csWellPathOverlay implements csISeisOverlay {
  final private csSeisView mySeisView;
  final private ArrayList<WellData> myWellDataList;
  private csWellPathOverlayAttr myOverlayAttr;
  
  private boolean myIsDrawing = false;

  private int myHdrIndexX;
  private int myHdrIndexY;

  public csWellPathOverlay( csSeisView seisView ) {
    mySeisView   = seisView;
    myWellDataList = new ArrayList<WellData>();
    myHdrIndexX = -1;
    myHdrIndexY = -1;
    myOverlayAttr = new csWellPathOverlayAttr();

    mySeisView.addSeisViewListener( new csISeisViewListener() {
      @Override
      public void changedSettings( csSeisDispSettings settings ) {}
      @Override
      public void vertScrollChanged( int scrollValue ) {}
      @Override
      public void horzScrollChanged( int scrollValue ) {}
      @Override
      public void sizeChanged( Dimension size ) {}
      @Override
      public void traceBufferChanged( csISeismicTraceBuffer traceBuffer ) {
        for( WellData wellData : myWellDataList ) {
          updateSeismicBuffer( wellData, traceBuffer );
        }
      }
      public void mouseExited() {}
    } );

  }
  public void updateHeaders( int hdrID_x, int hdrID_y ) {
    myHdrIndexX = hdrID_x;
    myHdrIndexY = hdrID_y;
  }
  public void update( csWellPathOverlayAttr overlayAttr, ArrayList<csWellPathAttr> attrList ) {
    myOverlayAttr = new csWellPathOverlayAttr( overlayAttr );
    for( csWellPathAttr attr_in : attrList ) {
      boolean found = false;
      for( WellData wellData : myWellDataList ) {
        if( attr_in.id == wellData.attr.id ) {
          found = true;
          if( wellData.attr.isDifferent(attr_in) ) {
            wellData.attr.setAttributes( attr_in );
          }
          break;
        }
      } // END: for current
      if( !found ) {
        WellData wellData = new WellData( attr_in );
        updateSeismicBuffer( wellData, mySeisView.getTraceBuffer() );
        myWellDataList.add( wellData );
      }
    } // END for attr_in
    // Remove well paths if needed
    for( WellData wellData : myWellDataList ) {
      boolean found = false;
      for( csWellPathAttr attr_in : attrList ) {
        if( attr_in.id == wellData.attr.id ) {
          found = true;
          break;
        }
      }
      if( !found ) {
        myWellDataList.remove( wellData );
        break;
      }
    }

    mySeisView.invalidate();
    mySeisView.repaint();
  }
  public csWellPathOverlayAttr getOverlayAttr() {
    return myOverlayAttr;
  }
  //--------------------------------------------------------------------------------------------------------------
  //
  private void updateSeismicBuffer( WellData wellData, csISeismicTraceBuffer buffer ) {
    while( myIsDrawing ) {
      try {
        Thread.sleep(100);
      } catch( Exception e ) {}
    }
    wellData.pointList.clear();

    // Set XY positions for each seismic trace
    Point2D.Double tracePositions[] = new Point2D.Double[buffer.numTraces()];
    for( int itrc = 0; itrc < buffer.numTraces(); itrc++ ) {
      csHeader[] headerValues = buffer.headerValues(itrc);
      tracePositions[itrc] = new Point2D.Double( headerValues[myHdrIndexX].doubleValue(), headerValues[myHdrIndexY].doubleValue() );
    }

    // Convert well positions into trace & sample positions of seismic section
    double sampleInt = mySeisView.getSampleInt();
    for( csPoint3D wellPos : wellData.attr.posXYZList ) {
      WellPointInSeismic wellPoint = new WellPointInSeismic();
      wellPoint.sampleDouble = wellPos.z / sampleInt;
      computeNearestTrace( wellPos, tracePositions, wellPoint );
      wellData.pointList.add( wellPoint );
    }

    // Decimate points: Maximum 2 points per sample and/or trace (purpose: speed up drawing)
//    System.out.println(" Number of points before decimation: " + wellData.pointList.size() );
    int counter = wellData.pointList.size()-1;
    WellPointInSeismic wellPoint_prev = wellData.pointList.get( counter );
    counter -= 1;
    while( counter > 0 ) {
      WellPointInSeismic wellPoint = wellData.pointList.get(counter);
      double traceDiff  = Math.abs( wellPoint.traceDouble  - wellPoint_prev.traceDouble );
      double sampleDiff = Math.abs( wellPoint.sampleDouble - wellPoint_prev.sampleDouble );
      if( traceDiff < 0.5 && sampleDiff < 0.5 ) {
        wellData.pointList.remove( counter );
      }
      else {
        wellPoint_prev = wellPoint;
      }
      counter -= 1;
    }
//    System.out.println(" Number of points AFTER decimation: " + wellData.pointList.size() + "\n");
    
  }
  //---------------------------------------------------------------------------------------------------------------
  //
  private void computeNearestTrace( csPoint3D wellPos, Point2D.Double[] pSeismic, WellPointInSeismic wellPointOut ) {

    // 1) Find trace that is located closest to well point
    double minDistance = Math.sqrt( Math.pow(wellPos.x-pSeismic[0].x,2.0) + Math.pow(wellPos.y-pSeismic[0].y,2.0) );
    int trcIndexMinDist = 0;
    for( int ipos = 1; ipos < pSeismic.length; ipos++ ) {
      double distance = Math.sqrt( Math.pow(wellPos.x-pSeismic[ipos].x,2.0) + Math.pow(wellPos.y-pSeismic[ipos].y,2.0) );
      if( distance < minDistance ) {
        minDistance  = distance;
        trcIndexMinDist = ipos;
      }
    }
    wellPointOut.distance = minDistance;
    if( trcIndexMinDist == 0 || trcIndexMinDist == pSeismic.length-1 ) {
      // Well point is at edge of seismic
      wellPointOut.traceDouble = (double)trcIndexMinDist;
      wellPointOut.isInFront = false;
      return;
    }

    Point2D.Double ps1 = new Point2D.Double( pSeismic[trcIndexMinDist-1].x, pSeismic[trcIndexMinDist-1].y );
    Point2D.Double ps2 = new Point2D.Double( pSeismic[trcIndexMinDist+1].x, pSeismic[trcIndexMinDist+1].y );
    Point2D.Double pw  = new Point2D.Double( wellPos.x, wellPos.y );
    Point2D.Double as  = new Point2D.Double(); // Gradient
    Point2D.Double aw  = new Point2D.Double(); // Gradient

    as.x = ps2.x - ps1.x;
    as.y = ps2.y - ps1.y;
    double norm = Math.sqrt( as.x*as.x + as.y*as.y );
    if( Math.abs(norm) < 1e-20 ) return;
    as.x /= norm;
    as.y /= norm;

    aw.x = -as.y;
    aw.y = as.x;

    double delta_x = ps1.x - pw.x;
    double delta_y = ps1.y - pw.y;

    double temp1 = delta_y*aw.x - delta_x*aw.y;
    double temp2 = as.x*aw.y - as.y*aw.x;

    if( Math.abs(temp2) < 1e-20 ) {
      return;
    }
   
    double ls = temp1 / temp2;
    
    Point2D.Double intersection = new Point2D.Double( ps1.x + ls * as.x, ps1.y + ls * as.y );

    double distIntersect = Math.sqrt( Math.pow( intersection.x-ps1.x, 2.0 ) + Math.pow( intersection.y-ps1.y, 2.0 ) );
    double dist2Traces   = Math.sqrt( Math.pow( ps2.x-ps1.x, 2.0 ) + Math.pow( ps2.y-ps1.y, 2.0 ) );

    double weight = distIntersect / dist2Traces;

    wellPointOut.traceDouble = (double)(trcIndexMinDist-1) + 2.0*weight;

    // 3) Compute if wellPoint is in front or behind seismic section
    Point2D.Double vectorWell = new Point2D.Double( intersection.x - pw.x, intersection.y - pw.y );
    double crossProduct = vectorWell.x * as.y - vectorWell.y * as.x;
    wellPointOut.isInFront   = (crossProduct <= 0);
  }

  //------------------------------------------------------------------------------------
  @Override
  public void draw( csSeisView seisview, Graphics2D g ) {
    if( seisview.getTraceBuffer() == null || myWellDataList.isEmpty() ) return;
    if( myIsDrawing ) return;
    myIsDrawing = true;
    RenderingHints rhints_save = g.getRenderingHints();
    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

    // 2) Draw lines between well positions
    for( WellData wellData : myWellDataList ) {
      ArrayList<WellPointInSeismic> wellPointList = wellData.pointList;
      if( !wellData.attr.showWell || wellPointList.isEmpty() ) continue;
      GeneralPath pathSolid  = new GeneralPath();
      GeneralPath pathDashed = new GeneralPath();
      GeneralPath pathShadow = new GeneralPath();
      boolean isInFront = wellPointList.get(0).isInFront;
      float xView0 = seisview.xModel2View( (float)wellPointList.get(0).traceDouble );
      float yView0 = seisview.yModel2View( (float)wellPointList.get(0).sampleDouble );
      pathDashed.moveTo( xView0, yView0 );
      pathSolid.moveTo( xView0, yView0 );
      if( myOverlayAttr.showShadow ) pathShadow.moveTo( xView0, yView0 );
  
      for( WellPointInSeismic pointModel : wellPointList ) {
        float xView = seisview.xModel2View( (float)pointModel.traceDouble );
        float yView = seisview.yModel2View( (float)pointModel.sampleDouble );
        boolean isWithinDistance = ( pointModel.distance <= myOverlayAttr.maxDisplayDistance );
        if( isWithinDistance ) {
          if( isInFront ) {
            pathSolid.lineTo( xView, yView );
            if( myOverlayAttr.showShadow ) {
              pathShadow.lineTo( xView, yView + myOverlayAttr.shadowDepth*(pointModel.distance/myOverlayAttr.maxDisplayDistance) );
            }
          }
          else {
            pathDashed.lineTo( xView, yView );
          }
          if( pointModel.isInFront != isInFront ) {
            if( pointModel.isInFront ) {
              pathSolid.moveTo( xView, yView );
              if( myOverlayAttr.showShadow ) {
                pathShadow.moveTo( xView, yView + myOverlayAttr.shadowDepth*(pointModel.distance/myOverlayAttr.maxDisplayDistance) );
              }
            }
            else {
              pathDashed.moveTo( xView, yView );
            }
            isInFront = pointModel.isInFront;
          }
        }
        else { // not within distance
          pathSolid.moveTo( xView, yView );
          if( myOverlayAttr.showShadow ) {
            pathShadow.moveTo( xView, yView + myOverlayAttr.shadowDepth*(pointModel.distance/myOverlayAttr.maxDisplayDistance) );
          }
          pathDashed.moveTo( xView, yView );
        }
  		}
  
      csWellPathAttr attr = wellData.attr;
      int red = attr.color.getRed();
      int blue = attr.color.getBlue();
      int green = attr.color.getGreen();
      g.setStroke( new BasicStroke(myOverlayAttr.lineSize) );
  
      g.setColor( new Color( 0, 0, 0, myOverlayAttr.transparency/2 ) );
      g.draw( pathShadow );
  
      g.setColor( new Color( red, green, blue, myOverlayAttr.transparency ) );
      g.draw( pathSolid );
      float[] dashPattern = { 10, 15 };
      g.setStroke( new BasicStroke(myOverlayAttr.lineSize, BasicStroke.CAP_ROUND, BasicStroke.JOIN_MITER, 10, dashPattern, 0 ) );
      g.draw( pathDashed );
    }
    g.setRenderingHints(rhints_save);
    myIsDrawing = false;
  }
/*
  private void draw2( csSeisView seisview, Graphics2D g ) {
    float[] dashPattern = { 10, 10 };

    // 2) Draw lines between well positions
    for( WellData wellData : myWellDataList ) {
      ArrayList<WellPointInSeismic> wellPointList = wellData.pointList;
      if( !wellData.attr.showWell || wellPointList.isEmpty() ) continue;
      float xViewPrev = seisview.xModel2View( (float)wellPointList.get(0).traceDouble );
      float yViewPrev = seisview.yModel2View( (float)wellPointList.get(0).sampleDouble );
      float yViewShadowPrev = yViewPrev + (float)(myOverlayAttr.shadowDepth*(wellPointList.get(0).distance/myOverlayAttr.maxDisplayDistance));

      csWellPathAttr attr = wellData.attr;
      int red = attr.color.getRed();
      int blue = attr.color.getBlue();
      int green = attr.color.getGreen();

      for( WellPointInSeismic pointModel : wellPointList ) {
        float xView = seisview.xModel2View( (float)pointModel.traceDouble );
        float yView = seisview.yModel2View( (float)pointModel.sampleDouble );
        boolean isWithinDistance = ( pointModel.distance <= myOverlayAttr.maxDisplayDistance ) ? true : false;
        if( isWithinDistance ) {
          Shape line = new Line2D.Double( xViewPrev, yViewPrev, xView, yView );
          int transparency = (int)( (double)myOverlayAttr.transparency * ( 1.0 - 0.9 * (pointModel.distance / myOverlayAttr.maxDisplayDistance) ) );
          System.out.println("Transparency: " + transparency + " " + pointModel.distance);
          if( pointModel.isInFront ) {
            g.setStroke( new BasicStroke(myOverlayAttr.lineSize) );
            if( myOverlayAttr.showShadow ) {
              float yViewShadow = yView + (float)(myOverlayAttr.shadowDepth*(pointModel.distance/myOverlayAttr.maxDisplayDistance));
              Shape lineShadow = new Line2D.Double( xViewPrev, yViewShadowPrev, xView, yViewShadow );
              g.setColor( new Color( 0, 0, 0, transparency/2 ) );
              g.draw( lineShadow );
              yViewShadowPrev = yViewShadow;
            }
            g.setColor( new Color( red, green, blue, transparency ) );
            g.draw( line );
          }
          else {
            g.setStroke( new BasicStroke(myOverlayAttr.lineSize, BasicStroke.CAP_ROUND, BasicStroke.JOIN_MITER, 10, dashPattern, 0 ) );
            g.setColor( new Color( red, green, blue, transparency ) );
            g.draw( line );
          }
        }
        else { // not within distance
          // Nothing to do
        }
        xViewPrev = xView;
        yViewPrev = yView;
      }
    }  
  }
*/
//*****************************************************************************************************
//
  public class WellPointInSeismic {
    double  traceDouble;
    double  sampleDouble;
    double  distance;
    boolean isInFront;
  }
  public class WellData {
    csWellPathAttr attr;
    ArrayList<csWellPathOverlay.WellPointInSeismic> pointList;
    public WellData( csWellPathAttr attr_in ) {
      attr      = new csWellPathAttr( attr_in );
      pointList = new ArrayList<csWellPathOverlay.WellPointInSeismic>();
    }
  }
}
