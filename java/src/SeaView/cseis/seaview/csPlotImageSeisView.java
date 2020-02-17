/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.seis.csTraceBuffer;
import cseis.seisdisp.csSeisView;
import cseis.seisdisp.csISeisOverlay;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import javax.swing.JFrame;

/**
 * Seismic display used for image plot.
 * @author 2009 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csPlotImageSeisView extends csSeisView {
  private int myNewHeight;
  private int myNewWidth;
  private int myNewMarginLeftRight;
  private int myNewMarginTopBottom;
  private int myMinSample;
//  int myMaxSample;
  private int myViewMin;

  public csPlotImageSeisView( JFrame parentFrame, csTraceBuffer buffer, double sampleInt ) {
    super( parentFrame, buffer, sampleInt );
    myMinSample = 0;
//    myMaxSample = myNumSamples - 1;
  }
  public int set( int width, int height, int marginTopBottom, int marginLeftRight, boolean evenTraceSpacing,
      float minTime, float maxTime, int colorBarWidth, int colorBarHeight ) {
    myNewHeight = height;
    myNewWidth  = width;
    myNewMarginLeftRight = marginLeftRight;  // 
    myNewMarginTopBottom = marginTopBottom;  //
//    myColorBarOption = colorBarOption;
//    myColorBarStep   = colorBarStep;
    myMinSample = (int)( minTime / mySampleInt + 0.5f );
    if( myMinSample < 0 ) myMinSample = 0;
    int maxSample = (int)( maxTime / mySampleInt + 0.5f );
    if( maxSample > myNumSamples-1 ) maxSample = myNumSamples-1;

    if( mySettings.isLogScale ) {
      // CHANGE: Log scale not working yet...
      resetScreenBuffer();
    }

    myMarginLeftRight = 0;
    myMarginTopBottom = 0;
    int xViewMax = (int)xModel2View( myNumTraces-1 ) + myViewPositionHorz;
    int yViewMin = (int)yModel2View( myMinSample ) + myViewPositionVert;
    int yViewMax = (int)yModel2View( maxSample ) + myViewPositionVert;

    mySettings.zoomHorz *= ((float)(myNewWidth-2*myNewMarginLeftRight)/(float)xViewMax);
    mySettings.zoomVert *= ((float)(myNewHeight-2*myNewMarginTopBottom)/(float)(yViewMax-yViewMin));

//    System.out.println("max view: " + xViewMax + " " + yViewMax);
//    System.out.println("zoom horz: " + mySettings.zoomHorz + " " + xViewMax + " vert: " + mySettings.zoomVert);
    if( evenTraceSpacing ) {
      if( mySettings.zoomHorz < 1 ) {
        int evenZoom = (int)( 1.0f/mySettings.zoomHorz - 0.01f ) + 1;
        if( evenZoom % 2 != 0 ) evenZoom += 1;
        mySettings.zoomHorz = 1.0f/( (float)( evenZoom ) );
      }
      else {
        mySettings.zoomHorz = (float)( (int)(mySettings.zoomHorz - 0.01f) + 1 );
      }
    }
    myNewWidth = (int)xModel2View( myNumTraces-1 ) + 2*myNewMarginLeftRight + myViewPositionHorz;
    myViewMin = (int)yModel2View( myMinSample ) + myViewPositionVert;

    return myNewWidth;
  }
  public void paintImage( Graphics2D g2 ) {
    myMarginLeftRight = myNewMarginLeftRight;
    myMarginTopBottom = myNewMarginTopBottom;
    myViewPositionHorz = 0;
    myViewPositionVert = (int)yModel2View( myMinSample );
    myViewPositionVert = myViewMin;
    super.resetScreenBuffer();
    myIsFullRepaint = true;
    super.repaintStep1( false );
    g2.drawImage( myBitmap, 0, 0, this );
    g2.setColor(  Color.white );
    g2.fillRect( 0, 0, myNewWidth, myNewMarginTopBottom );
    g2.fillRect( 0, myNewHeight-myNewMarginTopBottom, myNewWidth, myNewMarginTopBottom );

    for( int ioverlay = 0; ioverlay < mySeisOverlays.size(); ioverlay++ ) {
      csISeisOverlay ol = mySeisOverlays.get(ioverlay);
      ol.draw( this, g2 );
    }
 
  }
  /**
   * Override in order to control image size
   */
  public Rectangle getVisibleRect() {
    Rectangle rect = new Rectangle( myNewWidth, myNewHeight );
    return rect;
  }
// Does not work:
//  public Graphics2D paintIt() {
//    super.resetScreenBuffer();
//    myIsFullRepaint = true;
//    super.repaintStep1( false );
//    return myBitmap.createGraphics();
//  }
}


