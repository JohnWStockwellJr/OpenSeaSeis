/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.geom.GeneralPath;
import java.awt.*;

/**
 * Rubber band overlay.
 * @author Bjorn Olofsson
 */
public class csRubberBandOverlay implements csISeisOverlay {
  private Color myColor;
  private BasicStroke myStroke;
  private int myNumGridLines;
  private Point myPosStart;
  private Point myPosEnd;
  private boolean myIsInverted;

  public csRubberBandOverlay() {
    myNumGridLines = 4;
    myStroke = new BasicStroke(1.0f);
    myColor = Color.red;
    myPosStart = new Point(-1,-1);
    myPosEnd   = new Point(-1,-1);
    myIsInverted = true;
  }
  public java.awt.Point startPos() {
    return myPosStart;
  }
  public java.awt.Point endPos() {
    return myPosEnd;
  }
  public void setColor( Color color ) {
    setColor( color, 255 );
  }
  public void setColor( Color color, int transparency ) {
    if( transparency == 255 ) {
      myColor = color;
    }
    else {
      myColor = new Color( color.getRed(), color.getGreen(), color.getBlue(), transparency );
    }
    myIsInverted = false;
  }
  public void setLineWidth( float lineWidth ) {
    myStroke = new BasicStroke(lineWidth,java.awt.BasicStroke.CAP_BUTT,java.awt.BasicStroke.JOIN_ROUND);
  }
  public void updateStartPoint( Point posStart ) {
    myPosStart = posStart;
    myPosEnd   = posStart;
  }
  public void updateEndPoint( Point pos ) {
    myPosEnd   = pos;
  }
  @Override
  public void draw( csSeisView seisview, Graphics2D g ) {
    RenderingHints rhints = g.getRenderingHints();
    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

    // XOR current color with 0x00aaaaff instead of 0x00ffffff to shift color towards yellow and avoid blank in grey areas.
    if( myIsInverted ) {
      g.setXORMode( new Color( 0x00aaaaff ) );
    }
    else {
      g.setColor( myColor );
    }
    g.setStroke( myStroke );

    GeneralPath path = new GeneralPath();
    float xrange = myPosEnd.x - myPosStart.x;
    for( int ix = 0; ix < myNumGridLines; ix++ ) {
      float xp = (float)myPosStart.x + (float)ix / (float)(myNumGridLines-1) * xrange;
      path.moveTo( xp, myPosStart.y );
      path.lineTo( xp, myPosEnd.y );
    }
    float yrange = myPosEnd.y - myPosStart.y;
    for( int iy = 0; iy < myNumGridLines; iy++ ) {
      float yp = (float)myPosStart.y + (float)iy / (float)(myNumGridLines-1) * yrange;
      path.moveTo( myPosStart.x, yp );
      path.lineTo( myPosEnd.x, yp );
    }
    g.draw(path);

    if( myIsInverted ) {
      g.setPaintMode();  // Revert from XOR mode
    }
    g.setRenderingHints(rhints);
  }
}


