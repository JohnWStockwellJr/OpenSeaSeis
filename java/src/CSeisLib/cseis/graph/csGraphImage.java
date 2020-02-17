/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

import java.awt.*;
import java.awt.image.*;
import java.awt.geom.*;

/**
 * ..under development, not fully functional yet, not used by any class.
 */
public class csGraphImage extends BufferedImage {
  float myXScalar;
  float myXMin;
  float myXInset;
  float myYScalar;
  float myYMin;
  float myYInset;
  float myMaxValue;

  public csGraphImage(int width, int height, int imageType) {
    super(width, height, imageType);
    
    myXScalar = 1.0f;
    myYScalar = 1.0f;
    myXInset = 0;
    myYInset = 0;
    myXMin = 0;
    myYMin = 0;
  }
  public void redraw( float[] xValues, float[] yValues, int startSamp, int nValues ) {
    Graphics2D g = (Graphics2D)getGraphics();
    int circleSize = 11;

    g.setRenderingHint( RenderingHints.KEY_ANTIALIASING,  RenderingHints.VALUE_ANTIALIAS_ON );
    // Clear screen
    g.setColor( Color.white );
    g.fillRect( 0, 0, getWidth(), getHeight() );

    // Coordinate system
    g.setColor( Color.gray );
    GeneralPath pathAxes = new GeneralPath();
    pathAxes.moveTo( xData2Image(-myMaxValue), yData2Image(0) );
    pathAxes.lineTo( xData2Image(+myMaxValue), yData2Image(0) );
    pathAxes.moveTo( xData2Image(0), yData2Image(-myMaxValue) );
    pathAxes.lineTo( xData2Image(0), yData2Image(+myMaxValue) );
    g.draw( pathAxes );
//    System.out.println("Samples " + xData2Image(-myMaxValue) + ", " + yData2Image(-myMaxValue) + ", " + myXMin + ", " + myMaxValue );
    
    int endSamp = startSamp+nValues-1;
    g.setColor( Color.black );
    g.setStroke( new BasicStroke(2) );
    GeneralPath path = new GeneralPath();
    path.moveTo( xData2Image(xValues[startSamp]), yData2Image(yValues[startSamp]) );
    Shape shape = new Ellipse2D.Float( xData2Image(xValues[startSamp])-circleSize/2,
        yData2Image(yValues[startSamp])-circleSize/2, circleSize, circleSize );
    g.draw( shape );
    for( int i = startSamp+1; i <= endSamp; i++ ) {
      path.lineTo( xData2Image(xValues[i]), yData2Image(yValues[i]) );
      shape = new Ellipse2D.Float( xData2Image(xValues[i])-circleSize/2,
          yData2Image(yValues[i])-circleSize/2, circleSize, circleSize );
      g.draw( shape );
//      System.out.println("Data: " + xValues[i] + ", " + yValues[i] + "   xpos ypos: " + xData2Image(xValues[i]) + ", " + yData2Image(yValues[i]));
//      System.out.println("Samples " + i + ", " + xData2Image(xValues[i]) + ", " + yData2Image(yValues[i]) );
    }
    g.draw( path );
  
  }
  public float xData2Image( float data ) {
    return( (data-myXMin)*myXScalar + myXInset );
  }
  public float yData2Image( float data ) {
    return( (data-myYMin)*myYScalar + myYInset );
  }
  public void setXScaling( int inset, float min, float max ) {
    float diff = 2* Math.max( Math.abs(max), Math.abs(min) );
    if( diff == 0.0f ) diff = 1.0f;
    myMaxValue = diff/2;
    myXScalar = (getWidth()-2*inset) / diff;
    myXMin    = min;
    myXInset  = inset;
  }
  public void setYScaling( int inset, float min, float max ) {
    float diff = 2* Math.max( Math.abs(max), Math.abs(min) );
    if( diff == 0.0f ) diff = 1.0f;
    myMaxValue = diff/2;
    myYScalar = (getHeight()-2*inset) / diff;
    myYMin    = min;
    myYInset  = inset;
  }
}


