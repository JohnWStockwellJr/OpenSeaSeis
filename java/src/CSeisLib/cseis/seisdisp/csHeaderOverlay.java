/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import cseis.seis.csHeader;
import cseis.seis.csISeismicTraceBuffer;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;

/**
 * Trace header overlay.<br>
 * This class defines a graphical overlay which can be plotted on top of a seismic view (csSeisView).
 * This is done by implementing the standard interface csISeisOverlay.
 * @author 2011 Bjorn Olofsson
 */
public class csHeaderOverlay implements csISeisOverlay {
  public static int TYPE_LINE   = 0;
  public static int TYPE_CROSS  = 1;
  public static int TYPE_CIRCLE = 2;

  public static final String TEXT_LINE         = "line";
  public static final String TEXT_CROSS        = "cross";
  public static final String TEXT_CIRCLE       = "circle";

  public static final String[] TYPE_TEXT_FIELDS = {
    TEXT_LINE,
    TEXT_CROSS,
    TEXT_CIRCLE
  };

  private BasicStroke myStroke;
  private int myHeaderIndex;
  private String myHeaderName;
  private csHeaderOverlay.Attribute myAttr;
  public boolean myShowLabel;

  public csHeaderOverlay() {
    this("_none_", -1);
  }
  public csHeaderOverlay( String hdrName, int hdrIndex ) {
    setHeader( hdrName, hdrIndex );
    myAttr   = new Attribute();
    myShowLabel = true;
    setLineWidth(myAttr.lineWidth);
  }
  public void setHeader( String hdrName, int hdrIndex ) {
    myHeaderName  = hdrName;
    myHeaderIndex = hdrIndex;
  }
  public csHeaderOverlay.Attribute getAttributes() {
    return myAttr;
  }
  public void setColor( Color color ) {
    setColor( color, myAttr.transparency );
  }
  public void setColor( Color color, int transparency ) {
    myAttr.color = color;
    if( transparency == 255 ) {
      myAttr.color = color;
    }
    else {
      myAttr.color = new Color( color.getRed(), color.getGreen(), color.getBlue(), transparency );
    }
  }
  public void setLayoutAttributes( csHeaderOverlay.Attribute attr ) {
    myAttr.set(attr);
  }
  public void setLayoutAttributes( float multValue, float addValue ) {
    myAttr.multValue = multValue;
    myAttr.addValue  = addValue;
  }
  public void setType( int type ) {
    if( type >= 0 && type < TYPE_TEXT_FIELDS.length ) {
      myAttr.type = type;
    }
  }
  public void setTraceInterval( int trcInterval ) {
    myAttr.trcInterval = trcInterval;
  }
  public void setTransparency( int transparency ) {
    setColor( myAttr.color, transparency );
  }
  public void setShowLabel( boolean doShow ) {
    myShowLabel = doShow;
  }
  public void setSymbolSize( int size ) {
    myAttr.size      = size;
  }
  public void setLineWidth( float lineWidth ) {
    myAttr.lineWidth = lineWidth;
    myStroke = new BasicStroke(myAttr.lineWidth,java.awt.BasicStroke.CAP_BUTT,java.awt.BasicStroke.JOIN_ROUND);
  }
  public void draw( csSeisView seisview, Graphics2D g ) {
    RenderingHints rhints = g.getRenderingHints();
    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
    csISeismicTraceBuffer traceBuffer = seisview.getTraceBuffer();
    if( traceBuffer == null || myHeaderIndex < 0 || myHeaderIndex >= traceBuffer.numHeaders() ) return;
    Rectangle rectVisible = seisview.getVisibleRect();

    int minX = rectVisible.x;
    int maxX = rectVisible.width + minX;
    int minY = rectVisible.y;
    int maxY = rectVisible.height + minY;

    g.setColor( myAttr.color );
    g.setStroke( myStroke );

    int x0 = (int)seisview.xView2Trace(minX);
    int minTrace = Math.max( x0 - (x0 % myAttr.trcInterval), 0 );
    int maxTrace = Math.min((int)seisview.xView2Trace(maxX) + 1,traceBuffer.numTraces()-1);

    if( myAttr.type == TYPE_LINE ) {
      GeneralPath path = new GeneralPath();
      for( int itrc = minTrace; itrc <= maxTrace; itrc += myAttr.trcInterval ) {
        csHeader[] headers = traceBuffer.headerValues(itrc);
        float valueTime = headers[myHeaderIndex].floatValue()*myAttr.multValue + myAttr.addValue;
        float valueSampleIndex = (float)(valueTime / seisview.getSampleInt());
        float yp = seisview.yModel2View( valueSampleIndex );
        float xp = seisview.xModel2View( itrc );
        if( itrc == minTrace ) path.moveTo(xp, yp);
        else path.lineTo(xp, yp);
      }
      g.draw(path);
    }
    else {
      for( int itrc = minTrace; itrc <= maxTrace; itrc += myAttr.trcInterval ) {
        csHeader[] headers = traceBuffer.headerValues(itrc);
        float valueTime = headers[myHeaderIndex].floatValue()*myAttr.multValue + myAttr.addValue;
        float valueSampleIndex = (float)(valueTime / seisview.getSampleInt());
        float yp = seisview.yModel2View( valueSampleIndex );
        if( yp >= minY && yp <= maxY ) {
          float xp = seisview.xModel2View( itrc );
          if( myAttr.type == TYPE_CROSS ) {
            Shape shape = csHeaderOverlay.createDiagonalCross( xp,yp,myAttr.size );
            g.draw(shape);
          }
          else if( myAttr.type == TYPE_CIRCLE ) {
            Shape shape = csHeaderOverlay.createCircleDot( xp,yp,myAttr.size );
            g.draw(shape);
          }
        }
      }
    }
    g.setColor(myAttr.color);

    csHeader[] headers = traceBuffer.headerValues(minTrace);
    float valueTime = headers[myHeaderIndex].floatValue()*myAttr.multValue + myAttr.addValue;
    float yp = seisview.yModel2View( (float)(valueTime / seisview.getSampleInt()) );
    float dist = 10.0f;
    yp = Math.max(yp-dist,minY+dist);
    yp = Math.min(yp,maxY);

    if( myShowLabel ) {
      Font fontSave = g.getFont();
      g.setFont( new Font(fontSave.getFontName(), Font.BOLD, fontSave.getSize()));
      g.drawString(myHeaderName,minX+dist,yp);
      g.setFont(fontSave);
    }
    g.setRenderingHints(rhints);
  }
  
  public static Shape createDiagonalCross( float xp, float yp, final int symbolSize ) {
    final GeneralPath p0 = new GeneralPath();
    float length = (float)symbolSize/2.0f;
    p0.moveTo(xp-length, yp-length);
    p0.lineTo(xp+length, yp+length);
    p0.moveTo(xp-length, yp+length);
    p0.lineTo(xp+length, yp-length);
    return p0;
  }
  public static Shape createCircleDot( float xp, float yp, final int symbolSize ) {
    float length = (float)symbolSize/2.0f;
    int x0 = (int)(xp-length+0.5);
    int y0 = (int)(yp-length+0.5);
    Shape shape = new Ellipse2D.Float(x0,y0,symbolSize,symbolSize);
    final GeneralPath p0 = new GeneralPath(shape);
    Shape shape2 = new Rectangle.Float(xp,yp,1,1);
    p0.append(shape2, false);
    return p0;
  }
  public class Attribute {
    public int size;
    public Color color;
    public float multValue;
    public float addValue;
    public int trcInterval;
    public int transparency;
    public float lineWidth;
    public int type;
    public Attribute() {
      size      = 9;
      color     = Color.yellow;
      multValue = 1.0f;
      addValue  = 0.0f;
      trcInterval = 1;
      transparency = 255;
      lineWidth   = 2.0f;
      type = csHeaderOverlay.TYPE_LINE;
    }
    public void set( Attribute attr ) {
      size      = attr.size;
      color     = attr.color;
      multValue = attr.multValue;
      addValue  = attr.addValue;
      trcInterval = attr.trcInterval;
      transparency = attr.transparency;
      lineWidth = attr.lineWidth;
      type = attr.type;
    }
  }
}

