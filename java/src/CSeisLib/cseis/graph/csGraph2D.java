/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

import javax.swing.*;

import java.awt.event.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;
import java.util.*;

/**
 * Prototype for 2D graph panel.<br>
 * Will currently only display a single 'curve' at once.
 * @author Bjorn Olofsson
 * 2015 - Added support for arbitrary number of graphs
 */
@SuppressWarnings("serial")
public class csGraph2D extends JPanel implements MouseListener, MouseMotionListener, csCurveAttrChangeListener {
  private java.text.DecimalFormat myAnnotationFormat;
  private java.text.DecimalFormat myExponentFormat;
  public static final float[] DASH_PATTERN = { 3.0f, 3.0f };
  public static final float DASH_PHASE = 3.0f;
  
  protected int my2TopXAxis;
  protected int my2TopOuterBorder;
  protected int my2HeightOuterBorder;
  protected int my2TopInnerBorder;
  protected int my2HeightInnerBorder;
  protected int my2TopInnerGraph;
  protected int my2HeightInnerGraph;

  protected int my2LeftYAxis;
  protected int my2LeftOuterBorder;
  protected int my2WidthOuterBorder;
  protected int my2LeftInnerBorder;
  protected int my2WidthInnerBorder;
  protected int my2LeftInnerGraph;
  protected int my2WidthInnerGraph;

  protected Point2D.Float my2InnerTitleCenterPos = null;
  protected Rectangle2D.Float my2InnerViewPortRect;
  protected Rectangle2D.Float my2OuterViewPortRect;
  protected Rectangle2D.Float my2OuterGraphRect;
  protected Rectangle2D.Float my2OuterBorderRect;
  protected Rectangle2D.Float my2InnerBorderRect;

  // Current view position (= scroll position) of graph
  protected int myViewPositionVert = 0;
  protected int myViewPositionHorz = 0;
  protected int myPrevViewPositionVert = 0;
  protected int myPrevViewPositionHorz = 0;
  protected int myNextViewPositionVert = -1;
  protected int myNextViewPositionHorz = -1;
  protected boolean myIsScrolling;

  /// Total width/height of area available for graph painting, borders etc
  /// Part of the total area may be hidden (outside viewport)
  protected int myTotalWidth;
  protected int myTotalHeight;
  /// Width/height of currently allocated bitmap image
  /// This covers only the currently visible area
  protected int myImageWidth;
  protected int myImageHeight;

  public static final int DEFAULT_PADDING_LABEL  = 4;
  public static final int DEFAULT_SIZE_LABEL     = 20;
  public static final int DEFAULT_SIZE_TITLE     = 30;
  public static final int DEFAULT_BORDER_PADDING = 10;
  public static final int DEFAULT_MAJOR_TIC_LENGTH  = 9;
  public static final int DEFAULT_MINOR_TIC_LENGTH  = 7;

  protected float myTotalCurveMinX;
  protected float myTotalCurveMaxX;
  protected float myTotalCurveMinY;
  protected float myTotalCurveMaxY;
  protected float myTotalCurveMinY2;
  protected float myTotalCurveMaxY2;

  // Helper fields for quicker painting
  // Lower left corner of labels/title, centre text position
  // ...needs to be adjusted depending on length of text string
  protected Point2D.Float myTitleCenterPos;
  protected Point2D.Float myXLabelCenterPos;
  protected Point2D.Float myYLabelCenterPos;

  /// Translation of actual data plot in X/Y direction of visible area, in pixels
  protected csGraphAttributes myGraphAttr;

  protected HashMap<Integer,csCurve> myCurves;

  private float myXModel2ViewScalar;
  private float myYModel2ViewScalar;
  private float myY2Model2ViewScalar;

  protected ArrayList<float[]> myDataX;
  protected ArrayList<float[]> myDataY;
  protected csGraphDispDialog myDialogSettings;
  protected csGraphPopupMenu myPopupMenu;
  protected int myGridViewInc;

  protected float myMinGridX;
  protected float myMinGridY;
  protected int myNumGridLinesX;
  protected int myNumGridLinesY;

  protected ArrayList<csCurveAttributes> myAttr;
  protected float myMaxTotal;
  protected BufferedImage myBitmap;
  protected boolean myDoDoubleBuffering = true;

  /// Power of ten for X axis tic increment
  protected float myPowerX;
  /// Power of ten for Y axis tic increment
  protected float myPowerY;
  protected boolean myIsFullRepaint;

  ArrayList<csIGraph2DListener> myListeners;

  public csGraph2D() {
    myGraphAttr = new csGraphAttributes();

    myIsFullRepaint = true;
    myAnnotationFormat = new java.text.DecimalFormat("0");
    myExponentFormat = new java.text.DecimalFormat("0");
    myCurves = new HashMap<Integer,csCurve>(1);
    myPowerX = 1.0f;
    myPowerY = 1.0f;
    myMaxTotal = 0.0f;
    myXModel2ViewScalar  = 1.0f;
    myYModel2ViewScalar  = 1.0f;
    myY2Model2ViewScalar = 1.0f;
    my2OuterBorderRect = null;
    my2OuterGraphRect  = null;
    myTotalCurveMinX = -1;
    myTotalCurveMaxX = 1;
    myTotalCurveMinY = -1;
    myTotalCurveMaxX = 1;
    myTotalCurveMinY2 = -1;
    myTotalCurveMaxY2 = 1;

    myBitmap = null;

    myTotalWidth  = 200;
    myTotalHeight = 200;
    myImageWidth  = 0;
    myImageHeight = 0;

    myGridViewInc = 40;

    myViewPositionVert = 0;
    myViewPositionHorz = 0;
    myPrevViewPositionVert = 0;
    myPrevViewPositionHorz = 0;
    myNextViewPositionVert = -1;
    myNextViewPositionHorz = -1;

//    myDialogSettings = new csGraphDispDialog( this );
    myPopupMenu      = new csGraphPopupMenu( this );
    myListeners = new ArrayList<csIGraph2DListener>();

    addMouseListener( this );
    addMouseMotionListener( this );
  }
  public void setNewVariableDim( int totalWidth, int totalHeight ) {
    if( !myGraphAttr.isFixedXDim ) {
      myTotalWidth = totalWidth;
    }
    if( !myGraphAttr.isFixedYDim ) {
      myTotalHeight = totalHeight;
    }
    setNewSize( new Dimension(myTotalWidth,myTotalHeight) );
    setGraphLayout();
    setCurveLayout();
    myIsFullRepaint = true;
    repaint();
  }
  public void setFixedDim( boolean fixedWidth, boolean fixedHeight, int totalWidth, int totalHeight ) {
    myGraphAttr.isFixedXDim = fixedWidth;
    myGraphAttr.isFixedYDim = fixedHeight;
    if( myGraphAttr.isFixedXDim ) {
      myTotalWidth = totalWidth;
    }
    if( myGraphAttr.isFixedYDim ) {
      myTotalHeight = totalHeight;
    }
    setNewSize( new Dimension(myTotalWidth,myTotalHeight) );
    setGraphLayout();
    setCurveLayout();
    myIsFullRepaint = true;
    repaint();
  }
  public void addGraph2DListener( csIGraph2DListener listener ) {
    myListeners.add(listener);
  }
  public void removeGraph2DListener( csIGraph2DListener listener ) {
    myListeners.remove(listener);
  }
  //-----------------------------------------------------------------
  //
  public void removeCurve( int curveID ) {
    csCurve curve = myCurves.get(curveID);
    if( curve == null ) return;
    curve.attr.removeListener( this );
    myCurves.remove(curveID);
    resetCurveMinMax();
    setCurveLayout();
    repaint();
    fireGraphChangedEvent();
  }
  public void updateCurve( int curveID, float[] valuesX, float[] valuesY, boolean refresh ) {
    csCurve curve = myCurves.get(curveID);
    if( curve == null ) {
      curve = new csCurve(valuesX,valuesY,new csCurveAttributes(),curveID);
    }
    else {
      curve.dataX = valuesX;
      curve.dataY = valuesY;
    }
    addCurve(curve,refresh);
  }
  public void addCurve( int curveID, float[] valuesX, float[] valuesY, csCurveAttributes attr_in ) {
    if( attr_in == null ) attr_in = new csCurveAttributes();
    addCurve( new csCurve( valuesX, valuesY, attr_in, curveID ) );
  }
  public void addCurve( csCurve curve ) {
    addCurve( curve, true );
  }
  public void addCurve( csCurve curve, boolean refresh ) {
    if( curve == null ) return;
    curve.attr.addListener( this );
    myCurves.put( curve.curveID, curve );
    setCurve( curve.curveID );
    if( refresh ) {
      resetCurveMinMax();
      setCurveLayout();
      repaint();
    }
    fireGraphChangedEvent();
  }
  private void setCurve( int curveID ) {
    csCurveAttributes attr = myCurves.get(curveID).attr;
    float[] valuesX = myCurves.get(curveID).dataX;
    float[] valuesY = myCurves.get(curveID).dataY;
    // Compute min/max values
    float minX = Float.MAX_VALUE;
    float maxX = -Float.MAX_VALUE;
    float minY = Float.MAX_VALUE;
    float maxY = -Float.MAX_VALUE;
    int nValues = valuesX.length;
    for( int ivalue = 0; ivalue < nValues; ivalue++ ) {
      if( valuesX[ivalue] < minX ) minX = valuesX[ivalue];
      if( valuesX[ivalue] > maxX ) maxX = valuesX[ivalue];
      if( valuesY[ivalue] < minY ) minY = valuesY[ivalue];
      if( valuesY[ivalue] > maxY ) maxY = valuesY[ivalue];
    }
    attr.setMinMax(minX, maxX, minY, maxY);
  }
  //---------------------------------------------------------------------------
  //
  private void resetCurveMinMax() {
    float minX = Float.MAX_VALUE;
    float maxX = Float.MIN_VALUE;
    float minY = Float.MAX_VALUE;
    float maxY = Float.MIN_VALUE;
    float minY2 = Float.MAX_VALUE;
    float maxY2 = Float.MIN_VALUE;
    for( Iterator<csCurve> iter = myCurves.values().iterator(); iter.hasNext(); ) {
      csCurveAttributes attr = iter.next().attr;
      if( attr.axisXType == csCurveAttributes.X_AXIS_BOTTOM ) {
        if( attr.minX() < minX ) minX = attr.minX();
        if( attr.maxX() > maxX ) maxX = attr.maxX();
      }
      if( attr.axisYType == csCurveAttributes.Y_AXIS_LEFT ) {
        if( attr.minY() < minY ) minY = attr.minY();
        if( attr.maxY() > maxY ) maxY = attr.maxY();
      }
      else {
        if( attr.minY() < minY2 ) minY2 = attr.minY();
        if( attr.maxY() > maxY2 ) maxY2 = attr.maxY();
      }
    }
    myTotalCurveMinX = minX;
    myTotalCurveMaxX = maxX;
    myTotalCurveMinY = minY;
    myTotalCurveMaxY = maxY;
    myTotalCurveMinY2 = minY2;
    myTotalCurveMaxY2 = maxY2;

    if( myTotalCurveMaxX == myTotalCurveMinX ) {
      myTotalCurveMinX -= 1;
      myTotalCurveMaxX += 1;
    }
    if( myTotalCurveMaxY == myTotalCurveMinY ) {
      myTotalCurveMinY -= 1;
      myTotalCurveMaxY += 1;
    }
    if( myTotalCurveMaxY2 == myTotalCurveMinY2 ) {
      myTotalCurveMinY2 -= 1;
      myTotalCurveMaxY2 += 1;
    }
  }
  public void repaintAll() {
    myIsFullRepaint = true;
    repaint();
  }
  //-----------------------------------------------------------------
  public void setNewSize( Dimension size ) {
    setPreferredSize( new Dimension(size) );
    setSize( new Dimension(size) );   // This will lead to a call to repaint()
  }
  public void paintComponent( Graphics g ) {
    boolean snapShot_isScrolling = myIsScrolling;
    Graphics2D g2 = (Graphics2D)g;
    paintGraph(myIsFullRepaint);
    if( myBitmap != null ) {
      g2.drawImage( myBitmap, 0, 0, null );
    }
    myIsFullRepaint = false;

    if( snapShot_isScrolling ) {
      myPrevViewPositionVert = myViewPositionVert;
      myPrevViewPositionHorz = myViewPositionHorz;
      myIsScrolling = false;
      boolean isVertScroll = myNextViewPositionVert >= 0 && myNextViewPositionVert != myViewPositionVert;
      boolean isHorzScroll = myNextViewPositionHorz >= 0 && myNextViewPositionHorz != myViewPositionHorz;
      if( isVertScroll ) {
        resetViewPositionVert( myNextViewPositionVert );
      }
      else if( isHorzScroll ) {
        resetViewPositionHorz( myNextViewPositionHorz );
      }
    }
  }
  //-----------------------------------------------------------------
  public int getNumCurves() {
    return myCurves.size();
  }
  //-----------------------------------------------------------------
  public java.util.Set<Integer> getCurveIDs() {
    return myCurves.keySet();
  }
  //-----------------------------------------------------------------
  public java.util.Collection<csCurve> getCurves() {
    return myCurves.values();
  }
  //-----------------------------------------------------------------
  //
  public void paintGraph( boolean doReset ) {
    Rectangle rect = getVisibleRect();
    if( rect.height == 0 || rect.width == 0 ) return;
    if( myGraphAttr.isFixedXDim && rect.width > myTotalWidth ) rect.width = myTotalWidth;  // Reduce painted area to actual maximum graph width in case enclosing frame/window is wider than graph
    if( doReset || myImageWidth != rect.width || myImageHeight != rect.height || (myBitmap == null && myDoDoubleBuffering)) {
      myImageWidth  = rect.width;
      myImageHeight = rect.height;
      if( myDoDoubleBuffering ) {
        myBitmap = new BufferedImage( myImageWidth, myImageHeight, BufferedImage.TYPE_INT_RGB );
      }
      if( !myGraphAttr.isFixedXDim ) myTotalWidth  = myImageWidth;
      if( !myGraphAttr.isFixedYDim ) myTotalHeight = myImageHeight;
      setGraphLayout();
      if( myCurves.size() > 0 ) setCurveLayout();
    }
//      resetCurveDependencies();
//    }
//    else if( doReset ) {
//      resetGraphLayout();
//    }
    if( myBitmap == null ) return;
    Graphics2D g2 = (Graphics2D)myBitmap.getGraphics();
    g2.setRenderingHint( RenderingHints.KEY_ANTIALIASING,  RenderingHints.VALUE_ANTIALIAS_ON );

    // If graph is inside view port, paint border & labels relative to current visible area
    // Otherwise, border etc is painted relative to left-most edge of total area (translation before painting border)
    if( myGraphAttr.isGraphViewPort ) {
      paintLabelsBorder( g2 );
      if( myGraphAttr.backgroundColor != null ) {
        g2.setColor( myGraphAttr.backgroundColor );
        g2.fill( rect );
      }
      else {
        g2.clearRect( rect.x, rect.y, rect.width, rect.height );
      }
      paintLabelsBorder( g2 );
      g2.setClip(my2OuterViewPortRect);
      g2.translate(-myViewPositionHorz, -myViewPositionVert);
      paintCurves( g2 );
      if( myCurves.size() > 0 ) {
        g2.setClip((int)my2OuterViewPortRect.x+myViewPositionHorz, (int)my2OuterViewPortRect.y+myViewPositionVert,
                (int)my2OuterViewPortRect.width, (int)my2OuterViewPortRect.height);
        paintGrid( g2 );
        g2.setClip(null);
        float xposMinAxis = my2OuterGraphRect.x + myViewPositionHorz;
        float xposMaxAxis = my2OuterGraphRect.x + my2OuterGraphRect.width + myViewPositionHorz;
        paintXAxis( g2, xposMinAxis, xposMaxAxis );
        g2.translate(myViewPositionHorz, myViewPositionVert);
        g2.setClip(null);
        xposMinAxis = my2OuterGraphRect.x;
        xposMaxAxis = my2OuterGraphRect.x + my2OuterGraphRect.width;
        paintYAxis( g2, xposMinAxis, xposMaxAxis );
      }
    }
    else {
      g2.translate(-myViewPositionHorz, -myViewPositionVert);
      if( myGraphAttr.backgroundColor != null ) {
        g2.setColor( myGraphAttr.backgroundColor );
        g2.fillRect( 0, 0, myTotalWidth, myTotalHeight );
      }
      else {
        g2.clearRect( 0, 0, myTotalWidth, myTotalHeight );
      }
      if( myCurves.size() > 0 ) {
        paintCurves( g2 );
        float xposMinAxis = my2OuterGraphRect.x;
        float xposMaxAxis = my2OuterGraphRect.x + my2OuterGraphRect.width;
        paintGrid( g2 );
        paintXAxis( g2, xposMinAxis, xposMaxAxis );
        paintYAxis( g2, xposMinAxis, xposMaxAxis );
      }
      paintLabelsBorder( g2 );
      g2.translate(myViewPositionHorz, myViewPositionVert);
    }
  }

  //-----------------------------------------------
  public float xModel2View( float xModel ) {
    return( (xModel-myTotalCurveMinX)*myXModel2ViewScalar + my2LeftInnerGraph );
  }
  //-----------------------------------------------
  public float yModel2View( float yModel ) {
    return( (myTotalCurveMaxY-yModel)*myYModel2ViewScalar + my2TopInnerGraph );
/*    if( myGraphAttr.axisScaleY == csGraphAttributes.AXIS_SCALE_LINEAR ) {
    }
    else {
      float diffFull = myTotalCurveMaxY-myTotalCurveMinY;
      float diffFullLog = Math.log10( myTotalCurveMaxY-myTotalCurveMinY );
      float diff = myTotalCurveMaxY-yModel;
      float diffLog = (float)Math.log10( Math.abs(diff) );
      this.myTotalCurveMinY
      return( ()*myYModel2ViewScalar + my2TopInnerGraph );
    }
*/
  }
  //-----------------------------------------------
  public float xView2Model( float xView ) {
    if( myXModel2ViewScalar == 0 ) return 0;
    return( (xView - my2LeftInnerGraph) / myXModel2ViewScalar + myTotalCurveMinX );
  }
  //-----------------------------------------------
  public float yView2Model( float yView ) {
    if( myYModel2ViewScalar == 0 ) return 0;
    return( ( my2TopInnerGraph - yView ) / myYModel2ViewScalar + myTotalCurveMaxY );
  }
  //-----------------------------------------------
  public float y2Model2View( float yModel ) {
    return( (myTotalCurveMaxY2-yModel)*myY2Model2ViewScalar + my2TopInnerGraph );
  }
  //-----------------------------------------------
  public float y2View2Model( float yView ) {
    if( myY2Model2ViewScalar == 0 ) return 0;
    return( ( my2TopInnerGraph - yView ) / myY2Model2ViewScalar + myTotalCurveMaxY2 );
  }
  private void setGraphLayout() {
    // Second Y axis/label not supported yet!

//-------------------------------------------------------------------------------------
    //...from top to bottom...
    // insetTop
    // Title
    // Border (outer)
    // Border padding
    // Inner title
    // Inner border
    // Inner border padding
    // Y axis 10^exponent
    // Graph padding (top)
    // Top of inner graph area (curves)
    //
    // Inner graph area (for drawing of curves)
    //
    // Graph padding
    // Axis & tics
    // X annotation
    // Inner border padding
    // Inner border
    // X label
    // Border padding
    // Border
    // insetBottom
    my2TopOuterBorder = myGraphAttr.insetTop;
    my2HeightOuterBorder = myTotalHeight - myGraphAttr.insetTop - myGraphAttr.insetBottom;
    if( myGraphAttr.title.length() > 0 ) {
      my2TopOuterBorder += DEFAULT_SIZE_TITLE;
      my2HeightOuterBorder -= DEFAULT_SIZE_TITLE;
    }
    my2TopInnerBorder = my2TopOuterBorder;
    my2HeightInnerBorder = my2HeightOuterBorder;
    if( myGraphAttr.showBorder ) {
      my2TopInnerBorder += myGraphAttr.borderPadding;
      my2HeightInnerBorder -= 2*myGraphAttr.borderPadding;
    }
    if( myGraphAttr.innerTitle.length() > 0 ) {
      my2TopInnerBorder += DEFAULT_SIZE_LABEL;
      my2HeightInnerBorder -= DEFAULT_SIZE_LABEL;
    }
    if( myGraphAttr.xLabel.length() > 0 ) {
      my2HeightInnerBorder -= DEFAULT_SIZE_LABEL;
    }
    int topGraph = my2TopInnerBorder;
    int heightGraph = my2HeightInnerBorder;
    if( myGraphAttr.showInnerBorder ) {
      topGraph += myGraphAttr.innerBorderPadding;
      heightGraph -= 2*myGraphAttr.innerBorderPadding;
    }
    if( myGraphAttr.showAxisY ) {
//      if( myGraphAttr.showYTics ) {
        if( myGraphAttr.showAxisYAnnotation ) {
          topGraph    += DEFAULT_SIZE_LABEL + DEFAULT_PADDING_LABEL; // additional space for 10^exponent on Y axis
          heightGraph -= (DEFAULT_SIZE_LABEL + DEFAULT_PADDING_LABEL); // additional space for 10^exponent on Y axis
        }
//      }
    }
    if( myGraphAttr.showAxisX ) {
//      if( myGraphAttr.showXTics ) {
        if( myGraphAttr.showAxisXAnnotation ) {
          heightGraph -= (DEFAULT_SIZE_LABEL + DEFAULT_PADDING_LABEL);
        }
//      }
    }
    my2TopXAxis = topGraph + heightGraph;
    topGraph    += myGraphAttr.graphPadding;
    heightGraph -= 2*myGraphAttr.graphPadding;
    my2TopInnerGraph    = topGraph;
    my2HeightInnerGraph = heightGraph;

    //...from left to right...
    // insetLeft
    // Border (outer)
    // Border padding
    // Y label
    // Inner border
    // Inner border padding
    // Y annotation
    // Y Axis & tics
    // Graph padding (left)
    // Left of inner graph area (curves)
    //
    // Inner graph area (for drawing of curves)
    //
    // Graph padding
    // X axis 10^exponent
    // Inner border padding
    // Inner border
    // Border padding
    // Border
    // insetBottom
    my2LeftOuterBorder = myGraphAttr.insetLeft;
    my2WidthOuterBorder = myTotalWidth - myGraphAttr.insetLeft - myGraphAttr.insetRight;

    my2LeftInnerBorder = my2LeftOuterBorder;
    my2WidthInnerBorder = my2WidthOuterBorder;
    if( myGraphAttr.showBorder ) {
      my2LeftInnerBorder += myGraphAttr.borderPadding;
      my2WidthInnerBorder -= 2*myGraphAttr.borderPadding;
    }
    if( myGraphAttr.yLabel.length() > 0 ) {
      my2LeftInnerBorder  += DEFAULT_SIZE_LABEL;
      my2WidthInnerBorder -= DEFAULT_SIZE_LABEL;
    }
    int leftGraph = my2LeftInnerBorder;
    int widthGraph = my2WidthInnerBorder;
    if( myGraphAttr.showInnerBorder ) {
      leftGraph += myGraphAttr.innerBorderPadding;
      widthGraph -= 2*myGraphAttr.innerBorderPadding;
    }
    if( myGraphAttr.showAxisX ) {
      if( myGraphAttr.showAxisXAnnotation ) {
        widthGraph -= (myGraphAttr.annotationSizeYAxis + DEFAULT_PADDING_LABEL); // additional space for 10^exponent on X axis
      }
    }
    if( myGraphAttr.showAxisY ) {
      if( myGraphAttr.showAxisYAnnotation ) {
        leftGraph  += (myGraphAttr.annotationSizeYAxis + DEFAULT_PADDING_LABEL);
        widthGraph -= (myGraphAttr.annotationSizeYAxis + DEFAULT_PADDING_LABEL);
      }
    }
    my2LeftYAxis = leftGraph;
    int leftOuterGraph = leftGraph;
    int widthOuterGraph = widthGraph;
    leftGraph    += myGraphAttr.graphPadding;
    widthGraph -= 2*myGraphAttr.graphPadding;
    my2LeftInnerGraph  = leftGraph;
    my2WidthInnerGraph = widthGraph;

    if( myGraphAttr.fixHorzGraphPosition ) {
      int diffLeft  = myGraphAttr.horzGraphPositionMin - my2LeftInnerGraph;
      int diffWidth = (myGraphAttr.horzGraphPositionMax - myGraphAttr.horzGraphPositionMin) - my2WidthInnerGraph;
      my2LeftInnerGraph  += diffLeft;
      my2WidthInnerGraph += diffWidth;
    }

    my2OuterGraphRect = new Rectangle2D.Float(
              leftOuterGraph, my2TopInnerGraph-myGraphAttr.graphPadding,
              widthOuterGraph, my2HeightInnerGraph + 2*myGraphAttr.graphPadding );
    my2InnerViewPortRect = new Rectangle2D.Float( my2LeftInnerGraph, my2TopInnerGraph,
            my2WidthInnerGraph - myTotalWidth + myImageWidth, my2HeightInnerGraph - myTotalHeight + myImageHeight );

    my2OuterViewPortRect = new Rectangle2D.Float( my2OuterGraphRect.x, my2OuterGraphRect.y,
            my2OuterGraphRect.width - myTotalWidth + myImageWidth, my2OuterGraphRect.height - myTotalHeight + myImageHeight );

    if( myGraphAttr.isGraphViewPort ) {
      my2WidthOuterBorder += myImageWidth - myTotalWidth;
      my2WidthInnerBorder += myImageWidth - myTotalWidth;
      my2HeightOuterBorder += myImageHeight - myTotalHeight;
      my2HeightInnerBorder += myImageHeight - myTotalHeight;
    }

    my2OuterBorderRect = new Rectangle2D.Float( my2LeftOuterBorder, my2TopOuterBorder, my2WidthOuterBorder, my2HeightOuterBorder );
    my2InnerBorderRect = new Rectangle2D.Float( my2LeftInnerBorder, my2TopInnerBorder, my2WidthInnerBorder, my2HeightInnerBorder );

    // Helper fields
    if( myGraphAttr.xLabel.length() > 0 ) {
      myXLabelCenterPos = new Point2D.Float(
              my2OuterGraphRect.x + my2OuterGraphRect.width/2,
              my2TopInnerBorder + my2HeightInnerBorder + DEFAULT_SIZE_LABEL-DEFAULT_PADDING_LABEL );
    }
    if( myGraphAttr.yLabel.length() > 0 ) {
      myYLabelCenterPos = new Point2D.Float(
              my2LeftInnerBorder - DEFAULT_PADDING_LABEL,
              my2OuterGraphRect.y + my2OuterGraphRect.height/2);
    }
    if( myGraphAttr.title.length() > 0 ) {
      myTitleCenterPos  = new Point2D.Float(
              my2LeftOuterBorder + my2WidthOuterBorder/2,
              my2TopOuterBorder - DEFAULT_SIZE_TITLE/2 );
    }
    if( myGraphAttr.innerTitle.length() > 0 ) {
      my2InnerTitleCenterPos  = new Point2D.Float(
              my2LeftInnerBorder + my2WidthInnerBorder/2,
              my2TopInnerBorder - DEFAULT_SIZE_LABEL/2 );
    }

    //-------------------------------------------------------------------------------
  }
  public boolean isFixedXDim() {
    return myGraphAttr.isFixedXDim;
  }
  public boolean isFixedYDim() {
    return myGraphAttr.isFixedYDim;
  }
  public void setCurveLayout() {
    if( !myGraphAttr.isFixedXDim || !myGraphAttr.isFixedYDim ) {
      myXModel2ViewScalar = my2WidthInnerGraph / (myTotalCurveMaxX-myTotalCurveMinX);
      myYModel2ViewScalar = my2HeightInnerGraph / (myTotalCurveMaxY-myTotalCurveMinY);

      float incx_model = Math.abs(xView2Model(myGraphAttr.gridXPixelInc)-xView2Model(0.0f));
      float incy_model = Math.abs(yView2Model(myGraphAttr.gridYPixelInc)-yView2Model(0.0f));

      int power1;
      double log1 = Math.log10(incx_model);
      if( log1 > 0 ) power1 = (int)log1;
      else power1 = -(int)(-log1) - 1;
      int scalar = (int)Math.round( incx_model / (float)Math.pow(10, power1) );
      myPowerX = (float)Math.pow(10, power1);
      float incx_model_red = myPowerX;
      if( scalar == 0 ) {
        incx_model_red = (float)Math.pow(10, power1-1);
      }
      else if( scalar < 3 ) incx_model_red *= scalar;
      else if( scalar < 7 ) incx_model_red *= 5.0;
      else { incx_model_red *= 10.0; myPowerX *= 10.0; }

      myGraphAttr.gridIncX = Math.abs(xModel2View(incx_model_red) - xModel2View(0.0f));

      log1 = Math.log10(incy_model);
      if( log1 > 0 ) power1 = (int)log1;
      else power1 = -(int)(-log1+0.001) - 1;
      scalar = (int)Math.round( incy_model / (float)Math.pow(10, power1) );
      myPowerY = (float)Math.pow(10, power1);
      float incy_model_red = myPowerY;
      if( scalar == 0 ) {
        incy_model_red = (float)Math.pow(10, power1-1);
      }
      else if( scalar < 3 ) incy_model_red *= scalar;
      else if( scalar < 7 ) incy_model_red *= 5.0;
      else { incy_model_red *= 10.0; myPowerY *= 10.0; }

      myGraphAttr.gridIncY = Math.abs(yModel2View(incy_model_red) - yModel2View(0.0f));
    }
  }
  public void paintCurves( Graphics2D g2 ) {
    for( Iterator<csCurve> iter = myCurves.values().iterator(); iter.hasNext(); ) {
      csCurve curve = iter.next();
      float[] dataX = curve.dataX;
      float[] dataY = curve.dataY;
      csCurveAttributes attr = curve.attr;
      int nValues = dataX.length;

      float xpos[] = new float[nValues];
      float ypos[] = new float[nValues];
      //------------------------------------------------------------------------
/*      if( myGraphAttr.axisScaleY == csGraphAttributes.AXIS_SCALE_LOG ) {
        float[] dataYOrig = dataY;
        double maxY = attr.maxY();
        if( myGraphAttr.useRefValueY ) {
          maxY = myGraphAttr.maxRefValueY;
        }
        dataY = new float[nValues];
        for( int ivalue = 0; ivalue < nValues; ivalue++ ) {
          dataY[ivalue] = (float)( 10 * Math.log( dataYOrig[ivalue]/maxY ) );
        }
      }
 */
      //------------------------------------------------------------------------
      int counter = 0;
      for( int ivalue = 0; ivalue < nValues; ivalue++ ) {
        if( dataX[ivalue] >= myTotalCurveMinX && dataX[ivalue] <= myTotalCurveMaxX &&
           dataY[ivalue] >= myTotalCurveMinY && dataY[ivalue] <= myTotalCurveMaxY ) {
          xpos[counter] = xModel2View(dataX[ivalue]);
          ypos[counter] = yModel2View(dataY[ivalue]);
          counter += 1;
        }
      }
      int nPositions = counter;
      if( nPositions == 0 ) return;

      //---------------------------------------------------------------------
      // Fill curves
      //
      if( attr.filledType != csCurveAttributes.FILLED_TYPE_NONE ) {
        GeneralPath pathFill = new GeneralPath();
        pathFill.moveTo( xpos[0], ypos[0] );
        for( int ipos = 1; ipos < nPositions; ipos++ ) {
          pathFill.lineTo( xpos[ipos], ypos[ipos] );
        }
        if( attr.filledType == csCurveAttributes.FILLED_TYPE_ZERO ) {
          pathFill.lineTo( xpos[nPositions-1], yModel2View(0) );
          pathFill.lineTo( xpos[0], yModel2View(0) );
        }
        else if( attr.filledType == csCurveAttributes.FILLED_TYPE_MIN ) {
          pathFill.lineTo( xpos[nPositions-1], yModel2View(myTotalCurveMinY) );
          pathFill.lineTo( xpos[0], yModel2View(myTotalCurveMinY) );
        }
        g2.setColor(attr.fillColor);
        g2.fill( pathFill );
      }
      //---------------------------------------------------------------------
      // Paint lines
      //
      if( attr.lineType != csCurveAttributes.LINE_TYPE_NONE ) {
        GeneralPath path = new GeneralPath();
        path.moveTo( xpos[0], ypos[0] );
        for( int ipos = 1; ipos < nPositions; ipos++ ) {
          path.lineTo( xpos[ipos], ypos[ipos] );
        }
        Stroke strokeSave = g2.getStroke();
        if( attr.lineType == csCurveAttributes.LINE_TYPE_SOLID ) {
          g2.setStroke( new BasicStroke(attr.lineSize) );
        }
        else {
          float[] dashPattern = new float[DASH_PATTERN.length];
          if( attr.lineType == csCurveAttributes.LINE_TYPE_DASH ) {
            dashPattern[0] = 3*attr.lineSize;
            dashPattern[1] = 4*attr.lineSize;
            g2.setStroke( new BasicStroke( 
                attr.lineSize, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_ROUND, 10.0f, dashPattern, 0.0f ) );
          }
          else if( attr.lineType == csCurveAttributes.LINE_TYPE_DOT ) {
            dashPattern[0] = 1;
            dashPattern[1] = 1 + (int)(1.3*attr.lineSize);
            g2.setStroke( new BasicStroke( 
                attr.lineSize, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10.0f, dashPattern, 0.0f ) );
          }
        }
        g2.setColor(attr.lineColor);
        g2.draw( path );
        g2.setStroke( strokeSave );
      }
      //---------------------------------------------------------------------
      // Paint points
      //
      if( attr.pointType != csCurveAttributes.POINT_TYPE_NONE ) {
        g2.setStroke( new BasicStroke(attr.lineSize) );
        g2.setColor(attr.pointColor);
        if( attr.pointType == csCurveAttributes.POINT_TYPE_CIRCLE || attr.pointSize == 1 ) {
          for( int ipos = 0; ipos < nPositions; ipos++ ) {
            Shape shape = new Ellipse2D.Float( xpos[ipos]-attr.pointSize/2, ypos[ipos]-attr.pointSize/2, attr.pointSize, attr.pointSize );
            g2.fill( shape );
          }
        }
        else if( attr.pointType == csCurveAttributes.POINT_TYPE_SQUARE ) {
          Rectangle2D.Double shape = new Rectangle2D.Double(0,0,attr.pointSize,attr.pointSize);
          int pointSizeHalf = attr.pointSize/2;
          for( int ipos = 0; ipos < nPositions; ipos++ ) {
            shape.x = xpos[ipos] - pointSizeHalf;
            shape.y = ypos[ipos] - pointSizeHalf;
            g2.fill( shape );
          }
        }
        else if( attr.pointType == csCurveAttributes.POINT_TYPE_CROSS ) {
          int pointSizeHalf = attr.pointSize/2;
          for( int ipos = 0; ipos < nPositions; ipos++ ) {
            GeneralPath path = new GeneralPath();
            double x0 = xpos[ipos];
            double y0 = ypos[ipos];
            path.moveTo(x0-pointSizeHalf,y0-pointSizeHalf);
            path.lineTo(x0+pointSizeHalf,y0+pointSizeHalf);
            path.moveTo(x0-pointSizeHalf,y0+pointSizeHalf);
            path.lineTo(x0+pointSizeHalf,y0-pointSizeHalf);
            g2.draw(path);
          }
        }
      }
    }
  }
//---------------------------------------------------------------------------
  public void paintGrid( Graphics2D g2 ) {
    float xpos0Axis = xModel2View( 0 );
    float ypos0Axis = yModel2View( 0 );

    float xposMinAxis = my2OuterGraphRect.x;
    float xposMaxAxis = my2OuterGraphRect.x + my2OuterGraphRect.width;
    float yposMaxAxis = my2OuterGraphRect.y;
    float yposMinAxis = my2OuterGraphRect.y + my2OuterGraphRect.height;

    float xpos0 = (float)(int)((Math.abs(xposMinAxis-xpos0Axis)/myGraphAttr.gridIncX)+2)*myGraphAttr.gridIncX *
            Math.signum(xposMinAxis-xpos0Axis) + xpos0Axis;
    float ypos0 = (float)(int)((Math.abs(yposMaxAxis-ypos0Axis)/myGraphAttr.gridIncY)+2)*myGraphAttr.gridIncY *
            Math.signum(yposMaxAxis-ypos0Axis) + ypos0Axis;

    int numGridX = (int)Math.round( (xposMaxAxis-xposMinAxis)/myGraphAttr.gridIncX ) + 4;
    int numGridY = (int)Math.round( Math.abs(yposMaxAxis-yposMinAxis)/myGraphAttr.gridIncY ) + 4;

    if( myGraphAttr.showGrid ) {
      GeneralPath path = new GeneralPath();
      for( int ix = 0; ix < numGridX; ix++ ) {
        float xpos = xpos0 + ix*myGraphAttr.gridIncX;
        if( xpos < my2OuterGraphRect.x || xpos > my2OuterGraphRect.x+my2OuterGraphRect.width ) continue;
        path.moveTo( xpos, yposMinAxis );
        path.lineTo( xpos, yposMaxAxis );
      }
      for( int iy = 0; iy < numGridY; iy++ ) {
        float ypos = ypos0 + iy*myGraphAttr.gridIncY;
        if( ypos < my2OuterGraphRect.y || ypos > my2OuterGraphRect.y+my2OuterGraphRect.height ) continue;
        path.moveTo( xposMinAxis, ypos );
        path.lineTo( xposMaxAxis, ypos );
      }
      g2.setColor( myGraphAttr.gridColor );
      g2.setStroke( new BasicStroke(1) );
      g2.draw( path );
    }
  }
//---------------------------------------------------------------------------
  public void paintXAxis( Graphics2D g2, float xposMinAxis, float xposMaxAxis ) {
    float xpos0Axis = xModel2View( 0 );
    float ypos0Axis = yModel2View( 0 );

    float yposMaxAxis = my2OuterGraphRect.y;
    float yposMinAxis = my2OuterGraphRect.y + my2OuterGraphRect.height;

    float xpos0;
    if( xposMinAxis > xpos0Axis ) {
      xpos0 = (float)(int)((xposMinAxis-xpos0Axis)/myGraphAttr.gridIncX - 1)*myGraphAttr.gridIncX + xpos0Axis;
    }
    else {
      xpos0 = (float)(int)((xpos0Axis - xposMinAxis)/myGraphAttr.gridIncX + 1)*-myGraphAttr.gridIncX + xpos0Axis;
    }

    int numGridX = (int)Math.round( (xposMaxAxis-xposMinAxis)/myGraphAttr.gridIncX ) + 4;

    // Paint tics
    int majorTicLength = DEFAULT_MAJOR_TIC_LENGTH/2;
    int minorTicLength = DEFAULT_MINOR_TIC_LENGTH/2;
    GeneralPath pathMajor = new GeneralPath();
    GeneralPath pathMinor = new GeneralPath();
    FontMetrics metrics = g2.getFontMetrics( myGraphAttr.annotationFont );
    g2.setFont(myGraphAttr.annotationFont);

    // Paint axes
    if( myGraphAttr.showAxisX ) {
      GeneralPath path = new GeneralPath();
      if( myGraphAttr.showZeroAxis ) {
        if( ypos0Axis >= yposMaxAxis && ypos0Axis <= yposMinAxis ) {
          path.moveTo( xposMinAxis, ypos0Axis );
          path.lineTo( xposMaxAxis, ypos0Axis );
        }
      }
      path.moveTo( xposMinAxis, yposMinAxis );
      path.lineTo( xposMaxAxis, yposMinAxis );
      g2.setColor( myGraphAttr.axisColor );
      g2.setStroke( new BasicStroke(2) );
      g2.draw( path );
      int numVal = ( myGraphAttr.minorXTicRatio == 0 ) ? numGridX : numGridX*myGraphAttr.minorXTicRatio;

      if( myGraphAttr.showXTics ) {
        for( int ix = 0; ix < numVal; ix++ ) {
          float xpos = xpos0 + (float)ix*myGraphAttr.gridIncX/(float)myGraphAttr.minorXTicRatio;
          if( xpos < xposMinAxis || xpos > xposMaxAxis ) continue;
          if( ix % myGraphAttr.minorXTicRatio == 0 ) {
            pathMajor.moveTo( xpos, yposMinAxis-majorTicLength );
            pathMajor.lineTo( xpos, yposMinAxis+majorTicLength );
          }
          else {
            pathMinor.moveTo( xpos, yposMinAxis-minorTicLength );
            pathMinor.lineTo( xpos, yposMinAxis+minorTicLength );
          }
        }
        g2.setColor( Color.black );
        g2.setStroke( new BasicStroke(2) );
        g2.draw( pathMajor );
        g2.setStroke( new BasicStroke(1) );
        g2.draw( pathMinor );
        g2.setClip(null);
      }
      if( myGraphAttr.showAxisXAnnotation ) {
        for( int ix = 0; ix < numVal; ix++ ) {
          float xpos = xpos0 + (float)ix*myGraphAttr.gridIncX/(float)myGraphAttr.minorXTicRatio;
          if( xpos < xposMinAxis || xpos > xposMaxAxis ) continue;
          if( ix % myGraphAttr.minorXTicRatio == 0 ) {
            String text = myAnnotationFormat.format(xView2Model(xpos)/myPowerX);
            int width = metrics.stringWidth( text );
            g2.drawString( text, (int)(xpos-width/2), (int)yposMinAxis+majorTicLength+metrics.getHeight() );
          }
        }
        String exponent = myExponentFormat.format(Math.log10(myPowerX));
        java.text.AttributedString as = new java.text.AttributedString("10"+exponent);
        as.addAttribute(java.awt.font.TextAttribute.FONT, myGraphAttr.annotationFont, 0, 2 );
        as.addAttribute(java.awt.font.TextAttribute.SUPERSCRIPT, java.awt.font.TextAttribute.SUPERSCRIPT_SUPER, 2, 2+exponent.length());
        g2.drawString( as.getIterator(), xposMaxAxis+DEFAULT_PADDING_LABEL, yposMinAxis+metrics.getHeight()/2 );
      }    
    }
  }
//---------------------------------------------------------------------------
  public void paintYAxis( Graphics2D g2, float xposMinAxis, float xposMaxAxis ) {
    float xpos0Axis = xModel2View( 0 );
    float ypos0Axis = yModel2View( 0 );

    float yposMaxAxis = my2OuterGraphRect.y;
    float yposMinAxis = my2OuterGraphRect.y + my2OuterGraphRect.height;

    float ypos0 = (float)(int)((Math.abs(yposMaxAxis-ypos0Axis)/myGraphAttr.gridIncY)+2)*myGraphAttr.gridIncY *
            Math.signum(yposMaxAxis-ypos0Axis) + ypos0Axis;

    int numGridY = (int)Math.round( Math.abs(yposMaxAxis-yposMinAxis)/myGraphAttr.gridIncY ) + 4;

    // Paint tics
    int majorTicLength = DEFAULT_MAJOR_TIC_LENGTH/2;
    int minorTicLength = DEFAULT_MINOR_TIC_LENGTH/2;
    GeneralPath pathMajor = new GeneralPath();
    GeneralPath pathMinor = new GeneralPath();
    FontMetrics metrics = g2.getFontMetrics( myGraphAttr.annotationFont );
    g2.setFont(myGraphAttr.annotationFont);
    g2.setColor( Color.black );

    if( myGraphAttr.showAxisY ) {
      GeneralPath path = new GeneralPath();
      path.moveTo( xposMinAxis, yposMinAxis );
      path.lineTo( xposMinAxis, yposMaxAxis );
      if( myGraphAttr.showZeroAxis ) {
        if( xpos0Axis >= xposMinAxis && xpos0Axis <= xposMaxAxis ) {
          path.moveTo( xpos0Axis, yposMinAxis );
          path.lineTo( xpos0Axis, yposMaxAxis );
        }
      }
      g2.setStroke( new BasicStroke(2) );
      g2.draw( path );
      
      int numVal = ( myGraphAttr.minorYTicRatio == 0 ) ? numGridY : numGridY*myGraphAttr.minorYTicRatio;
      if( myGraphAttr.showYTics ) {
        for( int iy = 0; iy < numVal; iy++ ) {
          float ypos = ypos0 + (float)iy*myGraphAttr.gridIncY/(float)myGraphAttr.minorYTicRatio;
          if( ypos < my2OuterGraphRect.y || ypos > my2OuterGraphRect.y+my2OuterGraphRect.height ) continue;
          if( iy % myGraphAttr.minorYTicRatio == 0 ) {
            pathMajor.moveTo(xposMinAxis-majorTicLength, ypos);
            pathMajor.lineTo(xposMinAxis+majorTicLength, ypos );
          }
          else {
            pathMinor.moveTo(xposMinAxis-minorTicLength, ypos);
            pathMinor.lineTo(xposMinAxis+minorTicLength, ypos );
          }
        } // END for iy
        
        g2.setColor( Color.black );
        g2.setStroke( new BasicStroke(2) );
        g2.draw( pathMajor );
        g2.setStroke( new BasicStroke(1) );
        g2.draw( pathMinor );
      } // END if showYTics
      if( myGraphAttr.showAxisYAnnotation ) {
        for( int iy = 0; iy < numVal; iy++ ) {
          float ypos = ypos0 + (float)iy*myGraphAttr.gridIncY/(float)myGraphAttr.minorYTicRatio;
          if( ypos < my2OuterGraphRect.y || ypos > my2OuterGraphRect.y+my2OuterGraphRect.height ) continue;
          if( iy % myGraphAttr.minorYTicRatio == 0 ) {
            String text = myAnnotationFormat.format(yView2Model(ypos)/myPowerY);
            int width = metrics.stringWidth( text );
            g2.drawString( text, (int)(xposMinAxis-majorTicLength*2-width), (int)ypos+metrics.getHeight()/2 );
          }
        }
        String text = "10";
        String exponent = myExponentFormat.format(Math.log10(myPowerY));
        java.text.AttributedString as = new java.text.AttributedString("10"+exponent);
        int width = metrics.stringWidth( text );
        as.addAttribute(java.awt.font.TextAttribute.FONT, myGraphAttr.annotationFont, 0, 2 );
        as.addAttribute(java.awt.font.TextAttribute.SUPERSCRIPT, java.awt.font.TextAttribute.SUPERSCRIPT_SUPER, 2, 2+exponent.length());
        g2.drawString( as.getIterator(), xposMinAxis-width/2, yposMaxAxis-DEFAULT_PADDING_LABEL );
      } // END if showAxisYAnnotation
    } // END: if showAxisY
  }
  //------------------------------------------------------------------------------------
  //
  public void paintLabelsBorder( Graphics2D g2 ) {
    // Paint Border
    if( myGraphAttr.showBorder ) {
      g2.setStroke( new BasicStroke(myGraphAttr.borderSize) );
      g2.setColor(myGraphAttr.borderColor);
      g2.draw( my2OuterBorderRect );
    }
    if( myGraphAttr.showInnerBorder ) {
      g2.setStroke( new BasicStroke(myGraphAttr.innerBorderSize) );
      g2.setColor(myGraphAttr.innerBorderColor);
      g2.draw( my2InnerBorderRect );
    }

    FontMetrics metrics = g2.getFontMetrics( g2.getFont() );

    g2.setColor( myGraphAttr.labelColor );
    g2.setFont( myGraphAttr.labelFont );
    if( myGraphAttr.xLabel.length() > 0 ) {
      g2.drawString( myGraphAttr.xLabel, myXLabelCenterPos.x - metrics.stringWidth( myGraphAttr.xLabel )/2, myXLabelCenterPos.y );
    }
    if( myGraphAttr.yLabel.length() > 0 ) {
      double x = myYLabelCenterPos.x;
      double y = myYLabelCenterPos.y + metrics.stringWidth( myGraphAttr.yLabel )/2;
      g2.rotate( -Math.PI/2, x, y );
      g2.drawString( myGraphAttr.yLabel, (int)x, (int)y );
      g2.rotate( Math.PI/2, x, y );
    }
    if( myGraphAttr.title.length() > 0 ) {
      g2.setFont( myGraphAttr.titleFont );
      g2.drawString( myGraphAttr.title, myTitleCenterPos.x - metrics.stringWidth( myGraphAttr.title )/2, myTitleCenterPos.y );
    }
    if( myGraphAttr.innerTitle.length() > 0 && my2InnerTitleCenterPos != null ) {
      g2.setFont( myGraphAttr.labelFont );
      g2.drawString( myGraphAttr.innerTitle, my2InnerTitleCenterPos.x - metrics.stringWidth( myGraphAttr.innerTitle )/2, my2InnerTitleCenterPos.y );
      g2.drawString( myGraphAttr.innerTitle, my2InnerTitleCenterPos.x - metrics.stringWidth( myGraphAttr.innerTitle )/2, my2InnerTitleCenterPos.y );
    }
    if( myGraphAttr.showInnerBorder ) {
      g2.setStroke( new BasicStroke(myGraphAttr.innerBorderSize) );
      g2.setColor(myGraphAttr.innerBorderColor);
        g2.draw( my2InnerBorderRect );
    }
  }

  public void setMaximumAxisValue( float maxValue ) {
    myMaxTotal = maxValue;
  }
  public csGraphAttributes getGraphAttributes()  {
    return myGraphAttr;
  }
  public csCurveAttributes getCurveAttributes( int curveID )  {
    if( myCurves.size() <= 0 ) {
      return null;
    }
    csCurve curve = myCurves.get( curveID );
    if( curve == null ) return null;
    return curve.attr;
  }
  public void setCurveAttributes( int curveID, csCurveAttributes attr )  {
    if( curveID >= myCurves.size() ) return;
    myCurves.get( curveID ).attr = attr;
  }
  public void setGraphAttributes( csGraphAttributes attr )  {
    myGraphAttr = attr;
    repaintAll();
  }
  public void showDialog() {
    if( myDialogSettings == null ) {
      myDialogSettings = new csGraphDispDialog( this );
    }
    else {
      myDialogSettings.changedSettings(myGraphAttr);
    }
    if( !myDialogSettings.isVisible() ) {
      myDialogSettings.setVisible(true);
    }
    else {
      myDialogSettings.toFront();
    }
  }


  //******************************************************************************
  //
  // Repainting and bitmap operations
  //
  //******************************************************************************
/*
  private void resetCurveDependencies() {
    if( myGraphAttr.autoScaleAxes ) {
      if( myGraphAttr.centreAxisX ) {
        float ratio = (float)myGraphHeight / (float)myGraphWidth;
        if( myGraphHeight > myGraphWidth ) {
//          myGraphAttr.minAxisX = -myMaxTotal;
//          myGraphAttr.maxAxisX = myMaxTotal;
//          myGraphAttr.minAxisY = ratio * myGraphAttr.minAxisX;
//          myGraphAttr.maxAxisY = ratio * myGraphAttr.maxAxisX;
        }
        else {
//          myGraphAttr.minAxisY = -myMaxTotal;
//          myGraphAttr.maxAxisY = myMaxTotal;
//          myGraphAttr.minAxisX = myGraphAttr.minAxisY / ratio;
//          myGraphAttr.maxAxisX = myGraphAttr.maxAxisY / ratio;
        }
      }
      else {
        //...
      }
    }
    // Compute scaling
    myModel2ViewScalarX = (float)myGraphWidth  / (myGraphAttr.maxAxisX-myGraphAttr.minAxisX);
    myModel2ViewScalarY = (float)myGraphHeight / (myGraphAttr.maxAxisY-myGraphAttr.minAxisY);

    //-----------------------------------------------------------
    // Grid lines
    // Compute auto-grid AFTER setting model scalars etc.
    //
    if( myGraphAttr.autoScaleAxes ) {
      myGraphAttr.gridIncX = view2ModelX( myGridViewInc ) - view2ModelX( 0 );
      myGraphAttr.gridIncY = view2ModelY( myGridViewInc ) - view2ModelY( 0 );
    }
    int nGridLinesLeft  = (int)Math.abs( myGraphAttr.minAxisX / myGraphAttr.gridIncX );
    int nGridLinesRight = (int)Math.abs( myGraphAttr.maxAxisX / myGraphAttr.gridIncX );
    myMinGridX = -(nGridLinesLeft * myGraphAttr.gridIncX);
    myNumGridLinesX = nGridLinesLeft + nGridLinesRight + 1;
    int nGridLinesTop    = (int)Math.abs( myGraphAttr.maxAxisY / myGraphAttr.gridIncY );
    int nGridLinesBottom = (int)Math.abs( myGraphAttr.minAxisY / myGraphAttr.gridIncY );
    myMinGridY = -(nGridLinesBottom * myGraphAttr.gridIncY);
    myNumGridLinesY = nGridLinesBottom + nGridLinesTop + 1;
  }
*/
  //----------------------------------------------
  //
  //

  public void mouseClicked( MouseEvent arg0 ) {
  }
  public void mouseEntered( MouseEvent arg0 ) {
  }
  public void mouseExited( MouseEvent arg0 ) {
  }
  public void mousePressed( MouseEvent e ) {
    int xpos = e.getX();
    int ypos = e.getY();

    if (SwingUtilities.isRightMouseButton(e)) {
      myPopupMenu.show( this, xpos, ypos );
    }
  }
  public void mouseReleased( MouseEvent arg0 ) {
  }
  @Override
  public void mouseDragged(MouseEvent e) {
  }
  @Override
  public void mouseMoved(MouseEvent e) {
    int xpos = e.getX();
    int ypos = e.getY();
    float xmodel = this.xView2Model(xpos);
    float ymodel = this.yView2Model(ypos);
    fireGraph2DEvent(xmodel,ymodel);
  }
  //-----------------------------------------------------
  public int getViewPositionVert() {
    return myViewPositionVert;
  }
  public int getViewPositionHorz() {
    return myViewPositionHorz;
  }
//  private void resetViewPositions() {
//    Dimension size = getPreferredSize();
//    Rectangle rect = getVisibleRect();
//    int maxViewPosVert = size.height - rect.height;
//    int maxViewPosHorz = size.width  - rect.width;
//    if( myViewPositionVert > maxViewPosVert ) myViewPositionVert = maxViewPosVert;
//    if( myViewPositionHorz > maxViewPosHorz ) myViewPositionHorz = maxViewPosHorz;
//
//    myPrevViewPositionVert = myViewPositionVert;
//    myPrevViewPositionHorz = myViewPositionHorz;
//  }
  public synchronized void resetViewPositionVert( int scrollValue ) {
    if( scrollValue != myViewPositionVert && !myIsScrolling ) {
      myNextViewPositionVert = -1;
      myViewPositionVert = scrollValue;
      myIsScrolling     = true;
      repaint();
    }
    else if( myIsScrolling ) {
      myNextViewPositionVert = scrollValue;
    }
    else {
      myNextViewPositionVert = -1;
    }
  }
  public synchronized void resetViewPositionHorz( int scrollValue ) {
    if( scrollValue != myViewPositionHorz && !myIsScrolling ) {
      myNextViewPositionHorz = -1;
      myViewPositionHorz    = scrollValue;
      myIsScrolling   = true;
      repaint();
    }
    else if( myIsScrolling ) {
      myNextViewPositionHorz = scrollValue;
    }
    else {
      myNextViewPositionHorz = -1;
    }
  }
  private void fireGraph2DEvent( float xModel, float yModel ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).graph2DValues( xModel, yModel );
    }
  }
  private void fireGraphChangedEvent() {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).graphChanged();
    }
  }
  @Override
  public void curveAttrChanged() {
    repaint();
  }
}

