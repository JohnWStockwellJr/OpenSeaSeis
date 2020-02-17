/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

import java.awt.Color;
import java.awt.Font;

/**
 * Attributes of 2D graph csGraph2D.
 * @author Bjorn Olofsson
 */
public class csGraphAttributes {
  public static final int DEFAULT_GRID_PIXEL_INC = 40;
  public static final int DEFAULT_BORDER_PADDING = 8;
  public static final int DEFAULT_INNER_BORDER_PADDING = 2;
  public static final int DEFAULT_BORDER_LABEL   = 20;
  public static final int DEFAULT_BORDER_TITLE   = 30;

  public static final int AXIS_SCALE_LINEAR = 95;
  public static final int AXIS_SCALE_LOG    = 96;
//  public static final int AXIS_SCALE_LOG_POS = 96;

  /// True if graph has a second Y axis
  public boolean hasY2Axis;
  public boolean centreAxisX;
  public boolean centreAxisY;
  /// Approximate increment between grid/major tic lines, in screen pixels
  public int gridXPixelInc;
  public int gridYPixelInc;
  /// Grid/main tic increment in main model values
  public float gridIncX;
  public float gridIncY;
  /// Ratio between major and minor X tics
  public int minorXTicRatio;
  /// Ratio between major and minor Y tics
  public int minorYTicRatio;
  ///
  public boolean isTrueScale;

  public Color backgroundColor;
  public Color gridColor;
  public Color axisColor;
  public Color borderColor;
  public Color innerBorderColor;

  public boolean fixHorzGraphPosition;
  public int horzGraphPositionMin;
  public int horzGraphPositionMax;

  public float minAxisX;
  public float maxAxisX;
  public float minAxisY;
  public float maxAxisY;
  public float minAxisY2;
  public float maxAxisY2;
  /// True if graph has fixed X dimension. Virtual graph width is constant, so model2view X scalar is constant
  public boolean isFixedXDim;
  /// True if graph has fixed Y dimension. Virtual graph width is constant, so model2view Y scalar is constant
  public boolean isFixedYDim;

  /// True to show border
  public boolean showBorder;
  /// True to show inner graph border
  public boolean showInnerBorder;
  /// Show grid
  public boolean showGrid;
  /// Show xtics
  public boolean showXTics;
  public boolean showYTics;
  /// Show X axis
  public boolean showAxisX;
  /// Show Y axis
  public boolean showAxisY;
  /// Show axis annotation. Is only shown for axes which have showAxisX/showAxisY set
  public boolean showAxisXAnnotation;
  public boolean showAxisYAnnotation;
  /// Show zero axis
  public boolean showZeroAxis;

  //--------------------------------------------------
  // Layout sizes
  //
  /// Insets around whole graph view, outside of labels, title and border
  public int insetLeft;
  public int insetRight;
  public int insetTop;
  public int insetBottom;
  /// Padding outside of actual graph
  public int graphPadding;
  /// Size of border line
  public int borderSize;
  /// Size of border line
  public int innerBorderSize;
  /// Padding on outside of border
  public int borderPadding;
  public int innerBorderPadding;
  /// Area allocated for XY axis annotation, if any
  public int annotationSizeXAxis;
  public int annotationSizeYAxis;

  public Color labelColor;
  public String xLabel;
  public String yLabel;
  public String yLabel2;
  public String title;
  public String innerTitle;
  public Font titleFont;
  public Font labelFont;
  public Font annotationFont;
  // true if graph curves are plotted within 'view port', and border etc stays visible at same position
  public boolean isGraphViewPort;

  /// True if Y reference value shall be used:
  public boolean useRefValueY;
  /// Maximum reference value for Y axis, used for dB display
  public float maxRefValueY;
  
  //------------------------------------------------------
  // Attributes which are not really used anywhere yet...
  //
  /// Linear/log scale
  public int axisScaleX;
  /// Linear/log scale
  public int axisScaleY;
  public boolean autoScaleAxes;

  public csGraphAttributes() {
    gridXPixelInc  = DEFAULT_GRID_PIXEL_INC;
    gridYPixelInc  = DEFAULT_GRID_PIXEL_INC;
    insetLeft     = 10;
    insetRight    = 10;
    insetTop      = 10;
    insetBottom   = 10;
    autoScaleAxes = true;
    axisScaleX  = AXIS_SCALE_LINEAR;
    axisScaleY  = AXIS_SCALE_LINEAR;
    hasY2Axis = false;
    showZeroAxis = true;
    showAxisX   = true;
    showAxisY   = true;
    showAxisXAnnotation = true;
    showAxisYAnnotation = true;
    centreAxisX = true;
    centreAxisY = true;
    showBorder  = true;
    showInnerBorder = false;
    showGrid    = true;
    isTrueScale = true;
    gridIncX    = 0;
    gridIncY    = 0;
    backgroundColor = Color.white;
    gridColor = Color.lightGray;
    axisColor = Color.black;
    borderColor = Color.black;
    innerBorderColor = Color.gray;
    isGraphViewPort = false;
    showXTics = true;
    showYTics = true;
    minorXTicRatio = 5;
    minorYTicRatio = 5;
    xLabel = "X axis";
    yLabel = "Y axis";
    yLabel2 = "Y2 axis";
    title  = "2D graph";
    innerTitle  = "";
    isFixedXDim = false;
    isFixedYDim = false;
    graphPadding = 10;
    fixHorzGraphPosition = false;
    horzGraphPositionMin = 0;
    horzGraphPositionMax = 0;
    borderSize  = 2;
    innerBorderSize = 1;
    borderPadding = DEFAULT_BORDER_PADDING;
    innerBorderPadding = DEFAULT_INNER_BORDER_PADDING;
    labelColor = Color.darkGray;
    labelFont = new Font("SansSerif", Font.PLAIN, 12 );
    titleFont = new Font("SansSerif", Font.BOLD, 14 );
    annotationFont = new Font("SansSerif", Font.PLAIN, 11 );
    annotationSizeXAxis = 30;
    annotationSizeYAxis = 30;

    useRefValueY = false;
    maxRefValueY = 0;

    minAxisX    = -1;
    maxAxisX    = 1;
    minAxisY    = -1;
    maxAxisY    = 1;
    minAxisY2   = -1;
    maxAxisY2   = 1;
  }
  public csGraphAttributes( csGraphAttributes a ) {
    gridXPixelInc  = a.gridXPixelInc;
    gridYPixelInc  = a.gridYPixelInc;
    insetLeft     = a.insetLeft;
    insetTop      = a.insetTop;
    insetRight    = a.insetRight;
    insetBottom   = a.insetBottom;
    showZeroAxis  = a.showZeroAxis;
    showAxisX     = a.showAxisX;
    showAxisY     = a.showAxisY;
    showAxisXAnnotation = a.showAxisXAnnotation;
    showAxisYAnnotation = a.showAxisYAnnotation;
    showXTics       = a.showXTics;
    showYTics       = a.showYTics;
    minorXTicRatio  = a.minorXTicRatio;
    minorYTicRatio  = a.minorYTicRatio;
    hasY2Axis       = a.hasY2Axis;
    backgroundColor = a.backgroundColor;
    gridColor       = a.gridColor;
    axisColor       = a.axisColor;
    borderColor     = a.borderColor;
    innerBorderColor= a.innerBorderColor;

    isGraphViewPort = a.isGraphViewPort;

    centreAxisX = a.centreAxisX;
    centreAxisY = a.centreAxisY;
    showBorder  = a.showBorder;
    showInnerBorder  = a.showInnerBorder;

    annotationFont = a.annotationFont;
    annotationSizeXAxis = a.annotationSizeXAxis;
    annotationSizeYAxis = a.annotationSizeYAxis;

    showGrid    = a.showGrid;
    gridIncX    = a.gridIncX;
    gridIncY    = a.gridIncY;

    isTrueScale = a.isTrueScale;
    xLabel      = a.xLabel;
    yLabel      = a.yLabel;
    yLabel2     = a.yLabel2;
    title       = a.title;
    innerTitle  = a.innerTitle;
    graphPadding = a.graphPadding;
    borderSize  = a.borderSize;
    innerBorderSize  = a.innerBorderSize;
    borderPadding = a.borderPadding;
    innerBorderPadding = a.innerBorderPadding;
    isFixedXDim = a.isFixedXDim;
    isFixedYDim = a.isFixedYDim;
    labelColor  = a.labelColor;
    labelFont   = a.labelFont;
    titleFont   = a.titleFont;
    useRefValueY = a.useRefValueY;
    maxRefValueY = a.maxRefValueY;

    fixHorzGraphPosition = a.fixHorzGraphPosition;
    horzGraphPositionMin = a.horzGraphPositionMin;
    horzGraphPositionMax = a.horzGraphPositionMax;

    autoScaleAxes = a.autoScaleAxes;
    axisScaleX  = a.axisScaleX;
    axisScaleY  = a.axisScaleY;

    minAxisX    = a.minAxisX;
    maxAxisX    = a.maxAxisX;
    minAxisY    = a.minAxisY;
    maxAxisY    = a.maxAxisY;
    minAxisY2   = a.minAxisY2;
    maxAxisY2   = a.maxAxisY2;
  }
}


