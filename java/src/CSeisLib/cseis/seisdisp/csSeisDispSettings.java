/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.Color;
import cseis.general.csColorMap;
import cseis.general.csCustomColorMap;

/**
 * Seismic display settings.<br>
 * Set of user settings that define the way seismic traces are displayed.
 * @author Bjorn Olofsson
 */
public class csSeisDispSettings {

  public static final int WIGGLE_FILL_NONE   = 0;
  public static final int WIGGLE_FILL_POS    = 1;
  public static final int WIGGLE_FILL_NEG    = 2;
  public static final int WIGGLE_FILL_POSNEG = 3;
  
  public static final int WIGGLE_COLOR_FIXED    = 41;
  public static final int WIGGLE_COLOR_VARIABLE = 42;

  public static final int SCALE_TYPE_SCALAR = 31;
  public static final int SCALE_TYPE_RANGE  = 32;
  public static final int SCALE_TYPE_TRACE  = 33;

  public static final int TRACE_SCALING_MAXIMUM = 51;
  public static final int TRACE_SCALING_AVERAGE = 52;
  
  public static final int WIGGLE_TYPE_LINEAR = 11;
  public static final int WIGGLE_TYPE_CUBIC  = 12;

  public static final float POLARITY_NORMAL    = 1.0f;
  public static final float POLARITY_REVERSED  = -1.0f;

  public static final int VA_TYPE_2DSPLINE = 21;
  public static final int VA_TYPE_VERTICAL = 22;
  public static final int VA_TYPE_DISCRETE = 23;

  public static final int PLOT_DIR_VERTICAL   = 0;
  public static final int PLOT_DIR_HORIZONTAL = 1;

  /// Vertical zoom factor = Screen pixels per sample
  public float zoomVert;
  /// Horizontal zoom factor = Screen pixels per trace
  public float zoomHorz;
  /// Scalar to be applied to amplitudes
  public float dispScalar;
  /// Scalar to be applied to amplitudes, for full trace scaling
  public float fullTraceScalar;
  /// Min value when range is specified
  public float minValue;
  /// Max value when range is specified
  public float maxValue;
  
  /// Highlight color
  public Color highlightColor;
  /// Highlight flag
  public boolean isHightlightOn;

  // Attributes
//  public int fillWiggle;
  /// Wiggle type (linear,cubic)
  public int wiggleType;
  /// Variable intensity type (discrete,vertical,spline)
  public int viType;
  /// 'Scale' type: How is seismic trace scaled? By scalar, range, or trace-by-trace
  public int scaleType;
  /// Further definition of trace-by-trace scaling: By trace maximum, or average
  public int traceScaling;
  /// Show or hide wiggle trace
  public boolean showWiggle;
  /// true if positive wiggle is filled
  public boolean isPosFill;
  /// true if negative wiggle is filled
  public boolean isNegFill;
  /// true if wiggle area shall be filled by variable color
  public boolean isVariableColor;
  /// Show or hide variable intensity display
  public boolean isVIDisplay;
  /// Apply trace clipping yes/no
  public boolean doTraceClipping;
  /// Show/hide zero line
  public boolean showZeroLines;
  /// Number of traces where seismic wiggle is clipped
  public float traceClip;
  /// Polarity of display (+/-1.0)
  public float polarity;
  /// true if vertical axis shall be plotted in log scale
  public boolean isLogScale;
  /// Title of seismic display
  public String  title;
  /// Title of vertical axis
  public String  titleVertAxis;

  /// Color of positive wiggle area
  public Color wiggleColorPos;
  /// Color of negative wiggle area
  public Color wiggleColorNeg;
  /// Color map of positive/negative wiggle area (in case of variable color display)
  public csColorMap wiggleColorMap;
  /// Color map of variable intensity display
  public csColorMap viColorMap;
  
  /// Show/hide horizontal time lines
  public boolean showTimeLines;
  /// true if horizontal lines shall be computed automatically
  public boolean isTimeLinesAuto;
  /// in vertical units (ms or Hz..)
  public double timeLineMinorInc;
  /// in vertical units (ms or Hz..)
  public double timeLineMajorInc;

  /// Plot direction: Vertical or horizontal
  public int plotDirection;

  public csSeisDispSettings( csSeisDispSettings ds ) {
    set( ds );
  }
  public void set( csSeisDispSettings ds ) {
    if( ds == this ) return;
    zoomVert   = ds.zoomVert;
    zoomHorz   = ds.zoomHorz;
    traceClip  = ds.traceClip;
    dispScalar = ds.dispScalar;
    fullTraceScalar = ds.fullTraceScalar;
    scaleType       = ds.scaleType;
    traceScaling    = ds.traceScaling;
    polarity        = ds.polarity;
    minValue     = ds.minValue;
    maxValue     = ds.maxValue;
    title        = ds.title;
    
    highlightColor = ds.highlightColor;
    isHightlightOn = ds.isHightlightOn;

    // Wiggle settings
//    fillWiggle = WIGGLE_FILL_POS;
    wiggleType      = ds.wiggleType;
    isVariableColor = ds.isVariableColor;
    doTraceClipping = ds.doTraceClipping;
    showWiggle      = ds.showWiggle;
    isPosFill       = ds.isPosFill;
    isNegFill       = ds.isNegFill;
    isVIDisplay     = ds.isVIDisplay;
    viType          = ds.viType;
    showZeroLines   = ds.showZeroLines;

    // Time lines
    showTimeLines = ds.showTimeLines;
    timeLineMinorInc = ds.timeLineMinorInc;
    timeLineMajorInc = ds.timeLineMajorInc;
    isTimeLinesAuto  = ds.isTimeLinesAuto;
    titleVertAxis = ds.titleVertAxis;

    // Colors
    wiggleColorPos = ds.wiggleColorPos;
    wiggleColorNeg = ds.wiggleColorNeg;
    // Setting color maps as follows is problematic! Should be changed to something more stable.
    // Problems occur with min/max values and scalars
    if( wiggleColorMap == null ) {
      wiggleColorMap = new csCustomColorMap( ds.wiggleColorMap );
    }
    else if( wiggleColorMap != ds.wiggleColorMap ) {
      wiggleColorMap.resetColors( ds.wiggleColorMap );
    }
    if( viColorMap == null ) {
      viColorMap     = new csCustomColorMap( ds.viColorMap );
    }
    else if( viColorMap != ds.viColorMap ) {
      viColorMap.resetColors( ds.viColorMap );
    }    
    plotDirection = ds.plotDirection;
    isLogScale    = ds.isLogScale;
  }
  public csSeisDispSettings() {
    setDefaults();
  }
  public void setDefaults() {
    zoomVert   = 1.0f;
    zoomHorz   = 4.0f;
    traceClip  = 2.0f;
    dispScalar = 1.0f;
    fullTraceScalar = 1.0f;
    scaleType    = SCALE_TYPE_SCALAR;
    traceScaling = TRACE_SCALING_AVERAGE;
    polarity     = POLARITY_NORMAL;
    minValue     = 0.0f;
    maxValue     = 0.0f;
    title        = "";
    titleVertAxis = "DEFAULT";
    
    highlightColor = Color.yellow;
    isHightlightOn = false;

    // Wiggle settings
//    fillWiggle = WIGGLE_FILL_POS;
    wiggleType      = WIGGLE_TYPE_LINEAR;
    isVariableColor = false;
    doTraceClipping = true;
    showWiggle      = true;
    isPosFill       = true;
    isNegFill       = false;
    isVIDisplay     = false;
    viType          = VA_TYPE_DISCRETE;
    showZeroLines   = false;

    // Time lines
    showTimeLines = true;
    timeLineMinorInc = 100;
    timeLineMajorInc = 500;
    isTimeLinesAuto  = true;

    // Colors
    wiggleColorPos = Color.black;
    wiggleColorNeg = Color.lightGray;
    wiggleColorMap = new csColorMap( csColorMap.GRAY_WB, csColorMap.COLOR_MAP_TYPE_32BIT );
    viColorMap     = new csColorMap( csColorMap.DEFAULT, csColorMap.COLOR_MAP_TYPE_32BIT );

    plotDirection = PLOT_DIR_VERTICAL;
    isLogScale    = false;
  }
  public void dump() {
    System.out.println("zoomVert   " + zoomVert);
    System.out.println("zoomHorz   " + zoomHorz);
    System.out.println("traceClip  " + traceClip);
    System.out.println("dispScalar " + dispScalar);
    System.out.println("fullTraceScalar " + fullTraceScalar);
    System.out.println("scaleType       " + scaleType);
    System.out.println("traceScaling    " + traceScaling);
    System.out.println("polarity        " + polarity);
    System.out.println("minValue     " + minValue);
    System.out.println("maxValue     " + maxValue);

    System.out.println("highlightColor " + highlightColor.getRGB());
    System.out.println("isHightlightOn " + isHightlightOn);

    // Wiggle settings
//    fillWiggle = WIGGLE_FILL_POS);
    System.out.println("wiggleType      " + wiggleType);
    System.out.println("isVariableColor " + isVariableColor);
    System.out.println(" doTraceClipping " + doTraceClipping);
    System.out.println("showWiggle      " + showWiggle);
    System.out.println("isPosFill       " + isPosFill);
    System.out.println("isNegFill       " + isNegFill);
    System.out.println("isVADisplay     " + isVIDisplay);
    System.out.println("vaType          " + viType);
    System.out.println("showZeroLines   " + showZeroLines);

    // Time lines
    System.out.println("showTimeLines " + showTimeLines);
    System.out.println("timeLineMinorInc " + timeLineMinorInc);
    System.out.println(" timeLineMajorInc " + timeLineMajorInc);
    System.out.println(" isTimeLinesAuto  " + isTimeLinesAuto);

    // Colors
    System.out.println(" wiggleColorPos " + wiggleColorPos.getRGB());
    System.out.println("  wiggleColorNeg " + wiggleColorNeg.getRGB());
    System.out.println(" wiggleColorMap " + wiggleColorMap );
    System.out.println("  vaColorMap     " + viColorMap );
    System.out.println(" plotDirection " + plotDirection);
    System.out.println("  isLogScale    " + isLogScale);

    System.out.println("  title    " + title);
    System.out.println("  titleVertAxis    " + titleVertAxis);
  }
}


