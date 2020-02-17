/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import cseis.general.csColorMap;
import cseis.general.csCustomColorMap;
import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Properties;

/**
 * Seismic display properties.<br>
 * Store and retrieve seismic display properties from external file.
 * @author Bjorn Olofsson
 */
public class csSeisDispProperties {
  private File myPropertiesFile;
  private final java.util.List<csCustomColorMap> myCustomColorMapList;

  public static final String PROPERTY_zoomVert         = "zoomVert";
  public static final String PROPERTY_zoomHorz         = "zoomHorz";
  public static final String PROPERTY_dispScalar       = "dispScalar";

  public static final String PROPERTY_fullTraceScalar  = "fullTraceScalar";
  public static final String PROPERTY_minValue         = "minValue";
  public static final String PROPERTY_maxValue         = "maxValue";
  public static final String PROPERTY_highlightColor   = "highlightColor";
  public static final String PROPERTY_isHightlightOn   = "isHightlightOn";
  public static final String PROPERTY_wiggleType       = "wiggleType";
  public static final String PROPERTY_viType           = "viType";
  public static final String PROPERTY_scaleType        = "scaleType";
  public static final String PROPERTY_traceScaling     = "traceScaling";
  public static final String PROPERTY_showWiggle       = "showWiggle";
  public static final String PROPERTY_isPosFill        = "isPosFill";
  public static final String PROPERTY_isNegFill        = "isNegFill";
  public static final String PROPERTY_isVariableColor  = "isVariableColor";
  public static final String PROPERTY_isVIDisplay      = "isVIDisplay";
  public static final String PROPERTY_doTraceClipping  = "doTraceClipping";
  public static final String PROPERTY_showZeroLines    = "showZeroLines";
  public static final String PROPERTY_traceClip        = "traceClip";
  public static final String PROPERTY_polarity         = "polarity";
  public static final String PROPERTY_isLogScale       = "isLogScale";
  public static final String PROPERTY_showTimeLines    = "showTimeLines";
  public static final String PROPERTY_isTimeLinesAuto  = "isTimeLinesAuto";
  public static final String PROPERTY_timeLineMajorInc = "timeLineMajorInc";
  public static final String PROPERTY_timeLineMinorInc = "timeLineMinorInc";
  public static final String PROPERTY_plotDirection    = "plotDirection";
  public static final String PROPERTY_wiggleColorPos   = "wiggleColorPos";
  public static final String PROPERTY_wiggleColorNeg   = "wiggleColorNeg";
  public static final String PROPERTY_wiggleColorMap   = "wiggleColorMap";
  public static final String PROPERTY_viColorMap       = "viColorMap";
  public static final String PROPERTY_title            = "title";
  public static final String PROPERTY_titleVertAxis    = "titleVertAxis";
  
  public csSeisDispProperties( String filenameProperties ) {
    this( filenameProperties, null );
  }
  public csSeisDispProperties( String filenameProperties, java.util.List<csCustomColorMap> cmapList ) {
    myPropertiesFile = new File( filenameProperties );

    myCustomColorMapList = cmapList;
  }
  public String filename() {
    return myPropertiesFile.getAbsolutePath();
  }
  /**
   * Retrieve seismic display settings from properties file
   * @param s Seismic display settings
   */
  public void load( csSeisDispSettings s ) throws IOException {
    Properties p = new Properties();
    p.load(new FileInputStream( myPropertiesFile ));
    String value;
    value = p.getProperty( PROPERTY_zoomVert );
    if( value != null ) s.zoomVert        = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_zoomHorz );
    if( value != null )  s.zoomHorz        = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_dispScalar );
    if( value != null )  s.dispScalar      = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_fullTraceScalar );
    if( value != null )  s.fullTraceScalar = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_minValue );
    if( value != null )  s.minValue        = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_maxValue );
    if( value != null )  s.maxValue        = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_highlightColor );
    if( value != null )  s.highlightColor  = new Color( Integer.parseInt( value ) );
    value = p.getProperty( PROPERTY_isHightlightOn );
    if( value != null )  s.isHightlightOn  = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_wiggleType );
    if( value != null )  s.wiggleType      = Integer.parseInt( value );
    value = p.getProperty( PROPERTY_viType );
    if( value != null )  s.viType          = Integer.parseInt( value );
    value = p.getProperty( PROPERTY_scaleType );
    if( value != null )  s.scaleType       = Integer.parseInt( value );
    value = p.getProperty( PROPERTY_traceScaling );
    if( value != null )  s.traceScaling    = Integer.parseInt( value );
    value = p.getProperty( PROPERTY_showWiggle );
    if( value != null )  s.showWiggle      = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_isPosFill );
    if( value != null )  s.isPosFill       = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_isNegFill );
    if( value != null )  s.isNegFill       = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_isVariableColor );
    if( value != null )  s.isVariableColor = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_isVIDisplay );
    if( value != null )  s.isVIDisplay     = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_doTraceClipping );
    if( value != null )  s.doTraceClipping = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_showZeroLines );
    if( value != null )  s.showZeroLines   = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_traceClip );
    if( value != null )  s.traceClip       = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_polarity );
    if( value != null )  s.polarity        = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_isLogScale );
    if( value != null )  s.isLogScale      = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_showTimeLines );
    if( value != null )  s.showTimeLines   = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_isTimeLinesAuto );
    if( value != null )  s.isTimeLinesAuto = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_timeLineMajorInc );
    if( value != null )  s.timeLineMajorInc = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_timeLineMinorInc );
    if( value != null )  s.timeLineMinorInc = Float.parseFloat( value );
    value = p.getProperty( PROPERTY_plotDirection );
    if( value != null )  s.plotDirection    = Integer.parseInt( value );
    value = p.getProperty( PROPERTY_wiggleColorPos );
    if( value != null )  s.wiggleColorPos  = new Color( Integer.parseInt( value ) );
    value = p.getProperty( PROPERTY_wiggleColorNeg );
    if( value != null )  s.wiggleColorNeg  = new Color( Integer.parseInt( value ) );
    value = p.getProperty( PROPERTY_wiggleColorMap );
    if( value != null ) {
      int colorMapIndex = csColorMap.getDefaultMapIndex( value );
      if( colorMapIndex >= 0 ) {
        s.wiggleColorMap  = new csColorMap( colorMapIndex, csColorMap.COLOR_MAP_TYPE_32BIT );
      }
      else {
        csCustomColorMap cmap = getCustomColorMap( value );
        if( cmap != null ) s.wiggleColorMap = new csCustomColorMap( cmap );
      }
    }
    value = p.getProperty( PROPERTY_viColorMap );
    if( value != null ) {
      int colorMapIndex = csColorMap.getDefaultMapIndex( value );
      if( colorMapIndex >= 0 ) {
        s.viColorMap  = new csColorMap( colorMapIndex, csColorMap.COLOR_MAP_TYPE_32BIT );
      }
      else {
        csCustomColorMap cmap = getCustomColorMap( value );
        if( cmap != null ) {
          s.viColorMap = new csCustomColorMap( cmap );
        }
      }
    }
    value = p.getProperty( PROPERTY_title );
    if( value != null ) s.title = value;
    value = p.getProperty( PROPERTY_titleVertAxis );
    if( value != null ) s.titleVertAxis = value;
  }
  /**
   * Write seismic display settings to properties file
   * @param s Seismic display settings
   */
  public void write( csSeisDispSettings s ) throws IOException {
    Properties p = new Properties();
    p.setProperty( PROPERTY_zoomVert, s.zoomVert+"" );
    p.setProperty( PROPERTY_zoomHorz, s.zoomHorz+"" );
    p.setProperty( PROPERTY_dispScalar, s.dispScalar+"" );
    p.setProperty( PROPERTY_fullTraceScalar, s.fullTraceScalar+"" );
    p.setProperty( PROPERTY_minValue, s.minValue+"" );
    p.setProperty( PROPERTY_maxValue, s.maxValue+"" );
    p.setProperty( PROPERTY_highlightColor, s.highlightColor.getRGB()+"" );
    p.setProperty( PROPERTY_isHightlightOn, s.isHightlightOn+"" );
    p.setProperty( PROPERTY_wiggleType, s.wiggleType+"" );
    p.setProperty( PROPERTY_viType, s.viType+"" );
    p.setProperty( PROPERTY_scaleType, s.scaleType+"" );
    p.setProperty( PROPERTY_traceScaling, s.traceScaling+"" );
    p.setProperty( PROPERTY_showWiggle, s.showWiggle+"" );
    p.setProperty( PROPERTY_isPosFill, s.isPosFill+"" );
    p.setProperty( PROPERTY_isNegFill, s.isNegFill+"" );
    p.setProperty( PROPERTY_isVariableColor, s.isVariableColor+"" );
    p.setProperty( PROPERTY_isVIDisplay, s.isVIDisplay+"" );
    p.setProperty( PROPERTY_doTraceClipping, s.doTraceClipping+"" );
    p.setProperty( PROPERTY_showZeroLines, s.showZeroLines+"" );
    p.setProperty( PROPERTY_traceClip, s.traceClip+"" );
    p.setProperty( PROPERTY_polarity, s.polarity+"" );
    p.setProperty( PROPERTY_isLogScale, s.isLogScale+"" );
    p.setProperty( PROPERTY_showTimeLines, s.showTimeLines+"" );
    p.setProperty( PROPERTY_isTimeLinesAuto, s.isTimeLinesAuto+"" );
    p.setProperty( PROPERTY_timeLineMajorInc, s.timeLineMajorInc+"" );
    p.setProperty( PROPERTY_timeLineMinorInc, s.timeLineMinorInc+"" );
    p.setProperty( PROPERTY_plotDirection, s.plotDirection+"" );
    p.setProperty( PROPERTY_wiggleColorPos, s.wiggleColorPos.getRGB()+"" );
    p.setProperty( PROPERTY_wiggleColorNeg, s.wiggleColorNeg.getRGB()+"" );
    p.setProperty( PROPERTY_wiggleColorMap, s.wiggleColorMap.toString() );
    p.setProperty( PROPERTY_viColorMap, s.viColorMap.toString() );

    p.store(new FileOutputStream( myPropertiesFile ), null);
  }
  private csCustomColorMap getCustomColorMap( String name ) {
    for( csCustomColorMap cmap : myCustomColorMapList ) {
      if( name.compareTo( cmap.toString() ) == 0 ) {
        return cmap;
      }
    }
    return null;
  }

}

