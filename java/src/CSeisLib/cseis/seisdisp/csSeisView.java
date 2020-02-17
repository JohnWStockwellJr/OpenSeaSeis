/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import cseis.cmapgen.csCustomColorMapModel;
import cseis.general.csColorMap;
import cseis.seis.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.*;
import java.util.*;

/**
 * Seismic display/view.<br>
 * See also csSeisPane...<br>
 * <br>
 * csSeisView provides:<br>
 * <ul>
 *  <li> rendering of seismic data in the most common ways: wiggle and variable density display (discrete and interpolated)
 *  <li> user settings for scaling and colouring,
 *  <li> scrolling and zooming,
 *  <li> time lines,
 *  <li> trace highlighting,
 *  <li> popup dialog where display settings can be changed,
 *  <li> trace and sample information at current mouse position is sent to listeners,
 *  <li> overlays: Arbitrary Java graphics objects can be painted on top of the seismic view by the client application, through listeners,
 *  <li> Direction of trace plotting: Vertical (normal) or horizontal.<br>
 *    When directon of plotting is horizontal, all parameters are kept the same, only when drawing is performed, x plots as y, and y plots as x.
 *    Horizontal plotting is not fully implemented yet!
 * </ul>
 * @author Bjorn Olofsson
 */

@SuppressWarnings("serial")
public class csSeisView extends JPanel {
  //
  public static final int    INSET_DEFAULT = 10;
  public static final double EPSILON       = 10.0e-30;
  
  public static final int DEFAULT_INC_PIXELS = 50;
  public static final double DEFAULT_INC_MS [] =
  { 1,  2, 2.5,  5,  10 };
  public static final int DEFAULT_MAJOR_MS [] =
  { 10, 10, 50, 50, 100 };
  
  protected String myTEMPName = "DEFAULT";
  protected csISeismicTraceBuffer myTraceBuffer;
  protected double mySampleInt;
  protected int myNumTraces;    // Number of trace in seismic buffer
  protected int myNumSamples;   // Number of samples in seismic buffer

  protected int myMarginLeftRight;  // Width of inset (margin on left and right hand side)
  protected int myMarginTopBottom;  // Height of inset (margin at top and bottom)
  
  // Scrolling = View position of seismic display
  /// Current view position ( = scroll position)
  protected int myViewPositionVert = 0;
  protected int myViewPositionHorz = 0;
  protected int myPrevViewPositionVert = 0;
  protected int myPrevViewPositionHorz = 0;
  protected int myNextViewPositionVert = -1;
  protected int myNextViewPositionHorz = -1;
  protected boolean myIsScrolling;

  protected int myImageWidth;  // Width (in pixels) of bitmap image
  protected int myImageHeight; // Height (in pixels) of bitmap image
  protected boolean myIsFullRepaint = false;

  // Highlighting
  private int myHighlightTrace = 0;
  private boolean myIsHighlightOn = false;
  private boolean myShowCrosshair = false;
  private boolean mySnapCrosshair = false;
  private boolean myIsCrosshairPainted = false;
  private Point myCrosshairPosition = new Point(-1,-1);

  protected ArrayList<csISeisOverlay> mySeisOverlays;
  private ArrayList<csISeisViewListener> mySeisViewListeners;
  private csIVisibleRectListener myRectListener;

  /// Constant scalar to apply to each trace. This is usually 1.
  /// When using 'trace' scaling, the trace scalar becomes the inverse of the mean or max of each individual trace.
  protected float[] myTraceScalar;
  /// The 'wiggle central amplitude' is usually 0; this is where the wiggle 'zeor' line is plotted.
  /// When using 'range' scaling, the wiggle central amplitude becomes the mean value of the specified min/max range values.
  protected float myWiggleCentreAmplitude;
  /// General display scalar. This local field is for convenience. The display scalar is also stored in the display settings (mySettings)
  protected float myDispScalar;
  /// Trace 'step'. Only every myStepTrace'th trace is drawn to the screen.
  /// This is to enhance performance when the display is squeezed below one pixel/trace.
  protected int myStepTrace = 1;

  // Data buffer
  protected BufferedImage myBitmap;
//  private VolatileImage myVolatile;
  protected Image myVolatileBitmap;

  
  protected csSeisDispSettings mySettings;
  protected csSeisDispDialog myDialogSettings;
  
  protected csSeisViewPopupMenu myPopupMenu;
  protected JFrame myParentFrame;
  protected csSeisViewEventHandler myEventHandler;
  
// For test purposes:
//private Color[] myTMPColors = { Color.white, Color.blue, Color.yellow, Color.red, Color.gray, Color.orange };
//private int myTMPColorIndex = 0;
//private int myNumColors = myTMPColors.length;

  private int myColorBitType;
  private csSeisViewRenderer myRenderer;
  protected float myLogScaleRatio = 1.0f;
  /// Switch between horizontal & verticall scrolling
  private boolean mySwitchScroll = true;
  /// Mouse mode
  private int myMouseMode;
  
  /// Custom color maps: "Model" that can be used in several "views" such as in JComboBox, JList...
  /// If not set, a default color map will be used
  private csCustomColorMapModel myCustomColorMapModel;

  //-----------------------------------------------------------------
  //
  
  /**
   * Constructor.
   * Use this constructor only if seismic traces are already available when constructing this seismic view.
   * 
   * @param parentFrame  Parent fram
   * @param buffer       Trace buffer
   * @param sampleInt    Sample interval [ms]
   */
  public csSeisView( JFrame parentFrame, csTraceBuffer buffer, double sampleInt ) {
    this( parentFrame, null );
    updateTraceBuffer( buffer, sampleInt, true );
  }
  public csSeisView( JFrame parentFrame ) {
    this( parentFrame, null );
  }
  public csSeisView( JFrame parentFrame, csIVisibleRectListener rectListener ) {
    this( parentFrame, null, new csCustomColorMapModel(csColorMap.COLOR_MAP_TYPE_32BIT), csColorMap.COLOR_MAP_TYPE_32BIT, rectListener );
  }
//  public csSeisView( JFrame parentFrame, csSeisDispSettings dispSettings ) {
//    this( parentFrame, dispSettings, new csCustomColorMapModel(csColorMap.COLOR_MAP_TYPE_32BIT) );
//  }
//  public csSeisView( JFrame parentFrame, csSeisDispSettings dispSettings, csCustomColorMapModel cmapModel ) {
//    this( parentFrame, dispSettings, cmapModel, csColorMap.COLOR_MAP_TYPE_32BIT );
//  }
  /**
   * Constructor.
   * 
   * @param parentFrame  Parent JFrame, used for internal dialog windows..
   * @param dispSettings Display settings. Pass 'null' to use default display settings.
   * @param cmapModel    Custom color map model
   * @param colorMapType Color map type
   */
  public csSeisView( JFrame parentFrame, csSeisDispSettings dispSettings, csCustomColorMapModel cmapModel, int colorMapType ) {
    this( parentFrame, dispSettings, cmapModel, colorMapType, null );
  }
  public csSeisView( JFrame parentFrame, csSeisDispSettings dispSettings, csCustomColorMapModel cmapModel, int colorMapType, csIVisibleRectListener rectListener ) {
    super( new BorderLayout() );
    myParentFrame = parentFrame;
    myCustomColorMapModel = cmapModel;
    
    myTraceBuffer           = null;
    myTraceScalar           = null;
    myWiggleCentreAmplitude = 0.0f;
    myDispScalar = 1.0f;
    myBitmap = null;
    
    myColorBitType = csColorMap.COLOR_MAP_TYPE_32BIT;
    mySettings = dispSettings;
    if( mySettings == null ) {
      mySettings = new csSeisDispSettings();
      mySettings.wiggleColorMap.setColorMapType( csColorMap.COLOR_MAP_TYPE_32BIT );
      mySettings.viColorMap.setColorMapType( csColorMap.COLOR_MAP_TYPE_32BIT );
    }
    myRenderer = new csSeisViewRenderer32bit( this, mySettings );

    if( !java.awt.GraphicsEnvironment.isHeadless() ) {
      myDialogSettings    = new csSeisDispDialog( this, getDispSettings(), myCustomColorMapModel );
    }

    mySeisOverlays      = new ArrayList<csISeisOverlay>(1);
    mySeisViewListeners = new ArrayList<csISeisViewListener>(1);
    
    myIsScrolling     = false;

    myMarginLeftRight = 0;
    myMarginTopBottom = 0;
    myViewPositionVert = 0;
    myViewPositionHorz = 0;
    myPrevViewPositionVert = 0;
    myPrevViewPositionHorz = 0;
    myNextViewPositionVert = -1;
    myNextViewPositionHorz = -1;
    myMouseMode = csMouseModes.NO_MODE;

    myPopupMenu    = new csSeisViewPopupMenu( this );
    myEventHandler = new csSeisViewEventHandler( this );
    KeyboardFocusManager manager = KeyboardFocusManager.getCurrentKeyboardFocusManager();
    manager.addKeyEventDispatcher(new SeisViewKeyDispatcher(myEventHandler));
    if( myDialogSettings != null ) addSeisViewListener( myDialogSettings );
    addSeisViewListener( myRenderer );
    myRectListener = rectListener;
    if( myRectListener == null ) {
      myRectListener = new csIVisibleRectListener() {
        @Override
        public Rectangle getVisibleRect() {
          return csSeisView.this.getVisibleRect();
        }
      };
    }
  }
//-----------------------------------------------------
  public int getMouseMode() {
    return myMouseMode;
  }
//-----------------------------------------------------
  public void setMouseMode( int mouseMode ) {
    myMouseMode = mouseMode;
//    if( mouseMode == csMouseModes.NO_MODE ) setCursor( java.awt.Cursor.DEFAULT_CURSOR );
  }
  /**
   * Show display setting dialog
   *
   */
  public void showDispSettingsDialog() {
    if( myDialogSettings != null ) {
    	myDialogSettings.dispose();
      if( !myDialogSettings.isVisible() ) myDialogSettings.setVisible(true);
      else myDialogSettings.toFront();
    }
  }
  //-----------------------------------------------------
  /**
   * Update seismic trace buffer
   * 
   * @param buffer  New trace buffer
   * @param refresh true if display shall be fully refreshed
   */
  public void updateTraceBuffer( csISeismicTraceBuffer buffer, boolean refresh ) {
    updateTraceBuffer( buffer, mySampleInt, refresh );
  }
  /**
   * Update seismic trace buffer
   * 
   * @param buffer    New trace buffer
   * @param sampleInt Sample interval [ms]
   * @param refresh   true if display shall be fully refreshed
   */
  public void updateTraceBuffer( csISeismicTraceBuffer buffer, double sampleInt, boolean refresh ) {
    boolean firstCall = (myTraceBuffer == null);
    myTraceBuffer = buffer;
    mySampleInt   = sampleInt;
    
    myNumSamples = myTraceBuffer.numSamples();
    myNumTraces  = myTraceBuffer.numTraces();
    myTraceScalar = new float[myNumTraces];
    for( int itrc = 0; itrc < myNumTraces; itrc++ ) {
      myTraceScalar[itrc] = 1.0f;      
    }
    if( myNumSamples <= 0 || myNumTraces <= 0 ) return;
    myHighlightTrace = 0;
    
    if( refresh ) {
      // Estimate best scalar: Analyse five random traces
      int step = (int)(myNumTraces/5) + 1;
      float meanAbs = 0.0f;
      int counter = 0;
      for( int itrc = 0; itrc < myNumTraces; itrc += step ) {
        float[] samples = myTraceBuffer.samples( itrc );
        for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
          if( samples[isamp] == 0.0f ) continue;
          meanAbs += Math.abs( samples[isamp] );
          counter += 1;
        }
      }
      if( meanAbs != 0 && counter != 0 ) {
        meanAbs /= (float)counter;
        mySettings.dispScalar = 0.2f/meanAbs;
      }
      else {
        mySettings.dispScalar = 1.0f;
      }
    }
    if( firstCall && mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_RANGE ) {
      mySettings.minValue = myTraceBuffer.minTotalAmplitude();
      mySettings.maxValue = myTraceBuffer.maxTotalAmplitude();
    }
    resetTraceScaling();
    resetLeftRightMargins();

    myRenderer.reset( buffer, sampleInt );
    
    resetScreenBuffer();
    resetViewPositions();
    resetTimeLines();

    fireEventChangedSettings();

    myIsFullRepaint  = true;
    myIsScrolling    = false;

    fireEventChangedTraceBuffer( buffer );
  }
  /**
   * Set color bit type
   * 
   * @param colorMapType (see @class csColorMap)
   */
  public void setColorBits( int colorMapType ) {
    myColorBitType = colorMapType;
    removeSeisViewListener( myDialogSettings );
    removeSeisViewListener( myRenderer );

    mySettings.wiggleColorMap.setColorMapType( colorMapType );
    mySettings.viColorMap.setColorMapType( colorMapType );
    if( colorMapType == csColorMap.COLOR_MAP_TYPE_8BIT ) {
      myRenderer = new csSeisViewRenderer8bit( this, mySettings );
    }
    else if( colorMapType == csColorMap.COLOR_MAP_TYPE_32BIT ) {
      myRenderer = new csSeisViewRenderer32bit( this, mySettings );
    }
//    myDialogSettings    = new csSeisDispDialog( this, colorMapType, myCustomColorMapModel );
//    addSeisViewListener( myDialogSettings );
    addSeisViewListener( myRenderer );
    if( myTraceBuffer != null ) {
      updateTraceBuffer( myTraceBuffer, mySampleInt, false );
    }
  }
  private void resetViewPositions() {
    Dimension size = getPreferredSize();
    Rectangle rect = myRectListener.getVisibleRect();
    int maxViewPosVert = size.height - rect.height;
    int maxViewPosHorz = size.width  - rect.width;
    if( myViewPositionVert > maxViewPosVert ) myViewPositionVert = maxViewPosVert;
    if( myViewPositionHorz > maxViewPosHorz ) myViewPositionHorz = maxViewPosHorz;

    myPrevViewPositionVert = myViewPositionVert;
    myPrevViewPositionHorz = myViewPositionHorz;
  }
  //-----------------------------------------------------
  //
  public double getSampleInt() {
    return mySampleInt;
  }
  public csSeisDispSettings getDispSettings() {
    return new csSeisDispSettings( mySettings );
  }
  public JFrame getParentFrame() {
    return myParentFrame;
  }
  public csSeisDispDialog getDispDialog() {
    return myDialogSettings;
  }
  public csSeisViewPopupMenu getPopupMenu() {
    return myPopupMenu;
  }
  public csISeismicTraceBuffer getTraceBuffer() {
    return myTraceBuffer;
  }
  //-----------------------------------------------------
  public int getViewPositionVert() {
    return myViewPositionVert;
  }
  public int getViewPositionHorz() {
    return myViewPositionHorz;
  }
  public int getMarginLeftRight() {
    return myMarginLeftRight;
  }
  /**
   * Set display setting dialog
   * This method enables overwriting the default display settings dialog with a customised one
   * @param dialog
   * @param useDialogDispSettings  true if display settings as set in dialog shall be used
   */
  public void setDialog( csSeisDispDialog dialog, boolean useDialogDispSettings ) {
    if( myDialogSettings != null ) {
      removeSeisViewListener( myDialogSettings );
    }
    myDialogSettings = dialog;
    if( useDialogDispSettings ) {
      dialog.apply();
    }
    else {
      myDialogSettings.changedSettings( mySettings );
    }
    addSeisViewListener( myDialogSettings );
  }
  /**
   * Set popup menu.
   * This method enables overwriting the default popup menu with a customised one
   * @param menu
   */
  public void setPopupMenu( csSeisViewPopupMenu menu ) {
    myPopupMenu = menu;
  }
  public csSampleInfo getSampleInfo( int xpos, int ypos ) {
    csSampleInfo info = new csSampleInfo();
    if( mySettings.plotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
      info.sampleDouble = this.yView2Model( ypos );
      info.traceDouble = this.xView2Trace( xpos );
    }
    else {
      info.sampleDouble = this.xView2ModelPlotDirHorz( xpos );
      info.traceDouble = this.yView2TracePlotDirHorz( ypos );
    }
    info.sample = (int)( info.sampleDouble + 0.5f );
    info.time = info.sampleDouble * mySampleInt * 0.001;
    info.trace = (int)( info.traceDouble + 0.5f );
    if( info.trace < 0 ) {
      info.trace = 0;
    }
    else if( info.trace > myNumTraces-1 ) {
      info.trace = myNumTraces-1;
    }
    if( info.sample >= 0 && info.sample <= myNumSamples-1 ) {
      float[] samples = myTraceBuffer.samples( info.trace );
      info.amplitude  = samples[info.sample];
    }
    return info;
  }
  /**
   * 
   * @return Seismic event handler of this view. See @class csSeismicEventHandler for more information
   */
  public csSeisViewEventHandler getEventHandler() {
    return myEventHandler;
  }
  //-----------------------------------------------------
  public void updateDispSettings( csSeisDispSettings ds ) {
    boolean resetSpline = (
        (ds.isVIDisplay != mySettings.isVIDisplay || ds.wiggleType != mySettings.wiggleType) &&
        (ds.isVIDisplay || ds.wiggleType == csSeisDispSettings.WIGGLE_TYPE_CUBIC) );
    boolean resetZoom = ( ds.zoomHorz != mySettings.zoomHorz || ds.zoomVert != mySettings.zoomVert || ds.isLogScale != mySettings.isLogScale );
    boolean resetScaling =
        ds.scaleType != mySettings.scaleType ||
        ds.traceScaling != mySettings.traceScaling ||
        ds.dispScalar != mySettings.dispScalar ||
        ds.minValue != mySettings.minValue ||
        ds.maxValue != mySettings.maxValue ||
        ds.fullTraceScalar != mySettings.fullTraceScalar;
//    boolean resetColor =
//        ds.wiggleColorPos != mySettings.wiggleColorPos ||
//        ds.wiggleColorNeg != mySettings.wiggleColorNeg ||
//        ds.highlightColor != mySettings.highlightColor;
    boolean resetTimeLines =
        ds.isTimeLinesAuto != mySettings.isTimeLinesAuto ||
        ds.showTimeLines != mySettings.showTimeLines ||
        ds.timeLineMinorInc != mySettings.timeLineMinorInc ||
        ds.timeLineMajorInc != mySettings.timeLineMajorInc ||
        ds.zoomVert != mySettings.zoomVert ||
        ds.isLogScale != mySettings.isLogScale;
    boolean resetScreenBuffer =
      ds.isVIDisplay != mySettings.isVIDisplay || ds.isLogScale != mySettings.isLogScale;

    if( resetSpline ) myRenderer.resetSplineFields( resetSpline );
    mySettings.set( ds );
    
    if( (resetZoom || resetSpline) && myTraceBuffer != null ) resetLeftRightMargins();
    if( (resetScreenBuffer || resetZoom) && myTraceBuffer != null ) {
      resetScreenBuffer();
      resetViewPositions();
    }
    if( resetScaling ) {
      resetTraceScaling();
    }
    if( resetTimeLines ) {
      resetTimeLines();
    }
//    if( resetZoom && ds.isLogScale ) {
//      zoom( ds.zoomVert, ds.zoomHorz );
//    }
    myIsFullRepaint = true;
    repaint();
    fireEventChangedSettings();
  }
  //-----------------------------------------------------
  //
  protected void resetTraceScaling() {
    myWiggleCentreAmplitude = 0.0f;
    if( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_SCALAR ) {
      for( int itrc = 0; itrc < myNumTraces; itrc++ ) {
        myTraceScalar[itrc] = 1.0f;
      }
      myDispScalar = mySettings.dispScalar;
      mySettings.viColorMap.setScalar( myDispScalar );
      mySettings.wiggleColorMap.setScalar( myDispScalar );
    }
    else if( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_TRACE ) {
      if( mySettings.traceScaling == csSeisDispSettings.TRACE_SCALING_AVERAGE ) {
        for( int itrc = 0; itrc < myNumTraces; itrc++ ) {
          float amp = myTraceBuffer.meanAmplitude(itrc);
          if( amp != 0.0f ) myTraceScalar[itrc] = 0.1f/amp;
          else myTraceScalar[itrc] = 1.0f;
        }
      }
      else {
        for( int itrc = 0; itrc < myNumTraces; itrc++ ) {
          float amp = Math.max( -myTraceBuffer.minAmplitude(itrc), myTraceBuffer.maxAmplitude(itrc) ) ;
          if( amp > EPSILON ) myTraceScalar[itrc] = 1.0f/amp;
          else myTraceScalar[itrc] = 1.0f;
        }
      }
      myDispScalar = mySettings.fullTraceScalar;
      mySettings.viColorMap.setScalar( myDispScalar );
      mySettings.wiggleColorMap.setScalar( myDispScalar );
    }
    else { // ( scaleType == csSeisDispSettings.SCALE_TYPE_RANGE ) {
      myWiggleCentreAmplitude = 0.5f * ( mySettings.minValue + mySettings.maxValue );
      float amp = Math.max( Math.abs( mySettings.minValue), Math.abs( mySettings.minValue ) );
      if( amp > EPSILON ) {
        myDispScalar = 1.0f / amp;
      }
      else {
        myDispScalar = 1.0f;
      }
      //System.out.println("RANGE " + mySettings.minValue + " " + mySettings.maxValue + " " + myWiggleCentreAmplitude + " Scalar:" + myDispScalar );
      mySettings.viColorMap.setMinMax( mySettings.minValue , mySettings.maxValue );
      mySettings.wiggleColorMap.setMinMax( mySettings.minValue , mySettings.maxValue );
    }
    myRenderer.resetTraceScaling( myDispScalar, myWiggleCentreAmplitude, myTraceScalar );
  }
  //--------------------------------------------------------------------------------------
  //
  private void resetTimeLines() {
    if( mySettings.isTimeLinesAuto ) {
      float pixelsPerMS = mySettings.zoomVert / (float)mySampleInt;
      float minorIncNew = DEFAULT_INC_PIXELS / pixelsPerMS;

      if( minorIncNew < getSampleInt() ) {
        mySettings.timeLineMinorInc = getSampleInt();
        mySettings.timeLineMajorInc = 10*mySettings.timeLineMinorInc;
      }
      else {
        double powerOf10  = Math.floor( Math.log10( minorIncNew ) );
        double incReduced = minorIncNew / Math.pow( 10, powerOf10 );

        for( int i = 0; i < DEFAULT_INC_MS.length; i++ ) {
          if( incReduced < DEFAULT_INC_MS[i] ) {
            mySettings.timeLineMinorInc = DEFAULT_INC_MS[i] * Math.pow( 10,  powerOf10 );
            mySettings.timeLineMajorInc = DEFAULT_MAJOR_MS[i] * Math.pow( 10,  powerOf10 );
            break;
          }
        }        
      }
    }
  }
  //-----------------------------------------------------------------
  //
  public float getScalar() {
    if( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_SCALAR ) {
      return mySettings.dispScalar;
    }
    else if( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_TRACE ) {
      return mySettings.fullTraceScalar;
    }
    else {
      return 1.0f;
    }
  }
  public void setScalar( float scalar) {
    boolean reset = false;
    if( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_SCALAR && scalar != mySettings.dispScalar ) {
      mySettings.dispScalar = scalar;
      reset = true;
    }
    else if( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_TRACE && scalar != mySettings.fullTraceScalar ) {
      mySettings.fullTraceScalar = scalar;
      reset = true;
    }
    if( reset ) {
      resetTraceScaling();
      if( myTraceBuffer != null ) {
        if( mySettings.isVIDisplay || mySettings.wiggleType == csSeisDispSettings.WIGGLE_TYPE_CUBIC ) myRenderer.resetSplineFields();
        myIsFullRepaint = true;
        repaint();
      }
      fireEventChangedSettings();
    }
  }
  public void zoom( float zoomVert, float zoomHorz ) {
    Rectangle rect = myRectListener.getVisibleRect();
    float sampleIndexCentre = yView2Model( rect.height/2 );
    float traceIndexCentre = xView2Trace( rect.width/2 );
    zoom( zoomVert, zoomHorz, traceIndexCentre, sampleIndexCentre );
  }
  public void zoom( float zoomVert, float zoomHorz, float traceIndexCentre, float sampleIndexCentre ) {
    Rectangle rect = myRectListener.getVisibleRect();
    if( zoomVert != mySettings.zoomVert ) {
      mySettings.zoomVert = zoomVert;

      int yViewCentreNew  = (int)(yModel2View( sampleIndexCentre ) + 0.5f) + myViewPositionVert;
//      if( mySettings.plotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        myViewPositionVert  = yViewCentreNew - rect.height/2;
//      }
//      else {
//        myViewPositionVert  = yViewCentreNew - rect.width/2;
//      }
      if( myViewPositionVert < 0 ) myViewPositionVert = 0;
      resetTimeLines();
    }
    if( zoomHorz != mySettings.zoomHorz ) {
      mySettings.zoomHorz = zoomHorz;
      resetLeftRightMargins();
      int xViewCentreNew = (int)(xModel2View( traceIndexCentre ) + 0.5f) + myViewPositionHorz;
//      if( mySettings.plotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        myViewPositionHorz = xViewCentreNew - rect.width/2;
//      }
//      else {
//        myViewPositionHorz = xViewCentreNew - rect.height/2;
//      }
      if( myViewPositionHorz < 0 ) myViewPositionHorz = 0;
    }
    if( myTraceBuffer != null ) {
      resetScreenBuffer();
      resetViewPositions();
      myIsFullRepaint = true;
      if( mySettings.isVIDisplay || mySettings.wiggleType == csSeisDispSettings.WIGGLE_TYPE_CUBIC ) myRenderer.resetSplineFields();
      repaint();  // Only need repaint when using DataBufferInt method
    }
    fireEventChangedSettings();
    fireEventSizeChanged();
  }
  public void zoomHorz( float zoomLevel ) {
    if( zoomLevel <= 0.0 ) return;
    zoom( mySettings.zoomVert, zoomLevel );
  }
  public void zoomVert( float zoomLevel ) {
    if( zoomLevel <= 0.0 ) return;
    zoom( zoomLevel, mySettings.zoomHorz );
  }
  public float getZoomVert() {
    return mySettings.zoomVert;
  }
  public float getZoomHorz() {
    return mySettings.zoomHorz;
  }
  public boolean isLogScale() {
    return mySettings.isLogScale;
  }
  public float getDispScalar() {
    return mySettings.dispScalar;
  }
//-------------------------------------------------
  protected void resetLeftRightMargins() {
    if( mySettings.isVIDisplay ) {
      myMarginLeftRight  = (int)( xModel2View( 0.5f ) - xModel2View( 0, 0 ) )-1;
    }
    else if( mySettings.doTraceClipping ) {
      myMarginLeftRight  = (int)( xModel2View( mySettings.traceClip ) - xModel2View( 0 ) );
    }
    else {
      myMarginLeftRight  = (int)( xModel2View( 0.5f ) - xModel2View( 0, 0 ) )-1;
    }
    myMarginLeftRight  += INSET_DEFAULT;
  }
  /**
   * Create screen buffer for seismic display
   * Call whenever zoom level changes, or traces are added or removed from seismic buffer
   * 
   * @return true if creation was successful, false if display is currently not visible (for example when panel is minimised).
   */
  protected boolean resetScreenBuffer() {
    //  Compute absolute max x and y view pixels. Add current scroll value and margin
    int xViewMax = (int)xModel2View( myNumTraces-1 ) + myMarginLeftRight + myViewPositionHorz;
    int yViewMax = (int)yModel2View( myNumSamples-1 ) + myMarginTopBottom + myViewPositionVert;

    if( mySettings.plotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
      int saveY = yViewMax;
      yViewMax = xViewMax;
      xViewMax = saveY;
    }
    setPreferredSize( new Dimension(xViewMax,yViewMax) );
    setMinimumSize( new Dimension(xViewMax,yViewMax) );
    setSize( new Dimension(xViewMax,yViewMax) );   // This will lead to a call to repaint()
    
    Rectangle rectVisible = myRectListener.getVisibleRect();
    myImageWidth  = rectVisible.width;
    myImageHeight = rectVisible.height;
    if( xViewMax < myImageWidth ) {
      myImageWidth = xViewMax;
    }

    if( myImageHeight == 0 || myImageWidth == 0 ) {
        //      System.out.println("!!!!!!!!!!!!!!!!!!!!!Image width NULL" + myImageWidth + " " +  myImageHeight);
      return false;
    }

    //    System.out.println("resetScreen: " + myTEMPName + " Creating volatile image " + myImageWidth + " " +  myImageHeight);
    myVolatileBitmap = createImage( myImageWidth, myImageHeight );
    if( myColorBitType == csColorMap.COLOR_MAP_TYPE_32BIT ) {
      myBitmap = new BufferedImage( myImageWidth, myImageHeight, BufferedImage.TYPE_INT_RGB );
    }
    else {
      myBitmap = new BufferedImage( myImageWidth, myImageHeight, BufferedImage.TYPE_BYTE_INDEXED );
    }
    myRenderer.resetScreenBuffer( myImageWidth, myImageHeight, myBitmap.getRaster().getDataBuffer() );

    if( mySettings.isLogScale ) {
      float maxView    = (myTraceBuffer.numSamples()*mySettings.zoomVert);
      float maxLogView = (float)Math.log10(((myTraceBuffer.numSamples())*mySettings.zoomVert));
      myLogScaleRatio = maxView/maxLogView;
    }
    return true;
  }
  /**
   * Convert model value/sample index to pixel value on screen
   * 
   * @param sampleIndex Sample index, starting with 0
   * @return View Y coordinate, relative to current position of view in terms of scrolling.
   *         For example, a negative yView means y coordinate is outside of visible view area ('above' view area)
   */
  public float yModel2View( float sampleIndex ) {
    if( !mySettings.isLogScale ) {
      float ypixView = (sampleIndex*mySettings.zoomVert) - myViewPositionVert + myMarginTopBottom;
      return ypixView;
    }
    else {
      float ypixView = myLogScaleRatio*(float)Math.log10(((sampleIndex+1)*mySettings.zoomVert)) - myViewPositionVert + myMarginTopBottom;
//      System.out.println("yModel2View: " + ypixView + " " + sampleIndex + "  ---  " + maxView + " " + maxLogView);
      return ypixView;
    }
  }
  /**
   * Convert model value/trace index to pixel value on screen
   * 
   * @param amplitude Amplitude of trace sample
   * @param traceIndex Trace index
   * @return View X coordinate, relative to current position of view in terms of scrolling.
   *         For example, a negative xView means x coordinate is outside of visible view area ('left' of view area)
   */
  public float xModel2View( float amplitude, int traceIndex ) {
    float xView = (int)(traceIndex*mySettings.zoomHorz+0.5f) - myViewPositionHorz + myMarginLeftRight +
          amplitude*myDispScalar*myTraceScalar[traceIndex]*mySettings.zoomHorz;
    return xView;
  }
  public float xModel2View( float traceIndex ) {
    float xView = (int)(traceIndex*mySettings.zoomHorz+0.5f) - myViewPositionHorz + myMarginLeftRight;
    return xView;
  }
  /**
   * Convert pixel value to model value. Convert y pixel to sample index.
   * 
   * @param yView  Absolute view coordinate, independent of scroll value (as if full seismic section was visible)
   * @return
   */
  public float yView2Model( int yView ) {
    if( !mySettings.isLogScale ) {
      float yModel = (float)(yView + myViewPositionVert - myMarginTopBottom)/(float)mySettings.zoomVert;
      return yModel;
    }
    else {
      float yModel = (float)(Math.pow(10,(yView + myViewPositionVert - myMarginTopBottom)/myLogScaleRatio)/(float)mySettings.zoomVert) - 1;
//      System.out.println("yView2Model: " + yView + " " + yModel  + "  ---  " + yModel2View(yModel) );
      return yModel;
    }
  }
  /**
   * Convert pixel value (view x coordinate) to trace index.
   * 
   * @param xView Absolute view x coordinate, which is independent of the current scroll value (as if full seismic section was visible).
   * @return trace index at this view coordinate.
   */
  public float xView2Trace( int xView ) {
    return ((float)(xView + myViewPositionHorz - myMarginLeftRight)/(float)mySettings.zoomHorz);
  }
  public float yView2TracePlotDirHorz( int yView ) {
    return ((float)(myImageWidth - yView + myViewPositionHorz - myMarginLeftRight)/(float)mySettings.zoomHorz);
  }
  public float xView2ModelPlotDirHorz( int xView ) {
    float yModel;
    yModel = (float)(xView + myViewPositionVert - myMarginTopBottom)/(float)mySettings.zoomVert;
    return yModel;
  }
  //************************************************************************************
  //************************************************************************************
  //************************************************************************************
  //
  // Paint methods
  //
  /**
   * System repaint method, entry point for rendering the JPanel.
   */
  public void paintComponent( Graphics g ) {
    if( myBitmap == null || myTraceBuffer == null ) return;
    // Take snapshot of 'myIsScrolling' attribute, in case it changes during the course of the current repaint process
    // Maybe the myIsScrolling attribute should be protected by synchronized get() and set() methods..?
    boolean snapShot_isScrolling = myIsScrolling;

//    System.out.printf("REPAINT  Width/Height: %6d/%6d   Visible width/height: %6d/%6d   Scroll pos vert/horz: %6d/%6d  Next scroll: %6d/%6d   %b/%b\n",
//        myImageWidth, myImageHeight, getVisibleRect().width, getVisibleRect().height, myViewPositionVert, myViewPositionHorz, myNextViewPositionVert, myNextViewPositionHorz, snapShot_isScrolling, myIsScrolling );
    
    if( myIsScrolling && myShowCrosshair && myIsCrosshairPainted ) {
      removeCrosshair();
    }
    repaintStep1( snapShot_isScrolling );

    Graphics2D g2 = (Graphics2D)g;
    if( mySettings.plotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
      g2.drawImage( myVolatileBitmap, 0, 0, this );
    }
    else {
      g2.rotate(-Math.PI/2);
      g2.translate(-myImageWidth, 0);
      g2.drawImage( myVolatileBitmap, 0, 0, this );
      g2.translate(myImageWidth, 0);
      g2.rotate(Math.PI/2);
    }

    for( int ioverlay = 0; ioverlay < mySeisOverlays.size(); ioverlay++ ) {
      csISeisOverlay ol = mySeisOverlays.get(ioverlay);
      ol.draw( this, g2 );
    }

    // If the view is currently scrolling, we need to update view positions, and check if there is another scroll action waiting ('next' view position)
    if( snapShot_isScrolling ) {
      myPrevViewPositionVert = myViewPositionVert;
      myPrevViewPositionHorz = myViewPositionHorz;
      myIsScrolling = false;
      boolean isVertScroll = myNextViewPositionVert >= 0 && myNextViewPositionVert != myViewPositionVert;
      boolean isHorzScroll = myNextViewPositionHorz >= 0 && myNextViewPositionHorz != myViewPositionHorz;
      if( mySwitchScroll ) {
        mySwitchScroll = false;
        if( isVertScroll ) {
          resetViewPositionVert( myNextViewPositionVert );
        }
        else if( isHorzScroll ) {
          resetViewPositionHorz( myNextViewPositionHorz );
        }
      }
      // Perform potential scroll operations in different order each time
      else {
        mySwitchScroll = true;
        if(isHorzScroll) {
          resetViewPositionHorz( myNextViewPositionHorz );
        }
        else if (isVertScroll) {
          resetViewPositionVert( myNextViewPositionVert );
        }
      }
    }
 //    System.out.printf("                       %6d/%6d                         %6d/%6d                         %6d/%6d               %6d/%6d  %b/%b\n",
//        myImageWidth, myImageHeight, getVisibleRect().width, getVisibleRect().height, myViewPositionVert, myViewPositionHorz, myNextViewPositionVert, myNextViewPositionHorz, snapShot_isScrolling, myIsScrolling );
  }
  /**
   * Repaint method, step 1
   * Computes the actual area that needs repainting, taking into account the current scroll
   * Passes on repainting to step2, with given min/max pixels values to repaint
   * 
   * @param isScrolling Pass true if scrolling is currently taking place.
   * @return true if repaint action was successful, false if not
   */
  protected boolean repaintStep1( boolean isScrolling ) {
    Rectangle rectVisible = myRectListener.getVisibleRect();
    //    if( myTEMPName == "NAME_BACKGROUND" )     System.out.println("-------------Paint step 1 " + myImageHeight + " =? " + rectVisible.height + ", " + myImageWidth + " =? " + rectVisible.width );
//    long time1 = System.currentTimeMillis();   // Comment out to compute time paint operation
    
    if( myBitmap == null || myImageWidth != rectVisible.width || myImageHeight != rectVisible.height ) {
      if( !resetScreenBuffer() ) return false;      
      myIsFullRepaint = true;
    }
    if( rectVisible.width == 0 && rectVisible.height == 0 ) {
      return false;
    }
    
    int xdiff = myViewPositionHorz - myPrevViewPositionHorz;
    int ydiff = myViewPositionVert - myPrevViewPositionVert;
    int xdiffAbs = Math.abs(xdiff);
    int ydiffAbs = Math.abs(ydiff);
    int copyWidth  = myImageWidth-xdiffAbs;
    int copyHeight = myImageHeight-ydiffAbs;
    if( ydiffAbs >= rectVisible.height || xdiffAbs >= rectVisible.width ) {
      myIsFullRepaint = true;
    }
    if( myIsFullRepaint || copyWidth < 10 || copyHeight < 10 ) {
      myIsFullRepaint = false;
      myRenderer.repaintStep2( myMarginTopBottom, rectVisible.height-1, 0, rectVisible.width-1 );
      if( myVolatileBitmap != null ) {
        Graphics2D g = (Graphics2D)myVolatileBitmap.getGraphics();
        // The following line of code transfers all pixels to the managed buffer myVolatile
        // This can be extremely slow when working in a client/server mode, even on a fast network
        g.drawImage( myBitmap, 0, 0, this );
        g.dispose();
      }
    }
    else if( isScrolling ) {
      if( xdiff != 0 ) {
        int xfrom = xdiff;
        int xto   = 0;
        int xnew  = copyWidth;
        if( xdiff < 0 ) {
          xfrom = 0;
          xto   = -xdiff;
          xnew  = 0;
        }
        myRenderer.copyBufferHorz( copyWidth, myImageHeight, xfrom, xto );
        myRenderer.repaintStep2( myMarginTopBottom, rectVisible.height-1, xnew, xnew+xdiffAbs-1 );
        int xmin = Math.max( xnew-1, 0 );
        int xmax = xnew+xdiffAbs;
        /*
        if( this.mySettings.isVADisplay ) {
          xmin = Math.max( xmin-3*(int)mySettings.zoomHorz, 0 );
          xmax = Math.min( xmax+3*(int)mySettings.zoomHorz, myImageWidth );
        }
         */
        Graphics2D g = (Graphics2D)myVolatileBitmap.getGraphics();
        g.copyArea( xfrom, 0, copyWidth, rectVisible.height, xto-xfrom, 0 );
        g.drawImage( myBitmap, xmin, 0, xmax, rectVisible.height, xmin, 0, xmax, rectVisible.height, this );
        g.dispose();
      }
      if( ydiff != 0 ) {
        int yfrom = ydiff;
        int yto   = 0;
        int ynew  = copyHeight;
        if( ydiff < 0 ) {
          yfrom = 0;
          yto   = -ydiff;
          ynew  = 0;
        }
        myRenderer.copyBufferVert( myImageWidth, copyHeight, yfrom, yto );
        myRenderer.repaintStep2( ynew, ynew+ydiffAbs-1, 0, rectVisible.width-1 );

        // Only repaint pixels that have changed
        int ymin = Math.max( ynew-1, 0 );
        int ymax = ynew+ydiffAbs;
        Graphics2D g = (Graphics2D)myVolatileBitmap.getGraphics();
        g.copyArea( 0, yfrom, rectVisible.width, copyHeight, 0, yto-yfrom );
        g.drawImage( myBitmap, 0, ymin, rectVisible.width, ymax, 0, ymin, rectVisible.width, ymax, this );
        g.dispose();
      }
      // Full repaint:
      //      Graphics2D g = (Graphics2D)myVolatile.getGraphics();
      //      g.drawImage( myBitmap, 0, 0, this );

    }
    else {
      return false;
    }

    //System.out.println("Total time in milliseconds: " + (System.currentTimeMillis()-time1) );
    return true;
  }
  //************************************************************************************
  //************************************************************************************
  //************************************************************************************
  //
  // Methods related to scrolling
  //
  /**
   * Call when vertical scroll value has been reset.
   * Only repaint if scroll value really has changed, and when there is currently no repainting in progress.
   * 
   * @param scrollValue New vertical scroll value
   */
  protected synchronized void resetViewPositionVert( int scrollValue ) {
//    System.out.println("--- SCROLL VERT " + myViewPositionVert +  " " + scrollValue + " " + myIsScrolling);
    // myNextViewPositionVert = scrollValue;  // BUGFIX 26/12/08: Comment out, added else if/else clause. This lead to wrong scroll value when seismic was smaller than visible area
    if( scrollValue != myViewPositionVert && !myIsScrolling ) {
      myNextViewPositionVert = -1;
      myViewPositionVert = scrollValue;
      myIsScrolling     = true;
      repaint();
      fireEventVertScrollChanged(scrollValue);
    }
    else if( myIsScrolling ) {
      myNextViewPositionVert = scrollValue;
    }
    else {
      myNextViewPositionVert = -1;
    }
  }
  /**
   * Call when horizontal scroll value has been reset.
   * Only repaint if scroll value really has changed, and when there is currently no repainting in progress.
   * 
   * @param scrollValue New horizontal scroll value
   */
  protected synchronized void resetViewPositionHorz( int scrollValue ) {
//    System.out.println("--- SCROLL HORZ " + myViewPositionHorz +  " " + scrollValue + " " + myIsScrolling);
    // myNextViewPositionHorz = scrollValue;  // BUGFIX 26/12/08: Comment out, added else if/else clause. This lead to wrong scroll value when seismic was smaller than visible area
    if( scrollValue != myViewPositionHorz && !myIsScrolling ) {
      myNextViewPositionHorz = -1;
      myViewPositionHorz    = scrollValue;
      myIsScrolling   = true;
      // If there is no trace clipping, everything has to be repainted, because no prediction can be made regarding the number of traces to repaint
      if( !myIsFullRepaint ) myIsFullRepaint = !mySettings.doTraceClipping && mySettings.showWiggle;
      repaint();
      fireEventHorzScrollChanged(scrollValue);
    }
    else if( myIsScrolling ) {
      myNextViewPositionHorz = scrollValue;
    }
    else {
      myNextViewPositionHorz = -1;
    }
  }
  //************************************************************************************
  //************************************************************************************
  //************************************************************************************
  //
  // Methods related to special listeners
  //
  public void addOverlay( csISeisOverlay ol ) {
    for( int i = 0; i < mySeisOverlays.size(); i++ ) {
      if( mySeisOverlays.get(i).equals(ol) ) {
        return;
      }
    }
    mySeisOverlays.add( ol );
  }
  public void removeOverlay( csISeisOverlay ol ) {
    for( int i = 0; i < mySeisOverlays.size(); i++ ) {
      if( mySeisOverlays.get(i).equals(ol) ) {
        mySeisOverlays.remove(i);
        break;
      }
    }
  }
  public csISeisOverlay getOverlay( int index ) {
    if( index >= 0 && index < mySeisOverlays.size() ) {
      return mySeisOverlays.get(index);
    }
    return null;
  }
  public int getNumOverlays() {
    return mySeisOverlays.size();
  }
  public void removeAllOverlays() {
    mySeisOverlays.clear();
  }
//-----------------------------------------------------
  public void addSeisViewListener( csISeisViewListener listener ) {
    for( int i = 0; i < mySeisViewListeners.size(); i++ ) {
      if( mySeisViewListeners.get(i).equals(listener) ) {
        return;
      }
    }
    mySeisViewListeners.add( listener );
  }
  public void removeSeisViewListener( csISeisViewListener listener ) {
    for( int i = 0; i < mySeisViewListeners.size(); i++ ) {
      if( mySeisViewListeners.get(i).equals(listener) ) {
        mySeisViewListeners.remove(i);
        break;
      }
    }
  }
  public void fireEventChangedTraceBuffer( csISeismicTraceBuffer traceBuffer ) {
    for( int i = 0; i < mySeisViewListeners.size(); i++ ) {
      mySeisViewListeners.get(i).traceBufferChanged(traceBuffer);
    }
  }
  public void fireEventChangedSettings() {
    for( int i = 0; i < mySeisViewListeners.size(); i++ ) {
      mySeisViewListeners.get(i).changedSettings( mySettings );
    }
  }
  public void fireEventVertScrollChanged( int scrollValue ) {
    for( int i = 0; i < mySeisViewListeners.size(); i++ ) {
      mySeisViewListeners.get(i).vertScrollChanged( scrollValue );
    }
  }
  public void fireEventHorzScrollChanged( int scrollValue ) {
    for( int i = 0; i < mySeisViewListeners.size(); i++ ) {
      mySeisViewListeners.get(i).horzScrollChanged( scrollValue );
    }
  }
  public void fireEventSizeChanged() {
    for( int i = 0; i < mySeisViewListeners.size(); i++ ) {
      mySeisViewListeners.get(i).sizeChanged( getPreferredSize() );
    }
  }
  public void fireEventMouseExited() {
    for( int i = 0; i < mySeisViewListeners.size(); i++ ) {
      mySeisViewListeners.get(i).mouseExited();
    }
  }
  //************************************************************************************
  //************************************************************************************
  //************************************************************************************
  //
  // Methods related to trace highlighting and crosshair
  //
  public void mouseMoved( csSampleInfo info ) {
    if( myIsHighlightOn ) {
      if( info.trace != myHighlightTrace ) {
        if( myShowCrosshair ) removeCrosshair();
        removeHighlightTrace();
        highlightTrace( info.trace );
      }
    }
    if( myShowCrosshair ) {
      paintCrosshair( info );
    }
  }
  public void mouseExited() {
    if( myShowCrosshair ) {
      removeCrosshair();
      fireEventMouseExited();
    }
    myIsCrosshairPainted = false;
  }
  public void setCrosshairPosition( csSampleInfo info ) {
    paintCrosshair( info );
  }
  public void removeCrosshair() {
    if( myCrosshairPosition.x >= 0 || myCrosshairPosition.y >= 0 ) {
      paintCrosshair2( myCrosshairPosition.x, myCrosshairPosition.y );
      myIsCrosshairPainted = false;
      myCrosshairPosition.x = -1;
      myCrosshairPosition.y = -1;
      repaint();  // Only need repaint when using DataBufferInt method
    }
  }
  private void paintCrosshair( csSampleInfo info ) {
    boolean refresh = false;
    int xView;
    int yView;
    if( !mySnapCrosshair ) {
      xView = (int)( xModel2View((float)info.traceDouble) + 0.5f );
      yView = (int)( yModel2View((float)info.sampleDouble) + 0.5f );
    }
    else {
      xView = (int)( xModel2View((float)info.trace) + 0.5f );
      yView = (int)( yModel2View((float)info.sample) + 0.5f );
    }
    if( myCrosshairPosition.x != xView || myCrosshairPosition.y != yView ) {
      if( myIsCrosshairPainted ) {
        paintCrosshair2( myCrosshairPosition.x, myCrosshairPosition.y );
        refresh = true;
      }
      myIsCrosshairPainted = false;
    }
    myCrosshairPosition.x = xView;
    myCrosshairPosition.y = yView;
    if( !myIsCrosshairPainted ) {
      paintCrosshair2( myCrosshairPosition.x, myCrosshairPosition.y );
      refresh = true;
      myIsCrosshairPainted = true;
    }
    if( refresh ) {
      repaint();  // Only need repaint when using DataBufferInt method
    }
  }
  /**
   * Paint crosshair on top of seismic section.
   * 
   * @param xPos X position of crosshair
   * @param yPos Y position of crosshair
   */
  private void paintCrosshair2( int xPos, int yPos ) {
    Graphics2D g2 = (Graphics2D)myVolatileBitmap.getGraphics();
    // XOR current color with 0x00aaaaff instead of 0x00ffffff to shift crosshair color towards yellow.
    // This makes sure crosshair is seen even with a gray background
    g2.setXORMode( new Color( 0x00aaaaff ) );
    g2.setStroke( new BasicStroke(1) );
    g2.drawLine( 0, yPos, myImageWidth, yPos );
    g2.drawLine( xPos, 0, xPos, myImageHeight );
    g2.dispose();
  }
  public void setAutoHighlightTrace( boolean doHighlight ) {
    if( myIsHighlightOn && !doHighlight ) {
      removeHighlightTrace();
      myIsHighlightOn = doHighlight;
      myIsFullRepaint = true;
      repaint();
    }
    else {
      myIsHighlightOn = doHighlight;
    }
  }
  public void setCrosshair( boolean showCrosshair, boolean snapCrosshair ) {
    myShowCrosshair = showCrosshair;
    mySnapCrosshair = snapCrosshair;
  }
  /**
   * Highlight trace
   * @param traceIndex Trace index (starting with 0 for first trace)
   */
  public boolean isHighlightOn() {
    return myIsHighlightOn;
  }
  public int getHighlightTrace() {
    return myHighlightTrace;
  }
  public void highlightTrace( int traceIndex ) {
    highlightTrace( traceIndex, true );
  }
  private void highlightTrace( int traceIndex, boolean doHighlightTrace ) {
    if( myBitmap == null ) return;
    Rectangle rectVisible = myRectListener.getVisibleRect();

    if( traceIndex >= 0 && traceIndex < myNumTraces ) {
      myStepTrace = 1;
      if( mySettings.zoomHorz < 1.0 ) {
        myStepTrace = (int)(1.0 / mySettings.zoomHorz);
      }
      int xMin = 0;
      int xMax = 0;
      if( mySettings.isVIDisplay ) {
        xMin = (int)( xModel2View( traceIndex-myStepTrace*0.5f ) + 0.5f ) + 1;
        xMax = (int)( xModel2View( traceIndex+myStepTrace*0.5f ) + 0.5f );
        if( xMax >= myRectListener.getVisibleRect().width )  xMax = myRectListener.getVisibleRect().width-1;
        if( xMin < 0 ) xMin = 0;
        else if( xMax < xMin ) {
          xMin = xMax;
        }
        if( xMax >= 0 ) {
          myRenderer.repaintVADiscrete( 0, rectVisible.height-1, xMin, xMax, true );
          // Repaint time lines if necessary
          if( mySettings.showTimeLines ) myRenderer.repaintLines( 0, rectVisible.height-1, xMin, xMax );
        }
      }
      float traceClipPixelsPrev = xModel2View( myStepTrace*mySettings.traceClip ) - xModel2View( 0 );
      float xtracePrev = xModel2View( myHighlightTrace );
      myHighlightTrace = traceIndex;
      // Paint new highlighted trace
      if( mySettings.showWiggle ) {
        myRenderer.repaintWiggleTrace( 0, rectVisible.height-1, myHighlightTrace, myHighlightTrace, mySettings.highlightColor );
      }
      // Full repaint:
      //      Graphics2D g = (Graphics2D)myVolatile.getGraphics();
      //      g.drawImage( myBitmap, 0, 0, this );

      // Only repaint pixels that have changed
      Graphics2D g = (Graphics2D)myVolatileBitmap.getGraphics();
      float traceClipPixels = xModel2View( myStepTrace*mySettings.traceClip ) - xModel2View( 0 );
      float xtrace = xModel2View( myHighlightTrace );
      int xpos1Prev = Math.max( (int)( xtracePrev - traceClipPixelsPrev ) - 1, 0 );
      int xpos2Prev = Math.min( (int)( xtracePrev + traceClipPixelsPrev + 0.5f ) + 1, myImageWidth-1 );
      int xpos1 = Math.max( (int)( xtrace - traceClipPixels ) - 1, 0 );
      int xpos2 = Math.min( (int)( xtrace + traceClipPixels + 0.5f ) + 1, myImageWidth-1 );
      if( xpos2 < xpos1 ) xpos2 = xpos1;
      if( xpos2Prev < xpos1Prev ) xpos2Prev = xpos1Prev;
      if( xpos1 < xpos1Prev ) {
        if( xpos2 < xpos1Prev ) {
          g.drawImage( myBitmap, xpos1, 0, xpos2, myImageHeight-1, xpos1, 0, xpos2, myImageHeight-1, this );
          g.drawImage( myBitmap, xpos1Prev, 0, xpos2Prev, myImageHeight-1, xpos1Prev, 0, xpos2Prev, myImageHeight-1, this );
        }
        else {
          g.drawImage( myBitmap, xpos1, 0, xpos2Prev, myImageHeight-1, xpos1, 0, xpos2Prev, myImageHeight-1, this );
        }
      }
      else {
        if( xpos2Prev < xpos1 ) {
          g.drawImage( myBitmap, xpos1, 0, xpos2, myImageHeight-1, xpos1, 0, xpos2, myImageHeight-1, this );
          g.drawImage( myBitmap, xpos1Prev, 0, xpos2Prev, myImageHeight-1, xpos1Prev, 0, xpos2Prev, myImageHeight-1, this );
        }
        else {
          g.drawImage( myBitmap, xpos1Prev, 0, xpos2, myImageHeight-1, xpos1Prev, 0, xpos2, myImageHeight-1, this );
        }
      }
      g.dispose();
      
      repaint();  // Only need repaint when using DataBufferInt method
    }
  }
  public void removeHighlightTrace() {
    if( myBitmap == null ) return;
    Rectangle rectVisible = myRectListener.getVisibleRect();

    // Paint over previous highlighted trace
    int traceLeft  = myHighlightTrace;
    int traceRight = myHighlightTrace;
    myStepTrace = 1;
    if( mySettings.zoomHorz < 1.0 ) {
      myStepTrace = (int)(1.0 / mySettings.zoomHorz);
    }
    if( mySettings.isVIDisplay ) {
      int xMin = (int)( xModel2View( myHighlightTrace-myStepTrace*0.5f ) + 0.5f ) + 1;
      int xMax = (int)( xModel2View( myHighlightTrace+myStepTrace*0.5f ) + 0.5f );
      if( xMax >= myRectListener.getVisibleRect().width ) xMax = myRectListener.getVisibleRect().width-1;
      if( xMin < 0 ) {
        xMin = 0;
      }
      else if( xMax < xMin ) {
        xMin = xMax;
      }
      if( xMax >= 0 ) {
        switch( mySettings.viType ) {
          case csSeisDispSettings.VA_TYPE_DISCRETE:
            myRenderer.repaintVADiscrete( 0, rectVisible.height-1, xMin, xMax, false );
            break;
          case csSeisDispSettings.VA_TYPE_VERTICAL:
            myRenderer.repaintVAVertical( 0, rectVisible.height-1, xMin, xMax );
            break;
          case csSeisDispSettings.VA_TYPE_2DSPLINE:
            myRenderer.repaintVA2DSpline( 0, rectVisible.height-1, xMin, xMax );
            break;
        }
        // Repaint time lines if necessary
        if( mySettings.showTimeLines ) myRenderer.repaintLines( 0, rectVisible.height-1, xMin, xMax );
      }
      traceLeft = Math.max( (int)(traceLeft-mySettings.traceClip), 0 );
      traceRight = Math.min( (int)(traceRight+mySettings.traceClip), myNumTraces-1 );
    }
    if( mySettings.showWiggle || mySettings.isPosFill || mySettings.isNegFill ) {
      myRenderer.repaintWiggleTrace( 0, rectVisible.height-1, traceLeft, traceRight, Color.black );
    }
  }
  public void test( int iter ) {
  }
  /**
   * Override for case of horizontal plot direction
   */
  @Override
  public Rectangle getVisibleRect() {
    Rectangle rect = super.getVisibleRect();
    if( mySettings.plotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
      int saveHeight = rect.height;
      rect.height = rect.width;
      rect.width = saveHeight;
    }
    return rect;
  }

  //***********************************************************************************
  //***********************************************************************************
  //***********************************************************************************
  private class SeisViewKeyDispatcher implements KeyEventDispatcher {
    private KeyListener myListener;
    public SeisViewKeyDispatcher( KeyListener listener ) {
      super();
      myListener = listener;
    }
    @Override
    public boolean dispatchKeyEvent(KeyEvent e) {
      if( e.getID() == KeyEvent.KEY_PRESSED) {
        myListener.keyPressed( e );
      }
      else if( e.getID() == KeyEvent.KEY_RELEASED ) {
        myListener.keyReleased( e );
      }
      else if( e.getID() == KeyEvent.KEY_TYPED ) {
        myListener.keyTyped( e );
      }
      return false;
    }
  }

}


