/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import cseis.seis.*;
import cseis.general.csUnits;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import javax.swing.*;

import java.awt.event.*;
import java.awt.image.BufferedImage;

/**
 * Seismic scroll pane.<br>
 * <br>
 * A seismic display consists of mainly two Java objects:<br>
 * <ul>
 *  <li> The seismic display itself (osSeisView) and the scroll pane in which the display resides (osSeisPane).
 *  <li> csSeisPane encapsulates some of csSeisView's functionality and acts as the interface to the outside world.
 *  <li> However, the encapsulation is not perfect, so in some instances the methods of the csSeisView object have to be called directly.
 * </ul>
 * <br>
 * The csSeisPane scroll pane provides:<br>
 *  <ul>
 *  <li> vertical and horizontal scroll bars,
 *  <li>  vertical annotation area (side label), showing times,
 *  <li>  horizontal annotation area (side label), showing trace headers, and
 *  <li>  zoom functionality based on mouse movements inside side label areas.
 * </ul>
 * <br>
 * Note: JScrollPane was avoided because of some bugs/inconsistent features.
 * <br>
 * Usage:<br>
 *  To construct a full seismic display...<br>
 * <ol>
 *  <li> Create seismic view<br>
 *     csSeisView view = new csSeisView( parent_JFrame );
 *  <li> Construct seismic pane<br>
 *     csSeisPane pane = new csSeisPane( view );
 *  <li> Create trace header definitions<br>
 *     csHeader traceHeaders = new csHeader[numHeaders];<br>
 *     for( ...all trace headers... ) {<br>
 *       traceHeaders[i] = new csHeader( headerName, headerDesc, headerType );<br>
 *     }
 *  <li> Load seismic data, populate trace buffer<br>
 *     csTraceBuffer traceBuffer = new csTraceBuffer( numSamples, numHeaders );<br>
 *     for( ...all traces... ) {<br>
 *       float[] sampleValues  = new float[numSamples];<br>
 *       Number[] headerValues = new Number[numHeaders];<br>
 *       // ...Set sample and header values...<br>
 *       traceBuffer.addTrace( new csSeismicTrace( sampleValues, headerValues ) );<br>
 *     }
 * <li> Finally, update seismic pane with populated trace buffer<br>
 *     pane.updateSeismic( traceBuffer, sampleInt, traceHeaders );
 * </ol>
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csSeisPane extends JPanel implements csISeisViewListener, csIKeyListener, csIPanningListener {
  public static final int DEFAULT_WIDTH_SIDELABEL  = 64;
  public static final int DEFAULT_HEIGHT_SIDELABEL = 25;
  public static final int MIN_ZOOM_WIDTH = 3;
  
  public static final int ZOOM_VERT = 1;
  public static final int ZOOM_HORZ = 2;
  public static final int ZOOM_BOTH = ZOOM_VERT | ZOOM_HORZ;

  public static final int ZOOM_IN  = 11;
  public static final int ZOOM_OUT = 12;

  public static final int PRESSED_NONE = 0;
  public static final int PRESSED_START = 1;
  public static final int PRESSED_END = 2;
  public static final int PRESSED_BOTH = 3;

  public static final int DEFAULT_UNIT_INCREMENT = 4;
  public static final Color DEFAULT_COLOR = new JPanel().getBackground();
  
  private csHeaderDef[] myTraceHeaders;
  private double mySampleInt;
  private int myNumTraces;
  private int myNumSamples;
  private JScrollBar myScrollBarVert;
  private JScrollBar myScrollBarHorz;
  private csCornerHorzLabel myCornerUpperLeft;
  private csISeismicTraceBuffer myTraceBuffer = null;

  private csSeisView mySeisView;
  private int myPreviousScrollValueVert;
  private int myPreviousScrollValueHorz;

  private boolean myStopScrolling = false;
//  private int[] myTraceIndexList = null;
  // Side labels

  // Time lines & side labels
  /// Vertical side label: Time axis
  protected csSideLabelVert mySideLabelVert;
  /// Horizontal side label: Trace axis
  protected csSideLabelHorz mySideLabelHorz;
  protected JPanel myPanelTopLabel;
  /// Font size for side label annotation
  public static final int DEFAULT_LABEL_FONTSIZE = 10;
  /// Show sequential trace header value in annotation area
  protected boolean myShowSeqTraceNum;
  protected boolean myOmitRepeatedValues;
//  private boolean myPaintLines = true;
  private int myTraceLabelStep = 1;
  private boolean myIsTraceLabelStepAuto = true;
  private boolean myIsTraceStepInitialised = false;
  private Dimension mySeisViewSize = null;
  private int myPlotDirection;
  private int myFirstSeqTraceIndex = 0;
  
  private int myVerticalDomain;  // csUnit::DOMAIN_TIME, csUnit::DOMAIN_DEPTH or csUnit::DOMAIN_FREQ
  private int mySampleIntScalar;
  private String myTitleVerticalAxis;

  //--------------------
  private JFrame myParentFrame = null;
  private csAnnotationDialog myAnnotationDialog = null;
//  private boolean myPanSwitch;
  private JPanel myPanelViewPort;

//---------------------------------------------------------------------------
//
  /**
   * Constructor
   * Use this constructor seismic traces are already available before constructing the seismic pane.
   * 
   * @param seisView  Seismic view that shall be placed into seismic pane
   * @param buffer    Trace buffer containing all seismic traces
   * @param sampleInt Sample interval [ms]
   * @param headers   Header definitions for all trace headers
   *                   (size of array should be the same as number of trace headers in seismic trace buffer)
   */
  public csSeisPane( csSeisView seisView, csISeismicTraceBuffer buffer, double sampleInt, csHeaderDef[] headers ) {
    this( seisView );
    updateSeismic( buffer, sampleInt, headers );
  }
  /**
   * Constructor.
   * Use this constructor if seismic traces shall be added at a later stage.
   * 
   * @param seisView  Seismic view that shall be placed into seismic pane
   */
  public csSeisPane( csSeisView seisView ) {
    this( seisView, new csAnnotationAttributes() );
  }
  public csSeisPane( csSeisView seisView, csAnnotationAttributes attr ) {
    super( new BorderLayout() );

//    myPanSwitch = false;
    //    myTraceIndexList = null;
    myTitleVerticalAxis = "Time [seconds]";
    myVerticalDomain    = csUnits.DOMAIN_TIME;
    mySampleIntScalar   = 1000;
    
    myShowSeqTraceNum    = attr.showSequential;
    myOmitRepeatedValues = attr.omitRepeating;
    myTraceLabelStep     = attr.traceLabelStep;
    myIsTraceLabelStepAuto = !attr.fixedTraceLabelStep;
    myIsTraceStepInitialised = false;
    
    myPlotDirection = csSeisDispSettings.PLOT_DIR_VERTICAL;
    mySampleInt = 2.0;
    
    myPreviousScrollValueVert = 0;
    myPreviousScrollValueHorz = 0;

    myScrollBarVert = new JScrollBar( JScrollBar.VERTICAL );
    myScrollBarHorz = new JScrollBar( JScrollBar.HORIZONTAL );
    myScrollBarVert.setValues( 0, 50, 0, 0 );
    myScrollBarHorz.setValues( 0, 50, 0, 0 );
    myScrollBarVert.setBlockIncrement(200);
    myScrollBarHorz.setBlockIncrement(200);
    myScrollBarVert.setUnitIncrement(DEFAULT_UNIT_INCREMENT);
    myScrollBarHorz.setUnitIncrement(DEFAULT_UNIT_INCREMENT);
    myScrollBarVert.setVisibleAmount( 200 );
    myScrollBarHorz.setVisibleAmount( 200 );
    mySeisView = seisView;
    mySeisView.getEventHandler().addKeyListener(this);
    mySeisView.getEventHandler().addPanningListener(this);
    //---------------------------------------------
    // Side labels

    mySideLabelVert = new csSideLabelVert( new csISideLabelEventListener() {
      public void zoomEvent( Point p1, Point p2, int zoomMode ) {
        int diff = Math.abs( p2.y - p1.y );
        int mid  = ( p2.y + p1.y ) / 2;
        if( diff < MIN_ZOOM_WIDTH ) {
          zoom( zoomMode, csSeisPane.ZOOM_VERT );
        }
        else { // Zoom in
          float zoomVert = mySeisView.getZoomVert();
          int height = myPanelViewPort.getVisibleRect().height; //Math.min( mySeisView.getVisibleRect().height, mySeisView.getHeight() );
          if( zoomMode == csSeisPane.ZOOM_IN ) {
            zoomVert *= (float)height/(float)diff;
          }
          else {
            zoomVert *= (float)diff/(float)height;
          }
          mySeisView.zoom( zoomVert, mySeisView.getZoomHorz(), 0, mySeisView.yView2Model(mid) );
        }
      }
    });
    mySeisView.addSeisViewListener( this );
    mySideLabelHorz = new csSideLabelHorz( new csISideLabelEventListener() {
      public void zoomEvent( Point p1, Point p2, int zoomMode ) {
        mySideLabelHorz.resetTraceStep();
        int diff = Math.abs( p2.x - p1.x );
        int mid  = ( p2.x + p1.x ) / 2;
        if( diff < MIN_ZOOM_WIDTH ) {
          zoom( zoomMode, csSeisPane.ZOOM_HORZ );
        }
        else { // Zoom in
          float zoomHorz = mySeisView.getZoomHorz();
          int width  = myPanelViewPort.getVisibleRect().width; // Math.min( mySeisView.getVisibleRect().width, mySeisView.getWidth() );
          if( zoomMode == csSeisPane.ZOOM_IN ) {
            zoomHorz *= (float)width/(float)diff;
          }
          else {
            zoomHorz *= (float)diff/(float)width;
          }
          mySeisView.zoom( mySeisView.getZoomVert(), zoomHorz, mySeisView.xView2Trace(mid), 0 );
        }
      }
    });
    mySeisView.addSeisViewListener(this);

    //
    myCornerUpperLeft = new csCornerHorzLabel();
    // GridBagLayout necessary to place seismic view into upper left corner of scroll pane's view port.
    myPanelViewPort = new JPanel( new GridBagLayout() );
    myPanelViewPort.add( mySeisView, new GridBagConstraints(
        0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    myPanelViewPort.add( Box.createHorizontalGlue(), new GridBagConstraints(
        1, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    myPanelViewPort.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, 1, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    myPanelTopLabel = new JPanel( new BorderLayout() );
//    myCornerUpperLeft.setMinimumSize( new Dimension(0,0) );
//    myCornerUpperLeft.setPreferredSize( new Dimension( DEFAULT_WIDTH_SIDELABEL, DEFAULT_HEIGHT_SIDELABEL ) );
    myPanelTopLabel.add( myCornerUpperLeft, BorderLayout.WEST );
    myPanelTopLabel.add( mySideLabelHorz, BorderLayout.CENTER );
    
    //---------------------------------------------
    add( myScrollBarVert, BorderLayout.EAST );
    add( myScrollBarHorz, BorderLayout.SOUTH );
    add( mySideLabelVert, BorderLayout.WEST );
    add( myPanelTopLabel, BorderLayout.NORTH );
    add( myPanelViewPort, BorderLayout.CENTER );
    setMinimumSize( new Dimension(0,0) );
    
    myScrollBarVert.addAdjustmentListener( new AdjustmentListener() {
      public void adjustmentValueChanged( AdjustmentEvent e ) {
        if( myStopScrolling ) return;
        int valueNew = e.getValue();
        //System.out.println("---Scrollbar(vert): prev/now/max  " + myPreviousScrollValueVert + " " + myScrollBarVert.getValue() + ", " + myScrollBarVert.getMaximum() + " " + e.getValueIsAdjusting() );
        if( myPreviousScrollValueVert != valueNew ) {
          myPreviousScrollValueVert = valueNew;
          mySeisView.resetViewPositionVert( valueNew );
        }
      }
    });

    myScrollBarHorz.addAdjustmentListener( new AdjustmentListener() {
      public void adjustmentValueChanged( AdjustmentEvent e ) {
        if( myStopScrolling ) return;
        //System.out.println("Scrollbar(horz):  " + myScrollBarHorz.getValue() + ", " + myScrollBarHorz.getMaximum() + e.getValueIsAdjusting() );
        int valueNew = e.getValue();
        if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
          valueNew = myScrollBarHorz.getMaximum() - valueNew;
        }
        if( myPreviousScrollValueHorz != valueNew ) {
          myPreviousScrollValueHorz = valueNew;
          mySeisView.resetViewPositionHorz( valueNew );
        }
      }
    });
        
    mySeisView.addMouseWheelListener( new MouseWheelListener() {
      public void mouseWheelMoved( MouseWheelEvent e ) {
        int blockStep = myScrollBarVert.getBlockIncrement()/10 + 1;
        if( e.getScrollType() == MouseWheelEvent.WHEEL_UNIT_SCROLL ) {
          scrollVertical( e.getUnitsToScroll()*blockStep );
        }
        else {
          scrollVertical( blockStep );
        }
      }
    });

    addComponentListener( new ComponentAdapter() {
      public void componentResized(ComponentEvent e ) {
        setNewSize();
        mySeisView.resetViewPositionVert( myScrollBarVert.getValue() );
        mySeisView.resetViewPositionHorz( myScrollBarHorz.getValue() );
      }
    });
    
    repaint();
  }
  //-------------------------------------------------------------
  //
  /**
   * Update pane with new seismic trace buffer
   * 
   * @param buffer      Trace buffer containing all traces to be displayed
   * @param sampleInt   Sample interval [ms]
   * @param trcHeaders  Definition of all trace headers in seismic traces (see @class csHeader for more information)
   *                     (size of this array must be the same as number of trace headers in trace buffer)
   */
  public void updateSeismic( csISeismicTraceBuffer buffer, double sampleInt, csHeaderDef[] trcHeaders ) {
    updateSeismic( buffer, sampleInt, trcHeaders, 0, true );
  }
  /**
   * Update pane with new seismic trace buffer
   * 
   * @param buffer       Trace buffer containing all traces to be displayed
   * @param sampleInt    Sample interval [ms]
   * @param trcHeaders   Definition of all trace headers in seismic traces (see @class csHeader for more information)
   *                      (size of this array must be the same as number of trace headers in trace buffer)
   * @param firstTraceIndex  Index of first trace, usually 0. The trace index is used to annotate sequential trace numbers.
   *                          Set to >0 if the sequential trace annotation shall start at >1.
   * @param refresh          Set to true to make that seismic display is fully refreshed.
   */
  public void updateSeismic( csISeismicTraceBuffer buffer, double sampleInt, csHeaderDef[] trcHeaders, int firstTraceIndex, boolean refresh ) {
    myTraceBuffer  = buffer;
    myTraceHeaders = trcHeaders;
    myFirstSeqTraceIndex = firstTraceIndex;
    mySampleInt    = sampleInt;
    myNumTraces    = buffer.numTraces();
    myNumSamples   = buffer.numSamples();

    if( myNumTraces == 0 ) return;

    mySeisView.updateTraceBuffer( buffer, sampleInt, refresh );
    mySideLabelHorz.updateSeismic();

    setNewSize();
  }
  public void setMaxHorizontalScroll() {
//    System.out.println("Max value: " + myScrollBarHorz.getMaximum());
    myScrollBarHorz.setValue(myScrollBarHorz.getMaximum());
  }
  public void setMinHorizontalScroll() {
//    System.out.println("Min value: " + myScrollBarHorz.getMinimum());
    myScrollBarHorz.setValue(myScrollBarHorz.getMinimum());
  }
  public void keyPressed( KeyEvent event ) {
//    boolean isSHIFTPressed = event.isShiftDown();
//    boolean isCTRLPressed  = event.isControlDown();
    int stepVert = 1;
    int stepHorz = 1;
    if( event.isControlDown() ) {
      stepVert = myScrollBarVert.getBlockIncrement()/10 + 1;
      stepHorz = myScrollBarHorz.getBlockIncrement()/10 + 1;
    }
    else if( event.isShiftDown() ) {
      stepVert = myScrollBarVert.getBlockIncrement();
      stepHorz = myScrollBarHorz.getBlockIncrement();
    }
    int code = event.getKeyCode();
    switch( code ) {
    case KeyEvent.VK_UP:
      scrollVertical( -stepVert );
      break;
    case KeyEvent.VK_DOWN:
      scrollVertical( stepVert );
      break;
    case KeyEvent.VK_RIGHT:
      scrollHorizontal( stepHorz );
      break;
    case KeyEvent.VK_LEFT:
      scrollHorizontal( -stepHorz );
      break;
  }
    
  }
  public void keyReleased( KeyEvent event ) {
  }
  //-------------------------------------------------------------
  //
  public JScrollBar getVerticalScrollBar() {
    return myScrollBarVert;
  }
  public JScrollBar getHorizontalScrollBar() {
    return myScrollBarHorz;
  }
  /**
   * Scroll vertically, by step
   * @param scrollStep 
   */
  public void scrollVertical( int scrollStep ) {
    int newValue = myScrollBarVert.getValue()+scrollStep;
    if( newValue > myScrollBarVert.getMaximum() ) newValue = myScrollBarVert.getMaximum();
    else if( newValue < myScrollBarVert.getMinimum() ) newValue = myScrollBarVert.getMinimum();
    
    myScrollBarVert.setValue( newValue );
  }
  /**
   * Scroll horizontally, by step
   * @param scrollStep 
   */
  public void scrollHorizontal( int scrollStep ) {
    int newValue = myScrollBarHorz.getValue()+scrollStep;
    if( newValue > myScrollBarHorz.getMaximum() ) newValue = myScrollBarHorz.getMaximum();
    else if( newValue < myScrollBarHorz.getMinimum() ) newValue = myScrollBarHorz.getMinimum();
    
    myScrollBarHorz.setValue( newValue );
  }
  /**
   * Update header names for trace annotation
   * @param headers Header names
   * @param showSeqTraceNum true if sequential trace number shall be annotated
   */
  public void updateAnnotationTraceHeaders( csHeaderDef[] headers, csAnnotationAttributes attr ) {
    myShowSeqTraceNum     = attr.showSequential;
    myOmitRepeatedValues  = attr.omitRepeating;
    int newHeight = mySideLabelHorz.updateTraceHeaders(headers);
    myCornerUpperLeft.resetPreferredSize( newHeight );
    myCornerUpperLeft.updateTraceHeaders(headers);
    myIsTraceLabelStepAuto = !attr.fixedTraceLabelStep;
    if( attr.fixedTraceLabelStep ) myTraceLabelStep = attr.traceLabelStep;
    myIsTraceStepInitialised = false;
  }
  /**
   * 
   * @return Seismic display scalar
   */
  public float getScalar() {
    return mySeisView.getScalar();
  }
  /**
   * Set seismic display scalar
   * @param scalar Seismic display scalar
   */
  public void setScalar( float scalar) {
    mySeisView.setScalar(scalar);
  }
  public void fitToScreen() {
    csSampleInfo pos1 = new csSampleInfo();
    csSampleInfo pos2 = new csSampleInfo();
    pos1.traceDouble = -1; // mySeismicTraceBuffer.originalTraceNumber(0) - 1;
    pos2.traceDouble = this.myTraceBuffer.numTraces(); // mySeismicTraceBuffer.originalTraceNumber( getDisplayedNumTraces()-1 ) - 1;
    pos1.sampleDouble = 0;
    pos2.sampleDouble = myNumSamples - 1;
    // Can't figure out how to compute zoom factor in one step. Call zoom three times to get fitToScreen correct:
    zoomArea( pos1, pos2, true );
    zoomArea( pos1, pos2, true );
    zoomArea( pos1, pos2, true );
  }
  /**
   * Zoom
   * @param zoomVert Vertical zoom factor (= pixels per time unit)
   * @param zoomHorz Horizontal zoom factor (pixels per trace)
   */
  public void zoom( float zoomVert, float zoomHorz ) {
    mySeisView.zoom( zoomVert, zoomHorz );
    mySideLabelHorz.resetTraceStep();
  }
  /**
   * Zoom to a specific area
   * @param pos1  "Lower left" position of area to zoom to
   * @param pos2  "Upper right" position of area to zoom to
   */
  public void zoomArea( csSampleInfo pos1, csSampleInfo pos2 ) {
    zoomArea( pos1, pos2, false );
  }
  /**
   * Zoom to a specific area
   * @param pos1  "Lower left" position of area to zoom to
   * @param pos2  "Upper right" position of area to zoom to
   * @param includeMargin  Include margin in calculation of horizontal zoom factor
   */
  public void zoomArea( csSampleInfo pos1, csSampleInfo pos2, boolean includeMargin ) {
    int height = myPanelViewPort.getVisibleRect().height; //Math.min( mySeisView.getVisibleRect().height, mySeisView.getHeight() );
    int width  = myPanelViewPort.getVisibleRect().width; // Math.min( mySeisView.getVisibleRect().width, mySeisView.getWidth() );
    
    double traceDiffNew  = Math.abs( pos2.traceDouble - pos1.traceDouble );
    double sampleDiffNew = Math.abs( pos2.sampleDouble - pos1.sampleDouble );
    float traceMid  = 0.5f * (float)( pos2.traceDouble + pos1.traceDouble );
    float sampleMid = 0.5f * (float)( pos2.sampleDouble + pos1.sampleDouble );

    double traceDiffCurrent  = Math.abs( mySeisView.xView2Trace(0) - mySeisView.xView2Trace(width) );
    double sampleDiffCurrent = Math.abs( mySeisView.yView2Model(0) - mySeisView.yView2Model(height) );

    float zoomVert = mySeisView.getZoomVert();
    float zoomHorz = mySeisView.getZoomHorz();

    if( traceDiffNew < MIN_ZOOM_WIDTH && sampleDiffNew < MIN_ZOOM_WIDTH ) {
      zoom( csSeisPane.ZOOM_OUT, csSeisPane.ZOOM_BOTH );
    }
    else {
      traceDiffNew      = Math.max( MIN_ZOOM_WIDTH, traceDiffNew );
      sampleDiffNew     = Math.max( MIN_ZOOM_WIDTH, sampleDiffNew );
      traceDiffCurrent  = Math.max( MIN_ZOOM_WIDTH, traceDiffCurrent );
      sampleDiffCurrent = Math.max( MIN_ZOOM_WIDTH, sampleDiffCurrent );

      float zoomVertNew = zoomVert * (float)( sampleDiffCurrent / sampleDiffNew );
      float zoomHorzNew = zoomHorz * (float)( traceDiffCurrent / traceDiffNew );
      if( includeMargin ) {
        int widthNew = width - (int)( (float)mySeisView.getMarginLeftRight() * zoomHorzNew/zoomHorz ) - 10;
        traceDiffCurrent  = Math.abs( mySeisView.xView2Trace(0) - mySeisView.xView2Trace(widthNew) );
        traceDiffCurrent  = Math.max( MIN_ZOOM_WIDTH, traceDiffCurrent );
        zoomHorzNew = zoomHorz * (float)( traceDiffCurrent / traceDiffNew );
      }
      mySeisView.zoom( zoomVertNew, zoomHorzNew, traceMid, sampleMid );
    }
  }
  public void zoomHorz( float zoomLevel ) {
    if( zoomLevel <= 0.0 ) return;
    zoom( zoomLevel, mySeisView.getZoomHorz() );
  }
  public void zoomVert( float zoomLevel ) {
    if( zoomLevel <= 0.0 ) return;
    zoom( mySeisView.getZoomVert(), zoomLevel );
  }
  public void zoom( int inOut, int zoomType ) {
    float zoomVert = mySeisView.getZoomVert();
    float zoomHorz = mySeisView.getZoomHorz();
    if( (zoomType & csSeisPane.ZOOM_VERT) != 0 ) {
      if( inOut == csSeisPane.ZOOM_IN ) {
        zoomVert *= 2.0f;
      }
      else {
        zoomVert /= 2.0f;
      }
      if( zoomVert == 0 ) zoomVert = 1;
    }
    if( (zoomType & csSeisPane.ZOOM_HORZ) != 0 ) {
      if( inOut == csSeisPane.ZOOM_IN ) {
        zoomHorz *= 2.0f;
      }
      else {
        zoomHorz /= 2.0f;
      }
      if( zoomHorz == 0 ) zoomHorz = 4.0f;
    }
    zoom( zoomVert, zoomHorz );
  }
  public void paintComponent( Graphics g ) {
    super.paintComponent( g );
  }
  /**
   * Manually highlight trace
   * @param traceIndex Trace index starting at 0
   */
  public void manualHighlightTrace( int traceIndex ) {
    mySeisView.removeHighlightTrace();
    mySeisView.highlightTrace( traceIndex );
  }
  /**
   * Set automatic trace highlighting
   * @param doHighlight true to set automatic trace highlighting
   */
  public void setAutoHighlightTrace( boolean doHighlight ) {
    mySeisView.setAutoHighlightTrace( doHighlight );
  }
  /**
   * 
   * @param setOn true to set crosshair ON
   * @param setSnap true if crosshair shall be snapped to nearest full sample/trace
   */
  public void setCrosshair( boolean setOn, boolean setSnap ) {
    mySeisView.setCrosshair( setOn, setSnap );
  }
  public void setParentFrame( JFrame frame ) {
    myParentFrame = frame;
  }
  /**
   * Show annotation dialog
   */
  public void showAnnotationDialog() {
    if( myAnnotationDialog == null ) {
      if( myTraceHeaders == null ) return;
      myAnnotationDialog = new csAnnotationDialog( myParentFrame, this, myTraceHeaders );
    }
    myAnnotationDialog.setVisible( true );
  }
  /**
   * Show seismic display settings dialog
   */
  public void showSettingsDialog() {
    mySeisView.showDispSettingsDialog();
  }
  public void vertScrollChanged( int scrollValue ) {
    mySideLabelVert.repaint();
  }
  public void horzScrollChanged( int scrollValue ) {
    mySideLabelHorz.repaint();
  }
  public void sizeChanged( Dimension size ) {
  }
  private void setNewSize() {
    mySeisViewSize        = mySeisView.getPreferredSize();
    if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
      int saveHeight = mySeisViewSize.height;
      mySeisViewSize.height = mySeisViewSize.width;
      mySeisViewSize.width = saveHeight;
    }

    Rectangle rectVisible = mySeisView.getVisibleRect();
    int heightVisible     = rectVisible.height;
    int widthVisible      = rectVisible.width;
    int maxScrollValueVert = mySeisViewSize.height - heightVisible;
    int maxScrollValueHorz = mySeisViewSize.width  - widthVisible;
    int scrollValueVert    = mySeisView.getViewPositionVert();
    int scrollValueHorz    = mySeisView.getViewPositionHorz();

//    System.out.println("resetPane 1 - Scroll/seisview: " + myScrollBarVert.getValue() + " " + mySeisView.getViewPositionVert() + ", maxScroll/new:" + myScrollBarVert.getMaximum() + " " + maxScrollValueVert +
//        ", Seisview size(full/visible): " + mySeisViewSize.height + " " + heightVisible);
//    System.out.println("resetPane 1b- Scroll/seisview: " + myScrollBarHorz.getValue() + " " + mySeisView.getViewPositionHorz() + ", maxScroll/new:" + myScrollBarHorz.getMaximum() + " " + maxScrollValueHorz +
//        ", Seisview size(full/visible): " + mySeisViewSize.width + " " + widthVisible);

    myPreviousScrollValueVert = -1;
    myPreviousScrollValueHorz = -1;

    // Complex coding to avoid refreshing of scroll bars and seismic while setting new maximum values etc. ...not desirable.
    myStopScrolling = true;
    myScrollBarVert.setBlockIncrement( (heightVisible/2)+1 );
    myScrollBarVert.setMaximum( maxScrollValueVert );
    myScrollBarHorz.setBlockIncrement( (widthVisible/2)+1 );
    myScrollBarHorz.setMaximum( maxScrollValueHorz );
    myStopScrolling = false;

    myScrollBarVert.setValue( scrollValueVert <= maxScrollValueVert ? scrollValueVert : maxScrollValueVert );
    myScrollBarHorz.setValue( scrollValueHorz <= maxScrollValueHorz ? scrollValueHorz : maxScrollValueHorz );

//    System.out.println("resetPane 2 - Scroll/seisview: " + myScrollBarVert.getValue() + " " + ( scrollValueVert <= maxScrollValueVert ? scrollValueVert : maxScrollValueVert ) );
//    System.out.println("resetPane 2b- Scroll/seisview: " + myScrollBarHorz.getValue() + " " + ( scrollValueHorz <= maxScrollValueHorz ? scrollValueHorz : maxScrollValueHorz ) );
    mySideLabelVert.repaint();
    mySideLabelHorz.repaint();
//    if( myTraceBuffer != null && mySeisView != null ) {
//      System.out.println("SEISVIEW MIN/MAX: numtrc: " + myTraceBuffer.numTraces() + " pixels: " +
//            mySeisView.xModel2View(0) + " " + mySeisView.xModel2View(myTraceBuffer.numTraces()-1) + " diff: "
//            + (mySeisView.xModel2View(myTraceBuffer.numTraces()-1)-mySeisView.xModel2View(0)) + "  size: " +
//            mySeisViewSize.width + " pane: " + getPreferredSize().width);
//    }
  }

  public void traceBufferChanged( csISeismicTraceBuffer traceBuffer ){
    // Nothing to do
  }
  public void changedSettings( csSeisDispSettings settings ) {
    mySideLabelVert.changedSettings(settings);
    mySideLabelHorz.changedSettings(settings);
    Dimension sizeSeisView = mySeisView.getPreferredSize();
    boolean resetSize = false;
    if( myPlotDirection != settings.plotDirection ) {
      myPlotDirection = settings.plotDirection;
      // Reset scroll bars & side labels
      mySideLabelVert.resetPlotDirection();
      int newHeight = mySideLabelHorz.resetPlotDirection();
      myCornerUpperLeft.resetPreferredSize( newHeight );

      myPanelTopLabel.removeAll();
      //---------------------------------------------
      remove( myScrollBarVert );
      remove( myScrollBarHorz );
      if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
        remove( mySideLabelVert );
        add( mySideLabelHorz, BorderLayout.WEST );
        myPanelTopLabel.add( myCornerUpperLeft, BorderLayout.WEST );
        myPanelTopLabel.add( mySideLabelVert, BorderLayout.CENTER );
        myScrollBarVert.setOrientation(JScrollBar.HORIZONTAL);
        myScrollBarHorz.setOrientation(JScrollBar.VERTICAL);
        add( myScrollBarVert, BorderLayout.SOUTH );
        add( myScrollBarHorz, BorderLayout.EAST );

        remove( myPanelTopLabel );
        add( myScrollBarVert, BorderLayout.NORTH );
        add( myPanelTopLabel, BorderLayout.SOUTH );
      }
      else {
        remove( mySideLabelHorz );
        add( mySideLabelVert, BorderLayout.WEST );
        myPanelTopLabel.add( myCornerUpperLeft, BorderLayout.WEST );
        myPanelTopLabel.add( mySideLabelHorz, BorderLayout.CENTER );
        myScrollBarVert.setOrientation(JScrollBar.VERTICAL);
        myScrollBarHorz.setOrientation(JScrollBar.HORIZONTAL);
        remove( myPanelTopLabel );

        add( myPanelTopLabel, BorderLayout.NORTH );
        add( myScrollBarVert, BorderLayout.EAST );
        add( myScrollBarHorz, BorderLayout.SOUTH );
      }
      revalidate();
      // Trick to reset all settings in seisview:
      mySeisView.updateTraceBuffer( mySeisView.getTraceBuffer(), mySampleInt, true );
      resetSize = true;
    }
    if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
      int saveHeight = sizeSeisView.height;
      sizeSeisView.height = sizeSeisView.width;
      sizeSeisView.width  = saveHeight;
      setNewSize();
    }
    else if( mySeisViewSize != null && ( resetSize || sizeSeisView.width != mySeisViewSize.width || sizeSeisView.height != mySeisViewSize.height ) ) {
      setNewSize();
    }
  }

  //-----------------------------------------------------
  //*************************************************************
  // Side labels
  //
  public void setTraceLabelStep( int traceLabelStep, boolean isTraceLabelStepAuto ) {
    myTraceLabelStep = traceLabelStep;
    myIsTraceLabelStepAuto = isTraceLabelStepAuto;
    myIsTraceStepInitialised = false;
  }
  public void setVerticalDomain( int domainType ) {
    myVerticalDomain = domainType;
    if( myVerticalDomain == csUnits.DOMAIN_FREQ ) {
      myTitleVerticalAxis = "Frequency [Hz]";
      mySampleIntScalar   = 1;
    }
    else if( myVerticalDomain == csUnits.DOMAIN_TIME ) {
      myTitleVerticalAxis = "Time [seconds]";
      mySampleIntScalar   = 1000;
    }
    else if( myVerticalDomain == csUnits.DOMAIN_DEPTH ) {
      myTitleVerticalAxis = "Depth [kilometers]";
      mySampleIntScalar   = 1000;
    }
    mySideLabelVert.repaint();
  }
  public void setVerticalAxisTitle( String title ) {
    myTitleVerticalAxis = title;
  }
  @Override
  public void hasPanned(int dx, int dy) {
    scrollHorizontal(dx);
    scrollVertical(dy);
  }
  @Override
  public void mouseExited() {
  }

  public class csSideLabelVert extends JPanel {
    public static final int MIN_INTERVAL = 40;
    public static final int MAX_INTERVAL = 80;
    
    private DecimalFormat myDecimalFormat;
    private int myNumDecimals;
    private double myHorzLineMinorInc;
    private double myHorzLineMajorInc;
    private csSideLabelMouseHandler myMouseHandler;
//    private boolean myShowLines;

    // Bitmap
    private BufferedImage myBitmap;
    protected int myWidth;
    protected int myHeight;
    
    public csSideLabelVert( csISideLabelEventListener listener ) {
      super();
      myBitmap = null;
//      myShowLines   = false;
      myHorzLineMinorInc = 100;
      myHorzLineMajorInc = 500;
      setMinimumSize( new Dimension(0,0) );
      setPreferredSize( new Dimension(DEFAULT_WIDTH_SIDELABEL,0) );
      setBorder( BorderFactory.createCompoundBorder(
          BorderFactory.createRaisedBevelBorder(),
          BorderFactory.createLoweredBevelBorder() ) );
      setNumDecimals( 1 );
      myMouseHandler = new csSideLabelMouseHandler( this, listener );
      addMouseListener(myMouseHandler);
      addMouseMotionListener(myMouseHandler);
    }
    public void resetPlotDirection() {
      if( myPlotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        setPreferredSize( new Dimension(DEFAULT_WIDTH_SIDELABEL,0) );
      }
      else {
        setPreferredSize( new Dimension(0,DEFAULT_WIDTH_SIDELABEL) );
      }
    }
//  ---------------------------------------------------
    public void setNumDecimals( int numDecimals ) {
      myNumDecimals = numDecimals;
      String s = new String( "0" );
      if( myNumDecimals > 0 ) {
        s += ".0";
      }
      for( int i = 1; i < myNumDecimals; i++ ) {
        s += "0";
      }
      myDecimalFormat = new DecimalFormat( s );
    }
    public void paintComponent( Graphics g ) {
      super.paintComponent( g );
      Graphics2D g2 = (Graphics2D)g;
      Rectangle rect = getVisibleRect();
      if( myBitmap == null || myWidth != rect.width || myHeight != rect.height ) {
        myWidth  = rect.width;
        myHeight = rect.height;
        if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
          myWidth  = rect.height;
          myHeight = rect.width;
        }
        myBitmap = new BufferedImage( myWidth, myHeight, BufferedImage.TYPE_INT_RGB );
        if( myBitmap == null ) return;
      }
      repaintSideLabelVert( (Graphics2D)myBitmap.getGraphics() );
      if( myPlotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        g2.drawImage( myBitmap, 0, 0, null );
      }
      else { // if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
        java.awt.geom.AffineTransform affineTransform = new java.awt.geom.AffineTransform();
        affineTransform.translate( 0, myWidth );
        affineTransform.rotate( -Math.PI/2 );
        g2.drawImage( myBitmap, affineTransform, null );
      }
    }
    private void repaintSideLabelVert( Graphics2D g ) {
      Insets insets  = getInsets();
      FontMetrics metrics = g.getFontMetrics( g.getFont() );
      int labelHeight = metrics.getHeight();
      
      int sampMinorInc = (int)( myHorzLineMinorInc / mySampleInt + 0.5f );
      int sampMajorInc = (int)( myHorzLineMajorInc / mySampleInt + 0.5f );
      if( sampMinorInc <= 0 ) sampMinorInc = 1;
      if( sampMajorInc <= 0 ) sampMajorInc = 1;
      
      g.setColor( DEFAULT_COLOR );
      g.fillRect( 0, 0, myWidth, myHeight );
      g.setFont( new Font("Fixed", Font.PLAIN, DEFAULT_LABEL_FONTSIZE) );

      int minSamp = (int)mySeisView.yView2Model( -1 );
      int maxSamp = (int)mySeisView.yView2Model( myHeight+1 ) + 1;
//      System.out.println(" " + minSamp + " " + maxSamp + " " + myHeight);
      maxSamp = Math.min( myNumSamples-1, maxSamp );

      g.setStroke( new BasicStroke(1) );
      g.setColor( Color.black );

      int firstSamp = ( (int)( minSamp / sampMinorInc - 0.1f ) + 1 ) * sampMinorInc;
      int xmin, xmax;

      xmin = myWidth-4-insets.right;
      xmax = myWidth-1-insets.right;

      int xpixText = insets.left+1 + labelHeight-4;// Add labelHeight to accomodate for 'Time [seconds]' label

//      if( myShowLines ) {
        for( int isamp = firstSamp; isamp < maxSamp; isamp += sampMinorInc ) {
          int ypix = (int)mySeisView.yModel2View( isamp );
          String label = myDecimalFormat.format( (double)(isamp)*mySampleInt / (double)mySampleIntScalar );
          int labelWidth  = metrics.stringWidth( label );
          int xpix = Math.max( xmin-labelWidth-2, xpixText );
          if( isamp % sampMajorInc != 0 ) {
            g.drawString( label, xpix, ypix+DEFAULT_LABEL_FONTSIZE/2-1 );
            g.drawLine( xmin, ypix, xmax, ypix );
          }
          else {
            g.setStroke( new BasicStroke(2) );
            g.drawString( label, xpix, ypix+DEFAULT_LABEL_FONTSIZE/2-1 );
            g.drawLine( xmin-1, ypix, xmax, ypix );
            g.setStroke( new BasicStroke(1) );
          }
        }
//      }
      
      // Plot Axis label
      g.rotate( -Math.PI/2 );
      int labelWidth  = metrics.stringWidth( myTitleVerticalAxis );
      int startPos = -myHeight/2 - labelWidth/2;
//      g.setColor( Color.white );
//      g.fillRect( startPos, 0, labelWidth, labelHeight );
      g.setColor( Color.black );
      g.drawString( myTitleVerticalAxis, startPos, DEFAULT_LABEL_FONTSIZE+2 );
      g.rotate( Math.PI/2 );

      g.setColor( Color.red );
      g.setStroke( new BasicStroke(2) );
      
      if( myPlotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        switch( myMouseHandler.myMouseMode ) {
        case csSideLabelMouseHandler.MOUSE_DRAGGED:
          g.drawLine( 0, myMouseHandler.myMouseDraggedPos.y, myWidth, myMouseHandler.myMouseDraggedPos.y );
          // Fall through
        case csSideLabelMouseHandler.MOUSE_MOVED:
          g.drawLine( 0, myMouseHandler.myMouseMovedPos.y, myWidth, myMouseHandler.myMouseMovedPos.y );
        }
      }
      else {
        switch( myMouseHandler.myMouseMode ) {
        case csSideLabelMouseHandler.MOUSE_DRAGGED:
          g.drawLine( 0, myMouseHandler.myMouseDraggedPos.x, myHeight, myMouseHandler.myMouseDraggedPos.x );
          // Fall through
        case csSideLabelMouseHandler.MOUSE_MOVED:
          g.drawLine( 0, myMouseHandler.myMouseMovedPos.x, myHeight, myMouseHandler.myMouseMovedPos.x );
        }
      }
    }
    public void changedSettings( csSeisDispSettings settings ) {
      csSeisDispSettings ds = settings;
//      myShowLines = ds.showTimeLines;
//      if( !myShowLines ) return;
      if( (ds.timeLineMinorInc != myHorzLineMinorInc || ds.timeLineMajorInc != myHorzLineMajorInc) ) {
        long number = (long)( 1000*(double)ds.timeLineMinorInc );
        if( number == 0 ) return;
        myHorzLineMinorInc = ds.timeLineMinorInc;
        myHorzLineMajorInc = ds.timeLineMajorInc;
        recomputeDecimals();
        repaint();
      }
    }
    private void recomputeDecimals() {
      long number = (long)( 1000*(double)myHorzLineMinorInc );
      int minDecimal = 1;
      while( number != 0 && number % 10 == 0 ) {
        number /= 10;
        minDecimal += 1;
      }
      int numDecimals = 7 - minDecimal;
      if( numDecimals <= 0 ) numDecimals = 1;

      if( numDecimals != myNumDecimals ) setNumDecimals( numDecimals );
    }
  }

  //******************************************************************************
  //******************************************************************************
  //
  public class csSideLabelHorz extends JPanel {
    private csSideLabelMouseHandler myMouseHandler;
    private int[] myHeaderIndex;
    private String[] myDisplayedHeaderNames;
    private int myNumHeaders;
    private int myTraceStep;

    // Bitmap
    private BufferedImage myBitmap;
    protected int myWidth;
    protected int myHeight;

    public csSideLabelHorz( csISideLabelEventListener listener ) {
      super();
      myBitmap = null;
      myTraceStep = 1;
      myHeaderIndex = null;
      myDisplayedHeaderNames = null;
      myNumHeaders  = 0;
      setMinimumSize( new Dimension(0,0) );
      resetPlotDirection();
      setBorder( BorderFactory.createCompoundBorder(
          BorderFactory.createRaisedBevelBorder(),
          BorderFactory.createLoweredBevelBorder() ) );
      myMouseHandler = new csSideLabelMouseHandler( this, listener );
      addMouseListener(myMouseHandler);
      addMouseMotionListener(myMouseHandler);
    }
    public void resetPreferredSize( int newHeight ) {
      if( myPlotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        setPreferredSize( new Dimension(0,newHeight) );
      }
      else {
        setPreferredSize( new Dimension(newHeight,0) );
      }
    }
    protected int setNewSize() {
      int labelHeight = DEFAULT_LABEL_FONTSIZE+3;
      int newSize = DEFAULT_HEIGHT_SIDELABEL + myNumHeaders*(labelHeight);
      if( !myShowSeqTraceNum ) newSize -= (labelHeight);
      resetPreferredSize( newSize );
      return newSize;
    }
    public int resetPlotDirection() {
      return setNewSize();
    }
    public void updateSeismic() {
      if( myNumHeaders > 0 ) {
        for( int i = 0; i < myDisplayedHeaderNames.length; i++ ) {
          myHeaderIndex[i] = -1;
          for( int ih = 0; ih < myTraceHeaders.length; ih++ ) {
            if( myTraceHeaders[ih].name.compareTo(myDisplayedHeaderNames[i]) == 0 ) {
              myHeaderIndex[i] = ih;
            }
          }
        }
      }
    }
    public int updateTraceHeaders( csHeaderDef[] headers ) {
      myNumHeaders  = headers.length;
      int newSize = setNewSize();
      if( myNumHeaders == 0 ) {
        myHeaderIndex = null;
      }
      else {
        myDisplayedHeaderNames = new String[myNumHeaders];
        myHeaderIndex = new int[myNumHeaders];
        for( int i = 0; i < myDisplayedHeaderNames.length; i++ ) {
          myDisplayedHeaderNames[i] = headers[i].name;
        }
        updateSeismic();
      }
      myIsTraceStepInitialised = false;
      revalidate();
      repaint();
      return newSize;
    }
    public void paintComponent( Graphics g ) {
      if( myTraceBuffer == null ) return;
      Graphics2D g2 = (Graphics2D)g;
      Rectangle rect = getVisibleRect();
      if( myBitmap == null || myWidth != rect.width || myHeight != rect.height ) {
        myWidth  = rect.width;
        myHeight = rect.height;
        if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
          myWidth  = rect.height;
          myHeight = rect.width;
        }
        myBitmap = new BufferedImage( myWidth, myHeight, BufferedImage.TYPE_INT_RGB );
        if( myBitmap == null ) return;
      }
      repaintSideLabelHorz( (Graphics2D)myBitmap.getGraphics() );
      if( myPlotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        g2.drawImage( myBitmap, 0, 0, null );
      }
      else { // if( myPlotDirection = csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
        java.awt.geom.AffineTransform affineTransform = new java.awt.geom.AffineTransform();
        affineTransform.rotate( -Math.PI/2 );
        affineTransform.translate( -myWidth, 0 );
        g2.drawImage( myBitmap, affineTransform, null );
      }
    }
    private void repaintSideLabelHorz( Graphics2D g ) {
      g.setColor( DEFAULT_COLOR );
      g.fillRect( 0, 0, myWidth, myHeight );
      if( myNumTraces == 0 || myTraceBuffer.numTraces() != myNumTraces ) return;

      Insets insets  = getInsets();
      g.setFont( new Font("Fixed", Font.PLAIN, DEFAULT_LABEL_FONTSIZE) );
//        FontMetrics metrics = g.getFontMetrics( g.getFont() );
//        myLabelFontPixelHeight = metrics.getHeight();
//        int textWidth = metrics.stringHeight( "1000" );

      int minTrace = (int)mySeisView.xView2Trace( 0 );
      int maxTrace = (int)mySeisView.xView2Trace( myWidth+1 ) + 1;
      minTrace = Math.max( 0, minTrace );
      maxTrace = Math.min( myNumTraces-1, maxTrace );

      g.setStroke( new BasicStroke(1) );
      g.setColor( Color.black );
      FontMetrics metrics = g.getFontMetrics( g.getFont() );
      
      // Increase trace step automatically when size is too small
      if( myNumTraces > 1 && !myIsTraceStepInitialised ) {
        if( myIsTraceLabelStepAuto ) {
          String label = "";
          int labelWidth = 0;
          if( myShowSeqTraceNum ) {
            label = "" + myFirstSeqTraceIndex + myNumTraces;
            labelWidth = metrics.stringWidth( label );
          }
          int maxWidth = labelWidth;
          for( int itrc = 0; itrc < myNumTraces; itrc++ ) {
            csHeader[] values = myTraceBuffer.headerValues(itrc);
            for( int ih = 0; ih < myNumHeaders; ih++ ) {
              if( myHeaderIndex[ih] < 0 ) continue;
  // Temp fix to avoid problems when repaint happens while data is still being loaded. Needs more fundamental fix
              if( myHeaderIndex[ih] > values.length ) {
  //              System.out.println(" " + itrc + " " + myNumTraces + " " + myHeaderIndex[ih] + " " + values.length);
                return;
              }
              label = "" + values[myHeaderIndex[ih]];
              labelWidth = metrics.stringWidth( label );
              if( labelWidth > maxWidth ) maxWidth = labelWidth;
            }
          }
          int xpix1 = (int)mySeisView.xModel2View( 0, 0 );
          int xpix2 = (int)mySeisView.xModel2View( 0, myNumTraces-1 );
          double x_int = (double)(xpix2-xpix1)/(double)(myNumTraces-1);
          if( x_int <= 0 ) return;
//          myTraceStep = ((int)((double)maxWidth / (double)(myTraceLabelStep * x_int)) + 1 ) * myTraceLabelStep;
          myTraceStep = ((int)((double)maxWidth / (double)x_int) + 1 );
        }
        else {
          myTraceStep = myTraceLabelStep;
        }
        myIsTraceStepInitialised = true;
      }
      
      int firstTrace = ( (int)( minTrace / myTraceStep ) ) * myTraceStep;
      int ymin, ymax;

      int labelHeight = metrics.getHeight();
      ymin = myHeight-4-insets.bottom;
      ymax = myHeight-1-insets.bottom;
      int ypix = myHeight - (ymax - ymin) - insets.bottom - 2;
      int addOneSeqTrace = 0;
      if( myShowSeqTraceNum ) {
        addOneSeqTrace = 1;
        for( int itrc = firstTrace; itrc <= maxTrace; itrc += myTraceStep ) {
//          String label = "" + (itrc+1+myFirstSeqTraceIndex);
//          String label = "" + ( myTraceBuffer.originalTrace( itrc ).originalTraceNumber() );
          String label = "" + ( myTraceBuffer.originalTraceNumber(itrc) );
          int labelWidth  = metrics.stringWidth( label );
          int xpix = (int)mySeisView.xModel2View( 0, itrc );
          g.drawString( label, xpix-labelWidth/2, ypix );
          g.drawLine( xpix, ymin, xpix, ymax );
        }
      }
      DecimalFormat format = new DecimalFormat("0");
      csHeader[] valuesPrev = new csHeader[myNumHeaders];
      for( int ih = 0; ih < myNumHeaders; ih++ ) {
        valuesPrev[ih] = new csHeader();
      }
      for( int itrc = firstTrace; itrc <= maxTrace; itrc += myTraceStep ) {
        csHeader[] values = myTraceBuffer.headerValues(itrc);
        for( int ih = 0; ih < myNumHeaders; ih++ ) {
          if( myHeaderIndex[ih] < 0 ) continue;
          int ypixNew = ypix - (ih+addOneSeqTrace)*(labelHeight);
          double value = values[myHeaderIndex[ih]].doubleValue();
          if( myOmitRepeatedValues && value == valuesPrev[ih].doubleValue() && itrc > 0 ) continue;
          valuesPrev[ih].setValue( value );
          String label = "" + format.format( value );
          int labelWidth  = metrics.stringWidth( label );
          int xpix = (int)mySeisView.xModel2View( 0, itrc );
          g.drawString( label, xpix-labelWidth/2, ypixNew );
        }
      }
      g.setColor( Color.red );
      g.setStroke( new BasicStroke(2) );
      if( myPlotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        switch( myMouseHandler.myMouseMode ) {
          case csSideLabelMouseHandler.MOUSE_DRAGGED:
            g.drawLine( myMouseHandler.myMouseDraggedPos.x, 0, myMouseHandler.myMouseDraggedPos.x, myHeight );
            // Fall through
          case csSideLabelMouseHandler.MOUSE_MOVED:
            g.drawLine( myMouseHandler.myMouseMovedPos.x, 0, myMouseHandler.myMouseMovedPos.x, myHeight );
        }
      }
      else { //if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
        switch( myMouseHandler.myMouseMode ) {
        case csSideLabelMouseHandler.MOUSE_DRAGGED:
          g.drawLine( myWidth-myMouseHandler.myMouseDraggedPos.y, 0, myWidth-myMouseHandler.myMouseDraggedPos.y, myWidth );
          // Fall through
        case csSideLabelMouseHandler.MOUSE_MOVED:
          g.drawLine( myWidth-myMouseHandler.myMouseMovedPos.y, 0, myWidth-myMouseHandler.myMouseMovedPos.y, myWidth );
        }
      }
    }
    public void changedSettings( csSeisDispSettings settings ) {
      resetTraceStep();
      repaint();  // Always repaint
    }
    public void resetTraceStep() {
      myIsTraceStepInitialised = false;
    }
  }
  /**
   * Upper left-hand corner of trace annotation area
   * This box shows the names of all annotated trace headers
   * 
   * Author: Bjorn Olofsson
   * Date:   15 Aug 2007
   */
  public class csCornerHorzLabel extends JPanel {
    private String labelTrace = "Trace";
    private ArrayList<String> myTraceHeaderNames;
    // Bitmap
    private BufferedImage myBitmap;
    protected int myWidth;
    protected int myHeight;

    public csCornerHorzLabel() {
      super();
      myBitmap = null;
      myTraceHeaderNames = new ArrayList<String>();
      if( myShowSeqTraceNum ) myTraceHeaderNames.add( labelTrace );
      setMinimumSize( new Dimension(0,0) );
      setPreferredSize( new Dimension(DEFAULT_WIDTH_SIDELABEL,DEFAULT_HEIGHT_SIDELABEL) );
      resetPlotDirection();
      setBorder( BorderFactory.createLoweredBevelBorder() );
    }
    public void resetPreferredSize( int newHeight ) {
      if( myPlotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        setPreferredSize( new Dimension(DEFAULT_WIDTH_SIDELABEL,newHeight) );
      }
      else {
        setPreferredSize( new Dimension(newHeight,DEFAULT_WIDTH_SIDELABEL) );
      }
    }
    public void resetPlotDirection() {
      resetPreferredSize( DEFAULT_HEIGHT_SIDELABEL );
      revalidate();
    }
    public void updateTraceHeaders( csHeaderDef[] headers ) {
      myTraceHeaderNames.clear();
      if( myShowSeqTraceNum ) myTraceHeaderNames.add( labelTrace );
      for( int i = 0; i < headers.length; i++ ) {
        myTraceHeaderNames.add( headers[i].name );
      }
      revalidate();
      repaint();
    }
    public void paintComponent( Graphics g ) {
      Graphics2D g2 = (Graphics2D)g;
      Rectangle rect = getVisibleRect();
      if( myBitmap == null || myWidth != rect.width || myHeight != rect.height ) {
        myWidth  = rect.width;
        myHeight = rect.height;
        if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
          myWidth  = rect.height;
          myHeight = rect.width;
        }
        myBitmap = new BufferedImage( myWidth, myHeight, BufferedImage.TYPE_INT_RGB );
        if( myBitmap == null ) return;
      }
      repaintCornerLabelHorz( (Graphics2D)myBitmap.getGraphics() );
      if( myPlotDirection == csSeisDispSettings.PLOT_DIR_VERTICAL ) {
        g2.drawImage( myBitmap, 0, 0, null );
      }
      else { // if( myPlotDirection == csSeisDispSettings.PLOT_DIR_HORIZONTAL ) {
        java.awt.geom.AffineTransform affineTransform = new java.awt.geom.AffineTransform();
        affineTransform.rotate( -Math.PI/2 );
        affineTransform.translate( -myWidth, 0 );
        g2.drawImage( myBitmap, affineTransform, null );
      }
    }
    private void repaintCornerLabelHorz( Graphics2D g ) {
      g.setFont( new Font("SansSerif", Font.PLAIN, DEFAULT_LABEL_FONTSIZE) );
      g.setColor( DEFAULT_COLOR );
      g.fillRect( 0, 0, myWidth, myHeight );
      g.setStroke( new BasicStroke(1) );
      g.setColor( Color.black );
//      Dimension size = getPreferredSize();
      Insets insets  = getInsets();
//      int height = insets.top-insets.bottom;

      FontMetrics metrics = g.getFontMetrics( g.getFont() );
      int labelHeight = metrics.getHeight();
      for( int i = 0; i < myTraceHeaderNames.size(); i++ ) {
        String text = myTraceHeaderNames.get(i);
        int labelWidth  = metrics.stringWidth( text );
        int xpos = Math.max( insets.left, myWidth - insets.right - labelWidth - 2 );
        int ypos = myHeight - insets.bottom - 6 - i*(labelHeight);
        g.drawString( text, xpos, ypos );
      }
    }
  }
  public void test( int numIter ) {
    mySeisView.test( numIter );
  }
}



