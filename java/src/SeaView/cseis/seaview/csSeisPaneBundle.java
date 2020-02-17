/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.general.csColorBarPanel;
import cseis.general.csCustomColorMap;
import cseis.general.csFilename;
import cseis.graph.csCurve;
import cseis.graph.csCurveAttributes;
import cseis.graph.csGraph2D;
import cseis.graph.csGraphAttributes;
import cseis.graph.csSpectrumGraph;
import cseis.jni.csFFTObject;
import cseis.jni.csITraceHeaderScanNotifier;
import cseis.jni.csNativeFFTTransform;
import cseis.jni.csNativeRSFReader;
import cseis.jni.csNativeSegdReader;
import cseis.jni.csNativeSegyReader;
import cseis.segy.csSegyHeaderView;
import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicReader;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csSeismicTrace;
import cseis.seis.csTraceBuffer;
import cseis.jni.csVirtualSeismicReader;
import cseis.processing.csIProcessing;
import cseis.processing.csIProcessingSetupListener;
import cseis.processing.csProcessingAGC;
import cseis.processing.csProcessingFilter;
import cseis.processing.csProcessingInterpolation;
import cseis.processing.csProcessingDCRemoval;
import cseis.processing.csProcessingDialog;
import cseis.processing.csProcessingOverlay;
import cseis.processing.csProcessingTest;
import cseis.seisdisp.*;
import cseis.tools.csIStoppable;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import javax.imageio.ImageIO;
import javax.imageio.ImageReader;
import javax.imageio.stream.FileImageInputStream;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSplitPane;

/**
 * A csSeisPaneBundle object bundles all objects related to one seismic data set viewed in SeaView.<br>
 * <br>
 * This includes the seismic view (csSeisView) which displays the seismic data, its seismic pane (csSeisPane)
 * which adds side labels, a graph (csGraph2D) which plots trace header values and others.
 * The bundle also manages all dialog window objects which are related to this data set.<br>
 * <br>
 * The bundle implements a number of listener interfaces which are used to communicate entries in dialog windows,
 * mouse actions in the views and others.
 * 
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csSeisPaneBundle extends JPanel implements csIGraphPanelListener,
        csISampleInfoListener, csIRubberBandListener, csIStoppable,
        csITraceHeaderScanListener, csITraceHeaderScanNotifier, csIPaintDialogListener, csIPickDialogListener {
  private static int id_counter = 0;
  /// Number of traces to read before progress bar is refreshed:
  private static final int REFRESH_FREQUENCY = 20;
  
  public static int OPERATION_SUBTRACT = 0; // Do not change number values
  public static int OPERATION_ADD      = 1;
  public static int OPERATION_MULTIPLY = 2;
  public static int OPERATION_DIVIDE   = 3;
  public static String[] TEXT_OPERATION = { "-", "+", "*", "/" };

  private JFrame myParentFrame;
  private csSeisPaneManager myManager;

  private Thread myThread;
  private boolean myStopReadDataThread;
  private csIDataSetListener myListener;
  private csWellPathOverlay myWellPathOverlay;
  
  // Used to temporarily store currently applied selection parameters. May need to revert to these in case
  // new trace selection is incorrect
  private csTraceSelectionParam myTempSavedSelectParam;
  private csTraceSelectionParam mySelectParam;
  private int myMoveOption;
  private boolean myResetScalar;
  
  private csISeismicReader myReader;
  private csFilename myFilename;
  public JSplitPane splitPane;
  public csSeisPane seisPane;
  public csSeisView seisView;
  public csGraph2D graphPanel;
  public csColorBarPanel colorBarPanel;
  public csSeaViewPopupMenu popupMenu;
  public String title;
  private ArrayList<csIGraphPanelListener> myGraphPanelListeners;
  private csHeaderMonitor myHeaderMonitor;
  private csTraceMonitor myTraceMonitor;
  private csGraphSelectionDialog myGraphDialog = null;
  private csOverlayDialog myOverlayDialog = null;
  
  private csISeismicTraceBuffer mySeismicTraceBuffer;
  private csHeaderDef[] myTraceHeaderDef;
  private csHeaderDef[] mySortedTraceHeaderDef;
  private csAnnotationDialog myAnnotationDialog;
  private ArrayList<csIProcessing> myProcessingSteps;
  private csProcessingOverlay myProcessingOverlay;
  private csPaintOverlay myPaintOverlay;
  private csPaintDialog myPaintDialog;
  private csPickingDialog myPickingDialog;
  private csPickOverlay myPickOverlay;

  private int myHdrIndexTimeSamp1;
  private int myHdrIndexTimeSamp1_us;
  private int myHdrIndexWavenumber;
  private int myHdrIndexKNorm;
  private boolean myIsFrequencyDomain;
  /// File format: CSEIS, SEGY etc
  private int myFileFormat;
  
  private int myCurrentSelectedHdrIndex;
  private int myCurrentScanHdrInfoIndex;
  private boolean myDoRefresh;
  private boolean myIsSynced;
  private ArrayList<csScanHeaderInfo> myScanHdrInfoList;
  private csTraceHeaderScan myTraceHeaderScan;
  private ArrayList<csISeisPaneBundleListener> mySeisPaneBundleListeners;
  private ArrayList<csISeisPaneBundleMouseModeListener> myMouseModeListeners;
  private ArrayList<csISeisPaneBundleScrollListener> mySeisPaneBundleScrollListeners;
  public final int id;
  private JFileChooser myFileChooser;
  private int myLocalMouseMode;
  private int myGlobalMouseMode;
  
  csSeisPaneBundle( JFrame parentFrame, csSeisPaneManager manager, csSeisView seisViewIn,
          csTraceSelectionParam traceSelectionParam, csFilename filename, int fileFormat ) {
    super( new BorderLayout() );
    id = ++id_counter;
    myParentFrame = parentFrame;
    myManager = manager;
    myMoveOption = SeaView.MOVE_OPTION_SELECTION;
    myResetScalar = true;
    myReader = null;
    myFilename = filename;
    title = myFilename.filename;
    myFileFormat = fileFormat;
    seisView = seisViewIn;
    seisView.getEventHandler().addSampleInfoListener( this );
    seisView.getEventHandler().addRubberBandListener(this);
    myScanHdrInfoList = new ArrayList<csScanHeaderInfo>();
    myHdrIndexWavenumber = -1;
    myHdrIndexKNorm = -1;
    myLocalMouseMode = csMouseModes.NO_MODE;
    myGlobalMouseMode = csMouseModes.NO_MODE;

    myCurrentScanHdrInfoIndex = -1;
    myIsSynced = false;
    myProcessingSteps = new ArrayList<csIProcessing>();
    myProcessingOverlay = new csProcessingOverlay("");
    myPaintOverlay = null;
    myPickOverlay  = null;
    myPaintDialog = null;
    myPickingDialog = null;
    myWellPathOverlay = null;
    myTraceHeaderScan = null;

    myHdrIndexTimeSamp1 = -1;
    myHdrIndexTimeSamp1_us = -1;
    
    mySeisPaneBundleListeners = new ArrayList<csISeisPaneBundleListener>(1);
    myMouseModeListeners = new ArrayList<csISeisPaneBundleMouseModeListener>(1);
    mySeisPaneBundleScrollListeners = new ArrayList<csISeisPaneBundleScrollListener>(1);
    myGraphPanelListeners = new ArrayList<csIGraphPanelListener>(1);
    mySeismicTraceBuffer = null;
    
    mySelectParam = new csTraceSelectionParam();
    if( traceSelectionParam != null ) {
//      mySelectParam.selectOption    = traceSelectionParam.selectOption;
//      mySelectParam.firstTraceIndex = traceSelectionParam.firstTraceIndex;
      mySelectParam.numTraces       = traceSelectionParam.numTraces;
      mySelectParam.traceStep       = traceSelectionParam.traceStep;
//      mySelectParam.firstEnsIndex   = traceSelectionParam.firstEnsIndex;
      mySelectParam.selectedHdrName = traceSelectionParam.selectedHdrName;
      mySelectParam.numEns          = traceSelectionParam.numEns;
      mySelectParam.selectedHdrIndex = traceSelectionParam.selectedHdrIndex;
      mySelectParam.selectedHdrValue   = traceSelectionParam.selectedHdrValue;
      mySelectParam.numTracesHdr       = traceSelectionParam.numTracesHdr;
    }
    myTempSavedSelectParam = mySelectParam;
    
    seisView.addSeisViewListener( new csISeisViewListener() {
      @Override
      public void traceBufferChanged( csISeismicTraceBuffer traceBuffer ) {
        if( myPickOverlay != null ) myPickOverlay.traceBufferChanged( traceBuffer );
      }
      @Override
      public void changedSettings( csSeisDispSettings settings ) {
        fireUpdateScalarEvent();
        fireUpdateDispSettingsEvent();
        setColorBar( settings );
      }
      @Override
      public void vertScrollChanged( int scrollValue ) {
        fireVerticalScrollEvent( scrollValue );
      }
      @Override
      public void horzScrollChanged( int scrollValue ) {
        if( graphPanel != null && graphPanel.isVisible() ) {
          graphPanel.resetViewPositionHorz(scrollValue);
        }
        fireHorizontalScrollEvent( scrollValue );
      }
      @Override
      public void sizeChanged( Dimension size ) {
        if( graphPanel != null && graphPanel.isVisible() ) {
          resetGraphSize();
          graphPanel.resetViewPositionHorz(seisPane.getHorizontalScrollBar().getValue());
        }
      }
      @Override
      public void mouseExited() {
        myManager.mouseExitedBundle( csSeisPaneBundle.this );
      }
    });
  }
  public void setProcessingStep( String name ) {
    csIProcessing proc = retrieveProcessingStep( name );
    if( proc == null ) {
      if( name.compareTo( csProcessingAGC.NAME ) == 0 ) {
        proc = new csProcessingAGC( 500, 0, getSampleInt() );
      }
      else if( name.compareTo( csProcessingDCRemoval.NAME ) == 0 ) {
        proc = new csProcessingDCRemoval();
      }
      else if( name.compareTo( csProcessingFilter.NAME ) == 0 ) {
        proc = new csProcessingFilter( 2, 12, 220, 200, getNumSamples(), getSampleInt() );
      }
      else if( name.compareTo( csProcessingInterpolation.NAME ) == 0 ) {
        proc = new csProcessingInterpolation( 2, getSampleInt() );
      }
      else if( name.compareTo( csProcessingTest.NAME ) == 0 ) {
        proc = new csProcessingTest();
      }
      else {
        JOptionPane.showMessageDialog( myParentFrame,
                "No processing step found with name " + name + "\n" +
                "This error is likely due to a program bug.",
                "Error", JOptionPane.ERROR_MESSAGE );
        return;
      }
      myProcessingSteps.add( proc );
    }
    proc.setActive( false );
    setupProcessing( proc );
  }
  private void setupProcessing( csIProcessing proc ) {
    final csProcessingDialog dialog = new csProcessingDialog( myParentFrame );
    dialog.setProcessingStep( proc, new csIProcessingSetupListener() {
      @Override
      public void setupComplete( csIProcessing proc ) {
        for( int i = 0; i < myProcessingSteps.size(); i++ ) {
          myProcessingSteps.get(i).setActive( false );
        }
        proc.setActive( true );
        myProcessingOverlay.setName( proc.getName() );
        Thread thread = new Thread( new Runnable() {
          @Override
          public void run() {
            try {
              ((csTraceBuffer)mySeismicTraceBuffer).clearProcessing();
              seisView.addOverlay( myProcessingOverlay );
              applyProcessingToData( true );
            }
            catch( Exception e ) {
              e.printStackTrace();
            }
            dialog.stopProcessing();
          }
        });
        dialog.startProcessing();
        thread.start();
      }
    });
    dialog.setVisible(true);
  }
  private csIProcessing retrieveProcessingStep( String name ) {
    for( int i = 0; i < myProcessingSteps.size(); i++ ) {
      csIProcessing proc = myProcessingSteps.get(i);
      if( proc.getName().compareTo( name ) == 0 ) {
        return proc;
      }
    }
    return null;
  }
  public void clearAllProcessingSteps() {
    if( myProcessingSteps.size() > 0 ) {
      if( ((csTraceBuffer)mySeismicTraceBuffer).isProcessed() ) {
        ((csTraceBuffer)mySeismicTraceBuffer).clearProcessing();
        seisView.removeOverlay( myProcessingOverlay );
        seisView.updateTraceBuffer( mySeismicTraceBuffer, true );
      }
      for( int i = 0; i < myProcessingSteps.size(); i++ ) {
        myProcessingSteps.get(i).setActive( false );
      }
    }
  }
  public boolean startTraceHeaderScan( csHeaderDef headerDef ) {
    myTraceHeaderScan = new csTraceHeaderScan( myFilename.filenamePath, headerDef, this, this );
    boolean success = myTraceHeaderScan.start( myReader );
    if( !success ) {
      myTraceHeaderScan.stopOperation();
      myListener.stopScan( this, false );
    }
    return success;
  }
  public void cancelTraceHeaderScan() {
    if( myTraceHeaderScan != null && myTraceHeaderScan.isRunning() ) {
      myTraceHeaderScan.stopOperation();
      myTraceHeaderScan = null;
      myListener.stopScan( this, false );
    }
  }
  public void setSyncState( boolean isSync ) {
    myIsSynced = isSync;
  }
  public boolean isSynced() {
    return myIsSynced;
  }
  @Override
  public void scanCompleted( csHeaderDef headerDef ) {
    for( int i = 0; i < myScanHdrInfoList.size(); i++ ) {
      if( headerDef.name.compareTo( myScanHdrInfoList.get(i).getHeaderName() ) == 0 ) {
        myScanHdrInfoList.remove(i);
      }
    }
    // Add latest scanned header at top of list
    if( myTraceHeaderScan != null ) {
      myScanHdrInfoList.add( 0, new csScanHeaderInfo( headerDef, myTraceHeaderScan.retrieveSelectedValues() ) );
      myListener.stopScan( this, true );
      myTraceHeaderScan = null;
    }
  }
  public boolean isHeaderScanned( String hdrName ) {
    for( int i = 0; i < myScanHdrInfoList.size(); i++ ) {
      if( hdrName.compareTo(myScanHdrInfoList.get(i).getHeaderName()) == 0 ) return true;
    }
    return false;
  }
  public int getNumScannedHeaders() {
    return myScanHdrInfoList.size();
  }
//  public csHeader getScanHeaderValue( int scanHeaderIndex, int scanTraceIndex ) {
//    return myScanHdrInfoList.get(scanHeaderIndex).hdrValues[scanTraceIndex];
//  }
  /*
  public csHeader getScanHeaderValue( String hdrName, int scanTraceIndex ) {
    csScanHeaderInfo info = getScanHeaderInfo( hdrName );
    if( info == null || info.bundle.hdrValues == null || scanTraceIndex >= info.bundle.hdrValues.length ) return null;
    return info.bundle.hdrValues[scanTraceIndex];
  }
   * 
   */
  public csScanHeaderInfo getScanHeaderInfo( String hdrName ) {
    for( int i = 0; i < myScanHdrInfoList.size(); i++ ) {
      if( hdrName.compareTo(myScanHdrInfoList.get(i).getHeaderName()) == 0 ) {
        return myScanHdrInfoList.get(i);
      }
    }
    return null;
  }
  public csScanHeaderInfo getLatestScanHeaderInfo() {
    return( myScanHdrInfoList.isEmpty() ? null : myScanHdrInfoList.get(0) );
  }
  public csScanHeaderInfo getScanHeaderInfo( int scanHeaderIndex ) {
    if( scanHeaderIndex >= getNumScannedHeaders() ) return null;
    return myScanHdrInfoList.get(scanHeaderIndex);
  }
  public void fitToScreen() {
    seisPane.fitToScreen();
  }
  public int getDisplayedNumTraces() {
    return mySeismicTraceBuffer.numTraces();
  }
  public int getTotalNumTraces() {
    if( myReader != null ) {
      return myReader.numTraces();
    }
    else {
      return 0;
    }
  }
  public boolean isScanning() {
    return( myTraceHeaderScan != null && myTraceHeaderScan.isRunning() );
  }
  public boolean isSameAs( csSeisPaneBundle bundle ) {
    return( bundle != null && getFilenamePath().compareTo(bundle.getFilenamePath()) == 0 );
  }
  public boolean isFrequencyDomain() {
    return myIsFrequencyDomain;
  }
  public int getFileFormat() {
    return myFileFormat;
  }
  public String getFilenamePath() {
    return myFilename.filenamePath;
  }
  public String getPath() {
    return myFilename.path;
  }
  public void showSegyCharHdr() {
    String charHdr = ((csNativeSegyReader)myReader).charHeader();
    String text = "";
    for( int i = 0; i < 3200; i += 80 ) {
      text += charHdr.substring(i, i+80) + "\n";
    }
    JDialog dialog = new csSegyHeaderView( myParentFrame, "SEG-Y char header", text );
    dialog.setVisible(true);
  }
  public void showSegyBinHdr() {
    JDialog dialog = new csSegyHeaderView( myParentFrame, "SEG-Y binary header", ((csNativeSegyReader)myReader).binHeader().headerDump() );
    dialog.setVisible(true);
  }
  public void addSeisPaneBundleListener( csISeisPaneBundleListener listener ) {
    mySeisPaneBundleListeners.add( listener );
  } 
  public void removeSeisPaneBundleListener( csISeisPaneBundleListener listener ) {
    mySeisPaneBundleListeners.remove( listener );
  } 
  public void addSeisPaneBundleMouseModeListener( csISeisPaneBundleMouseModeListener listener ) {
    myMouseModeListeners.add( listener );
  } 
  public void removeSeisPaneBundleMouseModeListener( csISeisPaneBundleMouseModeListener listener ) {
    myMouseModeListeners.remove( listener );
  } 
  public void addSeisPaneBundleScrollListener( csISeisPaneBundleScrollListener listener ) {
    mySeisPaneBundleScrollListeners.add( listener );
  } 
  public void removeSeisPaneScrollBundleListener( csISeisPaneBundleScrollListener listener ) {
    mySeisPaneBundleScrollListeners.remove( listener );
  } 
  public String getTitle() {
    return title;
  }
  @Override
  public String toString() {
    return getTitle();
  }
  public void setColorBar( csSeisDispSettings settings ) {
    if( colorBarPanel != null && colorBarPanel.isVisible() ) {
      if( settings.isVIDisplay ) {
        colorBarPanel.setColorMap( settings.viColorMap );
      }
      else if( settings.isVariableColor ) {
        colorBarPanel.setColorMap( settings.wiggleColorMap );
      }
      else {
        colorBarPanel.setColorMap( new csCustomColorMap( Color.black ) );
      }
    }
  }
  public void setTraceToZero( int traceIndex ) {
    if( traceIndex >= 0 && traceIndex < mySeismicTraceBuffer.numTraces() ) {
      float[] samples = mySeismicTraceBuffer.samples(traceIndex);
      for( int isamp = 0; isamp < samples.length; isamp++ ) {
        samples[isamp] = 0.0f;
      }
      boolean isHighlight = seisView.isHighlightOn();
      seisPane.setAutoHighlightTrace( !isHighlight );
      seisPane.setAutoHighlightTrace( isHighlight );
    }
  }
  public void updateScalar( float scalar ) {
    seisPane.setScalar( scalar );
  }
  public void incScaling() {
    seisPane.setScalar( (float)(1.5*(double)seisPane.getScalar()) );
  }
  public void decScaling() {
    seisPane.setScalar( (float)((double)seisPane.getScalar()/1.5) );
  }
  public void showSettingsDialog() {
    seisPane.showSettingsDialog();
  }
  public void showOverlayDialog() {
    if( myOverlayDialog == null ) {
      myOverlayDialog = new csOverlayDialog( myParentFrame, seisView, myTraceHeaderDef );
    }
    myOverlayDialog.setVisible( true );
  }
  @Override
  public void closePaintDialog() {
//    if( myPaintDialog != null ) myPaintDialog.setVisible(false);
//    if( myPaintOverlay != null ) myPaintOverlay.deactivate();
//    if( seisView.getMouseMode() == csMouseModes.PAINT_MODE ) {
//      seisView.setMouseMode( csMouseModes.NO_MODE );
//    }
    setLocalMouseMode( csMouseModes.NO_MODE );
//    popupMenu.setPaintMode( false );
//    setPaintMode( false );
  }
  @Override
  public void clearPaintOverlay() {
    if( myPaintOverlay != null ) myPaintOverlay.clearOverlay();
  }
  @Override
  public void updatePaintOverlay( int lineSize, Color lineColor, int pointSize, Color pointColor ) {
    if( myPaintOverlay != null ) {
      myPaintOverlay.update( lineSize, lineColor, pointSize, pointColor );
    }
  }
  @Override
  public void closePickDialog() {
    setLocalMouseMode( csMouseModes.NO_MODE );
//    popupMenu.setPickMode( false );
//    setPickMode( false );
  }
  @Override
  public void updatePicks() {
    seisView.repaint();
  }
  @Override
  public void setPicksFromHeader( csHeaderDef hdrDef ) {
    int headerIndex = getTraceHeaderIndex( hdrDef.name );
    float[] picks = new float[seisView.getTraceBuffer().numTraces()];
    csISeismicTraceBuffer traceBuffer = seisView.getTraceBuffer();
    for( int itrc = 0; itrc < picks.length; itrc++ ) {
      picks[itrc] = traceBuffer.headerValues(itrc)[headerIndex].floatValue();
    }
    myPickOverlay.setActivePicks( picks );
    seisView.repaint();
  }
  @Override
  public void savePicks( csHorizon horizon, String name ) {
    if( myFileChooser == null ) myFileChooser = new JFileChooser();
    myFileChooser.setSelectedFile( new File(name + ".txt") );
    int option = myFileChooser.showSaveDialog( myPickingDialog );
    if( option == JFileChooser.APPROVE_OPTION ) {
      File file = myFileChooser.getSelectedFile();
      double sampleInt = myPickOverlay.getSampleInt();
      try {
        BufferedWriter writer = new BufferedWriter( new FileWriter(file) );
        String text = "";
        DecimalFormat floatFormat = new DecimalFormat("0.0");
        horizon.startIteration();
        while( horizon.hasNext() ) {
          csHorizonPickArray array = (csHorizonPickArray)horizon.next();
          float[] picks = array.picks;
          for( int itrc = 0; itrc < picks.length; itrc++ ) {
            float sampleIndex = picks[itrc];
            if( sampleIndex != csPickOverlay.NO_VALUE ) {
              text += " " + (itrc+array.trace1) + " " + floatFormat.format(sampleInt*sampleIndex) + "\n";
            }
          }
        }
        writer.write( text );
        writer.close();
        JOptionPane.showMessageDialog(myParentFrame, "Saved horizon picks to file\n" + file.getAbsolutePath(), "Info", JOptionPane.INFORMATION_MESSAGE);
      }
      catch (IOException ex) {
        JOptionPane.showMessageDialog(myParentFrame, "Error occurred when writing time picks:\n" + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
      }
    }
  }
  public void setMouseMode( int mouseMode ) {
    myGlobalMouseMode = mouseMode;
    if( myLocalMouseMode == csMouseModes.NO_MODE ) {
      seisView.setMouseMode(myGlobalMouseMode);
      fireUpdateMouseModeEvent( myGlobalMouseMode );
    }
    else {
      // Nothing to be done. Local mouse mode overrides global mouse mode
//      enableLocalMouseMode();
//      seisView.setMouseMode(myLocalMouseMode);
//      fireUpdateMouseModeEvent( myLocalMouseMode );
    }
  }
  private void disableLocalMouseMode() {
    if( myLocalMouseMode == csMouseModes.PAINT_MODE ) setPaintMode( false );
    else if( myLocalMouseMode == csMouseModes.PICK_MODE ) setPickMode( false );
  }
  private void enableLocalMouseMode() {
    if( myLocalMouseMode == csMouseModes.PAINT_MODE ) setPaintMode( true );
    else if( myLocalMouseMode == csMouseModes.PICK_MODE ) setPickMode( true);
  }
  public void setLocalMouseMode( int newLocalMouseMode ) {
    disableLocalMouseMode();
    myLocalMouseMode = newLocalMouseMode;
    enableLocalMouseMode();

    if( myLocalMouseMode == csMouseModes.NO_MODE ) {
      seisView.setMouseMode(myGlobalMouseMode);
      fireUpdateMouseModeEvent( myGlobalMouseMode );
    }
    else {
      seisView.setMouseMode(myLocalMouseMode);
      fireUpdateMouseModeEvent( myLocalMouseMode );
    }
  }
  private void setPaintMode( boolean doSet ) {
    if( doSet ) {
      if( myPaintOverlay == null ) myPaintOverlay = new csPaintOverlay( seisView );
      seisView.addOverlay( myPaintOverlay );
      myPaintOverlay.activate();
//      setLocalMouseMode( csMouseModes.PAINT_MODE );
      if( myPaintDialog == null ) myPaintDialog = new csPaintDialog( myParentFrame, this );
      myPaintDialog.setVisible(true);
    }
    else {
      if( myPaintDialog != null ) myPaintDialog.setVisible(false);
      if( myPaintOverlay != null ) myPaintOverlay.deactivate();
//      setLocalMouseMode( csMouseModes.NO_MODE );
    }
  }
  private void setPickMode( boolean doSet ) {
    if( doSet ) {
      if( myPickOverlay == null ) myPickOverlay = new csPickOverlay( seisView );
      seisView.addOverlay( myPickOverlay );
      myPickOverlay.activate();
//      setLocalMouseMode( csMouseModes.PICK_MODE );
      if( myPickingDialog == null ) myPickingDialog = new csPickingDialog( myParentFrame, myPickOverlay, title, this, mySortedTraceHeaderDef );
      myPickingDialog.setVisible(true);
    }
    else {
      if( myPickingDialog != null ) myPickingDialog.setVisible(false);
      if( myPickOverlay != null ) myPickOverlay.deactivate();
//      setLocalMouseMode( csMouseModes.NO_MODE );
    }
  }
  public void updateWellPaths( ArrayList<csWellPathAttr> wellPathAttrList, csWellPathOverlayAttr wellPathOverlayAttr ) {
    boolean resetHdrs = true;
    if( myWellPathOverlay == null ) {
      myWellPathOverlay = new csWellPathOverlay( seisView );
      seisView.addOverlay( myWellPathOverlay );
      resetHdrs = true;
    }
    else {
      csWellPathOverlayAttr wellPathOverlayAttr_old = myWellPathOverlay.getOverlayAttr();
      resetHdrs = ( wellPathOverlayAttr_old.hdrName_binx.compareTo(wellPathOverlayAttr.hdrName_binx) != 0 ) ||
        ( wellPathOverlayAttr_old.hdrName_biny.compareTo(wellPathOverlayAttr.hdrName_biny) != 0 );
    }
    if( resetHdrs ) {
      int hdrIndex_x = getTraceHeaderIndex( wellPathOverlayAttr.hdrName_binx );
      int hdrIndex_y = getTraceHeaderIndex( wellPathOverlayAttr.hdrName_biny );
      if( hdrIndex_x < 0 || hdrIndex_y < 0 ) {
        JOptionPane.showMessageDialog( myParentFrame,
            "Trace headers " + wellPathOverlayAttr.hdrName_binx + " and/or " + wellPathOverlayAttr.hdrName_biny + " not found in data set\n" +
            myFilename.filename + "\n" +
            "Well path(s) cannot be displayed.",
            "Trace header not found",
            JOptionPane.ERROR_MESSAGE );
        return;
      }
      myWellPathOverlay.updateHeaders( hdrIndex_x, hdrIndex_y );
    }
    myWellPathOverlay.update( wellPathOverlayAttr, wellPathAttrList );
  }
  public float getSampleInt() {
    if( myReader == null ) return 0;
    return myReader.sampleInt();
  }
  public int getNumSamples() {
    if( mySeismicTraceBuffer == null ) return 0;
    return mySeismicTraceBuffer.numSamples();
  }
  public csISeismicReader readImage( String filename ) {
    FileImageInputStream input = null;
    try {
      input = new FileImageInputStream( new File(filename) );
    } 
    catch (FileNotFoundException e1) {
      e1.printStackTrace();
      return null;
    }
    catch (IOException e2) {
      e2.printStackTrace();
      return null;
    }
    Iterator<ImageReader> iter = ImageIO.getImageReadersByFormatName("JPG");
    if( iter.hasNext() ) {
      ImageReader imageReader = iter.next();
      imageReader.setInput( input );
      BufferedImage image = null;
      try {
        image = imageReader.read(0);
      } catch (IOException e) { e.printStackTrace(); System.exit(-1); }
      int numTraces   = image.getWidth();
      int numSamples  = image.getHeight();
      float[][] imageDataArray = new float[numTraces][numSamples];
      for( int itrc = 0; itrc < numTraces; itrc++ ) {
        for( int isamp = 0; isamp < numSamples; isamp++ ) {
          imageDataArray[itrc][isamp] = (float)image.getRGB( itrc, isamp );
        }
      }
      int numHeaders  = 1;
      float sampleInt = 1.0f;
      myTraceHeaderDef = new csHeaderDef[numHeaders];
      myTraceHeaderDef[0] = new csHeaderDef("trace","Trace number",csHeaderDef.TYPE_INT);
      csVirtualSeismicReader reader = new csVirtualSeismicReader(
            numSamples, myTraceHeaderDef.length, sampleInt, myTraceHeaderDef, myReader.verticalDomain() );
      csTraceBuffer traceBuffer = reader.retrieveTraceBuffer();
      for( int itrc = 0; itrc < numTraces; itrc++ ) {
        csHeader[] headerValues = new csHeader[numHeaders];
        headerValues[0] = new csHeader(itrc+1);
        traceBuffer.addTrace( new csSeismicTrace( imageDataArray[itrc], headerValues, itrc+1 ) );
      }
      return reader;
    }
    else {
      try {
		input.close();
      } catch (IOException e) {
	    e.printStackTrace();
      }
      return null;
    }
  }
//----------------------------------------------------------
  public csISeismicReader combineDataSets(
          csSeisPaneBundle bundle2,
          int operation,
          float multiplier,
          float addValue )
  {
    float sampleInt = myReader.sampleInt();
    int numSamples = Math.min( myReader.numSamples(), bundle2.myReader.numSamples() );
    int numTraces = Math.min( getDisplayedNumTraces(), bundle2.getDisplayedNumTraces() );
    int numHeaders = myReader.numHeaders();
    csVirtualSeismicReader reader = new csVirtualSeismicReader(
            numSamples, myTraceHeaderDef.length, sampleInt, myTraceHeaderDef, myReader.verticalDomain() );
    csTraceBuffer traceBuffer = reader.retrieveTraceBuffer();
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      float[] samples1   = mySeismicTraceBuffer.samples(itrc);
      float[] samples2   = bundle2.mySeismicTraceBuffer.samples(itrc);
      float[] samplesOut = new float[numSamples];
      csHeader[] headerValuesOut = new csHeader[numHeaders];
      csHeader[] headerValuesIn = mySeismicTraceBuffer.headerValues(itrc);
      for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
        headerValuesOut[ihdr] = new csHeader( headerValuesIn[ihdr] );
      }
      if( operation == OPERATION_SUBTRACT ) {
        for( int isamp = 0; isamp < numSamples; isamp++ ) {
          samplesOut[isamp] = samples1[isamp] - samples2[isamp] * multiplier + addValue;
        }
      }
      else if( operation == OPERATION_ADD ) {
        for( int isamp = 0; isamp < numSamples; isamp++ ) {
          samplesOut[isamp] = samples1[isamp] + samples2[isamp] * multiplier + addValue;
        }
      }
      else if( operation == OPERATION_MULTIPLY ) {
        for( int isamp = 0; isamp < numSamples; isamp++ ) {
          samplesOut[isamp] = samples1[isamp] * samples2[isamp] * multiplier + addValue;
        }
      }
      else if( operation == OPERATION_DIVIDE ) {
        for( int isamp = 0; isamp < numSamples; isamp++ ) {
          float value2 = samples2[isamp] * multiplier + addValue;
          if( value2 != 0.0 ) {
            samplesOut[isamp] = samples1[isamp] / value2;
          }
          else {
            samplesOut[isamp] = 0.0f;
          }
        }
      }
      traceBuffer.addTrace( samplesOut, headerValuesOut );
    }
    return reader;
  }
  public int getVerticalDomain() {
    return myReader.verticalDomain();
  }
  public csHeaderDef[] getHeaderDef() {
    return myTraceHeaderDef;
  }
  public void zoom( int zoomType, int zoomMode ) {
    seisPane.zoom( zoomType, zoomMode );
  }
  public void testColorMap( csCustomColorMap cmap ) {
    csSeisDispSettings ds = seisView.getDispSettings();
    ds.wiggleColorMap.resetColors( cmap );
    ds.viColorMap.resetColors( cmap );

    seisView.updateDispSettings(ds);
  }
  //----------------------------------------------------------
  public void showTraceMonitor() {
    if( mySeismicTraceBuffer == null ) return;
    if( myTraceMonitor == null ) {
      myTraceMonitor = new csTraceMonitor( myParentFrame, mySeismicTraceBuffer, (float)seisView.getSampleInt(),
                          "Trace monitor " + getTitle() );
      myTraceMonitor.addWindowListener(new WindowAdapter() {
      @Override
        public void windowClosing( WindowEvent e ) {
          csTraceMonitor traceMonitor = (csTraceMonitor)e.getSource();
          seisView.getEventHandler().removeSampleInfoListener( traceMonitor );
        }
      });
      seisView.getEventHandler().addSampleInfoListener( myTraceMonitor );
    }
    else if( !myTraceMonitor.isVisible() ) {
      myTraceMonitor.updateBuffer( mySeismicTraceBuffer );
    }
    myTraceMonitor.setVisible( true );
  }
  public void showHeaderMonitor() {
    if( myTraceHeaderDef == null ) return;
    if( myHeaderMonitor == null ) {
      myHeaderMonitor = new csHeaderMonitor( myParentFrame, myTraceHeaderDef, 
            "Header monitor " + getTitle() );
    }
    myHeaderMonitor.setVisible( true );
  }
  public void resetGraphSize( ) {
    csGraphAttributes attr = graphPanel.getGraphAttributes();
    Dimension sizeSeisPane = seisPane.getPreferredSize();
    Dimension sizeSeisView = seisView.getPreferredSize();
    attr.fixHorzGraphPosition = true;
    attr.horzGraphPositionMin = sizeSeisPane.width - sizeSeisView.width + seisView.getMarginLeftRight() - 15;
    attr.horzGraphPositionMax = attr.horzGraphPositionMin + sizeSeisView.width - 2*seisView.getMarginLeftRight();
// should be:    attr.horzGraphPositionMax = sizeSeisPane.width - mySeisPane.getVerticalScrollBar().getPreferredSize().width;
    graphPanel.setFixedDim(true, false, sizeSeisPane.width, sizeSeisView.width);
// !CHANGE! review settings above. No -15 ...
  }
  public void updateGraph() {
    if( graphPanel == null ) return;
    java.util.Set<Integer> curveIDs = graphPanel.getCurveIDs();
    java.util.Iterator<Integer> iter = curveIDs.iterator();
    while( iter.hasNext() ) {
      int curveId = iter.next().intValue();
      if( curveId >= myTraceHeaderDef.length ) {
        graphPanel.removeCurve(curveId);
      }
      else {
        float[] valuesX = new float[mySeismicTraceBuffer.numTraces()];
        float[] valuesY = new float[mySeismicTraceBuffer.numTraces()];
        for( int itrc = 0; itrc < mySeismicTraceBuffer.numTraces(); itrc++ ) {
          valuesX[itrc] = itrc + 1 + mySelectParam.firstTraceIndex;  // !CHANGE! CHECK firstTraceIndex!
          valuesY[itrc] = (float)(mySeismicTraceBuffer.headerValues( itrc )[curveId].doubleValue());
        }
        graphPanel.updateCurve(curveId, valuesX, valuesY,true);
      }
    }
//    myGraphPanel.repaintAll();
 }
  @Override
  public void showGraph( boolean doShow ) {
    showGraphPanel( doShow );
  }
  public void showGraphPanel( boolean doShow ) {
    if( graphPanel == null ) {
      graphPanel = new csGraph2D();
      myGraphDialog = new csGraphSelectionDialog( myParentFrame, this, graphPanel, myTraceHeaderDef );
      myGraphPanelListeners = new ArrayList<csIGraphPanelListener>(1);
      myGraphPanelListeners.add( this );
      csGraphAttributes attr = graphPanel.getGraphAttributes();
      attr.isGraphViewPort = true;
      attr.insetBottom = 4;
      attr.insetRight  = seisPane.getVerticalScrollBar().getPreferredSize().width;
      attr.insetTop    = 0;
      attr.insetLeft   = 10;
      attr.borderPadding = 0;
      attr.showBorder = false;
      attr.showAxisXAnnotation = false;   // !CHANGE! Doesn't work with X axis annotation. Maybe due to 10^exponent.
      attr.showAxisYAnnotation = true;
      attr.showZeroAxis = false;
      attr.title = "";
      attr.xLabel = "";
      attr.showInnerBorder = false;
      attr.graphPadding    = 10;
      attr.fixHorzGraphPosition = true;
      attr.gridYPixelInc = 25;
      attr.gridXPixelInc = 40;
      resetGraphSize();
    }
    if( doShow ) {
      splitPane.setLeftComponent( graphPanel );
      splitPane.setDividerLocation(0.2);
      splitPane.setDividerSize( 8 );
    }
    else {
      splitPane.setLeftComponent( null );
      splitPane.setDividerSize( 0 );
    }
    splitPane.setOneTouchExpandable( doShow );
    myGraphDialog.setVisible( doShow );
  }
  public void addBackgroundImage() {
    if( seisView instanceof csBackgroundSeisView ) {
      // Test: Load same data into background:
      ((csBackgroundSeisView)seisView).updateBackground( mySeismicTraceBuffer, seisView.getSampleInt() );
    }
    else {
      JOptionPane.showMessageDialog( myParentFrame,
        "Background functionality is not supported for this type of data\n",
        "Unsupported option",
        JOptionPane.ERROR_MESSAGE );
    }
  }
  public void showAnnotationDialog() {
    createAnnotationDialog();
    myAnnotationDialog.setVisible( true );
  }
  private void createAnnotationDialog() {
    if( myAnnotationDialog == null ) {
      myAnnotationDialog = new csAnnotationDialog( myParentFrame, seisPane, myTraceHeaderDef );
      myAnnotationDialog.changeSettings( myManager.getAnnotationAttributes() );
    }
  }
  public void retrieveAnnotationSettings( csAnnotationAttributes annAttr ) {
    if( myAnnotationDialog != null ) myAnnotationDialog.getSettings( annAttr );
  }
  public csTraceSelectionParam getTraceSelectionParam() {
    return( new csTraceSelectionParam( mySelectParam ) );
  }
  public csCurve createCurve( String hdrName, csCurveAttributes attr ) {
    int hdrIndex = getTraceHeaderIndex( hdrName );
    if( hdrIndex < 0 ) return null;
    float[] valuesX = new float[mySeismicTraceBuffer.numTraces()];
    float[] valuesY = new float[mySeismicTraceBuffer.numTraces()];
    for( int itrc = 0; itrc < mySeismicTraceBuffer.numTraces(); itrc++ ) {
      valuesX[itrc] = itrc + 1 + mySelectParam.firstTraceIndex;    // !CHANGE! CHECK firstTraceIndex
      valuesY[itrc] = (float)(mySeismicTraceBuffer.headerValues( itrc )[hdrIndex].doubleValue());
    }
    return new csCurve(valuesX,valuesY,attr,hdrIndex);
  }
  private int getTraceHeaderIndex( String hdrName ) {
    if( mySeismicTraceBuffer.numTraces() == 0 ) return -1;
    for( int ihdr = 0; ihdr < myTraceHeaderDef.length; ihdr++ ) {
      if( myTraceHeaderDef[ihdr].name.compareTo( hdrName ) == 0 )  {
        return ihdr;
      }
    }
    return -1;
  }
  private void computeSpectrum( csSampleInfo posStart, csSampleInfo posEnd ) {
    int minSample = Math.min(posEnd.sample, posStart.sample);
    int maxSample = Math.max(posEnd.sample, posStart.sample);
    int minTrace = Math.min(posEnd.trace, posStart.trace);
    int maxTrace = Math.max(posEnd.trace, posStart.trace);
    int numTraces = maxTrace - minTrace + 1;
    int numSamples = maxSample - minSample + 1;
    if( numSamples < 10 ) {
      JOptionPane.showMessageDialog( myParentFrame,
        "Selected time window is too small,\n" +
        "no spectrum computed.\n" +
        "Please select a larger window (at least 10 samples).",
        "Window too small",
        JOptionPane.ERROR_MESSAGE );
      return;
    }
    csNativeFFTTransform transform = new csNativeFFTTransform(numSamples);

    float[] samples = mySeismicTraceBuffer.samples(minTrace);
    float[] samplesSelect = new float[numSamples];
    System.arraycopy(samples, minSample, samplesSelect, 0, numSamples);

    csFFTObject fftObject = transform.performFFT( samplesSelect, csNativeFFTTransform.FORWARD, (float)seisView.getSampleInt() );

    for( int itrc = minTrace+1; itrc <= maxTrace; itrc++ ) {
      samples = mySeismicTraceBuffer.samples(itrc);
      System.arraycopy(samples, minSample, samplesSelect, 0, numSamples);
      csFFTObject fftObjectTmp = transform.performFFT( samplesSelect, csNativeFFTTransform.FORWARD, (float)seisView.getSampleInt() );
      for( int i = 0; i < fftObject.numSamples; i++ ) {
        fftObject.amplitude[i] += fftObjectTmp.amplitude[i];
        fftObject.phase[i]     += fftObjectTmp.phase[i];
      }
    }
    for( int i = 0; i < fftObject.numSamples; i++ ) {
      fftObject.amplitude[i] /= numTraces;
      fftObject.phase[i]     /= numTraces;
    }
    double maxAmp = SeaView.convertAmp2DB( fftObject.amplitude );
    String text = "Amplitude spectrum, Time window " + (int)(minSample*seisView.getSampleInt()) + "-" + (int)(maxSample*seisView.getSampleInt()) +
            "ms, Traces " + (minTrace+1) + "-" + (maxTrace+1);
    
    JDialog dialog = new JDialog(myParentFrame,"Mean spectrum");
    dialog.getContentPane().setLayout(new BorderLayout());
    DecimalFormat format = new DecimalFormat("0.0000E0");
    csSpectrumGraph spectrum = new csSpectrumGraph(text,"Frequency [Hz]","Amplitude [dB relative to " + format.format(maxAmp) + "]");
    spectrum.update( fftObject );
    dialog.getContentPane().add( spectrum );
    dialog.pack();
    dialog.setVisible(true);
  }
  public void computeSpectrumOneTrace() {
    int numSamplesAll = mySeismicTraceBuffer.numTraces()*mySeismicTraceBuffer.numSamples();
    csNativeFFTTransform transform = new csNativeFFTTransform(numSamplesAll);
    float[] samples = new float[numSamplesAll];

    for( int itrc = 1; itrc < mySeismicTraceBuffer.numTraces(); itrc++ ) {
      float[] samples_tmp = mySeismicTraceBuffer.samples(itrc);
      int isamp = itrc*mySeismicTraceBuffer.numSamples();
      System.arraycopy(samples_tmp, 0, samples, isamp, mySeismicTraceBuffer.numSamples());
    }    
    csFFTObject fftObject = transform.performFFT( samples, csNativeFFTTransform.FORWARD, (float)seisView.getSampleInt() );
    double maxAmp = SeaView.convertAmp2DB( fftObject.amplitude );

    JDialog dialog = new JDialog(myParentFrame,"Spectrum - one long trace");
    DecimalFormat format = new DecimalFormat("0.0000E0");
    csSpectrumGraph spectrum = new csSpectrumGraph("Amplitude spectrum","Frequency [Hz]","Amplitude [dB relative to " + format.format(maxAmp) + "]");
    spectrum.update( fftObject );
    dialog.getContentPane().add( spectrum );
    dialog.pack();
    dialog.setVisible(true);
  }
  private void applyProcessingToData( boolean refresh ) {
    if( myProcessingSteps.size() > 0 ) {
      boolean hasProcessed = false;
      myProcessingOverlay.setName("");
      for( int i = 0; i < myProcessingSteps.size(); i++ ) {
        csIProcessing proc = myProcessingSteps.get(i);
        if( proc.isActive() ) {
          ((csTraceBuffer)mySeismicTraceBuffer).applyProcessing( proc );
          hasProcessed = true;
          myProcessingOverlay.addProcess( proc.getName() );
        }
      }
      if( refresh && hasProcessed ) {
        seisView.updateTraceBuffer( mySeismicTraceBuffer, true );
        seisView.repaint();
      }
    }
  }
  public float computeFKVelocity( float traceDouble, float frequency, float traceSpacing ) {
    float velocity = -1.0f;
    if( myHdrIndexKNorm >= 0 ) {
      int trace1 = Math.max( 0, (int)traceDouble );
      int trace2 = trace1 + 1;
      if( trace2 >= mySeismicTraceBuffer.numTraces() ) {
        trace2 = mySeismicTraceBuffer.numTraces()-1;
        trace1 = Math.max( 0, trace2 - 1 );
      }
      csHeader[] hdrValues1 = mySeismicTraceBuffer.headerValues( trace1 );
      csHeader[] hdrValues2 = mySeismicTraceBuffer.headerValues( trace2 );
      if( myHdrIndexKNorm < hdrValues1.length ) {
        float knorm1 = hdrValues1[myHdrIndexKNorm].floatValue();
        float knorm2 = hdrValues2[myHdrIndexKNorm].floatValue();
        float knormdiff = knorm2 - knorm1;
        float knorm = knorm1 + knormdiff * ( traceDouble - (float)trace1 );
        float k = Math.abs( knorm / (2.0f * traceSpacing) );
        velocity = frequency / Math.max( k, 1e-7f );
      }
    }
    else if( myHdrIndexWavenumber >= 0 ) {
      if( traceDouble < 0 ) traceDouble = 0.0f;
      csHeader[] hdrValues0 = mySeismicTraceBuffer.headerValues( 0 );
      csHeader[] hdrValues = mySeismicTraceBuffer.headerValues( (int)traceDouble );
      if( myHdrIndexWavenumber < hdrValues.length ) {
        float wavenumber = hdrValues[myHdrIndexWavenumber].floatValue() + ( traceDouble - (float)(int)traceDouble );
        velocity = Math.abs( 2.0f * traceSpacing * frequency * ( hdrValues0[myHdrIndexWavenumber].floatValue() / wavenumber ) );
      }
    }
    return velocity;
  }
 @Override
  public void rubberBandCompleted( csSampleInfo posStart, csSampleInfo posEnd ) {
    if( seisView.getMouseMode() == csMouseModes.ZOOM_MODE ) {
      seisPane.zoomArea( posStart, posEnd );
    }
    else if( seisView.getMouseMode() == csMouseModes.SPECTRUM_MODE ) {
      computeSpectrum(posStart,posEnd);
    }
  }
  @Override
  public void mouseMoved( Object source, csSampleInfo info_in ) {
    csSeisBundleSampleInfo info_out = new csSeisBundleSampleInfo( info_in );
    info_out.zoomVert = seisView.getZoomVert();
    info_out.originalTrace = mySeismicTraceBuffer.numTraces() <= info_in.trace ? 0 : mySeismicTraceBuffer.originalTraceNumber(info_in.trace);
    info_out.isFrequencyDomain = myIsFrequencyDomain;
    info_out.timeFullSample = info_in.sample*seisView.getSampleInt()/1000.0;
    if( myIsFrequencyDomain ) {
      info_out.info.time *= 1000;
      info_out.timeFullSample *= 1000.0;
    }
    csHeader[] hdrValues = mySeismicTraceBuffer.headerValues( info_in.trace );
    info_out.time_samp1_ms = 0;
    if( myHdrIndexTimeSamp1 >= 0 && myHdrIndexTimeSamp1 < hdrValues.length ) {
      info_out.time_samp1_ms = ( (long)hdrValues[myHdrIndexTimeSamp1].intValue() ) * 1000;
    }
    if( myHdrIndexTimeSamp1_us >= 0 && myHdrIndexTimeSamp1_us < hdrValues.length ) {
      info_out.time_samp1_ms += ( (long)hdrValues[myHdrIndexTimeSamp1_us].intValue() ) / 1000;
    }
    fireUpdateSampleEvent( info_out );
    if( info_in.trace >= 0 && info_in.trace < mySeismicTraceBuffer.numTraces() ) {
      if( myHeaderMonitor != null && myHeaderMonitor.isVisible() ) {
        myHeaderMonitor.updateValues( mySeismicTraceBuffer.headerValues( info_in.trace ) );
      }
    }
  }
  @Override
  public void mouseClicked( Object source, csSampleInfo info ) {
    if( seisView.getMouseMode() == csMouseModes.KILL_MODE ) {
      setTraceToZero( info.trace );
    }
    else if( seisView.getMouseMode() == csMouseModes.PICK_MODE ) {
      // ???
/*      csHeader[] hdrValues = mySeismicTraceBuffer.headerValues( info.trace );
      int time_samp1 = 0;
      for( int ihdr = 0; ihdr < myTraceHeaderDef.length; ihdr++ ) {
        if( myTraceHeaderDef[ihdr].name.compareTo("time_samp1") == 0 ) {
          time_samp1 = hdrValues[ihdr].intValue();
          break;
        }
      }
      System.out.println( " " + info.trace + " " + (double)((int)(1000000*info.time))/1000.0 + " " + time_samp1 );
*/
    }
  }
  private void fireVerticalScrollEvent( int value ) {
    for( int i = 0; i < mySeisPaneBundleScrollListeners.size(); i++ ) {
      mySeisPaneBundleScrollListeners.get(i).vertScrollChanged( this, value );
    }
  }
  private void fireHorizontalScrollEvent( int value ) {
    for( int i = 0; i < mySeisPaneBundleScrollListeners.size(); i++ ) {
      mySeisPaneBundleScrollListeners.get(i).horzScrollChanged( this, value );
    }
  }
//  private void fireMouseExitedEvent( int value ) {
//    for( int i = 0; i < mySeisPaneBundleScrollListeners.size(); i++ ) {
//      mySeisPaneBundleScrollListeners.get(i).horzScrollChanged( this, value );
//    }
//  }
  public void fireUpdateScalarEvent() {
    for( int i = 0; i < mySeisPaneBundleListeners.size(); i++ ) {
      mySeisPaneBundleListeners.get(i).updateBundleDisplayScalar( this );
    }
  }
   public void fireUpdateDispSettingsEvent() {
    for( int i = 0; i < mySeisPaneBundleListeners.size(); i++ ) {
      mySeisPaneBundleListeners.get(i).updateBundleDisplaySettings( this );
    }

  }
  public void fireUpdateSampleEvent( csSeisBundleSampleInfo info ) {
    for( int i = 0; i < mySeisPaneBundleListeners.size(); i++ ) {
      mySeisPaneBundleListeners.get(i).updateBundleSampleInfo( this, info );
    }
  }
  public void fireUpdateMouseModeEvent( int mouseMode ) {
    for( int i = 0; i < myMouseModeListeners.size(); i++ ) {
      myMouseModeListeners.get(i).updateBundleMouseMode( mouseMode );
    }
  }

  public void close() {
    if( myReader != null ) {
      myReader.closeFile();
      if( myReader instanceof csNativeSegyReader ) {
        ((csNativeSegyReader)myReader).disposeNative();
      }
      else if(myReader instanceof csNativeSegdReader) {
        ((csNativeSegdReader)myReader).disposeNative();
      }
      else if(myReader instanceof csNativeRSFReader) {
        ((csNativeRSFReader)myReader).disposeNative();
      }
      myReader = null;
    }
  }
  public int getOriginalTraceNumber( int traceIndex ) {
    return mySeismicTraceBuffer.originalTraceNumber( traceIndex );
  }
  public boolean displayedDataAtStart() {
    if( mySeismicTraceBuffer == null || getDisplayedNumTraces() <= 0 ) return( true );
    boolean isAtStart;
    if( mySelectParam.selectOption == csTraceSelectionParam.SELECT_TRACE ) {
      isAtStart = ( mySeismicTraceBuffer.originalTraceNumber(0) <= mySelectParam.traceStep );
    }
    else if( mySelectParam.selectOption == csTraceSelectionParam.SELECT_HEADER ) {
      if( myCurrentScanHdrInfoIndex < 0 || myCurrentScanHdrInfoIndex >= myScanHdrInfoList.size() ) return true;
      csScanHeaderInfo scanInfo = myScanHdrInfoList.get(myCurrentScanHdrInfoIndex );
      int traceIndexFirst = mySeismicTraceBuffer.originalTraceNumber(0) - 1;
      isAtStart = traceIndexFirst <= scanInfo.getCurrentHdrFirstTraceIndex();
    }
    else {
      isAtStart = ( mySeismicTraceBuffer.originalTraceNumber(0) == 1 );
    }
    return isAtStart;
  }
  public boolean displayedDataAtEnd() {
    if( mySeismicTraceBuffer == null || getDisplayedNumTraces() <= 0 ) return( true );
    boolean isAtEnd;
    if( mySelectParam.selectOption == csTraceSelectionParam.SELECT_TRACE ) {
      isAtEnd = ( mySeismicTraceBuffer.originalTraceNumber(getDisplayedNumTraces()-1) >= getTotalNumTraces()-(mySelectParam.traceStep-1) );
    }
    else if( mySelectParam.selectOption == csTraceSelectionParam.SELECT_HEADER ) {
      if( myCurrentScanHdrInfoIndex < 0 || myCurrentScanHdrInfoIndex >= myScanHdrInfoList.size() ) return true;
      csScanHeaderInfo scanInfo = myScanHdrInfoList.get(myCurrentScanHdrInfoIndex );
      int traceIndexLast  = mySeismicTraceBuffer.originalTraceNumber(mySeismicTraceBuffer.numTraces()-1) - 1;
      isAtEnd = traceIndexLast >= scanInfo.getCurrentHdrLastTraceIndex();
    }
    else {
      isAtEnd = ( mySeismicTraceBuffer.originalTraceNumber(getDisplayedNumTraces()-1) >= getTotalNumTraces() );
    }
    return isAtEnd;
  }
  public boolean hasRandomFileAccess() {
    if( myReader == null ) return false;
    return myReader.hasRandomFileAccess();
  }
  public void updateTraceSelection( csTraceSelectionParam param ) {
    myTempSavedSelectParam = mySelectParam;
    mySelectParam = new csTraceSelectionParam( param );
    if( mySelectParam.selectOption == csTraceSelectionParam.SELECT_ENS || mySelectParam.selectOption == csTraceSelectionParam.SELECT_HEADER ) {
      for( int ihdr = 0; ihdr < myTraceHeaderDef.length; ihdr++ ) {
        if( myTraceHeaderDef[ihdr].name.compareTo( mySelectParam.selectedHdrName ) == 0 ) {
          myCurrentSelectedHdrIndex = ihdr;
          break;
        }
      }
    }
    myCurrentScanHdrInfoIndex = -1;
    if( mySelectParam.selectOption != csTraceSelectionParam.SELECT_TRACE ) {
      for( int i = 0; i < myScanHdrInfoList.size(); i++ ) {
        if( mySelectParam.selectedHdrName.compareTo( myScanHdrInfoList.get(i).getHeaderName() ) == 0 ) {
          myCurrentScanHdrInfoIndex = i;
          break;
        }
      }
    }
  }
//  public void read() {
//    read( false, myListener, SeaView.MOVE_OPTION_SELECTION );
//  }
  public synchronized void read( boolean refresh, csIDataSetListener listener, int moveOption, boolean resetScalar ) {
    myListener    = listener;
    myMoveOption  = moveOption;
    myResetScalar = resetScalar;
    myDoRefresh   = refresh;
    myThread = new Thread( new Runnable() {
      @Override
      public void run() {
        int numHeaders = myReader.numHeaders();
        int numSamples = myReader.numSamples();
        int numTracesToBuffer = Math.max(mySelectParam.numTracesHdr/10,10);
        csScanHeaderInfo scanInfo = null;
        if( myCurrentScanHdrInfoIndex >= 0 ) {
          scanInfo = myScanHdrInfoList.get( myCurrentScanHdrInfoIndex );
        }
        csSeismicTrace trace = new csSeismicTrace( numSamples, numHeaders );
        boolean success = true;
        try {
          csTraceBuffer newTraceBuffer = new csTraceBuffer( myReader.numSamples(), numHeaders );
          myStopReadDataThread = false;
          int firstTraceIndex = 0;
          int lastTraceIndex = 0;
          if( mySeismicTraceBuffer != null && mySeismicTraceBuffer.numTraces() > 0 ) {
            firstTraceIndex = mySeismicTraceBuffer.originalTraceNumber(0) - 1;
            lastTraceIndex = mySeismicTraceBuffer.originalTraceNumber( mySeismicTraceBuffer.numTraces()-1 ) - 1;
          }
          //==============================================================================================
          if( mySelectParam.selectOption == csTraceSelectionParam.SELECT_TRACE ) {
            int traceIndex1 = 0;
            int minTraceIndex = firstTraceIndex % mySelectParam.traceStep;
            if( myMoveOption == SeaView.MOVE_OPTION_SELECTION ) {
              traceIndex1 = mySelectParam.firstTraceIndex;
            }
            else if( myMoveOption == SeaView.MOVE_OPTION_FORWARD ) {
              traceIndex1 = Math.min( getTotalNumTraces()-minTraceIndex-1, lastTraceIndex + mySelectParam.traceStep );
            }
            else if( myMoveOption == SeaView.MOVE_OPTION_BACKWARD ) {
              traceIndex1 = Math.max( minTraceIndex, firstTraceIndex - mySelectParam.numTraces * mySelectParam.traceStep );
            }
            else if( myMoveOption == SeaView.MOVE_OPTION_BEGIN ) {
              // traceIndex1 = 0; Nothing to be done
            }
            else if( myMoveOption == SeaView.MOVE_OPTION_END ) {
              traceIndex1 = Math.max( minTraceIndex, getTotalNumTraces() - mySelectParam.numTraces * mySelectParam.traceStep );
            }
            int counterDisplayedTraces = 0;  // Counter of traces that are actually going to be displayed
            int counterTraceSteps = mySelectParam.traceStep; // Make sure that first trace is read in
            if( hasRandomFileAccess() ) myReader.moveToTrace( traceIndex1, Math.max(mySelectParam.numTraces/10,10) );
            while( counterDisplayedTraces < mySelectParam.numTraces && myReader.getNextTrace( trace ) && !myStopReadDataThread ) {
              if( counterDisplayedTraces % REFRESH_FREQUENCY == 0 ) myListener.updateTrace( counterDisplayedTraces, traceIndex1 );
              traceIndex1 += 1;
              if( counterTraceSteps < mySelectParam.traceStep ) {
                counterTraceSteps += 1;
                continue;
              }
              trace.setOriginalTraceNumber( traceIndex1 );  // traceIndex1 has just been incremented by 1, making it the trace "number"
              newTraceBuffer.addTrace( trace );
              trace = new csSeismicTrace( numSamples, numHeaders );
              counterDisplayedTraces += 1;
              counterTraceSteps = 1;
            }
          } //
          //==============================================================================================
          else if( mySelectParam.selectOption == csTraceSelectionParam.SELECT_ENS ) {
            // 1) Determine if currently selected header has been scanned:
            int traceIndex1 = 0;
            // *** Ensemble selection ***
            // START if moveoption = backward or end AND data has not been scanned
            if( scanInfo == null && (myMoveOption == SeaView.MOVE_OPTION_BACKWARD || myMoveOption == SeaView.MOVE_OPTION_END) ) {
              if( myMoveOption == SeaView.MOVE_OPTION_END ) {
                traceIndex1 = getTotalNumTraces()-1;
              }
              else {
                traceIndex1 = mySeismicTraceBuffer.originalTraceNumber(0) - 2; // -2 to convert from trace number to index
              }
              int counterReadTraces = 0;
              int counterDisplayedEns = 0;
              csHeader peekValue = new csHeader();
              myReader.setHeaderToPeek( mySelectParam.selectedHdrName );
              myReader.peekHeaderValue( traceIndex1, peekValue );
              myListener.updateTrace( 0, traceIndex1 );
              while( counterDisplayedEns < mySelectParam.numEns && !myStopReadDataThread && traceIndex1 > 0 ) {
                // Step (1): Read in maximal numTraceToBuffer traces, backwards..
                int numTracesToRead = Math.min( numTracesToBuffer, traceIndex1 + 1 );
                traceIndex1 -= (numTracesToRead-1);
                csSeismicTrace[] tempTraces = new csSeismicTrace[numTracesToRead];
                myReader.moveToTrace( traceIndex1, numTracesToRead );
                for( int itrc = 0; itrc < numTracesToRead; itrc++ ) {
                  if( !myReader.getNextTrace( trace ) ) {
                    numTracesToRead = itrc;
                    break;
                  }
                  trace.setOriginalTraceNumber( traceIndex1+itrc+1 );
                  tempTraces[itrc] = trace;
                  trace = new csSeismicTrace( numSamples, numHeaders );
                }
                // Step (2): Check ensemble header value. Save to trace buffer. Stop when all ensembles have been read
                for( int itrc = numTracesToRead-1; itrc >= 0; itrc-- ) {
                  csHeader ensValue = tempTraces[itrc].headerValues()[myCurrentSelectedHdrIndex];
                  if( !ensValue.equals(peekValue) ) {
                    peekValue = ensValue;
                    counterDisplayedEns += 1;
                    if( counterDisplayedEns == mySelectParam.numEns ) {
                      break;
                    }
                  }
                  else { //if( tempTraces[itrc].headerValues()[myCurrentSelectedHdrIndex].longValue() == peekValue.value ) {
                    counterReadTraces += 1;
                    if( counterReadTraces % REFRESH_FREQUENCY == 0 ) myListener.updateTrace( counterDisplayedEns, traceIndex1 );
                    newTraceBuffer.addTraceAtStart( tempTraces[itrc] );
                  }
                }
                traceIndex1 -= 1;
                counterReadTraces += 1;
              }
            } // END if moveoption = end or backward
            //-------------------------------------------------------------------
            // *** Ensemble selection ***
            // START if moveoption = begin, forward or selection, or if scanInfo is set (=header has been scanned)
            else {
              if( myMoveOption == SeaView.MOVE_OPTION_SELECTION ) {
                if( scanInfo != null ) {
                  traceIndex1 = scanInfo.getEnsFirstTraceIndex( mySelectParam.firstEnsIndex );
                }
                else {
                  traceIndex1 = firstTraceIndex;
                }
              }
              else if( myMoveOption == SeaView.MOVE_OPTION_FORWARD ) {
                if( scanInfo != null ) {
                  mySelectParam.firstEnsIndex = Math.min( mySelectParam.firstEnsIndex+mySelectParam.numEns, scanInfo.getNumEnsembles()-1 );
                  traceIndex1 = scanInfo.getEnsFirstTraceIndex( mySelectParam.firstEnsIndex );
                }
                else {
                  traceIndex1 = lastTraceIndex + 1;
                }
              }
              else if( myMoveOption == SeaView.MOVE_OPTION_BEGIN ) {
                traceIndex1 = 0;
                if( scanInfo != null ) {
                  mySelectParam.firstEnsIndex = 0;
                }
              }
              // Next two options are possible only if trace header has been scanned:
              else if( myMoveOption == SeaView.MOVE_OPTION_BACKWARD ) {
                mySelectParam.firstEnsIndex = Math.max( mySelectParam.firstEnsIndex-mySelectParam.numEns, 0 );
                traceIndex1 = scanInfo.getEnsFirstTraceIndex( mySelectParam.firstEnsIndex );
              }
              else if( myMoveOption == SeaView.MOVE_OPTION_END ) {
                mySelectParam.firstEnsIndex = scanInfo.getNumEnsembles()-1;
                traceIndex1 = scanInfo.getEnsFirstTraceIndex( mySelectParam.firstEnsIndex );
              }
              int counterReadTraces = 0;
              int counterDisplayedEns = 0;
              myReader.moveToTrace( traceIndex1, numTracesToBuffer );
              myReader.setHeaderToPeek( mySelectParam.selectedHdrName );
              csHeader peekValue = new csHeader();
              myReader.peekHeaderValue( traceIndex1, peekValue );
              myListener.updateTrace( 0, traceIndex1 );
              while( counterDisplayedEns < mySelectParam.numEns && myReader.getNextTrace( trace ) && !myStopReadDataThread ) {
                if( counterReadTraces % REFRESH_FREQUENCY == 0 ) myListener.updateTrace( counterDisplayedEns, traceIndex1+counterReadTraces );
                csHeader ensValue = trace.headerValues()[myCurrentSelectedHdrIndex];
                counterReadTraces += 1;
                if( !ensValue.equals(peekValue) ) {
                  peekValue = ensValue;
                  counterDisplayedEns += 1;
                  if( counterDisplayedEns == mySelectParam.numEns ) {
                    break;
                  }
                  continue;
                }
                trace.setOriginalTraceNumber( traceIndex1+counterReadTraces );
                newTraceBuffer.addTrace( trace );
                trace = new csSeismicTrace( numSamples, numHeaders );
              }
//              myCurrentTraceIndex += counterReadTraces - 1;
            } // END if moveoption = begin, forward or selection
          } // END if ensemble trace selection
          //========================================================================================================
          else if( mySelectParam.selectOption == csTraceSelectionParam.SELECT_HEADER ) {
            int traceIndex1    = 0;
            boolean moveForwards = true;
            if( myMoveOption == SeaView.MOVE_OPTION_SELECTION ) {
              // Nothing to be done
            }
            else if( myMoveOption == SeaView.MOVE_OPTION_FORWARD ) {
              traceIndex1 = Math.min( lastTraceIndex+1, getTotalNumTraces()-1 );
            }
            else if( myMoveOption == SeaView.MOVE_OPTION_BACKWARD ) {
              traceIndex1 = Math.max( 0, firstTraceIndex - 1 );
              moveForwards = false;
            }
            else if( myMoveOption == SeaView.MOVE_OPTION_BEGIN ) {
              // Nothing to be done
            }
            else if( myMoveOption == SeaView.MOVE_OPTION_END ) {
              traceIndex1 = getTotalNumTraces() - 1;
              moveForwards = false;
            }
            if( scanInfo != null ) {
              scanInfo.reset_getNextHeaderTraceIndex( mySelectParam.selectedHdrValue, traceIndex1, moveForwards );
              if( moveForwards ) {
                myReader.moveToTrace( traceIndex1, numTracesToBuffer );
              }
              else {
                myReader.moveToTrace( traceIndex1, 1 );
              }
              traceIndex1 = -999;
              int nextTraceIndex = -999;
              int counterDisplayedTraces = 0;  // Counter of traces that are actually going to be displayed
              while( counterDisplayedTraces < mySelectParam.numTracesHdr &&
                     !myStopReadDataThread &&
                     ( nextTraceIndex = scanInfo.getNextHeaderTraceIndex() ) >= 0 ) {
                if( counterDisplayedTraces % REFRESH_FREQUENCY == 0 ) myListener.updateTrace( counterDisplayedTraces, traceIndex1 );
                if( nextTraceIndex-traceIndex1 != 1 ) myReader.moveToTrace( nextTraceIndex, numTracesToBuffer ); 
                traceIndex1 = nextTraceIndex;
                if( !myReader.getNextTrace( trace ) ) break;
                trace.setOriginalTraceNumber( traceIndex1+1 );
                if( moveForwards ) newTraceBuffer.addTrace( trace );
                else newTraceBuffer.addTraceAtFront( trace );
                trace = new csSeismicTrace( numSamples, numHeaders );
                counterDisplayedTraces += 1;
              }
            }
            //-------------------------------------------------------
            // *** Header selection, not scanned ***
            else {
              myReader.moveToTrace( traceIndex1, 1 );
              csHeader peekValue = new csHeader();
              peekValue.setValue( 1.0 );
              myReader.setHeaderToPeek( mySelectParam.selectedHdrName );
              int counterDisplayedTraces = 0;  // Counter of traces that are actually going to be displayed
              int counter = 0; // Simple counter;
              int prevTraceIndex = traceIndex1;
              boolean performMoveOperation = false;
              while( counterDisplayedTraces < mySelectParam.numTracesHdr && !myStopReadDataThread ) {
                if( counter % REFRESH_FREQUENCY == 0 ) myListener.updateTrace( counterDisplayedTraces, traceIndex1 );
                counter += 1;
                myReader.peekHeaderValue( traceIndex1, peekValue );
                if( peekValue.doubleValue() == mySelectParam.selectedHdrValue ) {
                  if( performMoveOperation ) {
                    if( !myReader.moveToTrace( traceIndex1, 1 ) ) break;
                    performMoveOperation = false;
                    prevTraceIndex = traceIndex1;
                  }
                  if( !myReader.getNextTrace( trace ) ) break;
                  trace.setOriginalTraceNumber( traceIndex1+1 );
                  if( moveForwards ) newTraceBuffer.addTrace( trace );
                  else newTraceBuffer.addTraceAtFront( trace );
                  trace = new csSeismicTrace( numSamples, numHeaders );
                  counterDisplayedTraces += 1;
                  if( counterDisplayedTraces == mySelectParam.numTracesHdr ) {
                    break;
                  }
                }
                if( moveForwards ) {
                  traceIndex1 += 1;
                  if( traceIndex1 == getTotalNumTraces() ) break;
                  if( !performMoveOperation && traceIndex1 > prevTraceIndex+1 ) performMoveOperation = true;
                }
                else {
                  traceIndex1 -= 1;
                  if( traceIndex1 < 0 ) break;
                  if( !performMoveOperation && traceIndex1 < prevTraceIndex-1 ) performMoveOperation = true;
                  myReader.moveToTrace( traceIndex1, 1 );
                }
              }
              
            }
          }
          //========================================================================================================

          if( newTraceBuffer.numTraces() > 0 ) { // if( counterDisplayedTraces > 0 ) {
            mySelectParam.firstTraceIndex = newTraceBuffer.originalTraceNumber(0) - 1;
            float sampleInt = myReader.sampleInt();
            if( mySeismicTraceBuffer != null ) mySeismicTraceBuffer.clear();
            mySeismicTraceBuffer = newTraceBuffer;
            applyProcessingToData( false );
            if( myTraceMonitor != null && myTraceMonitor.isVisible() ) myTraceMonitor.updateBuffer( mySeismicTraceBuffer );
              // Read trace header names, descriptions & types
            if( myDoRefresh ) {
              myTraceHeaderDef = new csHeaderDef[numHeaders];
              myHdrIndexTimeSamp1    = -1;
              myHdrIndexTimeSamp1_us = -1;
              myHdrIndexWavenumber   = -1;
              for( int i = 0; i < numHeaders; i++ ) {
                myTraceHeaderDef[i] = new csHeaderDef( myReader.headerName(i), myReader.headerDesc(i), myReader.headerType(i) );
                if( myReader.headerName(i).equals("time_samp1") ) {
                  myHdrIndexTimeSamp1 = i;
                }
                else if( myReader.headerName(i).equals("time_samp1_us") ) {
                  myHdrIndexTimeSamp1_us = i;
                }
                else if( myReader.headerName(i).equals("k-wavenumber") || myReader.headerName(i).equals("k_wavenumber") ) {
                  myHdrIndexWavenumber = i;
                }              
              }
              mySortedTraceHeaderDef = new csHeaderDef[myTraceHeaderDef.length];
              for( int ihdr = 0; ihdr < myTraceHeaderDef.length; ihdr++ ) {
                mySortedTraceHeaderDef[ihdr] = new csHeaderDef( myTraceHeaderDef[ihdr] );
              }
              java.util.Arrays.sort( mySortedTraceHeaderDef );
            } // END myDoRefresh
            // Update csSeisPane:
            seisPane.updateSeismic( mySeismicTraceBuffer, sampleInt, myTraceHeaderDef, mySelectParam.firstTraceIndex, myResetScalar );
            seisPane.setVerticalDomain( myReader.verticalDomain() );
            if( seisView.isLogScale() ) { // Trick to get log scale properly initialized:
              seisPane.zoom( seisView.getZoomVert(), seisView.getZoomHorz() );
            }
            seisView.repaint();  // !CHANGE! review if this call is always necessary for proper refresh

            if( myDoRefresh ) {
              if( myHeaderMonitor != null ) {
                myHeaderMonitor.updateHeaderNames( myTraceHeaderDef );
              }
              if( myAnnotationDialog != null ) {
                myAnnotationDialog.updateHeaderNames( myTraceHeaderDef );
              }
              if( myGraphDialog != null ) {
                myGraphDialog.removeAllCurves();
              }
              if( myOverlayDialog != null ) {
                myOverlayDialog.setTraceHeaders( myTraceHeaderDef );
              }
            }
            updateGraph();
          }
          else {
            JOptionPane.showMessageDialog( myParentFrame,
                "No traces read in.\n" +
                "Please check your trace selection settings under Tools->Trace Selection.\n",
                "Error message",
                JOptionPane.ERROR_MESSAGE );
            mySelectParam = myTempSavedSelectParam;
//            myManager.removeBundle( csSeisPaneBundle.this );
            success = false;
          }
        }
        catch( Exception e ) {
          JOptionPane.showMessageDialog( myParentFrame,
              "Error occurred when reading next trace from input file " + title + " . System message:\n\n" +
              e.getMessage(),
              "Error",
              JOptionPane.ERROR_MESSAGE );
            mySelectParam = myTempSavedSelectParam;
//          myManager.removeBundle( csSeisPaneBundle.this );
          success = false;
        }
        if( myGraphDialog != null ) {
          javax.swing.SwingUtilities.invokeLater( new Runnable() {
            @Override
            public void run() {
              myGraphDialog.setTraceHeaders( myTraceHeaderDef );
              resetGraphSize();
            }
          });
        }
        myListener.isComplete( success );
        myStopReadDataThread = true;
        myThread = null;
      }
    });
    myThread.start();
  }
  @Override
  public void traceHeaderScanNotify( int traceIndex ) {
    myListener.updateScan( this, traceIndex );
  }
  @Override
  public boolean traceHeaderScanContinue() {
    return( myTraceHeaderScan != null );
  }
  @Override
  public void stopOperation() {
    myStopReadDataThread = true;
  }
  public csHeaderDef[] getSortedTraceHeaderDef() {
    return mySortedTraceHeaderDef;
  }
  public void initializeReader( csISeismicReader reader, String filenamePath, int fileFormat ) {
    myReader = reader;
    myFilename.setFilenamePath( filenamePath );
    title = myFilename.filename;
    myIsFrequencyDomain = myReader.isFrequencyDomain();
    myFileFormat = fileFormat;
  }
  public void setShowFilenameInView( boolean showFilename ) {
    if( showFilename ) {
      seisView.addOverlay( new csSeisViewFilenameOverlay(getTitle()) );
    }
    else {
      int numOverlays = seisView.getNumOverlays();
      for( int ilay = 0; ilay < numOverlays; ilay++ ) {
        csISeisOverlay overlay = seisView.getOverlay(ilay);
        if( overlay instanceof csSeisViewFilenameOverlay ) {
          seisView.removeOverlay(overlay);
          break;
        }
      } // END for ilay
    } // END else
  }
  public void resetFilename( String filename ) {
    title = filename;
    int numOverlays = seisView.getNumOverlays();
    for( int ilay = 0; ilay < numOverlays; ilay++ ) {
      csISeisOverlay overlay = seisView.getOverlay(ilay);
      if( overlay instanceof csSeisViewFilenameOverlay ) {
        ((csSeisViewFilenameOverlay)overlay).setFilename( getTitle() );
        seisView.repaint();
        break;
      }
    } // END for ilay
  }
}
