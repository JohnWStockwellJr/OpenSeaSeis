/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.velanal.csVelocityAnalysisManager;
import cseis.cmapgen.csColorMapGenerator;
import cseis.cmapgen.csIColorMapGeneratorListener;
import cseis.cmapgen.csCustomColorMapModel;
import cseis.general.csCustomColorMap;
import cseis.tools.csProgressBarWindow;
import cseis.tools.csSnapShotFrame;
import cseis.general.csColorMap;
import cseis.general.csFilename;
import cseis.general.csStandard;
import cseis.jni.*;
import cseis.segy.csISegyAttrListener;
import cseis.segy.csSegyAttr;
import cseis.segy.csSegySetupDialog;
import cseis.seis.*;
import cseis.seisdisp.*;
import cseis.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.stream.FileImageOutputStream;
import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * SeaView - 2D Seismic Viewer<br>
 * <br>
 * This is the main class for the 2D seismic viewer.<br>
 * The 2D seismic viewer application uses Java CSeisLib and SeaSeis C++ libraries to load and display seismic data.
 * <br><br>
 * An image file with the same seismic display can be created in the background using class PlotImage.
 * PlotImage uses the same graphics engine and display settings as SeaView.
 * 
 * @author 2008 Bjorn Olofsson
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class SeaView extends JFrame implements csFileMenuListener {
  public static final String VERSION = "3.00";
  public static final String PROPERTIES_FILE_NAME = "seisdisp_settings_"+VERSION+".txt";
  public static final String PROPERTIES_FREQ_FILE_NAME = "seisdisp_settings_freq_"+VERSION+".txt";

  public static final int MOVE_OPTION_SELECTION = 40;
  public static final int MOVE_OPTION_FORWARD  = 41;
  public static final int MOVE_OPTION_BACKWARD = 42;
  public static final int MOVE_OPTION_BEGIN    = 43;
  public static final int MOVE_OPTION_END      = 44;

  public static final int FORMAT_ALL   = 0;
  public static final int FORMAT_SEGY  = 1;
  public static final int FORMAT_CSEIS = 2;
  public static final int FORMAT_SU    = 3;
  public static final int FORMAT_ASCII = 4;
  public static final int FORMAT_SEGD  = 5;
  public static final int FORMAT_RSF   = 6;
  public static final int FORMAT_JPG   = 7;
  
  /// Tool bar on top of seismic pane
  private csSeaViewToolBarTop  myToolBarTop;
  /// Tool bar to the left of seismic pane
  private csSeaViewToolBarLeft myToolBarLeft;
  private csMenuBar myMenuBar;
  /// Status bar displaying information messages
  private JLabel myStatusBar;
  private csSegyAttr mySegyAttr;
  /// Manages dock panes (dock panes contain seismic panes)
  private csDockPaneManager myDockPaneManager;
  private csDockPaneButtonSelection myDockPaneButtonSelection;
  /// Manages seismic bundles
  private csSeisPaneManager mySeisPaneManager;
  private csDockPane myCurrentNewDockPane;
  private csPreferences myPreferences;  // should really be deprecated
  private csSeaViewProperties myProperties;
  /// true if this instance of SeaView is a child process of another Java application
  private boolean myIsChildProcess;
  private csSnapShotFrame mySnapShotFrame;
  /// Manages velocity analysis
  private csVelocityAnalysisManager myVelocityAnalysisManager;

  // Progress bar which is shown when process is ongoing, such as reading in new file
  private csProgressBarWindow myProgressBar;
  private csTraceSelectionDialog myTraceSelectionDialog;
  private csWellPathDialog myWellPathDialog;
  private boolean myIsReadProcessOngoing = false;
  private boolean myIsScalarLocked = false;
  private boolean myOpenInNewPanel = true;
  
  private JFileChooser mySeismicFileChooser = null;
  private csFileChooser myFileChooserDump = null;
  private int myCurrentSelectedFileFormat;
  
  private csCustomColorMapModel myCustomColorMapModel;
  private ArrayList<String> myInputFilenameQueue;
//  private boolean myIsImageExport = false;
//  private int myExportHeight;
//  private int myExportWidth;
  
  private boolean myShowFKVelocity;
  private float myFKTraceSpacing;
  private boolean myShowAbsoluteTime;
  private boolean myIsRefreshFileOperation = false;
  private int myColorBitType;
  private int myMouseMode;

  protected static DecimalFormat floatFmt4 = new DecimalFormat("0.0000");
  protected static DecimalFormat floatFmt3 = new DecimalFormat("0.000");
  protected static DecimalFormat floatFmt2 = new DecimalFormat("0.00");
  protected static DecimalFormat floatFmt1 = new DecimalFormat("0.0");
  protected static DecimalFormat floatFmt0 = new DecimalFormat("0");
  protected static DecimalFormat formatAmplitude = new DecimalFormat( "0.00000000" );

  /**
   * Default constructor
   */
  public SeaView() {
    this(false);
  }
  /**
   * Constructor
   * @param isChildProcess Pass 'true' if SeaView is started as a Java child process and hence,
   *        when SeaView is closed the Java system should not shut down
   */
  public SeaView( boolean isChildProcess ) {
    super("SeaView, version " + VERSION);
    myIsChildProcess = isChildProcess;
    myColorBitType = csColorMap.COLOR_MAP_TYPE_32BIT;
    myMouseMode = csMouseModes.NO_MODE;
    myCurrentSelectedFileFormat = SeaView.FORMAT_ALL;

    myCustomColorMapModel = new csCustomColorMapModel( csColorMap.COLOR_MAP_TYPE_32BIT );

    myPreferences = new csPreferences();
    myPreferences.readPreferences();
    myProperties = new csSeaViewProperties();
    try {
      myProperties.load();
    }
    catch (IOException ex) {
    }
    csTraceSelectionParam.NUM_TRACES = myProperties.numTraces;
    myToolBarTop  = new csSeaViewToolBarTop( this, myProperties );
    myToolBarLeft = new csSeaViewToolBarLeft( this );
    myToolBarTop.setToolbarEnabled( false );
    myToolBarTop.setMoveButtonsEnabled( true, true );
    myToolBarLeft.setToolbarEnabled( false );
    
    myCustomColorMapModel.updateColorMaps( myProperties.retrieveColorMaps() );
    mySegyAttr = myProperties.segyAttr;
    myDockPaneManager = new csDockPaneManager( true );
    myDockPaneManager.resetWindowLayout( myProperties.windowLayout, false );
    myDockPaneManager.setActivePaneFeedback( myProperties.showActivePaneFeedback );
    myDockPaneManager.setMaxTabWidth( myProperties.maxTabWidth );
    myDockPaneButtonSelection = new csDockPaneButtonSelection();
    myDockPaneButtonSelection.scrollbars = false;
    mySeisPaneManager = new csSeisPaneManager( this, myProperties.annAttr );

    myStatusBar = new JLabel("----", JLabel.LEFT );
    mySnapShotFrame = null;
    myMenuBar = new csMenuBar( this, myProperties );
    java.util.List<String> list = myPreferences.getRecentFileList();
    if( list != null ) {
      for( int i = list.size()-1; i >= 0; i-- ) {  // Add files backwards. Last added file goes to top of list
        myMenuBar.addRecentFile( (String)list.get(i) );
      }
    }
    setJMenuBar( myMenuBar );
   
    myShowAbsoluteTime  = false;
    myShowFKVelocity      = false;
    myFKTraceSpacing = 1.0f;    
    getContentPane().add( myToolBarLeft, BorderLayout.WEST );
    getContentPane().add( myToolBarTop, BorderLayout.NORTH );
    getContentPane().add( myDockPaneManager, BorderLayout.CENTER );
    getContentPane().add( myStatusBar , BorderLayout.SOUTH );
    
    ImageIcon appIcon = cseis.resources.csResources.getIcon("seaview_icon_smooth.png");
    if( appIcon != null ) setIconImage( appIcon.getImage() );

    pack();
    Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
    int appWidth  = 700;
    int appHeight = 800;
    double screenWidth = screenSize.getWidth();
    double screenHeight = screenSize.getHeight();
    setSize( new Dimension(appWidth,appHeight) );
    setLocation( (int)((screenWidth-appWidth)/2), (int)((screenHeight-appHeight)/2) );

    myToolBarTop.addToolBarListener( new csISeaViewToolBarListener() {
      @Override
      public void updateDisplayScalar( float scalar ) {
        mySeisPaneManager.updateScalar( scalar, (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
      }
      @Override
      public void incDisplayScalar() {
        mySeisPaneManager.incScaling( (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
      }
      @Override
      public void decDisplayScalar() {
        mySeisPaneManager.decScaling( (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
      }
    });
    myDockPaneManager.addDockPaneSelectionListener( new csIDockPaneSelectionListener() {
      @Override
      public void dockPaneSelected( csDockPane pane ) {
        if( pane != null ) {
          csSeisPaneBundle bundle = (csSeisPaneBundle)pane.getMainPanel();
          if( bundle == null ) return;
          myMenuBar.updateSelectedBundle( bundle );
          setTitle( pane.attr().titleDescription );
          csSeisDispSettings settings = bundle.seisView.getDispSettings();
          if( settings.scaleType == csSeisDispSettings.SCALE_TYPE_SCALAR ) {
            myToolBarTop.updateScalar( settings.dispScalar, true );
          }
          else if( settings.scaleType == csSeisDispSettings.SCALE_TYPE_TRACE ) {
            myToolBarTop.updateScalar( settings.fullTraceScalar, true );
          }
          else {
            myToolBarTop.updateScalar( -1, false );
          }
          // If no current read process is taking place:
          if( !myIsReadProcessOngoing ) {
            if( myTraceSelectionDialog != null ) {
              myTraceSelectionDialog.updateBundle( bundle );
            }
            updateMoveButtons( bundle );
          }
        }
        else {
          myMenuBar.updateSelectedBundle( null );
          setTitle(" ");
        }
      }
      @Override
      public void dockPaneClosed( csDockPane pane ) {
        csSeisPaneBundle bundle = (csSeisPaneBundle)pane.getMainPanel();
        if( bundle != null ) {
          mySeisPaneManager.removeBundle( bundle );
        }
        if( (csSeisPaneBundle)myDockPaneManager.getActivePanel() == null ) {
          myMenuBar.updateSelectedBundle( null );
          setTitle(" ");
        }
      }
      @Override
      public void hideStateChanged( boolean hiddenPanesExist ) {
//        myToolBarTop.updateSelectPanes( hiddenPanesExist );
      }
    });
    super.addWindowListener(new WindowAdapter() {
      @Override
      public void windowClosing( WindowEvent e ) {
        exitApplication();
      }
    });
  }
  //********************************************************************************
  public void addProcessingStep( String name ) {
    csSeisPaneBundle bundleActive = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundleActive != null ) {
      bundleActive.setProcessingStep( name );
    }
  }
  public void resetProcessing() {
    csSeisPaneBundle bundleActive = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundleActive != null ) {
      bundleActive.clearAllProcessingSteps();
    }
  }
  public void setScalarLocked( boolean isLocked ) {
    myIsScalarLocked = isLocked;
  }
  public void setOpenInNewPanel( boolean doOpen ) {
    myOpenInNewPanel = doOpen;
  }
  /**
   * Update scalar value in GUI which is used to scale displayed data set
   * @param bundle The csSeisPaneBundle which initiated this call
   */
  private void updateSeismicDisplayScalar( csSeisPaneBundle bundle ) {
    csSeisPaneBundle bundleActive = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundleActive != null && bundle.equals(bundleActive) ) {
      csSeisDispSettings settings = bundle.seisView.getDispSettings();
      if( settings.scaleType == csSeisDispSettings.SCALE_TYPE_SCALAR ) {
        myToolBarTop.updateScalar( settings.dispScalar, true );
      }
      else if( settings.scaleType == csSeisDispSettings.SCALE_TYPE_TRACE ) {
        myToolBarTop.updateScalar( settings.fullTraceScalar, true );
      }
      else {
        myToolBarTop.updateScalar( -1, false );
      }
    }
  }
  /**
   * Show dialog window: In this dialog, windows are be selected to be viewed or hidden.
   */
  public void showSelectWindowsDialog() {
    final JDialog dialog = new JDialog( this, "Window/data selection", true );
    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );
    JPanel panelButtons = new JPanel( new FlowLayout(FlowLayout.RIGHT) );
    JButton buttonClose = new JButton("Close");
    panelButtons.add( buttonClose );
    panelAll.add( myDockPaneManager.createSelectionPanel(), BorderLayout.NORTH );
    panelAll.add( panelButtons, BorderLayout.SOUTH );
    dialog.getContentPane().add( panelAll );
    dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    dialog.pack();
    dialog.setLocationRelativeTo( this );
    buttonClose.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        dialog.dispose();
      }
    });
    dialog.setVisible(true);
  }
  /**
   * Show dialog window: This dialog displays the inbuilt SEGY header maps
   */
  public void showSegyTrcHdrMaps() {
    csSegySetupDialog.showBuiltInSegyTrcHdrMaps( this, 0 );
  }
  /**
   * Show dialog window: In this dialog, the user defines settings related to SEGY data loading.
   */
  public void showSegySetupDialog() {
    csSegySetupDialog dialog = new csSegySetupDialog( this, mySegyAttr );
    dialog.addSegyAttrListener( new csISegyAttrListener() {
      @Override
      public void updateSegyAttr(csSegyAttr attr) {
        mySegyAttr = new csSegyAttr( attr );
      }
    });
    dialog.setVisible(true);
  }
  /**
   * Show dialog window: Applies to active bundle/pane.
   * In this dialog, the seismic display settings are specified.
   */
  public void showSettingsDialog() {
    mySeisPaneManager.showSettingsDialog( (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
  }
  /**
   * Show dialog window: Applies to active bundle/pane.
   * In this dialog, trace header annotations are defined.
   */
  public void showAnnotationDialog() {
    mySeisPaneManager.showAnnotationDialog( (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
  }
  /**
   * Show dialog window: Applies to active bundle/pane.
   * In this dialog, trace selection is defined.
   */
  public void showSelectionDialog() {
    if( myTraceSelectionDialog == null ) {
      csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
      myTraceSelectionDialog = new csTraceSelectionDialog( this, bundle );
      myTraceSelectionDialog.addTraceSelectionListener( new csITraceSelectionListener() {
        @Override
        public void updateTraceSelection( csTraceSelectionParam param ) {
          csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
          if( bundle != null ) {
            bundle.updateTraceSelection( param );
            if( !myIsReadProcessOngoing ) {
              readSeismicBundle( bundle, SeaView.MOVE_OPTION_SELECTION );
            }
          }
        }
        @Override
        public void startScanTraceSelection( csHeaderDef headerDef ) {
          csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
          if( bundle != null ) {
            if( !bundle.startTraceHeaderScan( headerDef ) ) {
              //
            }
          }
        }
        @Override
        public void cancelScanTraceSelection() {
          csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
          if( bundle != null ) {
            bundle.cancelTraceHeaderScan();
          }
        }
      });
    }
    myTraceSelectionDialog.setVisible(true);
  }
  /**
   * Show dialog window: Applies to active bundle/pane.
   * In this dialog, trace header overlays shown on top of the seismic are defined.
   */
  public void showOverlayDialog() {
    mySeisPaneManager.showOverlayDialog( (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
  }
  /**
   * Show dialog window: Applies to active bundle/pane.
   * This dialog displays trace header values of the trace at the mouse location.
   */
  public void showHeaderMonitor() {
    mySeisPaneManager.showHeaderMonitor( (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
  }
  /**
   * Show dialog window: Applies to active bundle/pane.
   * This dialog displays trace sample values of the trace at the mouse location.
   */
  public void showTraceMonitor() {
    mySeisPaneManager.showTraceMonitor( (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
  }
  /**
   * Show dialog window: Applies to active bundle/pane.
   * In this dialog, a trace header graph shown on top of the seismic is defined.
   */
  public void showGraphPanel( boolean doShow ) {
    mySeisPaneManager.showGraphPanel( doShow, (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
  }
  /**
   * Show dialog window: Applies to active bundle/pane.
   */
  public void showPreferencesDialog() {
    csPreferencesDialog dialog = new csPreferencesDialog( this, myProperties );
    dialog.setVisible( true );
  }
  public void setDefaultNumTraces( int numTraces ) {
    if( numTraces != myProperties.numTraces ) {
      myProperties.numTraces = numTraces;
      csTraceSelectionParam.NUM_TRACES = numTraces;
    }
  }
  public void setMaxTabWidth( int maxTabWidth ) {
    if( maxTabWidth != myProperties.maxTabWidth ) {
      myProperties.maxTabWidth = maxTabWidth;
      myDockPaneManager.setMaxTabWidth(maxTabWidth);
    }
  }
  public void setShowFilename( boolean showFilename ) {
    if( showFilename != myProperties.showFilename ) {
      myProperties.showFilename = showFilename;
      mySeisPaneManager.setShowFilenameInView( showFilename, (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
    }
  }
  public void manageWellPaths() {
    if( myWellPathDialog == null ) {
      myWellPathDialog = new csWellPathDialog( this, myProperties.wellPathDirectoryPath );
      myWellPathDialog.addListener( new csIWellPathDialogListener() {
        @Override
        public void updateWellPaths( ArrayList<csWellPathAttr> wellPathAttrList, csWellPathOverlayAttr wellPathOverlayAttr ) {
          mySeisPaneManager.updateWellPaths( wellPathAttrList, wellPathOverlayAttr );
        }
        @Override
        public void updateWellPathDirectory( String directory ) {
          myProperties.wellPathDirectoryPath = directory;
        }
      });
    }
    myWellPathDialog.setVisible(true);
  }
  public void setShowActivePaneFeedback( boolean showActivePaneFeedback ) {
    myProperties.showActivePaneFeedback = showActivePaneFeedback;
    myDockPaneManager.setActivePaneFeedback(showActivePaneFeedback);
  }
  /**
   * Manage custom color maps, i.e. create, delete, edit custom color maps.
   */
  public void manageCustomColorMaps() {
    csColorMapGenerator generator = new csColorMapGenerator(
            this, myProperties.retrieveColorMaps(), myProperties.propertiesDirectoryPath );
    generator.setTestButton( true, "Test currently edited color map on active data set" );
    generator.addListener( myCustomColorMapModel );
    generator.addListener( new csIColorMapGeneratorListener() {
      @Override
      public void applyColorMap(csCustomColorMap cmap) {
        // Test new color map: Apply to active bundle
        csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
        if( bundle != null ) bundle.testColorMap( cmap );
      }
      @Override
      public void updateColorMaps( java.util.List<csCustomColorMap> list ) {
        myProperties.updateColorMaps( list );
      }
    });
    generator.setVisible(true);
  }
  /**
   * Show SEGY character header of SEGY data set in active seismic bundle.
   */
  public void showSegyCharHdr() {
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundle != null && bundle.getFileFormat() == SeaView.FORMAT_SEGY ) bundle.showSegyCharHdr();
  }
  /**
   * Show SEGY binary header of SEGY data set in active seismic bundle.
   */
  public void showSegyBinHdr() {
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundle != null && bundle.getFileFormat() == SeaView.FORMAT_SEGY ) bundle.showSegyBinHdr();
  }
  /**
   * Combine two seismic bundles. A new seismic data set is created by adding, subtracting etc.
   * two of currently displayed seismic data sets.
   */
  public void combineTwoSeismicBundles() {
    int numBundles = mySeisPaneManager.getNumBundles();
    if( numBundles < 1 ) return;
    csSeisPaneBundle[] bundleList = new csSeisPaneBundle[numBundles];
    for( int i = 0; i < numBundles; i++ ) {
      bundleList[i] = mySeisPaneManager.getBundle(i);
    }
    csCombineDataDialog dialog = new csCombineDataDialog( this, bundleList, (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
    dialog.setVisible(true);
  }
  /**
   * Retrieve general settings for trace header annotation.
   * @return Attributes defining trace header annotation
   */
  public csAnnotationAttributes getAnnotationAttributes() {
    return myProperties.annAttr;
  }
  /**
   * Trace highlighting. Applies to all visible seismic panes.
   * @param doHighlightTrace true if highlighting shall be turned on, false otherwise
   */
  public void highlightTrace( boolean doHighlightTrace ) {
    mySeisPaneManager.highlightTrace(doHighlightTrace);
  }
  /**
   * Show color bar yes/no. Applies to all visible seismic panes.
   * @param doShow true if color bar shall be displayed, false otherwise
   */
  public void showColorBar( boolean doShow ) {
    mySeisPaneManager.showColorBar(doShow);
  }
  /**
   * Zoom seismic view. Applies to active seismic bundle/pane.
   * @param zoomType The zoom 'type': csSeisPane.ZOOM_VERT or csSeisPane.ZOOM_HORZ
   * @param zoomMode The zoom 'mode': csSeisPane.ZOOM_IN or csSeisPane.ZOOM_OUT
   */
  public void zoom( int zoomType, int zoomMode ) {
    mySeisPaneManager.zoom( zoomType, zoomMode, (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
  }
  /**
   * Change number of bits used for view colors. Applies to all data views.
   * @param colorMapType 32bit or 8bit: csColorMap.COLOR_MAP_TYPE_32BIT or csColorMap.COLOR_MAP_TYPE_8BIT
   */
  public void changeColorBits( int colorBitType ) {
    myColorBitType = colorBitType;
    final Component glassPane = getGlassPane();
    glassPane.setCursor( Cursor.getPredefinedCursor( Cursor.WAIT_CURSOR ) );
    KeyAdapter keyAdapter = new KeyAdapter() {};
    MouseAdapter mouseAdapter = new MouseAdapter() {};
    glassPane.addKeyListener( keyAdapter );
    glassPane.addMouseListener( mouseAdapter );
    glassPane.setVisible( true );
    myCustomColorMapModel.updateColorMapType(colorBitType);
    mySeisPaneManager.setColorBits( colorBitType );
    glassPane.removeKeyListener( keyAdapter );
    glassPane.removeMouseListener( mouseAdapter );
    glassPane.setVisible( false );
  }
  public void updateSyncDispSettings( boolean doSync ) {
    mySeisPaneManager.updateSyncDispSettings( doSync );
  }
  /**
   * Load seismic display settings from external file. Apply to active bundle/pane.
   */
  public void loadDispSettings() {
    JFileChooser fc = new JFileChooser(myProperties.propertiesDirectoryPath);
    int returnVal = fc.showOpenDialog(SeaView.this);
    if (returnVal == JFileChooser.APPROVE_OPTION) {
      String filename = fc.getSelectedFile().getAbsolutePath();
      csSeisDispSettings settings = new csSeisDispSettings();           // initialized w/ defaults
      csSeisDispProperties properties = new csSeisDispProperties( filename, myProperties.retrieveColorMaps() );
      try {
        properties.load(settings);
        mySeisPaneManager.updateDisplaySettings( settings, (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
        myStatusBar.setText( "   ...loaded seismic display settings from file " + properties.filename());
      }
      catch(IOException ex) {
        JOptionPane.showMessageDialog( this, "Error occurred when reading display settings file:\n" +
          ex.getMessage(),
          "Error message", JOptionPane.ERROR_MESSAGE );
      }
    }
    myProperties.propertiesDirectoryPath = fc.getCurrentDirectory().getAbsolutePath();
  }
  //-----------------------------------------------------------------------------------------------------
  //
  /**
   * Save display settings of active bundle/pane to external file.
   */
  public void saveDispSettings() {
    JFileChooser fc = new JFileChooser(myProperties.propertiesDirectoryPath);
    int returnVal = fc.showSaveDialog(SeaView.this);
    if (returnVal == JFileChooser.APPROVE_OPTION) {
      String filename = fc.getSelectedFile().getAbsolutePath();
      csSeisDispProperties properties = new csSeisDispProperties( filename, myProperties.retrieveColorMaps() );
      try {
        csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
        if( bundle == null ) return;
        properties.write(bundle.seisView.getDispSettings());
        myProperties.propertiesDirectoryPath = fc.getCurrentDirectory().getAbsolutePath();
        myStatusBar.setText( "   ...saved seismic display settings to file " + properties.filename());
      }
      catch (IOException ex) {
        JOptionPane.showMessageDialog( this, "Error occurred when writing display settings file:\n" +
          ex.getMessage(),
          "Error message", JOptionPane.ERROR_MESSAGE );
      }
    }
  }
  /**
   * Set default display settings and apply to active bundle/pane.
   */
  public csSeisDispSettings getDefaultDispSettings( boolean isFrequencyDomain ) {
    csSeisDispSettings settings = new csSeisDispSettings();
    if( isFrequencyDomain ) {
      settings.isVIDisplay = true;
      settings.showWiggle  = false;
      settings.isPosFill   = false;
      settings.isNegFill   = false;
      settings.viType      = csSeisDispSettings.VA_TYPE_2DSPLINE;
      settings.viColorMap.resetColors( new csColorMap(csColorMap.RAINBOW_BLACK,settings.viColorMap.getColorMapType()) );
      settings.scaleType   = csSeisDispSettings.SCALE_TYPE_RANGE;
    }
    else { // Grey scale
      settings.isVIDisplay = true;
      settings.showWiggle  = false;
      settings.isPosFill   = false;
      settings.isNegFill   = false;
      settings.viType      = csSeisDispSettings.VA_TYPE_2DSPLINE;
      settings.viColorMap.resetColors( new csColorMap(csColorMap.GRAY_WB,settings.viColorMap.getColorMapType()) );
    }
    return settings;
  }
  public void setDefaultDispSettings() {
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundle == null ) return;
    csSeisDispSettings settings = getDefaultDispSettings( bundle.isFrequencyDomain() );
    settings.dispScalar = bundle.seisView.getDispScalar();
    if( bundle.isFrequencyDomain() ) {
      settings.minValue    = bundle.seisView.getTraceBuffer().minTotalAmplitude();
      settings.maxValue    = bundle.seisView.getTraceBuffer().maxTotalAmplitude();
    }
    bundle.seisView.updateDispSettings( settings );
  }
  /**
   * Load seismic display settings from default location and apply to active bundle/pane.
   * @param reportErrors Pass 'true' if errors shall be reported.
   * @return true if operation was successful.
   */
  public boolean applyDefaultDispSettings( boolean reportErrors ) {
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundle == null ) return false;
    float dispScalar = bundle.seisView.getDispScalar();
    csSeisDispSettings dispSettings = loadDefaultDispSettings( bundle.isFrequencyDomain(), reportErrors );
    dispSettings.dispScalar = dispScalar;  // Use original scalar, do not update
    dispSettings.minValue = bundle.seisView.getTraceBuffer().minTotalAmplitude();  // TEMP
    dispSettings.maxValue = bundle.seisView.getTraceBuffer().maxTotalAmplitude();  // TEMP
    bundle.seisView.updateDispSettings( dispSettings );
    if( dispSettings.isLogScale ) {
      myMenuBar.setIsLogScale( dispSettings.isLogScale );
    }
    return true;
  }
  /**
   * Load seismic display settings from default location.
   * @param isFrequencyDomain 'true' if frequency domain settings shall be loaded.
   * @param reportErrors Pass 'true' if errors shall be reported.
   * @return true if operation was successful.
   */
  public csSeisDispSettings loadDefaultDispSettings( boolean isFrequencyDomain, boolean reportErrors ) {
    csSeisDispProperties properties = null;
    java.util.List<csCustomColorMap> colorMapList = myProperties.retrieveColorMaps();
    if( isFrequencyDomain ) {
      properties = new csSeisDispProperties(csPreferences.getCseisDirectory().getAbsolutePath() + java.io.File.separatorChar+ PROPERTIES_FREQ_FILE_NAME,
              colorMapList );
    }
    else {
      properties = new csSeisDispProperties(csPreferences.getCseisDirectory().getAbsolutePath() + java.io.File.separatorChar + PROPERTIES_FILE_NAME,
              colorMapList );
    }
    csSeisDispSettings settings = new csSeisDispSettings();
    try {
      properties.load(settings);
      myStatusBar.setText( "   ...loaded seismic display settings from file " + properties.filename());
    }
    catch(IOException ex) {
      if( reportErrors ) {
        JOptionPane.showMessageDialog( this, "Error occurred when reading default display settings file:\n" +
          ex.getMessage(),
          "Error message", JOptionPane.ERROR_MESSAGE );
      }
      settings = getDefaultDispSettings( isFrequencyDomain );
    }
    return settings;
  }
  /**
   * Save seismic display settings to default location.
   * @param reportErrors Pass 'true' if error shall be reported.
   * @return true if operation was successful
   */
  public boolean saveDefaultDispSettings( boolean reportErrors ) {
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundle == null ) return false;
    csSeisDispProperties properties = null;
    if( bundle.isFrequencyDomain() ) {
      properties = new csSeisDispProperties(csPreferences.getCseisDirectory().getAbsolutePath() + java.io.File.separatorChar + PROPERTIES_FREQ_FILE_NAME,
              myProperties.retrieveColorMaps() );
    }
    else {
      properties = new csSeisDispProperties(csPreferences.getCseisDirectory().getAbsolutePath() + java.io.File.separatorChar + PROPERTIES_FILE_NAME,
              myProperties.retrieveColorMaps() );
    }
    
    try {
      properties.write(bundle.seisView.getDispSettings());
      myStatusBar.setText( "   ...saved seismic display settings to file " + properties.filename());
    }
    catch (IOException ex) {
      if( reportErrors ) {
        JOptionPane.showMessageDialog( this, "Error occurred when writing default display settings file:\n" +
          ex.getMessage(),
          "Error message", JOptionPane.ERROR_MESSAGE );
      }
      return false;
    }
    return true;
  }
  /**
   * Reset 'windows layout': Organizes visible windows (=seismic panes) in different ways.
   * @param layoutOption Defines how windows will be organized on screen, e.g. csDockPaneManager.LAYOUT_ONE_ROW etc.
   * @param refresh Pass 'true' if display shall be refreshed when done.
   */
  public void resetWindowsLayout( int layoutOption, boolean refresh ) {
    myDockPaneManager.resetWindowLayout(layoutOption, refresh);
  }
  /**
   * Display vertical/Y axis in log scale. Applies to active bundle/pane.
   * @param setLogScaleY Pass 'true' if Y axsi shall be displayed in log scale.
   */
  public void setLogScaleY( boolean setLogScaleY ) {
    mySeisPaneManager.setLogScaleY( setLogScaleY, (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
  }
  /**
   * Show crosshair on seismic views. Applies to all visible bundles/panes.
   * @param setOn Pass 'true' if crosshair shall be shown on top of seismic views.
   * @param setSnap Pass 'true' if crosshair shall be snapped to nearest trace & sample.
   */
  public void setCrosshair( boolean setOn, boolean setSnap ) {
    mySeisPaneManager.setCrosshair( setOn, setSnap );
  }
  /**
   * Fit view to 'screen': Zooms in seismic view to fit into its visible area. Applies to active bundle/pane
   */
  public void fitToScreen() {
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundle == null ) return;
    bundle.fitToScreen();
  }
  /**
   * Create snapshot.
   * @param includeSideLabels Pass 'true' if side labels shall be included in snapshot.
   */
  public void createSnapShot( boolean includeSideLabels ) {
    if( mySeisPaneManager.getNumBundles() == 0 ) return;
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundle == null ) {
      JOptionPane.showMessageDialog( this,
          "No seismic pane selected/active. ...this is likely a SeaView bug.\n" +
          "Click on top bar of seismic pane to select and repeat snapshot",
          "Select seismic pane", JOptionPane.INFORMATION_MESSAGE );
      return;
    }
    BufferedImage image = mySeisPaneManager.createSnapShot( includeSideLabels, bundle );
    if( mySnapShotFrame == null ) {
      mySnapShotFrame = new csSnapShotFrame(this,cseis.resources.csResources.getIcon("seaview_icon_smooth.png"));
    }
    mySnapShotFrame.addImage(image,bundle.getTitle());
    openSnapshotViewer();
  }
  /**
   * Open snapshow viewer application window.
   */
  public void openSnapshotViewer() {
    if( mySnapShotFrame == null ) {
      mySnapShotFrame = new csSnapShotFrame(this,cseis.resources.csResources.getIcon("seaview_icon_smooth.png"));
    }
    mySnapShotFrame.setVisible(true);
  }
  /**
   * Create screen dump (=image file) of current bundle/pane.
   * @param saveToDefault Pass 'true' if screen dump shall be saved to default location.
   */
  public void createScreenDump( boolean saveToDefault ) {
    if( myFileChooserDump == null ) {
      myFileChooserDump = new csFileChooser( "Export seismic display as jpeg image" );
      myFileChooserDump.setCurrentDirectory( new File(myProperties.screenDumpDirectoryPath) );
    }
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    boolean writeCurrentFile;
    File outFile = null;
    BufferedImage image = null;
    String filename = "";
    do {
      writeCurrentFile = false;
      if( saveToDefault ) {
        filename = bundle.getTitle();
        for( int i = filename.length()-1; i >= 0; i-- ) {
          if( filename.charAt(i) == '.' ) {
            filename = filename.substring( 0, i ) + ".jpg";
            break;
          }
          if( i == 0 ) {
            filename = filename + ".jpg";
          }
        }
        for( int i = filename.length()-1; i >= 0; i-- ) {
          if( filename.charAt(i) == '/' ) {
            filename = filename.substring( i+1 );
            break;
          }
        }
        filename = myFileChooserDump.getCurrentDirectory() + "/" + filename;
      }
      else {
        if( !myFileChooserDump.selectFileForWriting(this) ) return;
        filename = myFileChooserDump.getSelectedFile().getAbsolutePath();
        int length = filename.length();
        if( length < 5 || (filename.substring(length-4, length).compareTo(".jpg") != 0 &&
          filename.substring(length-5, length).compareTo(".jpeg") != 0 ) ) {
            filename += ".jpg";
        }
        else {
          writeCurrentFile = true;
        }
      }

      Rectangle rect = bundle.seisPane.getVisibleRect();
      image = new BufferedImage(rect.width, rect.height, BufferedImage.TYPE_INT_RGB);
      bundle.seisPane.paintAll( image.createGraphics() );
      outFile = new File(filename);
      if( !writeCurrentFile && outFile.exists() ) {
        int option = JOptionPane.showConfirmDialog(this,
          "File " + filename + " already exists.\nOK to overwrite?",
          "Confirm overwrite",
          JOptionPane.YES_NO_OPTION );
        if( option == JOptionPane.NO_OPTION ) {
          writeCurrentFile = false;
          return;
        }
        else if( option == JOptionPane.YES_OPTION ) {
          writeCurrentFile = true;
        }
        else {
          return;
        }
      }
      else {
        writeCurrentFile = true;
      }
    } while( !writeCurrentFile );

    try {
      FileImageOutputStream output = new FileImageOutputStream(outFile);

      Iterator<ImageWriter> iter = ImageIO.getImageWritersByFormatName("JPG");
      if (iter.hasNext()) {
         ImageWriter writer = (ImageWriter)iter.next();
         ImageWriteParam iwp = writer.getDefaultWriteParam();
         iwp.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
         iwp.setCompressionQuality(0.95f);
         
         writer.setOutput(output);
         IIOImage imageOut = new IIOImage(image, null, null);
         writer.write(null, imageOut, iwp);
      }
      myStatusBar.setText( "   ...export to jpeg file " + filename);
      output.close();
  	}
    catch( IOException e ) {
      JOptionPane.showMessageDialog( this,
              "Error occurred when creating screendump:\n" + e.getMessage(),
              "Error message",
              JOptionPane.ERROR_MESSAGE );
    }
  }
  /**
   * Refresh currently displayed file. Applies to active bundle/pane.
   *
   * @return false if problem occurred.
   */
  public boolean refreshFile() {
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundle == null ) return false;
    myIsRefreshFileOperation = true;
    boolean retValue = openFile( bundle.getFilenamePath(), bundle.getFileFormat() );
    myIsRefreshFileOperation = false;
    return retValue;
  }
  /**
   * Open several files.
   * @param filenameList List of input file names.
   */
  public void openSeveralFiles( java.util.List<String> filenameList ) {
    myInputFilenameQueue = new ArrayList<String>( filenameList );
    openNextFileInQueue();
  }
  private void openNextFileInQueue() {
    if( myInputFilenameQueue == null ) return;
    if( myInputFilenameQueue.isEmpty() ) {
      myInputFilenameQueue = null;
      return;
    }
//    try {
//      Thread.sleep(500);
//    } catch (InterruptedException ex) {
//    }
    SwingUtilities.invokeLater( new Runnable() {
      @Override
      public void run() {
        String filename = myInputFilenameQueue.remove(0);
        openFile( filename, SeaView.FORMAT_ALL );
      }
    });
  }
  /**
   * Start velocity analysis setup
   * @param reset True if velocity analysis shall be reset
   */
  public void openVelocityAnalysis( boolean reset ) {
    if( myVelocityAnalysisManager == null || reset ) {
      int numBundles = mySeisPaneManager.getNumBundles();
      if( numBundles < 1 ) return;
      csSeisPaneBundle[] bundleList = new csSeisPaneBundle[numBundles];
      for( int i = 0; i < numBundles; i++ ) {
       bundleList[i] = mySeisPaneManager.getBundle(i);
      }
      myVelocityAnalysisManager = new csVelocityAnalysisManager( this, bundleList, (csSeisPaneBundle)myDockPaneManager.getActivePanel(),
        myProperties.velFieldDirectoryPath);
    }
    myVelocityAnalysisManager.open();
  }
  /**
   * Open new input file.
   * @param filename The name of the input file.
   * @return true if operation was successful.
   */
  public boolean openFile( String filename ) {
    return openFile( filename, SeaView.FORMAT_ALL );
  }
  /**
   * Open new input file.
   * @return true if operationw as successful.
   */
  public boolean openFile() {
    return openFile(SeaView.FORMAT_ALL );
  }
  /**
   * Open new input file.
   * @param fileFormat The format of the file, e.g. SEGY.
   * @return true if operation was successful
   */
  public boolean openFile( int fileFormat ) {
    boolean isNew = false;
    if( mySeismicFileChooser == null ) {
      mySeismicFileChooser = new JFileChooser();
      mySeismicFileChooser.setCurrentDirectory( new File(myProperties.dataDirectoryPath) );
      isNew = true;
    }
    else {
      csSeisPaneBundle activeBundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
      if( activeBundle != null ) {
        mySeismicFileChooser.setCurrentDirectory( new File( activeBundle.getPath() ) );
      }
    }
    if( isNew || fileFormat != myCurrentSelectedFileFormat ) {
      javax.swing.filechooser.FileFilter[] filters = mySeismicFileChooser.getChoosableFileFilters();
      for( int i = 0; i < filters.length; i++ ) {
        mySeismicFileChooser.removeChoosableFileFilter(filters[i]);
      }
      mySeismicFileChooser.setAcceptAllFileFilterUsed(true);
      FileNameExtensionFilter filterSEGY = new FileNameExtensionFilter(
        "SEGY files (*.segy, *.sgy, *.SEGY, *.SGY)", "segy", "sgy", "SEGY", "SGY");
      FileNameExtensionFilter filterCSEIS = new FileNameExtensionFilter(
         "CSEIS files (*.cseis, *.oseis)", "cseis", "oseis");
      FileNameExtensionFilter filterSU = new FileNameExtensionFilter(
         "SU files (*.su, *.SU)", "su", "SU");
      FileNameExtensionFilter filterASCII = new FileNameExtensionFilter(
         "ASCII files (*.asc, *.txt, *.out)", "asc", "txt", "out");
      FileNameExtensionFilter filterSEGD = new FileNameExtensionFilter(
        "SEGD files (*.segd, *.sgd, *.gunlink, *.SEGD, *.SGD, *.raw)", "segd", "sgd", "gunlink", "SEGD", "SGD", "raw");
      FileNameExtensionFilter filterRSF = new FileNameExtensionFilter(
        "RSF files (*.rsf, *.RSF)", "rsf", "RSF");
      FileNameExtensionFilter filterJPG = new FileNameExtensionFilter(
        "JPEG files (*.jpg, *.jpeg, *.JPG, *.JPEG)", "jpg", "jpeg", "JPG", "JPEG");
      if( fileFormat == SeaView.FORMAT_ALL ) {
        mySeismicFileChooser.addChoosableFileFilter(filterCSEIS);
        mySeismicFileChooser.addChoosableFileFilter(filterSEGY);
        mySeismicFileChooser.addChoosableFileFilter(filterSU);
        mySeismicFileChooser.addChoosableFileFilter(filterASCII);
        mySeismicFileChooser.addChoosableFileFilter(filterSEGD);
        mySeismicFileChooser.addChoosableFileFilter(filterRSF);
        mySeismicFileChooser.addChoosableFileFilter(filterJPG);
      }
      else if( fileFormat == SeaView.FORMAT_CSEIS ) {
        mySeismicFileChooser.addChoosableFileFilter(filterCSEIS);
        mySeismicFileChooser.setFileFilter( filterCSEIS );
      }
      else if( fileFormat == SeaView.FORMAT_SEGY ) {
        mySeismicFileChooser.addChoosableFileFilter(filterSEGY);
        mySeismicFileChooser.setFileFilter( filterSEGY );
      }
      else if( fileFormat == SeaView.FORMAT_SU ) {
        mySeismicFileChooser.addChoosableFileFilter(filterSU);
      }
      else if( fileFormat == SeaView.FORMAT_ASCII ) {
        mySeismicFileChooser.addChoosableFileFilter(filterASCII);
        mySeismicFileChooser.setFileFilter( filterASCII );
      }
      else if( fileFormat == SeaView.FORMAT_SEGD ) {
        mySeismicFileChooser.addChoosableFileFilter(filterSEGD);
        mySeismicFileChooser.setFileFilter( filterSEGD );
      }
      else if( fileFormat == SeaView.FORMAT_RSF ) {
        mySeismicFileChooser.addChoosableFileFilter(filterRSF);
        mySeismicFileChooser.setFileFilter( filterRSF );
      }
      else if( fileFormat == SeaView.FORMAT_JPG ) {
        mySeismicFileChooser.addChoosableFileFilter(filterJPG);
        mySeismicFileChooser.setFileFilter( filterJPG );
      }
      myCurrentSelectedFileFormat = fileFormat;
    } // END isNew || fileformat != myCurrentSelectedFileFormat
    JFileChooser fc = mySeismicFileChooser;
    int option = fc.showOpenDialog( SeaView.this );
    if( option == JFileChooser.APPROVE_OPTION ) {
      String filename = fc.getSelectedFile().getAbsolutePath();
      return openFile( filename, fileFormat );
    }
    else {
      return false;
    }
  }
  
  /**
   * Open specific file.
   * @param filename The name of the file.
   * @param fileFormat The format of the file.
   * @return true if operation was successful.
   */
  public boolean openFile( String filename, int fileFormat ) {
    csISeismicReader newReader = null;
    try {
      File file = new File( filename );
      if( !file.exists() ) {
        throw( new Exception("File not found: " + filename) );
      }
      // If file format was not specifically selected, reset file format based on file extension:
      if( fileFormat == SeaView.FORMAT_ALL ) {
        if( filename.toLowerCase().endsWith("cseis") || filename.toLowerCase().endsWith("oseis") ) {
          fileFormat = SeaView.FORMAT_CSEIS;
        }
        else if( filename.toLowerCase().endsWith("segy") || filename.toLowerCase().endsWith("sgy") ) {
          fileFormat = SeaView.FORMAT_SEGY;
        }
        else if( filename.toLowerCase().endsWith("su") ) {
          fileFormat = SeaView.FORMAT_SU;
        }
        else if( filename.toLowerCase().endsWith("txt") || filename.toLowerCase().endsWith("asc") ) {
          fileFormat = SeaView.FORMAT_ASCII;
        }
        else if( filename.toLowerCase().endsWith("segd") || filename.toLowerCase().endsWith("sgd") || filename.toLowerCase().endsWith("gunlink") || filename.toLowerCase().endsWith("raw") ) {
          fileFormat = SeaView.FORMAT_SEGD;
        }
        else if( filename.toLowerCase().endsWith("rsf") ) {
          fileFormat = SeaView.FORMAT_RSF;
        }
        else if( filename.toLowerCase().endsWith("jpg") ) {
          fileFormat = SeaView.FORMAT_JPG;
        }
        else {   // Per default, assume it's a SEGY file:
          fileFormat = SeaView.FORMAT_SEGY;
        }
      }
      java.nio.ByteOrder endianOrder = java.nio.ByteOrder.nativeOrder();

      // Create new reader object based on file format
      switch( fileFormat ) {
        case SeaView.FORMAT_CSEIS:
          newReader = new csNativeSeismicReader( filename );
          break;
        case SeaView.FORMAT_SEGY:
          newReader = new csNativeSegyReader( filename, mySegyAttr.hdrMap, 200,
          mySegyAttr.endianSEGYData != java.nio.ByteOrder.BIG_ENDIAN, mySegyAttr.endianSEGYHdr != java.nio.ByteOrder.BIG_ENDIAN, mySegyAttr.autoScaleCoord );
          break;
        case SeaView.FORMAT_SU:
          newReader = new csNativeSegyReader( filename, csNativeSegyReader.HDR_MAP_SU_BOTH, 200,
            mySegyAttr.endianSUData != java.nio.ByteOrder.BIG_ENDIAN, mySegyAttr.endianSUHdr != java.nio.ByteOrder.BIG_ENDIAN, mySegyAttr.autoScaleCoord );
          break;
        case SeaView.FORMAT_SEGD:
          newReader = new csNativeSegdReader( filename );
          break;
        case SeaView.FORMAT_RSF:
          newReader = new csNativeRSFReader( filename );
          break;
        case SeaView.FORMAT_ASCII:
          newReader = new csASCIIReader( filename );
          break;
        case SeaView.FORMAT_JPG:
          newReader = new csImageSeismicReader( filename );
          break;
        default:  // Per default, assume it's a SEGY file:
          fileFormat = SeaView.FORMAT_SEGY;
          newReader = new csNativeSegyReader( filename, mySegyAttr.hdrMap, 200,
            mySegyAttr.endianSEGYData != endianOrder, mySegyAttr.endianSEGYHdr != endianOrder, mySegyAttr.autoScaleCoord );
      }
      if( newReader == null || newReader.numTraces() == 0 ) {
        throw( new Exception("Input file " + filename + "\ncontains no data") );
      }
      myMenuBar.addRecentFile( file.getAbsolutePath() );
      if( mySeismicFileChooser == null ) {
        mySeismicFileChooser = new JFileChooser();
      }
      mySeismicFileChooser.setCurrentDirectory(file);
      file = null;
    }
    catch( IOException exc ) {
      JOptionPane.showMessageDialog( this,
          "Unable to open file.\n" +
          "System message:\n" +
          "'" + exc.getMessage() + "'",
          "Error",
          JOptionPane.ERROR_MESSAGE );
      return false;
    }
    catch( Exception exc ) {
      JOptionPane.showMessageDialog( this,
          "Unable to open file.\n" +
          "System message:\n" +
          "'" + exc.getMessage() + "'",
          "Error",
          JOptionPane.ERROR_MESSAGE );
      return false;
    }

    myToolBarLeft.setToolbarEnabled( true );
    myToolBarTop.setToolbarEnabled( true );
    return readData( newReader, filename, fileFormat, myOpenInNewPanel );
  }
  /**
   * Read data set into SeaView.
   * @param reader The seismic reader object.
   * @param filenamePath The full path name to the file.
   * @param fileFormat   Flag defining file format
   * @param openInNewPane Open new data set in new pane (replace data in current if 'false')
   * @return true if operation was successful.
   */
  public synchronized boolean readData( csISeismicReader reader, String filenamePath, int fileFormat, boolean openInNewPane ) {
    csSeisPaneBundle bundle = null;
    boolean renameActivePane = true;
    for( int ibundle = 0; ibundle < mySeisPaneManager.getNumBundles(); ibundle++ ) {
      if( filenamePath.compareTo(mySeisPaneManager.getBundle(ibundle).getFilenamePath()) == 0 ) {
        int option = JOptionPane.YES_OPTION;
        if( !myIsRefreshFileOperation ) {
          option = JOptionPane.showConfirmDialog(this,
            "File\n" + filenamePath + ":\nRe-open/refresh?",
            "Confirm refresh",
            JOptionPane.YES_NO_OPTION );
        }
        if( option == JOptionPane.YES_OPTION ) {
          bundle = mySeisPaneManager.getBundle(ibundle);
          csSeisPaneBundle bundleActive = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
          if( !bundleActive.equals(bundle) ) {
            // Bundle which is to be updated is currently not active: Do not rename active bundle to avoid a name mixup
            renameActivePane = false;
          }
          if( myDockPaneManager.getState( bundle ) == csDockPaneManager.STATE_HIDDEN ) {
//            System.out.println("State of bundle is HIDDEN");
//            myDockPaneManager.changeDocking( bundle, true );
          }
          else {
//            System.out.println("State of bundle is DOCKED");
          }
        }
        else { //if( option == JOptionPane.NO_OPTION ) {
          return false;
        }
      }
    }
    if( bundle == null && !openInNewPane ) {
      bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    }
    readSeismicBundle( bundle, reader, filenamePath, fileFormat, SeaView.MOVE_OPTION_SELECTION, true, renameActivePane );
     
    return true;
  }
  
  /**
   * Read seismic bundle.
   * 
   * @param bundle The seismic bundle for which data traces shall be read.
   * @param moveOption The 'move' option defines which traces to read in.
   */
  private void readSeismicBundle( csSeisPaneBundle bundle, int moveOption ) {
    readSeismicBundle( bundle, null, null, bundle.getFileFormat(), moveOption, false, true );
  }
  /**
   * Read seismic bundle.
   * 
   * @param bundle The seismic bundle for which data traces shall be read.
   * @param reader The seismic reader object.
   * @param filenamePath The full pathname to the input file.
   * @param moveOption The 'move' option defines which traces to read in.
   * @param reinitializeReader Pass 'true' if the reader object shall be re-initialized before read operation.
   *                           This needs to be done if the file must be re-opened before reading (=refreshed).
   * @param renameActivePane   Pass 'true' if the active (dock) pane shall be renamed if necessary.
   */
   private void readSeismicBundle( csSeisPaneBundle bundle, csISeismicReader reader,
          String filenamePath, int fileFormat, int moveOption, boolean reinitializeReader, boolean renameActivePane ) {
    if( myIsReadProcessOngoing ) return;
    myIsReadProcessOngoing = true;
    boolean resetScalar = !myIsScalarLocked;
    if( bundle == null ) {
      csSeisDispSettings dispSettings = null;
      csSeisPaneBundle bundleActive = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
      if( bundleActive != null && reader.isFrequencyDomain() == bundleActive.isFrequencyDomain() ) {
        dispSettings = bundleActive.seisView.getDispSettings();
      }
      else { // Load default settings
        resetScalar = true;
        dispSettings = loadDefaultDispSettings( reader.isFrequencyDomain(), false );
        if( dispSettings.isLogScale ) {
          myMenuBar.setIsLogScale( dispSettings.isLogScale );
        }
      }
      csFilename filename = new csFilename( filenamePath );
      bundle = createSeismicBundle( dispSettings, filename, fileFormat );
      myCurrentNewDockPane = myDockPaneManager.addPanel( bundle, filename.filename, myDockPaneButtonSelection, filename.filenamePath );
      myCurrentNewDockPane.addSyncListener( new csIDockPaneSyncListener() {
        @Override
        public void syncStateChanged( csDockPane pane, boolean isSync ) {
          csSeisPaneBundle bundle = (csSeisPaneBundle)pane.getMainPanel();
          if( bundle != null ) {
            bundle.setSyncState( isSync );
          }
        }
      });
      reinitializeReader = true;
    }
    else {
      if( filenamePath != null ) {
        csFilename filename = new csFilename( filenamePath );
        // Reset bundle and dock pane titles if necessary (in case different data set is read into pane)
        bundle.resetFilename( filename.filename );
        setTitle( filename.filenamePath );
//        myDockPaneManager.resetActiveTitle( filename.filename, filename.filenamePath );
        if( renameActivePane ) myDockPaneManager.resetActiveTitle( filename.filename, filename.filenamePath );
      }
    }
    if( reinitializeReader ) bundle.initializeReader( reader, filenamePath, fileFormat );

    final Component glassPane = getGlassPane();
    glassPane.setCursor( Cursor.getPredefinedCursor( Cursor.WAIT_CURSOR ) );
    glassPane.addKeyListener( new KeyAdapter() {} );
    glassPane.addMouseListener( new MouseAdapter() {} );
    glassPane.setVisible( true );

    int numTotalTraces = bundle.getTotalNumTraces();
    csTraceSelectionParam param = bundle.getTraceSelectionParam();
    int numTracesToRead = param.numTraces;
    if( param.selectOption == csTraceSelectionParam.SELECT_ENS ) {
      numTracesToRead = param.numEns;
    }
    else if( param.selectOption == csTraceSelectionParam.SELECT_HEADER ) {
      numTracesToRead = param.numTracesHdr;
    }
    myProgressBar = new csProgressBarWindow( this, numTracesToRead, bundle, numTotalTraces );
    myProgressBar.setVisible(true);
    if( param.selectOption == csTraceSelectionParam.SELECT_ENS ) {
      myProgressBar.setTitle("Reading in ensembles...");
    }
    else {
      myProgressBar.setTitle("Reading in traces...");
    }
    csIDataSetListener listener = new csIDataSetListener() {
      @Override
      public void updateScan( csSeisPaneBundle bundle, int currentTraceIndex ) {
        if( myTraceSelectionDialog != null && bundle.isSameAs((csSeisPaneBundle)myDockPaneManager.getActivePanel()) ) {
          if( myTraceSelectionDialog.isVisible() ) myTraceSelectionDialog.updateScanProgress( currentTraceIndex );
        }
      }
      @Override
      public void stopScan( csSeisPaneBundle bundle, boolean success ) {
        if( myTraceSelectionDialog != null && bundle.isSameAs((csSeisPaneBundle)myDockPaneManager.getActivePanel()) ) {
          if( myTraceSelectionDialog.isVisible() ) myTraceSelectionDialog.stopScan( success );
        }
      }
      @Override
      public void updateTrace( int currentTraceIndex ) {
        myProgressBar.setValue( currentTraceIndex );
      }
      @Override
      public void updateTrace( int numDisplayedTraces, int currentTraceIndex ) {
        myProgressBar.setValues( numDisplayedTraces, currentTraceIndex );
      }
      @Override
      public void isComplete( boolean success ) {
        myProgressBar.dispose();
        if( glassPane != null ) {
          glassPane.setVisible( false );
        }
        if( !success && myCurrentNewDockPane != null ) {
          csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
          myDockPaneManager.closePane( myCurrentNewDockPane );
          if( bundle != null ) mySeisPaneManager.removeBundle(bundle);
        }
        else {
          csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
          updateMoveButtons( bundle );
          bundle.setMouseMode( myMouseMode );
          if( myTraceSelectionDialog != null ) {
            myTraceSelectionDialog.updateBundle( bundle );
          }
        }
        myIsReadProcessOngoing = false;
        myCurrentNewDockPane = null;
        if( myInputFilenameQueue != null ) openNextFileInQueue();
      }
    };
    bundle.read( true, listener, moveOption, resetScalar );
  }
  /**
   * Create a new seismic bundle object.
   * @param dispSettings Seismic display settings to be used for new object.
   * @return new seismic bundle object.
   */
  private csSeisPaneBundle createSeismicBundle( csSeisDispSettings dispSettings, csFilename filename, int fileFormat ) {
    csSeisView seisview = new csBackgroundSeisView( this, dispSettings, myCustomColorMapModel, myColorBitType );
    if( myProperties.showFilename ) {
      seisview.addOverlay( new csSeisViewFilenameOverlay( filename.filename ) );
    }
    csSeisPaneBundle bundleActive = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    csTraceSelectionParam traceSelectionParam = null;
    if( bundleActive != null ) traceSelectionParam = bundleActive.getTraceSelectionParam();
    csSeisPaneBundle bundle = mySeisPaneManager.add( seisview, traceSelectionParam, filename, fileFormat );
    bundle.addSeisPaneBundleListener( new csISeisPaneBundleListener() {
      @Override
      public void updateBundleDisplayScalar( csSeisPaneBundle bundle ) {
        updateSeismicDisplayScalar( bundle );
      }
      @Override
      public void updateBundleDisplaySettings( csSeisPaneBundle bundle ) {
        // Nothing to be done
      }
      @Override
      public void updateBundleSampleInfo( csSeisPaneBundle bundle, csSeisBundleSampleInfo sampleInfo ) {
        boolean printExtraInfo = sampleInfo.zoomVert > 1.0;
        String text = "";
        text +=    "Trace: " + intToString(sampleInfo.originalTrace);
        text += "   Sample: " + doubleToString1(sampleInfo.info.sampleDouble+1);
        if( printExtraInfo ) text += "(" + (sampleInfo.info.sample+1) + ")";
        if( !sampleInfo.isFrequencyDomain ) {
          text += "   Time: " + doubleToString4(sampleInfo.info.time);
          if( printExtraInfo ) text += "(" + doubleToString3(sampleInfo.timeFullSample) + ")";
        }
        else {
          text += "   Freq: " + doubleToString4(sampleInfo.info.time);
          if( printExtraInfo ) text += "(" + doubleToString3(sampleInfo.timeFullSample) + ")";
        }
        if( myShowAbsoluteTime ) {
          text += "   Date: " + computeAbsoluteTime( sampleInfo.time_samp1_ms, sampleInfo.info.time );
        }
        if( myShowFKVelocity ) {
          float fkVelocity = bundle.computeFKVelocity( (float)sampleInfo.info.traceDouble, (float)sampleInfo.info.time, myFKTraceSpacing );
          if( fkVelocity >= 0 ) text += "   FK velocity: " + fkVelocity;
        }
        text += "   Amplitude: " + amplitudeToString( sampleInfo.info.amplitude );
        myStatusBar.setText( text );
        if( bundle.isSynced() &&  mySeisPaneManager.isCrosshairSet() ) {
          mySeisPaneManager.setCrosshairPosition( bundle, sampleInfo.info );
        }
      }
    });
    return bundle;
  }
  /**
   * Set mouse mode. Applies to all seismic bundles/panes. The mouse mode determines what functionality is
   * available on mouse clicks and moves in the seismic views. Mouse modes are indicated by different mouse
   * cursors.
   * @param mouseMode Mouse modes include csMouseModes.NO_MODE, csMouseModes.ZOOM_MODE etc.
   * @param button The button which triggered the mouse mode change.
   */
  public void setMouseMode( int mouseMode, JToggleButton button ) {
    boolean isSelected = button.isSelected();
    if( isSelected ) {
      myToolBarLeft.resetToggleButtons();
      button.setSelected(isSelected);
      myMouseMode = mouseMode;
    }
    else {
      myMouseMode = csMouseModes.NO_MODE;
    }
    mySeisPaneManager.setMouseMode( myMouseMode );
  }
  /**
   * Implementation of csFileMenuListener interface.
   * @param e Object containing information about file menu event.
   */
  @Override
  public void fileSelected( csFileMenuEvent e ) {
    File file = e.file();
    boolean success = openFile( file.getAbsolutePath() );
    if( !success ) {
      myMenuBar.removeRecentFile( file.getAbsolutePath() );
    }
  }
  /**
   * Create String representation of absolute time & date.
   * @param time_samp1_ms  Millisecond fraction of current time.
   * @param time_s         Full UNIX seconds defining current time (=seconds since 1970-1-1).
   * @return 
   */
  public static String computeAbsoluteTime( long time_samp1_ms, double time_s ) {
    Calendar cal = Calendar.getInstance();
    java.util.TimeZone timeZone = java.util.TimeZone.getTimeZone("GMT");
    cal.setTimeZone( timeZone );
    cal.setTimeInMillis( time_samp1_ms + (long)(time_s*1000) );
    int year = cal.get( Calendar.YEAR );
    int month= cal.get( Calendar.MONTH ) + 1;  // Month 0 is January
    int day  = cal.get( Calendar.DAY_OF_MONTH );
    int hour = cal.get( Calendar.HOUR_OF_DAY );
    int min  = cal.get( Calendar.MINUTE );
    int sec  = cal.get( Calendar.SECOND );
    int msec = cal.get( Calendar.MILLISECOND );

    String text = "" + year + "-" + month + "-" + day + " " + hour + ":" + min + ":" + sec + "." + msec;

    return text;
  }
  //-----------------------------------------------------------------------------------------------------
  /**
   * Show absolute time at mouse location in status bar. Applies to all seismic bundles/panes.
   * @param doShow true if absolute time shall be shown in status bar.
   */
  public void showAbsoluteTime( boolean doShow ) {
    myShowAbsoluteTime = doShow;
  }
  //-----------------------------------------------------------------------------------------------------
  /**
   * Show computed velocity - Only works for Fk spectrum
   * @param doShow true if computed velocity shall be shown in status bar.
   */
  public void showVelocity( boolean doShow, float traceSpacing ) {
    myShowFKVelocity = doShow;
    myFKTraceSpacing = traceSpacing;
  }
//-----------------------------------------------------------------------------------------------------
  /**
   * Quit SeaView application. Save any unsaved changes to various SeaView properties and
   * default seismic display settings.
   */
  public void exitApplication() {
    if( myPreferences != null ) {
      myMenuBar.retrieveProperties(myProperties);
      if( mySeismicFileChooser != null && mySeismicFileChooser.getCurrentDirectory() != null ) {
        myProperties.dataDirectoryPath = mySeismicFileChooser.getCurrentDirectory().getAbsolutePath();
      }
      if( myFileChooserDump != null && myFileChooserDump.getCurrentDirectory() != null ) {
        myProperties.screenDumpDirectoryPath = myFileChooserDump.getCurrentDirectory().getAbsolutePath();
      }
      if( myVelocityAnalysisManager != null && myVelocityAnalysisManager.getVelFieldDir() != null ) {
        myProperties.velFieldDirectoryPath = myVelocityAnalysisManager.getVelFieldDir().getAbsolutePath();
      }
      myProperties.windowLayout = myDockPaneManager.getWindowLayout();
      mySeisPaneManager.retrieveAnnotationSettings( myProperties.annAttr, (csSeisPaneBundle)myDockPaneManager.getActivePanel() );
      myPreferences.setRecentFileList(myMenuBar.getRecentFileList());
      myPreferences.writePreferences();
      saveDefaultDispSettings( false );
      myProperties.segyAttr = mySegyAttr;
      try {
        myProperties.write();
      }
      catch( IOException e ) {
        // Nothing
      }
    }
    if( !myIsChildProcess ) System.exit(0);
  }
  /**
   * Convert amplitude to dB.
   * @param values  Amplitude values to convert.
   * @return Maximum amplitude.
   */
  public static double convertAmp2DB( double[] values ) {
    if( values.length == 0 ) return 0;
    double maxAmp = values[0];
    for( int isamp = 1; isamp < values.length; isamp++ ) {
      if( values[isamp] > maxAmp ) maxAmp = values[isamp];
    }
    for( int isamp = 0; isamp < values.length; isamp++ ) {
      values[isamp] = 20.0 * Math.log10(values[isamp]/maxAmp);
    }
    return maxAmp;
  }
  /**
   * Convert integer value to String.
   * @param value The value.
   * @return String representation of value.
   */
  public static String intToString( int value ) {
    if( value != csStandard.ABSENT_VALUE_INT ) {
      return String.valueOf(floatFmt0.format(value));
    }
    else {
      return "n/a";
    }
  }
  /**
   * Convert double value to String, showing 1 decimal place.
   * @param value The value.
   * @return String representation of value.
   */
  public static String doubleToString1( double value ) {
    if( value != csStandard.ABSENT_VALUE ) {
      return String.valueOf(floatFmt1.format(value));
    }
    else {
      return "n/a";
    }
  }
  /**
   * Convert double value to String, showing 3 decimal places.
   * @param value The value.
   * @return String representation of value.
   */
  public static String doubleToString3( double value ) {
    if( value != csStandard.ABSENT_VALUE ) {
      return String.valueOf(floatFmt3.format(value));
    }
    else {
      return "n/a";
    }
  }
  /**
   * Convert double value to String, showing 4 decimal places.
   * @param value The value.
   * @return String representation of value.
   */
  public static String doubleToString4( double value ) {
    if( value != csStandard.ABSENT_VALUE ) {
      return String.valueOf(floatFmt4.format(value));
    }
    else {
      return "n/a";
    }
  }
  /**
   * Convert amplitude value to String.
   * @param value The amplitude value.
   * @return String representation of amplitude value.
   */
  public static String amplitudeToString( double value ) {
    if( value != csStandard.ABSENT_VALUE ) {
      return formatAmplitude.format( value );
    }
    else {
      return "n/a";
    }
  }
  /**
   * Move forward in seismic data set. Load next traces or ensemble.
   */
  public void seismicForward() {
    seismicMove( SeaView.MOVE_OPTION_FORWARD );
  }
  /**
   * Move backward in seismic data set. Load previous traces or ensemble.
   */
  public void seismicBackward() {
    seismicMove( SeaView.MOVE_OPTION_BACKWARD );
  }
  /**
   * Move to beginning of seismic data set. Display first set of traces or ensemble.
   */
  public void seismicBegin() {
    seismicMove( SeaView.MOVE_OPTION_BEGIN );
  }
  /**
   * Move to end of seismic data set. Display last set of traces or ensemble.
   */
  public void seismicEnd() {
    seismicMove( SeaView.MOVE_OPTION_END );
  }
  /**
   * Move in seismic data set.
   * @param moveOption The 'move' option defines which traces to read in.
   */
  private void seismicMove( int moveOption ) {
    csSeisPaneBundle bundle = (csSeisPaneBundle)myDockPaneManager.getActivePanel();
    if( bundle == null || myIsReadProcessOngoing ) return;
    readSeismicBundle( bundle, moveOption );
  }
  /**
   * Update buttons which provide seismic move functionality (forward,backward,begin,end). Buttons are
   * enabled/disabled depending on where in the data set the current file pointer is located, and what
   * kind of trace selection is used for selecting traces.
   * @param bundle The current active seismic bundle object 
   */
  private void updateMoveButtons( csSeisPaneBundle bundle ) {
    boolean atStart = bundle.displayedDataAtStart();
    boolean atEnd   = bundle.displayedDataAtEnd();
    myToolBarTop.setMoveButtonsEnabled( atStart, atEnd );
  }

  //*********************************************************************************
  //*********************************************************************************
  //*********************************************************************************
  //
  /**
   * Main method, class SeaView.
   *<br>
   * Usage:   cseis/seaview/SeaView [filename]<br>
   *  <ul>
   *  <li> libraryPath  - Full path name of property java.library.path where libraries are located.
   *  <li> filename     - Name of file that shall be opened immediately
   *  </ul>
   * <br>
   * Requirements:<br>
   * <ul>
   * <li> Shared library 'csJNIlib' must be available within java.library.path.
   * </ul>
   */
  public static void main( String[] args ) {
    String path = System.getProperty( "java.library.path" );
    String libName = System.mapLibraryName( "csJNIlib" );
    try {
      System.load( path + java.io.File.separatorChar + libName );
    }
    catch( java.lang.UnsatisfiedLinkError e ) {
      JOptionPane.showMessageDialog( null,
          e.toString() + "\n" +
          "java.library.path = " + System.getProperty( "java.library.path" ) + "\n" +
          " - Seaview will not run.", "Error",
          JOptionPane.ERROR_MESSAGE );
      System.exit( -1 );
    }
    catch( Exception e2 ) {
      JOptionPane.showMessageDialog( null,
          e2.toString() + "\n" +
          "Unknown error occurred:" + e2.getMessage() + "\n" +
          " - Seaview will not run.", "Error",
          JOptionPane.ERROR_MESSAGE );
      System.exit( -1 );
    }
    
    try {
      UIManager.setLookAndFeel( "com.sun.java.swing.plaf.windows.WindowsLookAndFeel" );
    }
    catch( Exception e1 ) {
//      JOptionPane.showMessageDialog( null, e.toString(),
//          "Error", JOptionPane.ERROR_MESSAGE );
    }

    if( !csSeaViewProperties.fileExists() ) {
      String osName = System.getProperty("os.name");
      if( osName.startsWith("Windows") ) {
        csSeaViewLicenseDialog window = new csSeaViewLicenseDialog();
        window.setVisible(true);
        if( !window.acceptedLicense() ) {
          JOptionPane.showMessageDialog( window,
              "License terms not accepted. SeaView will not run.\n",
              "Info",
              JOptionPane.INFORMATION_MESSAGE );
          System.exit(0);
        }
//        if( window.acceptedSendEmail() ) {
//          JOptionPane.showMessageDialog( window,
//              "Email has been sent to " + window.getEmailAddress() + "\n",
//              "Info",
//              JOptionPane.INFORMATION_MESSAGE );
//        }
      } // END if osname
    }
    
    boolean createdNewDirectory = false;
    try {
      // Create CSEIS configuration folder, if it doesn't exist yet
      createdNewDirectory = cseis.general.csAbstractPreference.setCseisDirectory();
    }
    catch( Exception e ) {
      JOptionPane.showMessageDialog( null,
          e.toString() + "\nSeaView will not run.", "Error",
          JOptionPane.ERROR_MESSAGE );
      System.exit( -1 );
    }
    
    final SeaView seaview = new SeaView();
    seaview.setVisible(true);
    if( createdNewDirectory ) {
      JOptionPane.showMessageDialog( seaview,
          "The directory\n" + "'" + cseis.general.csAbstractPreference.getCseisDirectory() + "'\n" +
          "has been created on your system.\n\n" +
          "This directory is used to store configuration files.",
          "Info", JOptionPane.INFORMATION_MESSAGE );
    }
    final ArrayList<String> filenameList = new ArrayList<String>();
    for( int iarg = 0; iarg < args.length; iarg++ ) {
      if( args[iarg].length() > 0 ) {
        if( args[iarg].charAt(0) == '-' ) {
          String optionText = args[iarg].substring(1);
          if( optionText.compareToIgnoreCase("maximize") == 0 ) {
            seaview.setExtendedState( JFrame.MAXIMIZED_BOTH );
          }
        }
        else {
          filenameList.add(args[iarg]);
        }
      }
    }
    // Directly load specified file(s)
    if( filenameList.size() > 0 ) {
      SwingUtilities.invokeLater( new Runnable() {
        @Override
        public void run() {
          if( filenameList.size() == 1 ) {
            seaview.openFile( filenameList.get(0) );
          }
          else {
            seaview.openSeveralFiles(filenameList);
          }
        }
      });
    }
  }
}
