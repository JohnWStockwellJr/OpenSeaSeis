/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.seaview;

import javax.swing.ButtonGroup;
import javax.swing.JOptionPane;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.KeyEvent;
import java.awt.event.InputEvent;
import javax.swing.KeyStroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JMenuBar;
import javax.swing.JTextField;

import cseis.general.csColorMap;
import cseis.general.csStandard;
import cseis.processing.csProcessingAGC;
import cseis.processing.csProcessingDCRemoval;
import cseis.processing.csProcessingInterpolation;
import cseis.processing.csProcessingFilter;
import cseis.swing.csAboutDialog;
import cseis.swing.csRecentFileMenu;

/**
 * Menu bar of SeaView application.<br>
 * <br>
 * Defines all menu items and provides all functionality related to the menu bar.
 * @author 2007 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csMenuBar extends JMenuBar {
  private final static int MAX_RECENT_FILES = 10;
  
  private SeaView mySeaView;
  
  private JMenu menuFile;
  private JMenu menuView;
  private JMenu menuViewActive;
  private JMenu menuTools;
  private JMenu menuHelp;

  private JMenuItem menuOpen;
  private JMenu menuOpenFiletype;
  private JMenuItem menuOpenSegy;
  private JMenuItem menuOpenSegd;
  private JMenuItem menuOpenRSF;
  private JMenuItem menuOpenCseis;
  private JMenuItem menuOpenSU;
  private JMenuItem menuOpenASCII;
  private JMenuItem menuOpenJPG;
  private JMenuItem menuOpenBackground;

  private JMenuItem menuExit;
  private JMenuItem menuPreferences;
  private JMenuItem menuRefresh;

  private JMenu myMenuProcessing;
  private JMenuItem myMenuProcessingClear;
  private JMenuItem myMenuProcessingAGC;
  private JMenuItem myMenuProcessingFilter;
  private JMenuItem myMenuProcessingInterpolation;
  private JMenuItem myMenuProcessingDC;
  
  private JCheckBoxMenuItem menuHighlight;
  private JCheckBoxMenuItem menuCrosshairOn;
  private JCheckBoxMenuItem menuCrosshairSnap;

  private JMenuItem menuHdrMonitor;
  private JMenuItem menuTrcMonitor;
  private JMenuItem menuSettings;
  private JMenuItem menuAnnotation;
  private JMenuItem menuGraphSelection;
  private JMenuItem menuOverlay;
  private JMenuItem menuSnapshot;
  
  private JMenuItem menuCombineData;
  private JMenuItem menuColorMaps;
  private JMenuItem menuVelocityAnalysis;
  private JMenuItem menuWellPaths;

  private JCheckBoxMenuItem menuShowColorbar;
  private JCheckBoxMenuItem menuShowAbsoluteTime;
  private JCheckBoxMenuItem menuShowVelocity;
  private JCheckBoxMenuItem menuLogScaleY;

  private JMenuItem menuSelection;
  
  private JMenuItem menuAbout;  
  private JMenuItem menuImage;
  private JMenuItem menuImageAs;

  private JMenu menuDispSettings;
  private JMenuItem menuLoadDispSettings;
  private JMenuItem menuSaveDispSettings;
  private JMenuItem menuLoadDefaultDispSettings;
  private JMenuItem menuSaveDefaultDispSettings;
  private JMenuItem menuResetDispSettings;

  private csRecentFileMenu myMenuRecentFiles;
  private JRadioButtonMenuItem menuColors8bit;
  private JRadioButtonMenuItem menuColors32bit;

  private JMenu menuSEGY;
  private JMenuItem menuSEGYSetup;
  private JMenuItem menuSEGYCharHdr;
  private JMenuItem menuSEGYBinHdr;
  private JMenuItem menuSEGYTrcHdrMaps;

  public csMenuBar( SeaView seaview, csSeaViewProperties properties )
  {
    mySeaView = seaview;

    myMenuRecentFiles = new csRecentFileMenu("Recently opened files", MAX_RECENT_FILES);
    myMenuRecentFiles.addFileMenuListener(mySeaView);
    
    menuFile  = new JMenu("File");
    menuView  = new JMenu("View");
    menuViewActive = new JMenu("Active pane");
    menuTools = new JMenu("Tools");
    menuSEGY  = new JMenu("Segy/SU-Setup");
    menuHelp  = new JMenu("Help");

    menuOpen    = new JMenuItem("Open...");
    menuOpen.setAccelerator( KeyStroke.getKeyStroke(KeyEvent.VK_O, InputEvent.CTRL_DOWN_MASK) );
    menuImageAs   = new JMenuItem("Export as...");
    menuImageAs.setAccelerator( KeyStroke.getKeyStroke(KeyEvent.VK_E, InputEvent.SHIFT_DOWN_MASK | InputEvent.CTRL_DOWN_MASK) );
    menuImage   = new JMenuItem("Export");
    menuImage.setAccelerator( KeyStroke.getKeyStroke(KeyEvent.VK_E, InputEvent.CTRL_DOWN_MASK) );
    menuExit    = new JMenuItem("Exit");
    menuExit.setAccelerator( KeyStroke.getKeyStroke(KeyEvent.VK_Q, InputEvent.CTRL_DOWN_MASK) );
    menuPreferences = new JMenuItem("Preferences...");
    menuPreferences.setToolTipText("Set SeaView preferences");
    
    menuSelection = new JMenuItem("Trace selection...");
    menuSelection.setAccelerator( KeyStroke.getKeyStroke(KeyEvent.VK_T, InputEvent.CTRL_DOWN_MASK) );

    menuRefresh   = new JMenuItem("Refresh");
    menuRefresh.setAccelerator( KeyStroke.getKeyStroke(KeyEvent.VK_R, InputEvent.CTRL_DOWN_MASK) );
    menuRefresh.setToolTipText("Refresh/Re-read selected data set");

    menuCombineData = new JMenuItem("Combine data sets...");
    menuCombineData.setToolTipText("Perform operation on two data sets, create new data set");
    menuColorMaps = new JMenuItem("Manage color maps...");
    menuColorMaps.setToolTipText("Create/modify custom color maps...");
    menuWellPaths = new JMenuItem("Manage well paths...");
    menuWellPaths.setToolTipText("Display well paths on seismic sections");

    menuOpenFiletype  = new JMenu("Open file");
    menuOpenSegy    = new JMenuItem("Open SEG-Y file...check SEGY/SU-Setup options on menubar to make sure endian is correct!");
    menuOpenSegy.setToolTipText("Open SEG-Y file");
    menuOpenSegd    = new JMenuItem("Open SEG-D file...");
    menuOpenSegd.setToolTipText("Open SEG-D file");
    menuOpenRSF    = new JMenuItem("Open RSF file...");
    menuOpenRSF.setToolTipText("Open RSF (Seplib/Madagascar) file");
    menuOpenSU    = new JMenuItem("Open SU file...check SEGY/SU-Setup options on menubar to make sure endian is correct!");
    menuOpenSU.setToolTipText("Open SU file...");
    menuOpenCseis    = new JMenuItem("Open Seaseis file...");
    menuOpenCseis.setToolTipText("Open Seaseis (.cseis) file");
    menuOpenASCII    = new JMenuItem("Open ASCII file...");
    menuOpenASCII.setToolTipText("<html>Open ASCII text file containing a table of sample values.<br>" +
      "Each line contains a space-separated list of values.<br>" +
      "Each column is one trace, each line is one sample time.</html>");
    menuOpenJPG    = new JMenuItem("Open JPEG file...");
    menuOpenJPG.setToolTipText("Open JPEG image file as seismic");
    menuOpenBackground    = new JMenuItem("Open background...");
    menuOpenBackground.setToolTipText("Open background seismic/image/velocity for active pane");
    menuVelocityAnalysis    = new JMenuItem("Velocity analysis...");
    menuVelocityAnalysis.setToolTipText("<html>Velocity picking on semblance plot.<br>" +
            "<i>Note: A semblance plot must have been opened previously, consisting of one or several 1D/2D/3D CDP gathers.</i></html>");

    menuGraphSelection = new JMenuItem( csSeaViewActions.ACTION_TITLE[csSeaViewActions.ShowGraphAction], csSeaViewActions.getIcon(csSeaViewActions.ShowGraphAction) );
    menuGraphSelection.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.ShowGraphAction] );

    menuOverlay = new JMenuItem( csSeaViewActions.ACTION_TITLE[csSeaViewActions.ShowOverlayAction], csSeaViewActions.getIcon(csSeaViewActions.ShowOverlayAction) );
    menuOverlay.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.ShowOverlayAction] );

    menuSnapshot = new JMenuItem( "Open snap shot viewer...", csSeaViewActions.getIcon(csSeaViewActions.SnapShotAction) );
    menuSnapshot.setToolTipText( "Open snapshot viewer..." );

    menuDispSettings  = new JMenu("Display settings");

    menuLoadDispSettings = new JMenuItem("Load settings...");
    menuLoadDispSettings.setToolTipText("Load seismic display settings from external file");
    menuLoadDispSettings.setAccelerator( KeyStroke.getKeyStroke(KeyEvent.VK_L, InputEvent.CTRL_DOWN_MASK) );
    menuSaveDispSettings = new JMenuItem("Save settings...");
    menuSaveDispSettings.setToolTipText("Save seismic display settings to external file");
    menuSaveDispSettings.setAccelerator( KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_DOWN_MASK) );

    menuLoadDefaultDispSettings = new JMenuItem("Load settings (standard)");
    menuLoadDefaultDispSettings.setToolTipText("Load seismic display settings from standard file in home directory");
    menuSaveDefaultDispSettings = new JMenuItem("Save settings (standard)");
    menuSaveDefaultDispSettings.setToolTipText("Save seismic display settings to standard file in home directory");
    menuResetDispSettings = new JMenuItem("Reset");
    menuResetDispSettings.setToolTipText("Reset seismic display settings to fixed default");

    menuSEGYSetup = new JMenuItem("Segy/SU setup...");

    menuAbout   = new JMenuItem("About");   

    menuHighlight    = new JCheckBoxMenuItem("Highlight trace", false);
    menuCrosshairOn  = new JCheckBoxMenuItem("Show crosshair", false);
    menuCrosshairSnap= new JCheckBoxMenuItem("Snap crosshair", false);
    menuShowAbsoluteTime = new JCheckBoxMenuItem("Show absolute time", false);
    menuLogScaleY        = new JCheckBoxMenuItem("Log scale (Y axis)", false);
    menuShowColorbar     = new JCheckBoxMenuItem("Show colorbar", false);
    menuShowVelocity = new JCheckBoxMenuItem("Show velocity (FK only)", false);

    menuColors8bit  = new JRadioButtonMenuItem("8-bit colors",false);
    menuColors32bit = new JRadioButtonMenuItem("32-bit colors",false);
    ButtonGroup group2 = new ButtonGroup();
    group2.add(menuColors32bit);
    group2.add(menuColors8bit);
    menuColors32bit.setSelected( true );
    menuColors8bit.setToolTipText("8-bit color display: Use for example to speed up a slow remote connection");
    menuColors32bit.setToolTipText("Full 32-bit color display");
    
    menuHdrMonitor   = new JMenuItem("Monitor header values...");
    menuTrcMonitor   = new JMenuItem("Monitor sample values...");
    menuSettings     = new JMenuItem( "Display settings..." );
    menuAnnotation   = new JMenuItem( "Set trace annotation..." );

    menuSEGYCharHdr = new JMenuItem("Show SEG-Y char header");
    menuSEGYBinHdr  = new JMenuItem("Show SEG-Y binary header");
    menuSEGYTrcHdrMaps = new JMenuItem("Show SEG-Y trace header mapping");
    menuSEGYCharHdr.setEnabled(false);
    menuSEGYBinHdr.setEnabled(false);

    myMenuProcessing = new JMenu("Processing");
    myMenuProcessingClear = new JMenuItem("Clear processing");
    myMenuProcessingAGC = new JMenuItem("AGC");
    myMenuProcessingDC  = new JMenuItem("DC removal");
    myMenuProcessingFilter = new JMenuItem("Filter");
    myMenuProcessingInterpolation  = new JMenuItem("Interpolation");
    myMenuProcessing.add(myMenuProcessingClear);
    myMenuProcessing.addSeparator();
    myMenuProcessing.add(myMenuProcessingDC);
    myMenuProcessing.add(myMenuProcessingAGC);
    myMenuProcessing.add(myMenuProcessingFilter);
//    myMenuProcessing.add(myMenuProcessingInterpolation);
    
    menuDispSettings.add( menuLoadDispSettings );
    menuDispSettings.add( menuSaveDispSettings );
    menuDispSettings.addSeparator();
    menuDispSettings.add( menuLoadDefaultDispSettings );
    menuDispSettings.add( menuSaveDefaultDispSettings );
    menuDispSettings.addSeparator();
    menuDispSettings.add( menuResetDispSettings );

    menuFile.add(menuOpen);
    menuFile.add(menuOpenFiletype);
//    menuFile.add(menuOpenBackground);
    menuFile.addSeparator();
    menuFile.add(menuRefresh);
    menuFile.addSeparator();
    menuFile.add(menuImage);
    menuFile.add(menuImageAs);
    menuFile.addSeparator();
    menuFile.add(myMenuRecentFiles);
    menuFile.addSeparator();
    menuFile.add(menuDispSettings);
    menuFile.add(menuPreferences);
    menuFile.addSeparator();
    menuFile.add(menuExit);
    
    menuOpenFiletype.add(menuOpenCseis);
    menuOpenFiletype.add(menuOpenSegy);
    menuOpenFiletype.add(menuOpenSegd);
    menuOpenFiletype.add(menuOpenSU);
    menuOpenFiletype.add(menuOpenRSF);
    menuOpenFiletype.add(menuOpenASCII);
    menuOpenFiletype.add(menuOpenJPG);
    
    menuView.add(menuHighlight);
    menuView.add(menuCrosshairOn);
    menuView.add(menuCrosshairSnap);
    menuView.addSeparator();
    menuView.add(menuShowAbsoluteTime);
    menuView.add(menuShowColorbar);
    menuView.add(menuShowVelocity);
    menuView.addSeparator();
    menuView.add(menuSnapshot);
    menuView.addSeparator();
    menuView.add(menuColors8bit);
    menuView.add(menuColors32bit);

    menuViewActive.add(menuHdrMonitor);
    menuViewActive.add(menuTrcMonitor);
    menuViewActive.addSeparator();
    menuViewActive.add(menuSettings);
    menuViewActive.add(menuAnnotation);
    menuViewActive.addSeparator();
    menuViewActive.add(menuOverlay);
    menuViewActive.add(menuGraphSelection);
    menuViewActive.addSeparator();
    menuViewActive.add(menuLogScaleY);
    menuViewActive.addSeparator();
    menuViewActive.add(myMenuProcessing);
    
    menuTools.add( menuSelection );
    menuTools.addSeparator();
    menuTools.add(menuCombineData);
    menuTools.addSeparator();
    menuTools.add(menuColorMaps);
    menuTools.addSeparator();
    menuTools.add(menuWellPaths);
    menuTools.addSeparator();
    menuTools.add(menuVelocityAnalysis);

    menuSEGY.add(menuSEGYSetup);
    menuSEGY.add(menuSEGYTrcHdrMaps);
    menuSEGY.addSeparator();
    menuSEGY.add(menuSEGYCharHdr);
    menuSEGY.add(menuSEGYBinHdr);

    menuHelp.add(menuAbout);
    
    add(menuFile);
    add(menuView);
    add(menuViewActive);
    add(menuTools);
    add(menuSEGY);
    add(menuHelp);

    setToolTips();
    
    updateSelectedBundle( null );
    setMenuProperties( properties );

    menuExit.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
      String[] buttonLabels = {"Yes", "No"};
      int selectedOption = JOptionPane.showOptionDialog( mySeaView,
          "Are you sure you want to quit?\n",
          "Confirm exit",
          JOptionPane.YES_NO_OPTION,
          JOptionPane.QUESTION_MESSAGE,
          null, buttonLabels, buttonLabels[1] );
        if( selectedOption == JOptionPane.YES_OPTION ) {
          mySeaView.exitApplication();
        }
      }
    });
    menuPreferences.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showPreferencesDialog();
      }
    });
    myMenuProcessingDC.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.addProcessingStep( csProcessingDCRemoval.NAME );
      }
    });
    myMenuProcessingFilter.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.addProcessingStep( csProcessingFilter.NAME );
      }
    });
    myMenuProcessingClear.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.resetProcessing();
      }
    });
    myMenuProcessingAGC.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.addProcessingStep( csProcessingAGC.NAME );
      }
    });
    myMenuProcessingInterpolation.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.addProcessingStep( csProcessingInterpolation.NAME );
      }
    });
    menuSEGYSetup.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showSegySetupDialog();
      }
    });
    menuColors8bit.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.changeColorBits( csColorMap.COLOR_MAP_TYPE_8BIT);
      }
    });
    menuColors32bit.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.changeColorBits( csColorMap.COLOR_MAP_TYPE_32BIT );
      }
    });

    menuSelection.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showSelectionDialog();
      }
    });
    menuHighlight.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged( ItemEvent e ) {
        boolean isSelected = (e.getStateChange() == ItemEvent.SELECTED);
        mySeaView.highlightTrace( isSelected );
      }
    });
    menuShowAbsoluteTime.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged( ItemEvent e ) {
        boolean isSelected = (e.getStateChange() == ItemEvent.SELECTED);
        mySeaView.showAbsoluteTime( isSelected );
      }
    });
    menuShowVelocity.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged( ItemEvent e ) {
        final boolean isSelected = (e.getStateChange() == ItemEvent.SELECTED);
        if( !isSelected ) {
          mySeaView.showVelocity( isSelected, 1.0f );
        }
        else {
          final JDialog dialog = new JDialog( mySeaView, "Set trace spacing", true );
          final JTextField textfield = new JTextField("1.0");
  
          JPanel panelTop = new JPanel(new GridBagLayout());
          panelTop.setBorder( BorderFactory.createCompoundBorder(
              BorderFactory.createTitledBorder("Trace spacing"),
              csStandard.INNER_EMPTY_BORDER ) );
          
          panelTop.add( textfield, new GridBagConstraints(
              0, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
              GridBagConstraints.HORIZONTAL, new Insets( 0, 10, 0, 5 ), 0, 0 ) );
  
          JButton buttonClose = new JButton("Close");
          JPanel panelButtons = new JPanel( new GridBagLayout() );
          panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
              0, 0, 1, 1, 0.45, 1.0, GridBagConstraints.CENTER,
              GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
          panelButtons.add( buttonClose, new GridBagConstraints(
              1, 0, 1, 1, 0.45, 0.0, GridBagConstraints.SOUTHWEST,
              GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
          panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
              1, 0, 1, 1, 0.1, 1.0, GridBagConstraints.CENTER,
              GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
  
          JPanel panelAll = new JPanel(new BorderLayout());
          panelAll.setBorder( csStandard.DIALOG_BORDER );
  
          panelAll.add(panelTop,BorderLayout.CENTER);
          panelAll.add(panelButtons,BorderLayout.SOUTH);
          dialog.getContentPane().add(panelAll);
          dialog.pack();
          buttonClose.addActionListener( new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
              float traceSpacing = Float.parseFloat( textfield.getText() );
              mySeaView.showVelocity( true, traceSpacing );
              dialog.dispose();
            }
          });
          dialog.setLocationRelativeTo(mySeaView);
          dialog.setVisible(true);
        } // END if isSelected
      }
    });
    menuShowColorbar.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged( ItemEvent e ) {
        boolean isSelected = (e.getStateChange() == ItemEvent.SELECTED);
        mySeaView.showColorBar( isSelected );
      }
    });
    menuOverlay.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showOverlayDialog();
      }
    });
    menuGraphSelection.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showGraphPanel( true );
      }
    });
    
    menuSnapshot.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.openSnapshotViewer();
      }
    });
    menuCombineData.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.combineTwoSeismicBundles();
      }
    });
    menuWellPaths.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.manageWellPaths();
      }
     });
     menuColorMaps.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.manageCustomColorMaps();
      }
    });
    menuLogScaleY.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged( ItemEvent e ) {
        boolean isSelected = (e.getStateChange() == ItemEvent.SELECTED);
        mySeaView.setLogScaleY( isSelected );
      }
    });
    menuCrosshairOn.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged( ItemEvent e ) {
        mySeaView.setCrosshair( menuCrosshairOn.isSelected(), menuCrosshairSnap.isSelected() );
      }
    });
    menuCrosshairSnap.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged( ItemEvent e ) {
        mySeaView.setCrosshair( menuCrosshairOn.isSelected(), menuCrosshairSnap.isSelected() );
      }
    });
    menuSettings.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showSettingsDialog();
      }
    });
    menuHdrMonitor.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showHeaderMonitor();
      }
    });
    menuTrcMonitor.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showTraceMonitor();
      }
    });
    menuAnnotation.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showAnnotationDialog();
      }
    });

    menuOpen.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        if( mySeaView.openFile() ) {
//          setMenusEnabled(true);
        }
      }
    });
    menuOpenBackground.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
//        mySeaView.openBackgroundFile();
      }
    });
    menuVelocityAnalysis.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.openVelocityAnalysis(false);
      }
    });
    menuRefresh.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.refreshFile();
      }
    });
    menuOpenSegy.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.openFile(SeaView.FORMAT_SEGY);
      }
    });
    menuOpenSegd.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.openFile(SeaView.FORMAT_SEGD);
      }
    });
    menuOpenRSF.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.openFile(SeaView.FORMAT_RSF);
      }
    });
    menuOpenCseis.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.openFile(SeaView.FORMAT_CSEIS);
      }
    });
    menuOpenSU.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.openFile(SeaView.FORMAT_SU);
      }
    });
    menuOpenASCII.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.openFile(SeaView.FORMAT_ASCII);
      }
    });
    menuOpenJPG.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.openFile(SeaView.FORMAT_JPG);
      }
    });

    menuImage.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.createScreenDump( true );
      }
    });
    menuImageAs.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.createScreenDump( false );
      }
    });
    menuLoadDispSettings.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.loadDispSettings();
      }
    });
    menuSaveDispSettings.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.saveDispSettings();
      }
    });
    menuLoadDefaultDispSettings.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.applyDefaultDispSettings(true);
      }
    });
    menuSaveDefaultDispSettings.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.saveDefaultDispSettings(true);
      }
    });
    menuResetDispSettings.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.setDefaultDispSettings();
      }
    });
    menuSettings.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showSettingsDialog();
      }
    });
    menuAnnotation.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showAnnotationDialog();
      }
    });
    menuSEGYCharHdr.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showSegyCharHdr();
      }
    });
    menuSEGYBinHdr.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showSegyBinHdr();
      }
    });
    menuSEGYTrcHdrMaps.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showSegyTrcHdrMaps();
      }
    });

    menuAbout.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        csAboutDialog dialog = new csAboutDialog( mySeaView, "SeaView  (part of OpenSeaSeis)" );
        dialog.setVersionString( SeaView.VERSION );
        dialog.setContact("john@dix.mines.edu  ");
        dialog.setDate("2015  ");
        dialog.setLogo( cseis.resources.csResources.getIcon("seaview_icon_smooth.png") );
        dialog.setAdditionalComments("Copyright (c) Colorado School of Mines, 2015\n" +
                "All rights reserved.\n" +
                "Based on SeaSeis, developed by Bjorn Olofsson.");
        dialog.setVisible(true);
      }
    });
  }
//----------------------------------------------------------------
// Set menu items which were auto-saved in SeaView "properties" file
//
  public void setMenuProperties( csSeaViewProperties sp ) {
  }
  void setIsLogScale( boolean doSet ) {
    menuLogScaleY.setSelected( doSet );
  }
  public void retrieveProperties( csSeaViewProperties sp ) {
    if( menuColors8bit.isSelected() ) {
      sp.colorBitType = csColorMap.COLOR_MAP_TYPE_8BIT;
    }
    else {
      sp.colorBitType = csColorMap.COLOR_MAP_TYPE_32BIT;
    }
  }
  private void setToolTips() {
    menuOpen     .setToolTipText("Open seismic file");
    menuImage    .setToolTipText("Export to jpeg image file (default file name/location)");
    menuImageAs  .setToolTipText("Export to jpeg image file");
    menuExit     .setToolTipText("Exit");
    menuSelection.setToolTipText("Select traces to display");
    menuAbout    .setToolTipText("About");

    myMenuProcessingClear.setToolTipText( "Clear applied processing steps" );
    myMenuProcessingDC.setToolTipText( "Remove DC bias" );
    myMenuProcessingAGC.setToolTipText( "Apply AGC (automatic gain control)" );
    myMenuProcessingInterpolation.setToolTipText( "Apply interpolation (add traces)" );
    
    menuHighlight    .setToolTipText("Highlight trace at mouse position");
    menuCrosshairOn  .setToolTipText("Show crosshair");
    menuCrosshairSnap.setToolTipText("Snap crosshair to nearest sample/trace");
    menuHdrMonitor   .setToolTipText("Open header monitor window");
    menuTrcMonitor   .setToolTipText("Open trace sample monitor window");
    
    menuSettings.setToolTipText( "Set seismic display settings" );
    menuAnnotation.setToolTipText( "Select trace headers for trace annotation" );

    menuSEGYCharHdr   .setToolTipText("Show SEG-Y char header of selected data set");
    menuSEGYBinHdr    .setToolTipText("Show SEG-Y binary header of selected data set");
    menuSEGYTrcHdrMaps.setToolTipText("Show SEG-Y trace header mappings");
    menuSEGYSetup     .setToolTipText("Edit general SEG-Y/SU file settings");
    
    menuShowAbsoluteTime.setToolTipText("Compute absolute time at mouse position, show in status bar");
    menuShowColorbar.setToolTipText("<html>Show color bar.<br>" +
            "<i>Note: Only fully works for variable color maps, and range/trace scaling.</i></html>");
    menuLogScaleY.setToolTipText("<html>Plot Y axis in log scale<br>" +
            "<i>Minimum time/frequency is set automatically depending on zoom level.</i></html>");
  }
  public void removeRecentFile( String filenameFullPath ) {
    myMenuRecentFiles.removeFile( filenameFullPath );
  }
  public void addRecentFile( String filenameFullPath ) {
    myMenuRecentFiles.addFile( filenameFullPath );
  }
  public java.util.List<String> getRecentFileList() {
    return myMenuRecentFiles.getFileList();
  }

  public void updateSelectedBundle( csSeisPaneBundle bundle ) {
    boolean doSet = ( bundle != null );
    if( doSet ) {
      setIsLogScale( bundle.seisView.isLogScale() );
    }
    menuRefresh.setEnabled(doSet);
    menuViewActive.setEnabled(doSet);
    menuImage.setEnabled(doSet);
    menuImageAs.setEnabled(doSet);
    menuLoadDispSettings.setEnabled(doSet);
    menuSaveDispSettings.setEnabled(doSet);
    menuLoadDefaultDispSettings.setEnabled(doSet);
    menuSaveDefaultDispSettings.setEnabled(doSet);
    menuResetDispSettings.setEnabled(doSet);

    boolean isSEGY = ( bundle != null && bundle.getFileFormat() == SeaView.FORMAT_SEGY );
    menuSEGYCharHdr.setEnabled(isSEGY);
    menuSEGYBinHdr.setEnabled(isSEGY);
  }
}


