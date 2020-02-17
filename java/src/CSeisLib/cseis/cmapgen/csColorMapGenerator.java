/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.cmapgen;

import cseis.general.csTimerTaskIncrement;
import cseis.general.csITimerTaskIncrementListener;
import cseis.general.csITimerTaskDecrementListener;
import cseis.general.csTimerTaskDecrement;
import cseis.general.csColorBarPanel;
import cseis.general.csColorMap;
import cseis.general.csCustomColorMap;
import cseis.general.csStandard;
import cseis.seisdisp.csColorMapListItem;
import cseis.seisdisp.csComboColorMapRenderer;
import cseis.swing.csColorButton;
import cseis.swing.csColorChangeListener;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Timer;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JColorChooser;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListDataEvent;
import javax.swing.event.ListDataListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

/**
 * Color map generator.<br>
 * JFrame providing functionality to create, delete and edit custom color maps. Custom color maps
 * may be loaded from and stored to an external ASCII file.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csColorMapGenerator extends JFrame implements ChangeListener {
  public static double POSITION_INCREMENT = 0.01;
  private JPanel myPanelEditor;
  private csColorBarPanel myPanelColorBar;

  private int myMaxUndoElements = 50;
  private java.util.Stack<csCustomColorMap> myUndoStack;
  private JButton myButtonUndo;

  private JMenuBar myMenuBar;
  private JMenu myMenu;
  private JMenuItem myMenuOpen;
  private JMenuItem myMenuSave;
  private JMenuItem myMenuSaveAll;
  private JMenuItem myMenuTransfer;
  private JMenuItem myMenuClose;
  private JFileChooser myFileChooser;

  private JButton myButtonTestApply;
  private JButton myButtonCreate;
  private JButton myButtonUpdateSelected;
  private JButton myButtonRemoveSelected;
  private JButton myButtonAddColor;
  private JButton myButtonInterpColors;
  private JButton myButtonInterpPositions;
  private JTextField myTextName;
  
  private ArrayList<csColorMapGenerator.Bundle> myBundleList;
  private DecimalFormat myFormat0000;
  private JColorChooser myColorChooser;
  private Bundle mySelectedBundle;
  private boolean myIsUpdating;
  
  private csManagerCustomMapList myManagerCustomMaps;
  private JCheckBox myBoxManagerCustom;
  private JCheckBox myBoxManagerDefault;
  private JButton myButtonSelectCustom;
  private JButton myButtonSelectDefault;
  private JComboBox<csColorMapListItem> myComboDefaultMaps;
  private boolean myUpdatesWaiting = false;
  
  private ArrayList<csIColorMapGeneratorListener> myListeners;
  private String myDefaultDirectory;
  
  public csColorMapGenerator( Window window, java.util.List<csCustomColorMap> customColorMapList, String defaultDirectory ) {
    super("Color map generator");

    myListeners = new ArrayList<csIColorMapGeneratorListener>();
    myDefaultDirectory = defaultDirectory;

    myUndoStack  = new java.util.Stack<csCustomColorMap>();
    myButtonUndo = new JButton("Undo");
    myButtonUndo.setEnabled(false);
    myButtonUndo.setToolTipText("Undo last edit, recover previously edited color map");

    myButtonTestApply = new JButton("Test");
    myButtonCreate    = new JButton("Create");
    myButtonTestApply.setToolTipText("Test currently edited color map");
    myButtonCreate.setToolTipText("Create new custom color map");
    this.getRootPane().setDefaultButton(myButtonTestApply);

    myComboDefaultMaps = new JComboBox<csColorMapListItem>();
    myComboDefaultMaps.setRenderer( new csComboColorMapRenderer() );
    for( int imap = 0; imap < csColorMap.NUM_DEFAULT_MAPS; imap++ ) {
      csColorMap map = new csColorMap(imap,csColorMap.COLOR_MAP_TYPE_32BIT);
      myComboDefaultMaps.addItem( csColorMapListItem.createStandardItem( map ) );
    }

    myBoxManagerDefault = new JCheckBox("Auto-select",false);
    myBoxManagerDefault.setToolTipText("Replace edited color map with colors from selected default map");
    myButtonSelectDefault = new JButton("Select");
    myButtonSelectDefault.setToolTipText("Replace edited color map with colors from selected default map");

    csColorMapListModel listModel = new csColorMapListModel( customColorMapList );
    
    myManagerCustomMaps = new csManagerCustomMapList( listModel );
    myBoxManagerCustom = new JCheckBox("Auto-select",false);
    myBoxManagerCustom.setToolTipText("Replace edited color map with colors from selected custom map");
    myButtonSelectCustom = new JButton("Select");
    myButtonSelectCustom.setToolTipText("Replace edited color map with colors from selected custom map");
    
    myButtonRemoveSelected = new JButton("Remove");
    myButtonUpdateSelected  = new JButton("Update");
    myButtonUpdateSelected.setToolTipText("Update selected custom color map");
    myButtonRemoveSelected.setToolTipText("Remove selected custom map");

    myButtonRemoveSelected.setEnabled(false);
    myButtonUpdateSelected.setEnabled(false);
    
    myIsUpdating = false;
    mySelectedBundle = null;
    myColorChooser = new JColorChooser();
    myColorChooser.getSelectionModel().addChangeListener( this );
    myColorChooser.setPreviewPanel(new JPanel());  // Remove preview panel
    csCustomColorMap cmap = new csCustomColorMap();
    myPanelColorBar = new csColorBarPanel( cmap, csColorBarPanel.ORIENT_VERTICAL, csColorBarPanel.ANNOTATION_ADVANCED );
    myPanelColorBar.setMinMax(0.0, 1.0);
    myPanelColorBar.setPreferredSize( new Dimension(100, myPanelColorBar.getPreferredSize().height) );
    myPanelColorBar.setAutoSize(true);
    myPanelColorBar.setShowKneePoints(true);

    myMenuBar = new JMenuBar();
    myMenu = new JMenu("File");
    myMenuOpen = new JMenuItem("Open...");
    myMenuSave = new JMenuItem("Save...");
    myMenuSaveAll = new JMenuItem("Save all...");
    myMenuTransfer = new JMenuItem("Save changes");
    myMenuClose = new JMenuItem("Close");
    myMenu.add(myMenuOpen);
    myMenu.addSeparator();
    myMenu.add(myMenuTransfer);
    myMenu.addSeparator();
    myMenu.add(myMenuSave);
    myMenu.add(myMenuSaveAll);
    myMenu.addSeparator();
    myMenu.add(myMenuClose);
    myMenuBar.add(myMenu);
    super.setJMenuBar(myMenuBar);

    myMenuOpen.setToolTipText("Load color map from external ASCII file");
    myMenuSave.setToolTipText("Write selected custom color map to external ASCII file");
    myMenuSaveAll.setToolTipText("Write all custom color maps to external ASCII files");
    myMenuTransfer.setToolTipText("Save/transfer changes to main application");
    myMenuClose.setToolTipText("Close this dialog");

    myTextName = new JTextField( cmap.toString() );

    myFormat0000 = new DecimalFormat("0.0000");
    myBundleList = new ArrayList<csColorMapGenerator.Bundle>();
    csColorMapGenerator.Bundle bundle = new csColorMapGenerator.Bundle( Color.black, 1.0 );
    bundle.buttonRemove.setEnabled(false);
    myBundleList.add(bundle);
    bundle = new csColorMapGenerator.Bundle( Color.white, 0.0 );
    bundle.buttonRemove.setEnabled(false);
    myBundleList.add(bundle);
    
    JTextField textDummy = new JTextField("1.0");
    int preferredHeight = textDummy.getPreferredSize().height;
    
    myButtonAddColor = new JButton(cseis.resources.csResources.getIcon("icon_plus.png"));
    myButtonAddColor.setText("");
    myButtonAddColor.setMargin(new Insets(0,0,1,1));
    myButtonAddColor.setPreferredSize(new Dimension(preferredHeight,preferredHeight));

    myButtonInterpColors = new JButton("Color");
    myButtonInterpPositions = new JButton("Position");
    myButtonInterpColors.setToolTipText("Interpolate colors between first and last knee point");
    myButtonInterpPositions.setToolTipText("Distribute positions evenly between knee points");
    int height = myButtonInterpColors.getPreferredSize().height;
    int width = Math.max( myButtonInterpColors.getPreferredSize().width, myButtonInterpPositions.getPreferredSize().width );
    myButtonInterpColors.setPreferredSize( new Dimension(width,height) );
    myButtonInterpPositions.setPreferredSize( new Dimension(width,height) );
    
    myPanelEditor = new JPanel(new GridBagLayout());
    myPanelEditor.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Color map editor"),
        csStandard.INNER_EMPTY_BORDER ) );
    updateEditorPanel();
    JScrollPane paneEditor = new JScrollPane( myPanelEditor, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
    paneEditor.getVerticalScrollBar().setUnitIncrement( 50 );

    JPanel panelName = new JPanel(new GridBagLayout());
    panelName.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Color map name"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JPanel panelInterpolate = new JPanel(new GridBagLayout());
    panelInterpolate.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Editor interpolation"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JPanel panelFinalize = new JPanel(new GridBagLayout());
    panelFinalize.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Test current color map"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JPanel panelUndo = new JPanel(new GridBagLayout());
    panelUndo.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Undo operation"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JPanel panelSelect = new JPanel(new GridBagLayout());
    panelSelect.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Color map operation"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    panelName.add( myTextName, new GridBagConstraints(
        0, 0, 1, 1, 1.0, 0.0, GridBagConstraints.NORTHWEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    
    panelInterpolate.add( myButtonInterpColors, new GridBagConstraints(
        0, 0, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );
    panelInterpolate.add( myButtonInterpPositions, new GridBagConstraints(
        1, 0, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );

    panelFinalize.add( Box.createHorizontalGlue(), new GridBagConstraints(
        0, 0, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );
    panelFinalize.add( myButtonTestApply, new GridBagConstraints(
        1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 0, 5, 5 ), 0, 0 ) );
    panelFinalize.add( Box.createHorizontalGlue(), new GridBagConstraints(
        2, 0, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );

    panelUndo.add( Box.createHorizontalGlue(), new GridBagConstraints(
        0, 0, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );
    panelUndo.add( myButtonUndo, new GridBagConstraints(
        1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 0, 5, 5 ), 0, 0 ) );
    panelUndo.add( Box.createHorizontalGlue(), new GridBagConstraints(
        2, 0, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );

    panelSelect.add( myButtonCreate, new GridBagConstraints(
        0, 0, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );
    panelSelect.add( myButtonUpdateSelected, new GridBagConstraints(
        1, 0, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );
    panelSelect.add( Box.createHorizontalGlue(), new GridBagConstraints(
        0, 1, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );
    panelSelect.add( myButtonRemoveSelected, new GridBagConstraints(
        1, 1, 1, 1, 0.5, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 5 ), 0, 0 ) );

    JPanel panelRight = new JPanel( new GridBagLayout() );
    int yp = 0;
    panelRight.add( panelName, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.NORTH,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelRight.add( panelInterpolate, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelRight.add( panelSelect, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelRight.add( panelUndo, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelRight.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.5, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelRight.add( panelFinalize, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
//    panelRight.add( Box.createVerticalGlue(), new GridBagConstraints(
//        0, yp++, 1, 1, 1.0, 0.5, GridBagConstraints.CENTER,
//        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    
    JScrollPane paneManagerCustom  = new JScrollPane(myManagerCustomMaps,
            JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
    Dimension size = paneManagerCustom.getPreferredSize();
    paneManagerCustom.setPreferredSize(size);

    JPanel panelManagerDefault = new JPanel( new BorderLayout() );
    panelManagerDefault.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Default color maps"),
        csStandard.INNER_EMPTY_BORDER ) );
    JPanel panelNorthDefault = new JPanel( new FlowLayout( FlowLayout.LEFT ) );
    panelNorthDefault.add( myBoxManagerDefault );
    panelNorthDefault.add( myButtonSelectDefault );
    panelManagerDefault.add( panelNorthDefault, BorderLayout.NORTH );
    panelManagerDefault.add( myComboDefaultMaps, BorderLayout.CENTER );    

    JPanel panelManagerCustom = new JPanel( new BorderLayout() );
    panelManagerCustom.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Custom color maps"),
        csStandard.INNER_EMPTY_BORDER ) );
    JPanel panelNorthCustom = new JPanel( new FlowLayout( FlowLayout.LEFT ) );
    panelNorthCustom.add( myBoxManagerCustom );
    panelNorthCustom.add( myButtonSelectCustom );
    panelManagerCustom.add( panelNorthCustom, BorderLayout.NORTH );
    panelManagerCustom.add( paneManagerCustom, BorderLayout.CENTER );    

    //------------------------------------------------------------------
    
    JPanel panelManagers = new JPanel( new BorderLayout() );
    panelManagers.add( panelManagerCustom, BorderLayout.CENTER );
    panelManagers.add( panelManagerDefault, BorderLayout.SOUTH );

    JPanel panelEast = new JPanel(new BorderLayout());
    panelEast.add( panelRight,BorderLayout.WEST );
    panelEast.add( panelManagers, BorderLayout.CENTER );
    panelEast.add( myColorChooser, BorderLayout.SOUTH );
    
    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );
    panelAll.add( paneEditor,BorderLayout.CENTER );
    panelAll.add( panelEast,BorderLayout.EAST );
//    panelAll.add(panelButtons,BorderLayout.SOUTH);

    
    getContentPane().add( myPanelColorBar,BorderLayout.WEST );
    getContentPane().add( panelAll,BorderLayout.CENTER );
//    getContentPane().add( paneEditor,BorderLayout.EAST );
//    getContentPane().add( panelManagers,BorderLayout.EAST );

    pack();
    setLocationRelativeTo( window );

    //------------------------------------------------------------
    //
    myComboDefaultMaps.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if( myBoxManagerDefault.isSelected() ) {
          csColorMap cmap = ((csColorMapListItem)myComboDefaultMaps.getSelectedItem()).map;
          updateColorMap( cmap );
        }
      }
    });
    myManagerCustomMaps.addListDataListener( new ListDataListener() {
      @Override
      public void intervalAdded( ListDataEvent e ) {
        updateCustomListSelection();
      }
      @Override
      public void intervalRemoved(ListDataEvent e) {
        updateCustomListSelection();
      }
      @Override
      public void contentsChanged(ListDataEvent e) {
        updateCustomListSelection();
      }
    }); 
    myManagerCustomMaps.addListSelectionListener( new ListSelectionListener() {
      @Override
      public void valueChanged(ListSelectionEvent e) {
        updateCustomListSelection();
        if( e.getValueIsAdjusting() ) return;
        if( myBoxManagerCustom.isSelected() && !myManagerCustomMaps.isSelectionEmpty() ) {
          csColorMap cmap = myManagerCustomMaps.getSelectedMap();
          if( cmap != null ) updateColorMap( cmap );
        }
      }
    });
    myBoxManagerCustom.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        myButtonSelectCustom.setEnabled( !myBoxManagerCustom.isSelected() );
        if( myBoxManagerCustom.isSelected() ) {
          csColorMap cmap = myManagerCustomMaps.getSelectedMap();
          if( cmap != null ) updateColorMap( cmap );
        }
      }
    });
    myButtonSelectCustom.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        csColorMap cmap = myManagerCustomMaps.getSelectedMap();
        if( cmap != null ) updateColorMap( cmap );
      }
    });
    myBoxManagerDefault.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        myButtonSelectDefault.setEnabled( !myBoxManagerDefault.isSelected() );
        if( myBoxManagerDefault.isSelected() ) {
          csColorMap cmap = ((csColorMapListItem)myComboDefaultMaps.getSelectedItem()).map;
          updateColorMap( cmap );
        }
      }
    });
    myButtonSelectDefault.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        csColorMap cmap = ((csColorMapListItem)myComboDefaultMaps.getSelectedItem()).map;
        updateColorMap( cmap );
      }
    });

    myMenuOpen.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        openColorMap();
      }
    });
    myMenuSave.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        saveColorMap();
      }
    });
    myMenuSaveAll.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        saveAllColorMaps();
      }
    });
    myMenuTransfer.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        fireUpdateMapsEvent();
        JOptionPane.showMessageDialog( csColorMapGenerator.this,
                    "Custom maps were applied to main application.\n\n" +
                    "NOTE: To save color maps permanently, select 'Save all..' from menu.",
                    "Info", JOptionPane.INFORMATION_MESSAGE );
      }
    });
    myMenuClose.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if( myUpdatesWaiting ) {
          int option = JOptionPane.showConfirmDialog( csColorMapGenerator.this,
                      "Shall changes to custom color maps be applied to main application?\n\n" +
                      "NOTE: To save custom color maps permanently, cancel this dialog\n" +
                      "and select 'Save all..' from menu.",
                      "Confirm", JOptionPane.YES_NO_CANCEL_OPTION );
          if( option == JOptionPane.YES_OPTION ) {
            fireUpdateMapsEvent();
//            JOptionPane.showMessageDialog( csColorMapGenerator.this,
//                      "Custom maps were applied to main application.",
//                      "Info", JOptionPane.INFORMATION_MESSAGE );
          }
          else if( option == JOptionPane.CANCEL_OPTION ) {
            return;
          }
        }
        cancel();
      }
    });
    myTextName.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        myPanelColorBar.getColorMap().setName(myTextName.getText().trim());
      }
    });
    myTextName.addFocusListener( new FocusListener() {
      @Override
      public void focusGained(FocusEvent e) {
      }
      @Override
      public void focusLost(FocusEvent e) {
        myPanelColorBar.getColorMap().setName(myTextName.getText().trim());
      }
    });
    myButtonAddColor.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        addBundleAtEnd();
      }
    });
    myButtonUndo.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        performUndoOperation();
      }
    });
    myButtonInterpColors.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        distributeColorsEvenly();
      }
    });
    myButtonInterpPositions.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        distributePositionsEvenly();
      }
    });
    myButtonTestApply.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
        fireApplyEvent();
      }
    });
    myButtonCreate.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
        myManagerCustomMaps.addMap( myPanelColorBar.getColorMap() );
      }
    });
    myButtonUpdateSelected.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        csCustomColorMap cmap = myManagerCustomMaps.getSelectedMap();
        if( cmap == null ) return;
        int option = JOptionPane.showConfirmDialog( csColorMapGenerator.this,
                    "Are you sure you want to update custom color map\n" +
                    "  '" + cmap.toString() + "'\n" +
                    "with the specified colors & map name?",
                    "Confirm update", JOptionPane.YES_NO_OPTION );
        if( option == JOptionPane.YES_OPTION ) {
          apply();
          myManagerCustomMaps.updateColorMap( myManagerCustomMaps.getSelectedMapIndex(), myPanelColorBar.getColorMap() );
        }
      }
    });
    myButtonRemoveSelected.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        csCustomColorMap cmap = myManagerCustomMaps.getSelectedMap();
        if( cmap == null ) return;
        int option = JOptionPane.showConfirmDialog( csColorMapGenerator.this,
                    "Are you sure you want to remove custom color map " + cmap.toString() + "?",
                    "Confirm deletion", JOptionPane.YES_NO_OPTION );
        if( option == JOptionPane.YES_OPTION ) {
          myManagerCustomMaps.removeSelectedMap();
        }
      }
    });
    myUpdatesWaiting = false;
  }
  public void setTestButton( boolean enable, String toolTip ) {
    myButtonTestApply.setToolTipText(toolTip);
    myButtonTestApply.setEnabled(enable);
  }
  private void addToUndoStack( csCustomColorMap cmap ) {
    myUndoStack.push( new csCustomColorMap(cmap) );
    updateUndoOperation();
  }
  private void performUndoOperation() {
    if( !myUndoStack.isEmpty() ) {
      csCustomColorMap cmap = myUndoStack.pop();
      updateColorMap( cmap );
    }
    updateUndoOperation();
  }
//  private void emptyUndoStack() {
//    myUndoStack.clear();
//    updateUndoOperation();
//  }
  private void updateUndoOperation() {
    if( myUndoStack.size() > myMaxUndoElements ) myUndoStack.setSize(myMaxUndoElements);
    myButtonUndo.setEnabled( !myUndoStack.isEmpty() );
  }
  private void updateCustomListSelection() {
    boolean isSelection = myManagerCustomMaps.getSelectedMapIndex() >= 0;
    myButtonUpdateSelected.setEnabled( isSelection );
    myButtonRemoveSelected.setEnabled( isSelection );
    myUpdatesWaiting = true;
  }
  public void addListDataListener( ListDataListener listener ) {
    myManagerCustomMaps.addListDataListener( listener);
  }
  public void removeListDataListener( ListDataListener listener ) {
    myManagerCustomMaps.removeListDataListener( listener);
  }
  public void addListener( csIColorMapGeneratorListener listener ) {
    myListeners.add(listener);
  }
  public void removeListener( csIColorMapGeneratorListener listener ) {
    myListeners.remove(listener);
  }
  private void openColorMap() {
    if( myFileChooser == null ) myFileChooser = new JFileChooser( myDefaultDirectory );
    int option = myFileChooser.showOpenDialog( this );
    if( option == JFileChooser.APPROVE_OPTION ) {
      File file = myFileChooser.getSelectedFile();
      try {
        BufferedReader reader = new BufferedReader( new FileReader(file) );
        myPanelColorBar.getColorMap().load(reader);
        myPanelColorBar.repaint();
        resetBundles();
        myTextName.setText(myPanelColorBar.getColorMap().toString());
      }
      catch( FileNotFoundException ex ) {
        JOptionPane.showMessageDialog(
                this,
                "File " + file.getName() + ": \n" + ex.getMessage(),
                "Error", JOptionPane.ERROR_MESSAGE);
        return;
      }
      catch( IOException ex ) {
        JOptionPane.showMessageDialog(
                this,
                "File " + file.getName() + ": \n" + ex.getMessage(),
                "Error", JOptionPane.ERROR_MESSAGE);
        return;
      }
    }
  }
  private void resetBundles() {
    myBundleList.clear();
    csCustomColorMap cmap = myPanelColorBar.getColorMap();
    int numKneePoints = cmap.getNumKneePoints();
    for( int i = 0; i < numKneePoints; i++ ) {
      csColorMapGenerator.Bundle bundle = new csColorMapGenerator.Bundle(
              cmap.getColorKneePoint(numKneePoints-1-i), cmap.getWeightKneePoint(numKneePoints-1-i) );
      myBundleList.add(bundle);
    }
    myBundleList.get(0).buttonRemove.setEnabled(false);
    myBundleList.get(numKneePoints-1).buttonRemove.setEnabled(false);
    updateEditorPanel();
  }
  private void updateColorMap( csColorMap cmap ) {
    if( cmap instanceof csCustomColorMap ) myTextName.setText( cmap.toString() );
    myPanelColorBar.updateColors( cmap.getColorKneePoints(), cmap.getWeightKneePoints() );
    myPanelColorBar.repaint();
    resetBundles();
  }
  public static void saveColorMap( csCustomColorMap cmap, String filename ) throws IOException, FileNotFoundException {
    File file = new File( filename );
    BufferedWriter writer = new BufferedWriter( new FileWriter(file) );
    cmap.write(writer);
    writer.close();
  }
  protected void saveColorMap() {
    if( myFileChooser == null ) myFileChooser = new JFileChooser( myDefaultDirectory );
    int option = myFileChooser.showSaveDialog( this );
    if( option == JFileChooser.APPROVE_OPTION ) {
      File file = myFileChooser.getSelectedFile();
      try {
        csColorMapGenerator.saveColorMap( myPanelColorBar.getColorMap(), file.getAbsolutePath() );
      }
      catch( FileNotFoundException ex ) {
        JOptionPane.showMessageDialog(
                this,
                "File " + file.getName() + ": \n" + ex.getMessage(),
                "Error", JOptionPane.ERROR_MESSAGE);
        return;
      }
      catch( IOException ex ) {
        JOptionPane.showMessageDialog(
                this,
                "File " + file.getName() + ": \n" + ex.getMessage(),
                "Error", JOptionPane.ERROR_MESSAGE);
        return;
      }
      JOptionPane.showMessageDialog(
              this,
              "Saved color map to file...\n  " + file.getName(),
              "Info", JOptionPane.INFORMATION_MESSAGE);
    }
  }
  private void saveAllColorMaps() {
    if( myManagerCustomMaps.getNumMaps() == 0 ) return;
    ArrayList<csCustomColorMap> cmapList = new ArrayList<csCustomColorMap>();
    for( int imap = 0; imap < myManagerCustomMaps.getNumMaps(); imap++ ) {
      cmapList.add( myManagerCustomMaps.getMap(imap) );
    }
    final csColorMapSaveDialog dialog = new csColorMapSaveDialog( this, myDefaultDirectory, cmapList );
    dialog.setVisible( true );
  }
  private void distributePositionsEvenly() {
    myIsUpdating = true;
    int numBundles = myBundleList.size();
    double dpos = 1.0 / (double)( numBundles - 1 );
    for( int ibundle = 0; ibundle < numBundles; ibundle++ ) {
      double position = dpos * (double)( ibundle );
      if( ibundle == 0 ) position = 0.0;
      else if( ibundle == numBundles-1 ) position = 1.0;
      myBundleList.get(myBundleList.size()-1-ibundle).setPosition(position);
    }
    apply();
    myIsUpdating = false;
  }
  private void distributeColorsEvenly() {
    myIsUpdating = true;
    Color[] colors = new Color[2];
    double[] weights = new double[2];
    colors[0] = myPanelColorBar.getColorMap().getColor(0);
    colors[1] = myPanelColorBar.getColorMap().getColor(1);
    weights[0] = 0;
    weights[1] = 1;
    csColorMap cmap = new csColorMap( colors, weights, csColorMap.NUM_COLORS_DEFAULT );
    int numBundles = myBundleList.size();
    double dpos = 1.0 / (double)( numBundles - 1 );
    for( int ibundle = 0; ibundle < numBundles; ibundle++ ) {
      double position = dpos * (double)( ibundle );
      if( ibundle == 0 ) position = 0.0;
      else if( ibundle == numBundles-1 ) position = 1.0;
      Color color = cmap.getColor( (float) position );
      myBundleList.get(myBundleList.size()-1-ibundle).setColor(color);
    }
    apply();
    myIsUpdating = false;
  }
  private boolean apply() {
    Color[] colors = new Color[myBundleList.size()];
    double[] positions = new double[myBundleList.size()];
    for( int ibundle = 0; ibundle < myBundleList.size(); ibundle++ ) {
      colors[myBundleList.size()-1-ibundle]    = myBundleList.get(ibundle).getColor();
      positions[myBundleList.size()-1-ibundle] = myBundleList.get(ibundle).getPosition();
    }
    addToUndoStack( myPanelColorBar.getColorMap() );
    myPanelColorBar.getColorMap().setName(myTextName.getText().trim());
    myPanelColorBar.updateColors( colors, positions );
    invalidate();
    myIsUpdating = false;
    return true;
  }
  private void cancel() {
    dispose();
  }
  private void removeBundle( Bundle bundle ) {
    if( myIsUpdating ) return;
    myIsUpdating = true;
    myBundleList.remove( bundle );
    updateEditorPanel();
    apply();
    myIsUpdating = false;
  }
  private void addBundleAtEnd() {
    addBundle( null );
  }
  private void addBundle( Bundle bundle ) {
    int index = getBundleIndex( bundle );
    if( index <= 0 ) {
      addBundle( bundle, false );
    }
    else {
      addBundle( bundle, true );
    }
  }
  private void addBundle( Bundle bundle, boolean interpolate ) {
    if( myIsUpdating ) return;
    myIsUpdating = true;
    int index = getBundleIndex( bundle );
    if( bundle == null ) {
      index = myBundleList.size()-1;
    }
    double position = myBundleList.get(index).getPosition();
    Color color = myBundleList.get(index).getColor();
    if( index > 0 ) {
      position = 0.5 * ( position + myBundleList.get(index-1).getPosition() );
      if( bundle != null ) {
        color = myPanelColorBar.getColorMap().getColor( (float)position );
      }
    }
    else if( index == 0 ) {
      position = 0.5 * ( position + myBundleList.get(index+1).getPosition() );
      index = 1;
    }
    csColorMapGenerator.Bundle newBundle = new csColorMapGenerator.Bundle( color, position );
    myBundleList.add(index,newBundle);

    updateEditorPanel();
    apply();
    myIsUpdating = false;
  }
  private int getBundleIndex( Bundle bundle ) {
    if( bundle == null ) return -1;
    for( int ibundle = 0; ibundle < myBundleList.size(); ibundle++ ) {
      Bundle bundleCurrent = myBundleList.get(ibundle);
      if( bundle.equals(bundleCurrent) ) {
        return ibundle;
      }
    }
    return -1;
  }
  private void updateEditorPanel() {
    myPanelEditor.removeAll();
    int yp = 0;
    for( int ibundle = 0; ibundle < myBundleList.size(); ibundle++ ) {
      Bundle bundle = myBundleList.get(ibundle);
      myPanelEditor.add( bundle.buttonAddAbove, new GridBagConstraints(
          0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 10, 0, 5 ), 0, 0 ) );
      myPanelEditor.add( Box.createHorizontalGlue(), new GridBagConstraints(
          1, yp++, 4, 1, 1.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      myPanelEditor.add( bundle.buttonRemove, new GridBagConstraints(
          0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 10, 0, 5 ), 0, 0 ) );
      myPanelEditor.add( bundle.buttonColor, new GridBagConstraints(
          1, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
      myPanelEditor.add( bundle.buttonSelect, new GridBagConstraints(
          2, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
      myPanelEditor.add( bundle.textPosition, new GridBagConstraints(
          3, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
      myPanelEditor.add( bundle.buttonInc, new GridBagConstraints(
          4, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 5, 0, 0 ), 0, 0 ) );
      myPanelEditor.add( bundle.buttonDec, new GridBagConstraints(
          5, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    }

    myPanelEditor.add( myButtonAddColor, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 5 ), 0, 0 ) );
    myPanelEditor.add( Box.createHorizontalGlue(), new GridBagConstraints(
        1, yp++, 5, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    myPanelEditor.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, 6, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 5, 0, 5 ), 0, 0 ) );

    myPanelEditor.revalidate();
    invalidate();
  }
  private void setBundlePosition( Bundle bundle, double position ) {
    if( myIsUpdating ) return;
    myIsUpdating = true;
    int index = getBundleIndex( bundle );
    if( index < 0 ) return;
    int numBundles = myBundleList.size();
    if( index == 0 || index == numBundles-1 ) {
      myIsUpdating = false;
      return;  // no change allowed for first and last bundle
    }
    double positionSmaller = myBundleList.get(index+1).getPosition();
    double positionBigger  = myBundleList.get(index-1).getPosition();
    if( position > positionBigger ) {
      position = positionBigger;
    }
    else if( position < positionSmaller ) {
      position = positionSmaller;
    }
    bundle.setPosition(position);
    apply();
    myIsUpdating = false;
  }

  private void selectBundle( Bundle bundle ) {
    mySelectedBundle = bundle;
    int index = getBundleIndex( bundle );
    for( int ibundle = 0; ibundle < myBundleList.size(); ibundle++ ) {
      if( ibundle != index ) {
        Bundle bundleCurrent = myBundleList.get(ibundle);
        bundleCurrent.setSelected( false );
      }
    }
    if( index < 0 ) return;
    mySelectedBundle.setSelected( true );
  }
  
  @Override
  public void stateChanged(ChangeEvent e) {
    if( mySelectedBundle != null ) {
      Color color = myColorChooser.getColor();
      mySelectedBundle.setColor( color );
    }
  }
  private void fireApplyEvent() {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).applyColorMap( myPanelColorBar.getColorMap() );
    }
  }
  private void fireUpdateMapsEvent() {
    myUpdatesWaiting = false;
    ArrayList<csCustomColorMap> list = new ArrayList<csCustomColorMap>();
    for( int i = 0; i < myManagerCustomMaps.getNumMaps(); i++ ) {
      list.add( myManagerCustomMaps.getMap( i ) );
    }
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).updateColorMaps( list );
    }
  }

  //----------------------------------------------------------------
  //
  //
  class Bundle implements csITimerTaskIncrementListener, csITimerTaskDecrementListener,
          csColorChangeListener {
    protected JButton buttonAddAbove;
    protected JButton buttonRemove;
    protected csColorButton buttonColor;
    protected JTextField textPosition;
    protected JButton buttonInc;
    protected JButton buttonDec;
    protected JRadioButton buttonSelect;
    private Color myColor;
    private double myPosition;
    private csTimerTaskIncrement myTaskInc;
    private csTimerTaskDecrement myTaskDec;

    public Bundle( Color color, double position ) {
      myColor    = color;
      myPosition = position;

      textPosition = new JTextField( myFormat0000.format(myPosition) );
      int preferredHeight = textPosition.getPreferredSize().height;
      textPosition.setPreferredSize( new Dimension(80,preferredHeight) );
      buttonColor = new csColorButton( csColorMapGenerator.this, color, 0, false );
      buttonColor.addColorChangeListener( this );
      
      buttonSelect = new JRadioButton("");
      buttonAddAbove = new JButton(cseis.resources.csResources.getIcon("icon_plus.png"));
      buttonAddAbove.setText("");
      buttonRemove = new JButton(cseis.resources.csResources.getIcon("icon_minus.png"));
      buttonRemove.setText("");
      buttonInc = new JButton(cseis.resources.csResources.getIcon("icon_plus.png"));
      buttonInc.setText("");
      buttonDec = new JButton(cseis.resources.csResources.getIcon("icon_minus.png"));
      buttonDec.setText("");
      buttonAddAbove.setMargin(new Insets(0,0,1,1));
      buttonRemove.setMargin(new Insets(0,0,1,1));
      buttonInc.setMargin(new Insets(0,0,1,1));
      buttonDec.setMargin(new Insets(0,0,1,1));
      buttonAddAbove.setPreferredSize(new Dimension(preferredHeight,preferredHeight));
      buttonRemove.setPreferredSize(new Dimension(preferredHeight,preferredHeight));
      buttonInc.setPreferredSize(new Dimension(preferredHeight,preferredHeight));
      buttonDec.setPreferredSize(new Dimension(preferredHeight,preferredHeight));

      buttonColor.setToolTipText("Select color");
      buttonAddAbove.setToolTipText("Add color");
      buttonRemove.setToolTipText("Remove color");
      buttonInc.setToolTipText("Increment position");
      buttonDec.setToolTipText("Decrement position");
      buttonSelect.setToolTipText("Select to edit this color");
      textPosition.setToolTipText("Color 'position' in the range [0,1]");
      
      buttonInc.addMouseListener( new MouseAdapter() {
        @Override
        public void mousePressed( MouseEvent e ) {
          increment();
          myTaskInc = new csTimerTaskIncrement( Bundle.this );
          Timer timer = new Timer();
          timer.scheduleAtFixedRate( myTaskInc, 500, 100 );
        }
        @Override
        public void mouseReleased( MouseEvent e ) {
          myTaskInc.cancel();
          myTaskInc = null;
        }
      });
      buttonDec.addMouseListener( new MouseAdapter() {
        @Override
        public void mousePressed( MouseEvent e ) {
          decrement();
          myTaskDec = new csTimerTaskDecrement( Bundle.this );
          Timer timer = new Timer();
          timer.scheduleAtFixedRate( myTaskDec, 500, 100 );
        }
        @Override
        public void mouseReleased( MouseEvent e ) {
          myTaskDec.cancel();
          myTaskDec = null;
        }
      });
      buttonSelect.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          selectBundle( Bundle.this );
        }
      });
      buttonAddAbove.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          addBundle( Bundle.this );
        }
      });
      buttonColor.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          selectBundle( Bundle.this );
          myColorChooser.setColor( myColor );
        }
      });
      textPosition.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          updatePosition();
        }
      });
      textPosition.addFocusListener( new FocusListener() {
        @Override
        public void focusGained(FocusEvent e) {
        }
        @Override
        public void focusLost(FocusEvent e) {
          updatePosition();
        }
      });
      buttonRemove.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          removeBundle( Bundle.this );
        }
      });
    }
    public Color getColor() {
      return myColor;
    }
    public void setColor( Color color ) {
      myColor = color;
      buttonColor.setColor(color);
      apply();
    }
    public void updatePosition() {
      try {
        double position = Double.parseDouble( textPosition.getText() );
        setBundlePosition( Bundle.this, position );
      }
      catch( NumberFormatException e ) {
        textPosition.setText(myFormat0000.format(myPosition));
      }
    }
    public void setSelected( boolean set ) {
      buttonSelect.setSelected( set );
    }
    public double getPosition() {
      return myPosition;
    }
    @Override
    public void increment() {
      double position = getPosition();
      position += POSITION_INCREMENT;
      setBundlePosition(this,position);
    }
    @Override
    public void decrement() {
      double position = getPosition();
      position -= POSITION_INCREMENT;
      setBundlePosition(this,position);
    }
    public void setPosition( double position ) {
      myPosition = position;
      textPosition.setText( myFormat0000.format(myPosition) );
    }

    @Override
    public void colorChanged(Object obj, Color color) {
      myColor = color;
      apply();
    }
  }
  //------------------------------------------------------------------------
  //
  public static void main( String [] args ) {
    JFrame frame = new JFrame();
    ArrayList<csCustomColorMap> cmapList = new ArrayList<csCustomColorMap>(csColorMap.NUM_DEFAULT_MAPS);
    for( int imap = 0; imap < csColorMap.NUM_DEFAULT_MAPS; imap++ ) {
      csColorMap cmap = new csColorMap(imap,csColorMap.COLOR_MAP_TYPE_32BIT);
      cmap.setMinMax(0.0, 1.0);
      cmapList.add( new csCustomColorMap(cmap) );
    }
    csColorMapGenerator g = new csColorMapGenerator( frame, cmapList, null ) {
      @Override
      public void setVisible( boolean set ) {
        super.setVisible( set );
        if( !set ) System.exit(0);
      }
    };
    g.setVisible(true);
  }
}

