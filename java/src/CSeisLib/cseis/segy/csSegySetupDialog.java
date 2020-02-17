/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.segy;

import cseis.general.csStandard;
import cseis.general.csAbstractPreference;
import cseis.jni.csNativeSegyReader;
import cseis.jni.csSegyTrcHeaderDefinition;
import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * Dialog window where user can set SEGY attributes.<br>
 * SEGY attributes define how SEGY/SU files are read in, i.e. endian format, automatic coordinate scaling,
 * and the trace header 'map'. The trace header map defines byte locations etc in the SEGY/SU trace header.
 * @author 2013 Felipe Punto
 */
// TODO: Add functionality so that users can define their own custom SEGY trace header maps.
@SuppressWarnings("serial")
public class csSegySetupDialog extends JDialog {

  private JRadioButton myButtonSEGYDataBig;
  private JRadioButton myButtonSEGYDataLittle;
  private JRadioButton myButtonSEGYHdrBig;
  private JRadioButton myButtonSEGYHdrLittle;
  private JRadioButton myButtonSUDataBig;
  private JRadioButton myButtonSUDataLittle;
  private JRadioButton myButtonSUHdrBig;
  private JRadioButton myButtonSUHdrLittle;

  private JCheckBox myButtonSEGYScaleHdrs;
  private JRadioButton myButtonPresetMap;
  private JRadioButton myButtonCustomMap;
  private JComboBox<String> myComboPresetMap;
//  private JComboBox<String> myComboCustomMap;
  private JButton myButtonShowPresetMap;
  private JButton myButtonShowCustomMap;
  
  private JButton myButtonClose;
  private JButton myButtonDefault;
  
  private csSegyAttr mySegyAttr;
  private JFrame myParentFrame;
  private ArrayList<csISegyAttrListener> myListeners;

  public csSegySetupDialog( JFrame parentFrame, csSegyAttr attr ) {
    super( parentFrame, "SEGY/SU setup and analysis tool" );
    setModal(false);
    myParentFrame = parentFrame;
    mySegyAttr = new csSegyAttr( attr );
    myListeners = new ArrayList<csISegyAttrListener>();
    
    myButtonSEGYDataBig    = new JRadioButton("");
    myButtonSEGYDataLittle = new JRadioButton("");
    myButtonSEGYHdrBig     = new JRadioButton("");
    myButtonSEGYHdrLittle  = new JRadioButton("");
    myButtonSUDataBig      = new JRadioButton("");
    myButtonSUDataLittle   = new JRadioButton("");
    myButtonSUHdrBig       = new JRadioButton("");
    myButtonSUHdrLittle    = new JRadioButton("");

    ButtonGroup g1 = new ButtonGroup();
    ButtonGroup g2 = new ButtonGroup();
    ButtonGroup g3 = new ButtonGroup();
    ButtonGroup g4 = new ButtonGroup();
    g1.add(myButtonSEGYDataBig);
    g1.add(myButtonSEGYDataLittle);
    g2.add(myButtonSUDataBig);
    g2.add(myButtonSUDataLittle);
    g3.add(myButtonSEGYHdrBig);
    g3.add(myButtonSEGYHdrLittle);
    g4.add(myButtonSUHdrBig);
    g4.add(myButtonSUHdrLittle);
    
    myButtonSEGYScaleHdrs = new JCheckBox("Auto-scale coord/elev headers",attr.autoScaleCoord);

    myComboPresetMap = new JComboBox<String>( csNativeSegyReader.NAME_DEFAULT_MAPS );
//    myComboCustomMap  = new JComboBox<String>();
    myButtonPresetMap = new JRadioButton("Preset:");
    myButtonCustomMap  = new JRadioButton("Custom:");
    myButtonShowPresetMap = new JButton("Show");
    myButtonShowCustomMap  = new JButton("Show");
    
    ButtonGroup g5 = new ButtonGroup();
    g5.add(myButtonPresetMap);
    g5.add(myButtonCustomMap);

    myButtonClose   = new JButton("Close");
    myButtonDefault = new JButton("Default");
    
    myButtonShowPresetMap.setToolTipText("Show SEG-Y trace header mappings");
    myButtonSEGYScaleHdrs.setToolTipText("Automatically scale coordinate and elevation headers according to SEG-Y format description");
    myButtonDefault.setToolTipText("Reset to defaults");
    myButtonClose.setToolTipText("Close this dialog");
    
    JPanel panelSEGYEndian = new JPanel( new GridBagLayout() );
    panelSEGYEndian.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("SEGY endian format"),
        csStandard.INNER_EMPTY_BORDER ) );
    int xp = 0;
    int yp = 0;
    panelSEGYEndian.add( Box.createHorizontalGlue(), new GridBagConstraints(
      xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSEGYEndian.add( new JLabel("Big endian"), new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    panelSEGYEndian.add( new JLabel("Little endian"), new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    xp = 0;
    yp++;
    panelSEGYEndian.add( new JLabel("Header"), new GridBagConstraints(
      xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSEGYEndian.add( myButtonSEGYHdrBig, new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    panelSEGYEndian.add( myButtonSEGYHdrLittle, new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    xp = 0;
    yp++;
    panelSEGYEndian.add( new JLabel("Data"), new GridBagConstraints(
      xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSEGYEndian.add( myButtonSEGYDataBig, new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    panelSEGYEndian.add( myButtonSEGYDataLittle, new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    
    JPanel panelSUEndian = new JPanel( new GridBagLayout() );
    panelSUEndian.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("SU endian format"),
        csStandard.INNER_EMPTY_BORDER ) );
    xp = 0;
    yp = 0;
    panelSUEndian.add( Box.createHorizontalGlue(), new GridBagConstraints(
      xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSUEndian.add( new JLabel("Big endian"), new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    panelSUEndian.add( new JLabel("Little endian"), new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    xp = 0;
    yp++;
    panelSUEndian.add( new JLabel("Header"), new GridBagConstraints(
      xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSUEndian.add( myButtonSUHdrBig, new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    panelSUEndian.add( myButtonSUHdrLittle, new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    xp = 0;
    yp++;
    panelSUEndian.add( new JLabel("Data"), new GridBagConstraints(
      xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSUEndian.add( myButtonSUDataBig, new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    panelSUEndian.add( myButtonSUDataLittle, new GridBagConstraints(
      xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );

    JPanel panelParam = new JPanel( new GridBagLayout() );
    panelParam.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("SEGY/SU parameter setup"),
        csStandard.INNER_EMPTY_BORDER ) );

    xp = 0;
    yp = 0;
    panelParam.add( panelSEGYEndian, new GridBagConstraints(
        xp, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelParam.add( panelSUEndian, new GridBagConstraints(
        xp, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelParam.add( myButtonSEGYScaleHdrs, new GridBagConstraints(
        xp, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelParam.add( Box.createVerticalGlue(), new GridBagConstraints(
        xp, yp++, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    

    JPanel panelMap = new JPanel( new GridBagLayout() );
    panelMap.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Trace header mapping"),
        csStandard.INNER_EMPTY_BORDER ) );
    xp = 0;
    yp = 0;
    panelMap.add( myButtonPresetMap, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 5, 0, 0, 5 ), 0, 0 ) );
    panelMap.add( myComboPresetMap, new GridBagConstraints(
        xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 0, 0, 0 ), 0, 0 ) );
    panelMap.add( myButtonShowPresetMap, new GridBagConstraints(
        xp, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 5, 5, 0, 0 ), 0, 0 ) );
    xp = 0;
//    panelMap.add( myButtonCustomMap, new GridBagConstraints(
//        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
//        GridBagConstraints.NONE, new Insets( 5, 0, 0, 5 ), 0, 0 ) );
//    panelMap.add( myComboCustomMap, new GridBagConstraints(
//        xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
//        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
//    panelMap.add( myButtonShowCustomMap, new GridBagConstraints(
//        xp, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
//        GridBagConstraints.NONE, new Insets( 5, 0, 0, 0 ), 0, 0 ) );
     
    JPanel panelTop = new JPanel( new GridBagLayout() );

    xp = 0;
    yp = 0;
    panelTop.add( panelParam, new GridBagConstraints(
        xp, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTop.add( panelMap, new GridBagConstraints(
        xp, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTop.add( Box.createVerticalGlue(), new GridBagConstraints(
        xp, yp++, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( myButtonDefault, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.8, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    
    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );

    panelAll.add(panelTop,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);
    getContentPane().add(panelAll);
    pack();
    
    setButtons( mySegyAttr );

    myButtonDefault.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        setButtons( new csSegyAttr() );
        applyButtons(mySegyAttr);
      }
    });
    myButtonClose.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        applyButtons(mySegyAttr);
        setVisible(false);
        dispose();
      }
    });
    myButtonShowPresetMap.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        csSegySetupDialog.showBuiltInSegyTrcHdrMaps( myParentFrame, myComboPresetMap.getSelectedIndex() );
      }
    });
    myButtonShowCustomMap.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        showCustomSegyTrcHdrMaps();
      }
    });
    
    myButtonSEGYHdrBig.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        applyButtons(mySegyAttr);
      }
    });
    myButtonSEGYDataBig.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        applyButtons(mySegyAttr);
      }
    });
    myButtonSUHdrBig.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        applyButtons(mySegyAttr);
      }
    });
    myButtonSUDataBig.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        applyButtons(mySegyAttr);
      }
    });
    myButtonSEGYScaleHdrs.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        applyButtons(mySegyAttr);
      }
    });
  }
  public void addSegyAttrListener( csISegyAttrListener listener ) {
    myListeners.add( listener );
  }
  private void fireUpdateEvent() {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).updateSegyAttr( mySegyAttr );
    }
  }
  private void setButtons( csSegyAttr attr ) {
    myButtonSEGYDataLittle.setSelected( attr.endianSEGYData == java.nio.ByteOrder.LITTLE_ENDIAN );
    myButtonSEGYDataBig.setSelected( attr.endianSEGYData == java.nio.ByteOrder.BIG_ENDIAN );
    myButtonSEGYHdrLittle.setSelected( attr.endianSEGYHdr == java.nio.ByteOrder.LITTLE_ENDIAN );
    myButtonSEGYHdrBig.setSelected( attr.endianSEGYHdr == java.nio.ByteOrder.BIG_ENDIAN );

    myButtonSUDataLittle.setSelected( attr.endianSUData == java.nio.ByteOrder.LITTLE_ENDIAN );
    myButtonSUDataBig.setSelected( attr.endianSUData == java.nio.ByteOrder.BIG_ENDIAN );
    myButtonSUHdrLittle.setSelected( attr.endianSUHdr == java.nio.ByteOrder.LITTLE_ENDIAN );
    myButtonSUHdrBig.setSelected( attr.endianSUHdr == java.nio.ByteOrder.BIG_ENDIAN );

    myButtonSEGYScaleHdrs.setSelected(attr.autoScaleCoord);
    myButtonPresetMap.setSelected(!attr.isCustomMap);
    myButtonCustomMap.setSelected(attr.isCustomMap);

    int mapIndex = 0;
    while( mapIndex < csNativeSegyReader.NUM_DEFAULT_MAPS ) {
      if( csNativeSegyReader.DEFAULT_MAPS[mapIndex] == attr.hdrMap ) {
        break;
      }
      mapIndex += 1;
    }
    if( mapIndex < 0 || mapIndex >= csNativeSegyReader.NUM_DEFAULT_MAPS ) {
      mapIndex = 0;
    }
    myComboPresetMap.setSelectedIndex(mapIndex);
  }

  private void applyButtons( csSegyAttr attr ) {
    attr.endianSEGYData = myButtonSEGYDataLittle.isSelected() ? java.nio.ByteOrder.LITTLE_ENDIAN : java.nio.ByteOrder.BIG_ENDIAN;
    attr.endianSEGYHdr  = myButtonSEGYHdrLittle.isSelected() ? java.nio.ByteOrder.LITTLE_ENDIAN : java.nio.ByteOrder.BIG_ENDIAN;
    attr.endianSUData = myButtonSUDataLittle.isSelected() ? java.nio.ByteOrder.LITTLE_ENDIAN : java.nio.ByteOrder.BIG_ENDIAN;
    attr.endianSUHdr  = myButtonSUHdrLittle.isSelected() ? java.nio.ByteOrder.LITTLE_ENDIAN : java.nio.ByteOrder.BIG_ENDIAN;

    attr.autoScaleCoord = myButtonSEGYScaleHdrs.isSelected();

    attr.isCustomMap = myButtonCustomMap.isSelected();
    int mapIndex = myComboPresetMap.getSelectedIndex();
    if( mapIndex < 0 || mapIndex >= csNativeSegyReader.NUM_DEFAULT_MAPS ) {
      mapIndex = 0;
      myComboPresetMap.setSelectedIndex(mapIndex);
    }
    attr.hdrMap = csNativeSegyReader.DEFAULT_MAPS[mapIndex];
    fireUpdateEvent();
  }
  public csSegyAttr getSegyAttr() {
    return( new csSegyAttr(mySegyAttr) );
  }
  public static void showBuiltInSegyTrcHdrMaps( JFrame parentFrame, int selectedIndex ) {
    ArrayList<csSegyTrcHeaderDefinition> hdrMaps = new ArrayList<csSegyTrcHeaderDefinition>();

    for( int i = 0; i < csNativeSegyReader.NUM_DEFAULT_MAPS; i++ ) {
      csSegyTrcHeaderDefinition hdrMap = new csSegyTrcHeaderDefinition( csNativeSegyReader.NAME_DEFAULT_MAPS[i] );
      csNativeSegyReader.trcHdrMap( csNativeSegyReader.DEFAULT_MAPS[i], hdrMap );
      hdrMaps.add( hdrMap );
    }
    JDialog dialog = new csSegyHeaderView( parentFrame, "SEG-Y trace header mapping", hdrMaps, selectedIndex );
    dialog.setVisible(true);
  }
  public void showCustomSegyTrcHdrMaps() {
    
  }
  //----------------------------------------------------------
  public static void main( String[] args ) {
    try {
      csAbstractPreference.setCseisDirectory();
    }
    catch( Exception e ) {
      
    }
    String path = System.getProperty( "java.library.path" );
    String libName = System.mapLibraryName( "csJNIlib" );
    try {
      System.load( path + "/" + libName );
    }
    catch( java.lang.UnsatisfiedLinkError e ) {
      JOptionPane.showMessageDialog( null,
          e.toString() + "\n" +
          "java.library.path = " + System.getProperty( "java.library.path" ) + "\n" +
          " - SegySetupDialog will not open.", "Error",
          JOptionPane.ERROR_MESSAGE );
      System.exit( -1 );
    }

    JFrame frame = new JFrame("New frame");
    frame.setVisible(true);
    csSegySetupDialog d = new csSegySetupDialog( frame, new csSegyAttr() );
    d.setVisible(true);
  }
}

