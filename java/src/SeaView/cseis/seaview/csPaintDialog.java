/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import java.awt.*;
import javax.swing.*;
import cseis.general.csStandard;
import cseis.swing.csColorButton;
import cseis.swing.csColorChangeListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Dialog window providing painting functionality in seismic view.<br>
 */
@SuppressWarnings("serial")
public class csPaintDialog extends JDialog {
  final static public Integer[] SIZES_LINES  = { 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  final static public Integer[] SIZES_POINTS = { 3, 5, 7, 9, 11, 13, 15, 17 };

  final private csColorButton myColorButtonLine;
  final private csColorButton myColorButtonPoint;
  final private JComboBox<Integer> myComboLineSize;
  final private JComboBox<Integer> myComboPointSize;
  final private JButton myButtonClose;
  final private JButton myButtonClear;

  final private csIPaintDialogListener myListener;
  
  public csPaintDialog( JFrame frame, csIPaintDialogListener listener ) {
    super(frame,"Paint dialog",false);
    super.setResizable(false);
    myListener = listener;

    myColorButtonPoint = new csColorButton( this, Color.blue );
    myColorButtonLine  = new csColorButton( this, Color.red );
    myButtonClose = new JButton("Close");
    myButtonClear = new JButton("Clear");
    myComboLineSize = new JComboBox<Integer>( SIZES_LINES );
    myComboPointSize  = new JComboBox<Integer>( SIZES_POINTS );
    myComboLineSize.setSelectedIndex( 4 );
    myComboPointSize.setSelectedIndex( 4 );
    
    JPanel panelParam = new JPanel(new GridBagLayout());
    panelParam.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Paint attributes"),
        csStandard.INNER_EMPTY_BORDER ) );

    int xp = 0;
    int yp = 0;
    panelParam.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelParam.add( new JLabel("Size"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.8, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelParam.add( new JLabel("Color"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    xp = 0;
    yp++;
    panelParam.add( new JLabel("Lines:"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelParam.add( myComboLineSize, new GridBagConstraints(
        xp++, yp, 1, 1, 0.8, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelParam.add( myColorButtonLine, new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    xp = 0;
    yp++;
    panelParam.add( new JLabel("Points:"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelParam.add( myComboPointSize, new GridBagConstraints(
        xp++, yp, 1, 1, 0.8, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelParam.add( myColorButtonPoint, new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 0, 0, 0 ), 0, 0 ) );
    xp = 0;
    yp++;
    panelParam.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, yp, 3, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    getRootPane().setDefaultButton(myButtonClose);

    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( myButtonClear, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 20, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.9, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );
    panelAll.add(panelParam,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);

    myColorButtonPoint.addColorChangeListener( new csColorChangeListener() {
      @Override
      public void colorChanged(Object obj, Color color) {
        apply();
      }
    });
    myColorButtonLine.addColorChangeListener( new csColorChangeListener() {
      @Override
      public void colorChanged(Object obj, Color color) {
        apply();
      }
    });
    myComboLineSize.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    myComboPointSize.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    
    myButtonClose.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        if( myListener != null ) myListener.closePaintDialog();
        dispose();
      }
    });
    myButtonClear.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        if( myListener != null ) myListener.clearPaintOverlay();
      }
    });
    
    getContentPane().add(panelAll);
    pack();
    setLocationRelativeTo(frame);
  }
  public void apply() {
    Color colorLine  = myColorButtonLine.getColor();
    Color colorPoint = myColorButtonPoint.getColor();
    int sizeLine  = (Integer)myComboLineSize.getSelectedItem();
    int sizePoint = (Integer)myComboPointSize.getSelectedItem();
    myListener.updatePaintOverlay( sizeLine, colorLine, sizePoint, colorPoint );
  }
}


