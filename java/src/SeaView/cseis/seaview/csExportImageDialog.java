/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.*;
import cseis.general.csStandard;

/**
 * Dialog window providing functionality to save seismic view to external image file.<br>
 * <br>
 * This class is currently unused since the functionality to create an image with arbitrary pixels resolution
 * is not implemented yet in class csSeisView and related classes.
 * 
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csExportImageDialog extends JDialog {
  private JTextField myTextPixelWidth;
  private JTextField myTextPixelHeight;
  private JCheckBox  myBoxKeepRatio;
  private double myRatio;

  private JButton myButtonApply;
  private JButton myButtonClose;
  
  private csIExportImageListener myListener;

  public csExportImageDialog( JFrame parentFrame, int width, int height, csIExportImageListener listener ) {
    super( parentFrame, "Export image file", true );

    myRatio = (double)width / (double)height;
    myTextPixelWidth  = new JTextField( "" + width );
    myTextPixelHeight = new JTextField( "" + height );
    myBoxKeepRatio = new JCheckBox("Keep ratio");

    JPanel panelSet = new JPanel(new GridBagLayout());
    panelSet.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Set image resolution"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    int xp = 0;
    int yp = 0;
    panelSet.add( new JLabel("Image size (pixels):"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 0 ), 0, 0 ) );
    panelSet.add( myTextPixelHeight, new GridBagConstraints(
        xp++, yp, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelSet.add( new JLabel(" x ", JLabel.CENTER), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSet.add( myTextPixelWidth, new GridBagConstraints(
        xp++, yp, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    yp++;
    panelSet.add( myBoxKeepRatio, new GridBagConstraints(
        0, yp, xp, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 10, 0, 10 ), 0, 0 ) );

    myButtonApply = new JButton("Export as..");
    myButtonClose = new JButton("Close");
    this.getRootPane().setDefaultButton(myButtonApply);
    
    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonApply, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

//   Nuts and bolts
    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );

    yp = 0;

    panelAll.add(panelSet,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);
    getContentPane().add(panelAll);

    pack();
    setLocationRelativeTo( parentFrame );

    myTextPixelHeight.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if( myBoxKeepRatio.isSelected() ) {
          try {
            int height = Integer.parseInt( myTextPixelHeight.getText() );
            int width = (int)Math.round( (double)height * myRatio );
            myTextPixelWidth.setText("" + width);
          } catch( NumberFormatException exc ) { }
        }
      }
    });
    myTextPixelWidth.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if( myBoxKeepRatio.isSelected() ) {
          try {
            int width = Integer.parseInt( myTextPixelWidth.getText() );
            int height = (int)Math.round( (double)width / myRatio );
            myTextPixelHeight.setText("" + height);
          } catch( NumberFormatException exc ) { }
        }
      }
    });
    myButtonApply.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if( apply() ) cancel();
      }
    });
    myButtonClose.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        cancel();
      }
    });
    this.addWindowListener(new WindowAdapter() {
      @Override
      public void windowClosing(WindowEvent e) {
        cancel();
      }
    });
  }
  private void cancel() {
    dispose();
  }
  private boolean apply() {
    try {
      int width = Integer.parseInt( myTextPixelWidth.getText() );
      int height = Integer.parseInt( myTextPixelHeight.getText() );
      myListener.exportImage( width, height );
    }
    catch( NumberFormatException exc ) {
      return false;
    }
    return true;
  }
  public static void main( String[] args ) {
//    csCombineDataDialog dialog = new csCombineDataDialog(null);
//    dialog.setVisible(true);
  }
}


