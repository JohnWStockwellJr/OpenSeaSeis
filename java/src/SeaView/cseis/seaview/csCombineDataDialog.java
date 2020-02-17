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
import cseis.seis.csISeismicReader;

/**
 * Dialog window providing functionality to combine two data sets.<br>
 * <br>
 * The user selects two data sets (=seismic bundles) and a mathematical operation which shall 
 * be applied when combining the two data sets.<br>
 * The new data set (=seismic reader object) is created by a method defined in csSeismicBundle
 * and subsequently sent to SeaView for data retrieval and display in a new seismic bundle object.
 * 
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csCombineDataDialog extends JDialog {
  private SeaView mySeaView;
  
  private JComboBox<csSeisPaneBundle>  myComboBundle1;
  private JComboBox<csSeisPaneBundle>  myComboBundle2;
  private JComboBox<String>  myComboOperation;
  private JTextField myTextMultiplier;
  private JTextField myTextAddValue;

  private JButton myButtonApply;
  private JButton myButtonClose;

  public csCombineDataDialog( SeaView seaview, csSeisPaneBundle[] bundleList, csSeisPaneBundle selectBundle ) {
    super( seaview, "Combine data sets", true );  // true: Make dialog modal

    mySeaView = seaview;

    myComboBundle1 = new JComboBox<csSeisPaneBundle>( bundleList );
    myComboBundle2 = new JComboBox<csSeisPaneBundle>( bundleList );
    myComboOperation = new JComboBox<String>( csSeisPaneBundle.TEXT_OPERATION );
    if( selectBundle != null ) {
      myComboBundle1.setSelectedItem( selectBundle );
    }
    else {
      myComboBundle1.setSelectedIndex( bundleList.length-1 );
    }
    int index2 = myComboBundle1.getSelectedIndex() - 1;
    if( index2 < 0 ) index2 = Math.min( bundleList.length-1, index2+2 );
    myComboBundle2.setSelectedIndex(index2);
    
    myTextMultiplier = new JTextField( "1.0" );
    myTextAddValue   = new JTextField( "0.0" );

    JPanel panelOperation = new JPanel(new GridBagLayout());
    panelOperation.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Operation"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    int xp = 0;
    int yp = 0;
    panelOperation.add( myComboBundle1, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 0 ), 0, 0 ) );
    panelOperation.add( myComboOperation, new GridBagConstraints(
        xp++, yp, 1, 1, 0.33, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 10 ), 0, 0 ) );
    panelOperation.add( myComboBundle2, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelOperation.add( new JLabel("*",JLabel.CENTER), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 10 ), 0, 0 ) );
    panelOperation.add( myTextMultiplier, new GridBagConstraints(
        xp++, yp, 1, 1, 0.33, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelOperation.add( new JLabel("+",JLabel.CENTER), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    panelOperation.add( myTextAddValue, new GridBagConstraints(
        xp++, yp, 1, 1, 0.33, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    yp += 1;
    panelOperation.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, xp, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    myButtonApply = new JButton("Apply");
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

    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );

    yp = 0;

    panelAll.add(panelOperation,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);
    getContentPane().add(panelAll);

    pack();
    setLocationRelativeTo( mySeaView );

    myButtonApply.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
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
    String text = "";
    try {
      text = myTextMultiplier.getText();
      float multiplier = Float.parseFloat( text );
      text = myTextAddValue.getText();
      float addValue = Float.parseFloat( text );
      int operation = myComboOperation.getSelectedIndex();
      csSeisPaneBundle bundle1 = (csSeisPaneBundle)myComboBundle1.getSelectedItem();
      csSeisPaneBundle bundle2 = (csSeisPaneBundle)myComboBundle2.getSelectedItem();
      String message = "";
      if( bundle1.getSampleInt() != bundle2.getSampleInt() ) {
        message += "Data sets have different sample intervals:\n" +
              bundle1.getSampleInt() + " != " + bundle2.getSampleInt() + "\n";
      }
      if( bundle1.getNumSamples() != bundle2.getNumSamples() ) {
        message += "Data sets have different number of samples:\n" +
              bundle1.getNumSamples() + " != " + bundle2.getNumSamples() + "\n";
      }
      if( bundle1.getDisplayedNumTraces() != bundle2.getDisplayedNumTraces() ) {
        message += "Data sets have different number of displayed traces:\n" +
              bundle1.getDisplayedNumTraces() + " != " + bundle2.getDisplayedNumTraces() + "\n";
      }
      if( message.length() > 0 ) {
        message += "Combine data sets anyway?";
        int option = JOptionPane.showConfirmDialog( this, message, "Confirm", JOptionPane.YES_NO_OPTION);
        if( option != JOptionPane.YES_OPTION ) return false;
      }
      csISeismicReader reader = bundle1.combineDataSets( bundle2, operation, multiplier, addValue );
      String operationText = csSeisPaneBundle.TEXT_OPERATION[operation];
      String title = bundle1.getTitle() + " " + operationText + " " + bundle2.getTitle() +
              " * " + multiplier + " + " + addValue;
      mySeaView.readData( reader, title, SeaView.FORMAT_CSEIS, true );
    }
    catch( NumberFormatException exc ) {
      return false;
    }
    catch( Exception exc ) {
    }
    return true;
  }
  public static void main( String[] args ) {
//    csCombineDataDialog dialog = new csCombineDataDialog(null);
//    dialog.setVisible(true);
  }
}


