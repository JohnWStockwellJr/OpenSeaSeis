/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.tools;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import cseis.general.csStandard;

/**
 * Progress bar window.<br>
 * Shows progress while loading seismic traces
 * @author 2007 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csProgressBarWindow extends JWindow {
  private JProgressBar myBar;
  private JProgressBar myBar2;
  private JFrame myParentFrame;
  private csIStoppable myStoppable;
  private JLabel myLabel;

  public csProgressBarWindow( JFrame parentFrame, int numTraces, csIStoppable stoppable ) {
    this( parentFrame, numTraces, stoppable, -1 );
  }
  public csProgressBarWindow( JFrame parentFrame, int numTraces, csIStoppable stoppable, int numTotalTraces ) {
    setAlwaysOnTop(true);
    super.setModalExclusionType( Dialog.ModalExclusionType.TOOLKIT_EXCLUDE );
    myParentFrame = parentFrame;
    myStoppable = stoppable;

    myBar = new JProgressBar();
    myBar.setBorderPainted(true);
    myBar.setValue(0);
    if( numTraces > 0 ) {
      myBar.setStringPainted(true);
      myBar.setIndeterminate(false);
      myBar.setMaximum(numTraces);
    }
    else {
      myBar.setStringPainted(false);
      myBar.setIndeterminate(true);
    }

    if( numTotalTraces >= 0 ) {
      myBar2 = new JProgressBar();
      myBar2.setValue(0);
      if( numTraces > 0 ) {
        myBar2.setStringPainted(true);
        myBar2.setIndeterminate(false);
        myBar2.setMaximum(numTotalTraces);
      }
      else {
        myBar2.setStringPainted(false);
        myBar2.setIndeterminate(true);
      }
      myBar2.setBackground( Color.ORANGE );
      myBar2.setForeground( Color.GREEN );
      myBar2.setStringPainted( true );
    }

    //  ----------------------------
    //  Button panel
    JButton buttonStop = new JButton("Stop");
    JPanel panelButtons = new JPanel( new GridBagLayout() );
    int xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.5, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( buttonStop, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.5, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    buttonStop.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        if( myStoppable != null ) myStoppable.stopOperation();
      }
    });

    JPanel panelUpper = new JPanel( new BorderLayout() );
    panelUpper.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createEtchedBorder(Color.white, Color.black), csStandard.DIALOG_BORDER ) );
    myLabel = new JLabel( "Reading in traces...", JLabel.CENTER );
    panelUpper.add( myLabel, BorderLayout.NORTH );
    panelUpper.add( myBar, BorderLayout.CENTER );
    if( numTotalTraces >= 0 ) {
      panelUpper.add( myBar2, BorderLayout.SOUTH );
    }

    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder(
        BorderFactory.createCompoundBorder(
        BorderFactory.createCompoundBorder( BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder() ),
        csStandard.DIALOG_BORDER ) );

    panelAll.add(panelUpper,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);
    getContentPane().add( panelAll, BorderLayout.CENTER );
    pack();
    setLocationRelativeTo( parentFrame );
  }
  public void setTitle( String text ) {
    myLabel.setText( text );
    myLabel.repaint();
  }
  public void setValue( int value ) {
    myBar.setValue(value);
    setLocationRelativeTo( myParentFrame );
  }
  public void setValues( int value, int value2 ) {
    myBar2.setValue(value2);
    myBar2.setString( "" + value2 + "/" + myBar2.getMaximum() );
    setValue( value );
  }
}


