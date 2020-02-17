/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.processing;

import cseis.general.csStandard;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.MouseAdapter;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;

/**
 * Generic processing dialog.<br>
 * This dialog provides a generic frame for setting up parameters for processing steps.
 * Each specific processing step provides the inner "main panel" with the parameters.<br>
 * The csIProcessingSetupListener interface is used to communicate to the listener when the parameter setup
 * is complete.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csProcessingDialog extends JDialog {
  private JPanel myMainPanel;
  private JButton myButtonApply;
  private JButton myButtonClose;
  private csIProcessing myProcessingStep;
  private csIProcessingSetupListener myListener;
  private JProgressBar myProgressBar;
  private Component myGlassPane;
  
  public csProcessingDialog( Window parentWindow ) {
    super( parentWindow, "Processing step" );
    setModal(true);
    
    myProgressBar = new JProgressBar();
    myProcessingStep = null;
    myListener = null;
    
    myMainPanel = new JPanel( new BorderLayout() );
    myMainPanel.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Parameters"),
        csStandard.INNER_EMPTY_BORDER ) );

    myButtonApply = new JButton("Apply");
    myButtonClose = new JButton("Close");
    this.getRootPane().setDefaultButton(myButtonApply);
    
    JPanel panelButtons = new JPanel( new GridBagLayout() );
    int xp = 0;
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

    JPanel panelBottom = new JPanel(new BorderLayout());
    panelBottom.add( new JLabel("Progress: "), BorderLayout.WEST );
    panelBottom.add( myProgressBar, BorderLayout.CENTER );
    panelBottom.add( panelButtons, BorderLayout.SOUTH );
    
    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );

    
    panelAll.add(myMainPanel,BorderLayout.CENTER);
    panelAll.add(panelBottom,BorderLayout.SOUTH);
    getContentPane().add(panelAll);

    pack();
    setLocationRelativeTo( parentWindow );

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
    if( myProcessingStep != null ) {
      String errorMessage = myProcessingStep.retrieveParameters();
      if( errorMessage == null ) {
        myListener.setupComplete( myProcessingStep );
      }
      else {
        JOptionPane.showMessageDialog( this,
                "Error occurred when retrieving user input parameters:\n" +
                errorMessage,
                "Error",
                JOptionPane.ERROR_MESSAGE );
      }
    }
    return true;
  }
  public void startProcessing() {
    myGlassPane = getGlassPane();
    myGlassPane.setCursor( Cursor.getPredefinedCursor( Cursor.WAIT_CURSOR ) );
    myGlassPane.addKeyListener( new KeyAdapter() {} );
    myGlassPane.addMouseListener( new MouseAdapter() {} );
    myGlassPane.setVisible(true);
    myProgressBar.setIndeterminate(true);
  }
  public void stopProcessing() {
    myProgressBar.setIndeterminate(false);
    myGlassPane.setVisible(false);
    myGlassPane = null;
  }
  public void setProcessingStep( csIProcessing proc, csIProcessingSetupListener listener ) {
    myListener = listener;
    myProcessingStep = proc;
    setTitle( "Processing - " + proc.getName() );
    myMainPanel.removeAll();
    myMainPanel.add( proc.getParameterPanel() );
    invalidate();
    pack();
  }
  
}

