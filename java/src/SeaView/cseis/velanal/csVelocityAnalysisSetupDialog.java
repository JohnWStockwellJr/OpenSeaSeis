/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.velanal;

import cseis.general.csStandard;
import cseis.seaview.csSeisPaneBundle;
import cseis.seis.csHeaderDef;
import cseis.velocity.csVelEnsembleInfo;
import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

@SuppressWarnings("serial")
public class csVelocityAnalysisSetupDialog extends JDialog {
  private final JComboBox<csSeisPaneBundle>  myComboBundle;
  private JComboBox<csHeaderDef> myComboHdrVel;
  private JComboBox<csHeaderDef> myComboHdrEns1;
  private JComboBox<csHeaderDef> myComboHdrEns2;

  private final JCheckBox myBoxHdrEns1;
  private final JCheckBox myBoxHdrEns2;
  private final JButton myButtonOK;
  private final JButton myButtonCancel;
  
  private csIVelSetupDialogListener myListener;

  public csVelocityAnalysisSetupDialog( JFrame parentFrame, csSeisPaneBundle[] bundleList, csSeisPaneBundle semblanceBundle, csIVelSetupDialogListener listener ) {
    super( parentFrame, "Velocity analysis setup", true );
    myListener = listener;

    myComboBundle = new JComboBox<csSeisPaneBundle>( bundleList );
    if( semblanceBundle != null ) {
      myComboBundle.setSelectedItem( semblanceBundle );
    }
    else {
      myComboBundle.setSelectedIndex( bundleList.length-1 );
    }
    
    myComboHdrVel  = new JComboBox<csHeaderDef>();
    myComboHdrEns1 = new JComboBox<csHeaderDef>();
    myComboHdrEns2 = new JComboBox<csHeaderDef>();
    myBoxHdrEns1 = new JCheckBox("Ensemble header 1", true);
    myBoxHdrEns2 = new JCheckBox("Ensemble header 2", false);
    JLabel labelHdrVel = new JLabel("Velocity header: ");

    myBoxHdrEns2.setEnabled(true);
    myComboHdrEns2.setEnabled(false);
    
    if( semblanceBundle != null ) updateTraceHeaders( semblanceBundle );

    labelHdrVel.setToolTipText("<html>Select trace header containing RMS velocity<br><i>...this should be VEL_RMS if generated with module $SEMBLANCE</i></html>");
    myComboHdrVel.setToolTipText("<html>Select trace header containing RMS velocity<br><i>...this should be VEL_RMS if generated with module $SEMBLANCE</i></html>");
    myBoxHdrEns1.setToolTipText("<html>Select trace header identifying image gather<br><i>Examples: CMP(2D), COL(3D)</i></html>");
    myBoxHdrEns2.setToolTipText("<html>Select second trace header identifying 3D image gather<br><i>Example: ROW(</i></html>");
    myComboHdrEns1.setToolTipText("<html>Select trace header identifying image gather<br><i>Examples: CMP(2D), COL(3D)</i></html>");
    myComboHdrEns2.setToolTipText("<html>Select second trace header identifying 3D image gather<br><i>Example: ROW(3D)</i></html>");
        
    JPanel panelSemblance = new JPanel(new GridBagLayout());
    panelSemblance.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Semblance setup"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelCMP = new JPanel(new GridBagLayout());
    panelCMP.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("CMP gather setup"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelStack = new JPanel(new GridBagLayout());
    panelStack.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Stack section setup"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    int xp = 0;
    int yp = 0;
    panelSemblance.add( new JLabel("Semblance data set: "), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 5, 0, 0 ), 0, 0 ) );
    panelSemblance.add( myComboBundle, new GridBagConstraints(
        xp++, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 0 ), 0, 0 ) );
    xp = 0;
    panelSemblance.add( labelHdrVel, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 10, 5, 0, 0 ), 0, 0 ) );
    panelSemblance.add( myComboHdrVel, new GridBagConstraints(
        xp++, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 10, 10, 0, 0 ), 0, 0 ) );
    xp = 0;
    panelSemblance.add( myBoxHdrEns1, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 10, 5, 0, 0 ), 0, 0 ) );
    panelSemblance.add( myComboHdrEns1, new GridBagConstraints(
        xp++, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 10, 10, 0, 0 ), 0, 0 ) );
    xp = 0;
    panelSemblance.add( myBoxHdrEns2, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 10, 5, 0, 0 ), 0, 0 ) );
    panelSemblance.add( myComboHdrEns2, new GridBagConstraints(
        xp++, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 10, 10, 0, 0 ), 0, 0 ) );

    panelSemblance.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, xp, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JLabel cmplabel = new JLabel("...to be implemented");
    JLabel stacklabel = new JLabel("...to be implemented");
    int fontsize = stacklabel.getFont().getSize();
    cmplabel.setFont( new Font(Font.SANS_SERIF,Font.ITALIC,fontsize));
    cmplabel.setEnabled(false);
    stacklabel.setFont( new Font(Font.SANS_SERIF,Font.ITALIC,fontsize));
    stacklabel.setEnabled(false);
    panelCMP.setEnabled(false);
    panelStack.setEnabled(false);
    
    panelCMP.add( cmplabel, new GridBagConstraints(
        0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 10, 5, 0, 0 ), 0, 0 ) );
    
    panelStack.add( stacklabel, new GridBagConstraints(
        0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 10, 5, 0, 0 ), 0, 0 ) );
    
    myButtonOK = new JButton("OK");
    myButtonCancel = new JButton("Cancel");
    this.getRootPane().setDefaultButton(myButtonOK);
    
    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonOK, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonCancel, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelSetup = new JPanel(new GridBagLayout());
    yp = 0;
    panelSetup.add( panelSemblance, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSetup.add( panelCMP, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSetup.add( panelStack, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    
    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );

    yp = 0;

    panelAll.add(panelSetup,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);
    getContentPane().add(panelAll);

    pack();
    setLocationRelativeTo( parentFrame );

    myButtonOK.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
        cancel();
      }
    });
    myButtonCancel.addActionListener(new ActionListener() {
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

    myBoxHdrEns1.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        boolean isSelected = myBoxHdrEns1.isSelected();
        myComboHdrEns1.setEnabled( isSelected );
        myBoxHdrEns2.setEnabled( isSelected );
        if( !isSelected ) {
          myBoxHdrEns2.setSelected(false);
          myComboHdrEns2.setEnabled(false);
        }
      }
    });
    myBoxHdrEns2.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        myComboHdrEns2.setEnabled( myBoxHdrEns2.isSelected() );
      }
    });
    myComboBundle.addItemListener( new ItemListener() {
      @Override
      public void itemStateChanged( ItemEvent e ) {
        if( e.getStateChange() == ItemEvent.SELECTED && myComboBundle.getSelectedIndex() >= 0 ) {
          csSeisPaneBundle bundle = (csSeisPaneBundle)myComboBundle.getSelectedItem();
          updateTraceHeaders( bundle );
        }
      }
    });
  }
  public void addVelDialogListener( csIVelSetupDialogListener listener ) {
    myListener = listener;
  }
  private void cancel() {
    dispose();
  }
  private boolean apply() {
    csSeisPaneBundle bundle = (csSeisPaneBundle)myComboBundle.getSelectedItem();
    csVelEnsembleInfo info; // Note: Cannot set header index here since header names were sorted by name, not by position in input file
    String nameVel   = ((csHeaderDef)myComboHdrVel.getSelectedItem()).name;
    if( myBoxHdrEns1.isSelected() ) {
      String nameEns1  = ((csHeaderDef)myComboHdrEns1.getSelectedItem()).name;
      if( myBoxHdrEns2.isSelected() ) {
        String nameEns2  = ((csHeaderDef)myComboHdrEns2.getSelectedItem()).name;
        info = new csVelEnsembleInfo( nameVel, nameEns1, nameEns2 );
      }
      else {
        info = new csVelEnsembleInfo( nameVel, nameEns1 );
      }
    }
    else {
      info = new csVelEnsembleInfo( nameVel );
    }
    if( myListener != null ) myListener.setupVelocityAnalysis( bundle, info );
    return true;
  }  
  public void updateTraceHeaders( csSeisPaneBundle bundle ) {
    csHeaderDef[] headerDef = bundle.getSortedTraceHeaderDef();
    if( headerDef != null ) {
      myComboHdrVel.setModel( new DefaultComboBoxModel<csHeaderDef>(headerDef) );
      myComboHdrEns1.setModel( new DefaultComboBoxModel<csHeaderDef>(headerDef) );
      myComboHdrEns2.setModel( new DefaultComboBoxModel<csHeaderDef>(headerDef) );
      boolean velFound = false;
      boolean ens1Found = false;
      boolean ens2Found = false;
      myComboHdrVel.setSelectedIndex( 0 );
      myComboHdrEns1.setSelectedIndex( 0 );
      myComboHdrEns2.setSelectedIndex( 0 );
      for( int i = 0; i < headerDef.length; i++ ) {
        String name = headerDef[i].name;
        if( !velFound && name.startsWith("vel") ) {
          myComboHdrVel.setSelectedIndex(i);
          velFound = true;
        }
        else if( !ens1Found && name.equals("cmp") ) {
          myComboHdrEns1.setSelectedIndex(i);
          ens1Found = true;
        }
        else if( !ens2Found && name.equals("col") ) {
          myComboHdrEns2.setSelectedIndex(i);
          ens2Found = true;
        }
      }
    }    
  }
}
