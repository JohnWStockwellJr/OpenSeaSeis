/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.general.csStandard;
import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * SeaView preferences dialog. Set some specific default parameters for the application.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csPreferencesDialog extends JDialog {
  private SeaView mySeaView;
  
  private JTextField myTextNumTraces;
  private JCheckBox  myBoxShowFilename;
  private JCheckBox  myBoxShowActivePaneFeedback;
  private JTextField myTextMaxTabWidth;

  private JButton myButtonApply;
  private JButton myButtonClose;

  public csPreferencesDialog( SeaView seaview, csSeaViewProperties properties ) {
    super( seaview, "SeaView preferences", true );
    mySeaView = seaview;

    myTextMaxTabWidth = new JTextField( "" + properties.maxTabWidth );
    myTextMaxTabWidth.setToolTipText("Maximum width of tab label [in pixels]");
    myTextNumTraces  = new JTextField( "" + properties.numTraces );
    myTextNumTraces.setToolTipText("Number of traces to read in at once for new input files");
    myBoxShowFilename = new JCheckBox( "Show filename as label on seismic data", properties.showFilename );
    myBoxShowFilename.setToolTipText("Plot filename on top of seismic data");
    myBoxShowActivePaneFeedback = new JCheckBox( "Show graphical feedback", properties.showActivePaneFeedback );
    myBoxShowActivePaneFeedback.setToolTipText("<html>Highlight pane when selected as active<br><i>This may be slow when running SeaView on a remote computer</i></html>");
    
    JPanel panelSet = new JPanel(new GridBagLayout());
    panelSet.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Set preferences"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    int xp = 0;
    int yp = 0;
    panelSet.add( new JLabel("Number of traces to read in:"), new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 0 ), 0, 0 ) );
    panelSet.add( myTextNumTraces, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelSet.add( new JLabel("Maximum tab width:"), new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 0 ), 0, 0 ) );
    panelSet.add( myTextMaxTabWidth, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelSet.add( myBoxShowFilename, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelSet.add( myBoxShowActivePaneFeedback, new GridBagConstraints(
        0, yp, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

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

    panelAll.add(panelSet,BorderLayout.CENTER);
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
    try {
      int numTraces = Integer.parseInt( myTextNumTraces.getText() );
      mySeaView.setDefaultNumTraces( numTraces );
    }
    catch( NumberFormatException exc ) {
    }
    try {
      int maxTabWidth = Integer.parseInt( myTextMaxTabWidth.getText() );
      mySeaView.setMaxTabWidth( maxTabWidth );
    }
    catch( NumberFormatException exc ) {
    }
    boolean showFilename = myBoxShowFilename.isSelected();
    mySeaView.setShowFilename( showFilename );
    boolean showActivePaneFeedback = myBoxShowActivePaneFeedback.isSelected();
    mySeaView.setShowActivePaneFeedback( showActivePaneFeedback );
    return true;
  }
}

