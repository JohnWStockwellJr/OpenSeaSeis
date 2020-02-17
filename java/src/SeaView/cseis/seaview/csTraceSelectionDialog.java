/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Dictionary;
import java.util.Enumeration;

import cseis.general.csStandard;
import cseis.seis.csHeaderDef;
import cseis.swing.csColorRadioButton;

/**
 * Trace selection dialog for seismic viewer.<br>
 * <br>
 * In this dialog, the user selects which traces shall be displayed. Three main options are available:<br>
 * <ul>
 *  <li> Trace selection: Simply displays a specified maximum number of consecutive traces.
 *  <li> Ensemble selection: Displays all <i>consecutive</i> traces with a common trace header value = one ensemble.
 *  <li> Header value selection: Displays all traces with a common trace header value
 *       <i>independent of their location in the data file</i>.
 * </ul>
 * @author 2007 Bjorn Olofsson
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csTraceSelectionDialog extends JDialog {
  private static final int TEXT_MIN_WIDTH = 95;
  private static final String TEXT_BUTTON_SCAN   = "Scan";
  private static final String TEXT_BUTTON_CANCEL = "Cancel";
  private static final String TOOLTIP_BUTTON_SCAN = "<html>Scan all traces for value of selected trace header<br>" +
          "<i>This enables full functionality for ensemble and header selection modes</i></html>";
  private static final String TOOLTIP_BUTTON_CANCEL = "Cancel scanning process";
  private static final Color COLOR_SCANNED = Color.green;
  private static final Color COLOR_NOT_SCANNED = Color.red;
  
  private csSeisPaneBundle myCurrentBundle;
  
  private JTextField myTextNumTraces;
  private JTextField myTextFirstTrace;
  private JTextField myTextStepTrace;
  private JSlider    mySliderFirstTrace;
  
  private JTextField myTextNumEns;
  private JTextField myTextFirstEns;
  private JSlider    mySliderFirstEns;

  private JTextField myTextHeaderValue;
  private JTextField myTextHeaderNumTraces;
  private JSlider    mySliderHeaderValue;
  
  private JRadioButton myButtonDisplayTraces;
  private JRadioButton myButtonDisplayEns;
  private JRadioButton myButtonDisplayHeader;

  private JTextField myTextScanHdrName;
  private JComboBox<csHeaderDef> myComboScanHdr;
  private JButton myButtonScan;
  private JProgressBar myProgressScan;
  private csColorRadioButton myButtonScanFlag;
  
  private JButton myButtonApply;
  private JButton myButtonClose;
  private JButton myButtonOK;
  private ArrayList<csITraceSelectionListener> myListeners;

  private csScanHeaderInfo myScanHeaderInfo;
  private boolean myIsAdjustingSliderHeaderValue = false;
  
  public csTraceSelectionDialog( JFrame parentFrame, csSeisPaneBundle bundle ) {
    super( parentFrame, "Trace selection" );

    myScanHeaderInfo = null;
    csTraceSelectionParam param = null;
    if( bundle != null ) {
      param = bundle.getTraceSelectionParam();
    }
    if( param == null ) param = new csTraceSelectionParam();
    
    myCurrentBundle = null;
    myListeners = new ArrayList<csITraceSelectionListener>();

    myTextFirstTrace = new JTextField( "" + (param.firstTraceIndex+1) );
    myTextNumTraces  = new JTextField( "" + param.numTraces );
    myTextHeaderNumTraces  = new JTextField( "" + param.numTraces );
    myTextStepTrace  = new JTextField( "" + param.traceStep );
    myTextFirstEns   = new JTextField( "" + (param.firstEnsIndex+1) );
    myTextNumEns     = new JTextField( "" + param.numEns );
    myTextScanHdrName = new JTextField( "" + param.selectedHdrName );
    myTextHeaderValue = new JTextField( "" + param.selectedHdrValue );
    JLabel labelFirstTrace = new JLabel("First trace:");
    JLabel labelNumTraces  = new JLabel("#:");
    JLabel labelStepTrace  = new JLabel("Inc.:");
    JLabel labelFirstEns   = new JLabel("First ensemble:");
    JLabel labelNumEns     = new JLabel("#:");
    JLabel labelHeaderValue = new JLabel("Header value:");
    JLabel labelNumHdr     = new JLabel("#:");

    myTextFirstTrace.setToolTipText("First trace to display, starting at 1");
    myTextNumTraces.setToolTipText("Number of traces to display at once");
    myTextStepTrace.setToolTipText("Trace increment, display every n'th trace");
    myTextFirstEns.setToolTipText("First ensemble to display, starting at 1");  
    myTextNumEns.setToolTipText("Number of ensembles to display at once");
    myTextHeaderValue.setToolTipText("Value of trace header to display (only show traces with this header value)");
    myTextHeaderNumTraces.setToolTipText("Number of traces to display at once");

    labelFirstTrace.setToolTipText("First trace to display, starting at 1");
    labelNumTraces.setToolTipText("Number of traces to display at once");
    labelStepTrace.setToolTipText("Trace increment, display every n'th trace");
    labelFirstEns.setToolTipText("First ensemble to display, starting at 1");  
    labelNumEns.setToolTipText("Number of ensembles to display at once");
    labelHeaderValue.setToolTipText("Value of trace header to display (only show traces with this header value)");
    labelNumHdr.setToolTipText("Number of traces to display at once");
    
    myButtonDisplayTraces  = new JRadioButton("");
    myButtonDisplayEns     = new JRadioButton("");
    myButtonDisplayHeader  = new JRadioButton("");
    ButtonGroup group1 = new ButtonGroup();
    group1.add(myButtonDisplayTraces);
    group1.add(myButtonDisplayEns);
    group1.add(myButtonDisplayHeader);
    myButtonDisplayTraces.setSelected( true );
    
    mySliderFirstTrace = new JSlider(JSlider.HORIZONTAL);
    mySliderFirstTrace.setSnapToTicks( false );
    mySliderFirstTrace.setMinimum( 0 );
    mySliderFirstTrace.setMaximum( 10 );
    mySliderFirstTrace.setPaintTicks(true);
    mySliderFirstTrace.setPaintLabels(true);
    mySliderFirstTrace.setValue( 0 );
    
    mySliderFirstEns = new JSlider(JSlider.HORIZONTAL);
    mySliderFirstEns.setSnapToTicks( false );
    mySliderFirstEns.setMinimum( 0 );
    mySliderFirstEns.setMaximum( 0 );
    mySliderFirstEns.setPaintTicks(true);
    mySliderFirstEns.setPaintLabels(true);
    mySliderFirstEns.setValue( 0 );
    Dictionary<Integer,JComponent> table = mySliderFirstEns.createStandardLabels(1);
    mySliderFirstEns.setLabelTable( table );

    mySliderHeaderValue = new JSlider(JSlider.HORIZONTAL);
    mySliderHeaderValue.setSnapToTicks( false );
    mySliderHeaderValue.setMinimum( 0 );
    mySliderHeaderValue.setMaximum( 0 );
    mySliderHeaderValue.setPaintTicks(true);
    mySliderHeaderValue.setPaintLabels(true);
    mySliderHeaderValue.setValue( 0 );
    table = mySliderHeaderValue.createStandardLabels(1);
    mySliderHeaderValue.setLabelTable( table );

    myComboScanHdr = new JComboBox<csHeaderDef>();
    myButtonScan = new JButton( TEXT_BUTTON_SCAN );
    myButtonScan.setEnabled( false );
    myProgressScan = new JProgressBar();
    myButtonScanFlag = new csColorRadioButton( COLOR_NOT_SCANNED );
    myButtonScanFlag.setSelected(true);
    myButtonScanFlag.setToolTipText("Green if selected trace header has been scanned, red otherwise");
    myButtonScan.setToolTipText( TOOLTIP_BUTTON_SCAN );
    myProgressScan.setToolTipText("Scan progress");
    
    int height = myTextNumTraces.getPreferredSize().height;
    myTextFirstTrace.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextNumTraces.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextStepTrace.setPreferredSize( new Dimension( TEXT_MIN_WIDTH/2, height ) );
    myTextFirstEns.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextNumEns.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextHeaderValue.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextHeaderNumTraces.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );

    mySliderFirstTrace.setPreferredSize( new Dimension(TEXT_MIN_WIDTH, (int)(height*2.5) ) );
    mySliderFirstEns.setPreferredSize( new Dimension(TEXT_MIN_WIDTH, (int)(height*2.5)) );
    mySliderHeaderValue.setPreferredSize( new Dimension(TEXT_MIN_WIDTH, (int)(height*2.5)) );

    myButtonOK    = new JButton("OK");
    myButtonApply = new JButton("Apply");
    myButtonClose = new JButton("Close");

    if( bundle != null ) {
      myCurrentBundle = bundle;
//      updateBundle( bundle );
    }
    else {
      updateTraceSlider( 1 );
    }
    enableGUIElements();
    
    JPanel panelTrace = new JPanel(new GridBagLayout());
    panelTrace.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Trace selection"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JPanel panelEns = new JPanel(new GridBagLayout());
    panelEns.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Ensemble selection"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JPanel panelHdr = new JPanel(new GridBagLayout());
    panelHdr.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Header value selection"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    int yp = 0;
/*    panelTrace.add( Box.createHorizontalBox(), new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 0 ), 0, 0 ) );
    panelTrace.add( new JLabel("First",JLabel.CENTER), new GridBagConstraints(
        1, yp, 1, 1, 0.4, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTrace.add( new JLabel("Amount",JLabel.CENTER), new GridBagConstraints(
        2, yp, 1, 1, 0.4, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTrace.add( new JLabel("Increment",JLabel.CENTER), new GridBagConstraints(
        3, yp++, 1, 1, 0.2, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
*/

    JPanel panelTraceBottom  = new JPanel(new GridBagLayout());
    JPanel panelEnsBottom    = new JPanel(new GridBagLayout());
    JPanel panelHdrBottom = new JPanel(new GridBagLayout());

    yp = 0;
    int xp = 0;
    panelTraceBottom.add( labelFirstTrace, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelTraceBottom.add( myTextFirstTrace, new GridBagConstraints(
        xp++, yp, 1, 1, 0.4, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTraceBottom.add( labelNumTraces, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    panelTraceBottom.add( myTextNumTraces, new GridBagConstraints(
        xp++, yp, 1, 1, 0.4, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTraceBottom.add( labelStepTrace, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
    panelTraceBottom.add( myTextStepTrace, new GridBagConstraints(
        xp++, yp++, 1, 1, 0.2, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
//    panelTraceBottom.add( Box.createVerticalGlue(), new GridBagConstraints(
//        0, yp++, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
//        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    yp = 0;
    xp = 0;
    panelEnsBottom.add( labelFirstEns, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelEnsBottom.add( myTextFirstEns, new GridBagConstraints(
        xp++, yp, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelEnsBottom.add( labelNumEns, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 5 ), 0, 0 ) );
    panelEnsBottom.add( myTextNumEns, new GridBagConstraints(
        xp++, yp, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
//    panelEnsBottom.add( Box.createVerticalGlue(), new GridBagConstraints(
//        xp++, yp, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
//        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    xp = 0;
    yp = 0;
    panelHdrBottom.add( labelHeaderValue, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelHdrBottom.add( myTextHeaderValue, new GridBagConstraints(
        xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelHdrBottom.add( labelNumHdr, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 5 ), 0, 0 ) );
    panelHdrBottom.add( myTextHeaderNumTraces, new GridBagConstraints(
        xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
//    panelHdrBottom.add( Box.createVerticalGlue(), new GridBagConstraints(
//        0, yp++, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
//        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    
    yp = 0;
    panelTrace.add( myButtonDisplayTraces, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 5, 10, 0, 0 ), 0, 0 ) );
    panelTrace.add( mySliderFirstTrace, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 5, 10, 0, 0 ), 0, 0 ) );
    panelTrace.add( panelTraceBottom, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 10, 0, 0 ), 0, 0 ) );

    yp = 0;
    panelEns.add( myButtonDisplayEns, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 5, 10, 0, 0 ), 0, 0 ) );
    panelEns.add( mySliderFirstEns, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 5, 10, 0, 0 ), 0, 0 ) );
    panelEns.add( panelEnsBottom, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 10, 0, 0 ), 0, 0 ) );

    yp = 0;
    panelHdr.add( myButtonDisplayHeader, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 5, 10, 0, 0 ), 0, 0 ) );
    panelHdr.add( mySliderHeaderValue, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 5, 10, 0, 0 ), 0, 0 ) );
    panelHdr.add( panelHdrBottom, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 10, 0, 0 ), 0, 0 ) );


    JPanel panelScanLeft = new JPanel(new GridBagLayout());
    yp = 0;
    panelScanLeft.add( myComboScanHdr, new GridBagConstraints(
        0, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 10, 0, 5 ), 0, 0 ) );
    panelScanLeft.add( myTextScanHdrName, new GridBagConstraints(
        1, 0, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 10, 0, 5 ), 0, 0 ) );
    panelScanLeft.add( myButtonScanFlag, new GridBagConstraints(
        2, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 5, 10, 0, 5 ), 0, 0 ) );

    JPanel panelScanButtons = new JPanel( new GridBagLayout() );
    panelScanButtons.add( myButtonScan, new GridBagConstraints(
        0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 5, 10, 0, 5 ), 0, 0 ) );
    panelScanButtons.add( new JLabel("Progress:"), new GridBagConstraints(
        1, 1, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
        GridBagConstraints.NONE, new Insets( 5, 10, 0, 5 ), 0, 0 ) );
    panelScanButtons.add( myProgressScan, new GridBagConstraints(
        2, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 0, 0, 0 ), 0, 0 ) );
    
    JPanel panelScan = new JPanel( new BorderLayout() );
    panelScan.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Select trace header name/scan values"),
        csStandard.INNER_EMPTY_BORDER ) );

    panelScan.add( panelScanLeft, BorderLayout.CENTER );
    panelScan.add( panelScanButtons, BorderLayout.SOUTH );
    
    this.getRootPane().setDefaultButton(myButtonApply);

    // remove the binding for pressed
    getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW)
            .put(KeyStroke.getKeyStroke("ENTER"), "none");
    // retarget the binding for released
    getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW)
            .put(KeyStroke.getKeyStroke("released ENTER"), "press");

    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonOK, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonApply, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelTop = new JPanel( new GridBagLayout() );
    panelTop.add( panelTrace, new GridBagConstraints(
        0, 0, 1, 1, 1.0, 0.33, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTop.add( panelEns, new GridBagConstraints(
        0, 1, 1, 1, 1.0, 0.33, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTop.add( panelHdr, new GridBagConstraints(
        0, 2, 1, 1, 1.0, 0.33, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTop.add( panelScan, new GridBagConstraints(
        0, 3, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );

    yp = 0;

//    panelAll.add(panelTop,BorderLayout.NORTH);
    panelAll.add(panelTop,BorderLayout.CENTER);
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
    mySliderFirstTrace.addChangeListener( new ChangeListener() {
      @Override
      public void stateChanged( ChangeEvent e ) {
//        if( mySliderFirstTrace.getValueIsAdjusting() ) return;
        int value = mySliderFirstTrace.getValue();
        if( value < mySliderFirstTrace.getMinimum() || value > mySliderFirstTrace.getMaximum() ) return;
        myTextFirstTrace.setText( "" + value );
      }
    });
    mySliderHeaderValue.addChangeListener( new ChangeListener() {
      @Override
      public void stateChanged( ChangeEvent e ) {
//        if( mySliderHeaderValue.getValueIsAdjusting() ) return;
        if( myIsAdjustingSliderHeaderValue ) return;
        myIsAdjustingSliderHeaderValue = true;
        int sortedHdrIndex = mySliderHeaderValue.getValue();
        if( sortedHdrIndex < mySliderHeaderValue.getMinimum() || sortedHdrIndex > mySliderHeaderValue.getMaximum() || myCurrentBundle == null ) return;
        if( myScanHeaderInfo != null && sortedHdrIndex <= myScanHeaderInfo.getNumSortedHeaderValues() ) {
           Number number = myScanHeaderInfo.getSortedHeaderValue(sortedHdrIndex);
           myTextHeaderValue.setText( "" + number );
           mySliderHeaderValue.setToolTipText( "Trace #" + (sortedHdrIndex+1)+ "  Value " + number );
        }
        myIsAdjustingSliderHeaderValue = false;
      }
    });
    myTextHeaderValue.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        if( myCurrentBundle == null || myScanHeaderInfo == null ) return;
        boolean isScanned = ( myCurrentBundle != null && myScanHeaderInfo != null &&
            myCurrentBundle.isHeaderScanned( myScanHeaderInfo.getHeaderName() ) );
        if( !isScanned ) return;
        String value = myTextHeaderValue.getText();
        for( int i = 0; i < myScanHeaderInfo.getNumSortedHeaderValues(); i++ ) {
          if( value.compareTo(myScanHeaderInfo.getSortedHeaderValue(i).toString()) == 0 ) {
            mySliderHeaderValue.setValue(i);
            mySliderHeaderValue.setToolTipText( "Trace #" + (i+1) + "  Value " + value );
            break;
          }
        }
      }
    });
    mySliderFirstEns.addChangeListener( new ChangeListener() {
      @Override
      public void stateChanged( ChangeEvent e ) {
        int value = mySliderFirstEns.getValue();
        if( value < mySliderFirstEns.getMinimum() || value > mySliderFirstEns.getMaximum() || myCurrentBundle == null ) return;
        if( myScanHeaderInfo != null && value <= myScanHeaderInfo.getNumEnsembles() ) {
          Number number = myScanHeaderInfo.getEnsValue(value);
          mySliderFirstEns.setToolTipText( "Ensemble #" + (value+1) + "  Value " + number );
          myTextFirstEns.setText( "" + (value+1) ); // +1 to convert from index to ensemble number
        }
      }
    });
    myTextFirstEns.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        if( myCurrentBundle == null || myScanHeaderInfo == null ) return;
        String value = myTextFirstEns.getText();
        try {
          int valueInt = Integer.parseInt(value);
          if( valueInt < myScanHeaderInfo.getNumEnsembles() ) {
            mySliderFirstEns.setValue(valueInt-1);  //  // -1 to convert from trace number to trace index
            Number number = myScanHeaderInfo.getEnsValue(valueInt);
            mySliderFirstEns.setToolTipText( "Ensemble #" + (valueInt+1) + "  Value " + number );
          }
        } catch( NumberFormatException e2 ) {}
      }
    });
    
    myButtonDisplayTraces.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        enableGUIElements();
      }
    });
    myButtonDisplayEns.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        enableGUIElements();
      }
    });
    myButtonDisplayHeader.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        enableGUIElements();
      }
    });
    myComboScanHdr.addItemListener( new ItemListener() {
      @Override
      public void itemStateChanged( ItemEvent e ) {
        if( e.getStateChange() == ItemEvent.SELECTED && myComboScanHdr.getSelectedIndex() >= 0 ) {
          csHeaderDef headerDef = (csHeaderDef)myComboScanHdr.getSelectedItem();
//          String hdrName = myComboScanHdr.getSelectedItem().toString();
          myTextScanHdrName.setText( headerDef.name );
          setSelectedTraceHeader( headerDef );
        }
      }
    });
    myTextScanHdrName.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        String hdrName = myTextScanHdrName.getText();
        performScanComboAction( hdrName );
      }
    });
    myButtonScan.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if( myCurrentBundle == null ) return;
        updateScanAction( myButtonScan.getText().compareTo(TEXT_BUTTON_SCAN) == 0, true );
      }
    });
    this.addWindowListener(new WindowAdapter() {
      @Override
      public void windowClosing(WindowEvent e) {
        cancel();
      }
    });
  }
  @Override
  public void setVisible( boolean doSetVisible_in ) {
    boolean doSetVisible = doSetVisible_in;
    if( !isVisible() ) {
      updateBundle( myCurrentBundle, true );
    }
    super.setVisible( doSetVisible );
  }
  public void stopScan( boolean isSuccessful ) {
    updateScanAction( false, false );
    if( isSuccessful && myCurrentBundle != null && myCurrentBundle.getNumScannedHeaders() > 0 ) {
      myScanHeaderInfo = myCurrentBundle.getLatestScanHeaderInfo();
      myTextScanHdrName.setText( myScanHeaderInfo.getHeaderName() );
      if( performScanComboAction( myScanHeaderInfo.getHeaderName() ) ) {
        updateScanGUIComponents();
        enableGUIElements();
      }
    }
  }
  // Make sure header name exists in bundle data set before calling this function
  private void setSelectedTraceHeader( csHeaderDef headerDef ) {
    if( myCurrentBundle == null ) return;
    myScanHeaderInfo = myCurrentBundle.getScanHeaderInfo( headerDef.name );
    updateScanGUIComponents();
    enableGUIElements();
  }
  private void updateScanGUIComponents() {
    if( myCurrentBundle == null ) return;
    boolean isScanned = false;
    if( myScanHeaderInfo != null ) {
      isScanned = myCurrentBundle.isHeaderScanned( myScanHeaderInfo.getHeaderName() );
      int numHeaders = myScanHeaderInfo.getNumSortedHeaderValues();
      if( isScanned && numHeaders != 0 && myComboScanHdr.getSelectedIndex() >= 0 ) {
        myButtonScanFlag.changeColor( COLOR_SCANNED );
        updateHdrSlider( numHeaders );
        updateEnsSlider( myScanHeaderInfo.getNumEnsembles() );
      }
      else {
        myButtonScanFlag.changeColor( COLOR_NOT_SCANNED );
      }
    }
    else {
      myButtonScanFlag.changeColor( COLOR_NOT_SCANNED );
    }
    myButtonScan.setEnabled( myComboScanHdr.getSelectedIndex() >= 0 && myCurrentBundle.hasRandomFileAccess() );
    repaint();
  }
  
  private void updateScanAction( boolean startScan, boolean fireEvent ) {
    if( startScan ) {
      if( fireEvent ) fireStartScanEvent( (csHeaderDef)myComboScanHdr.getSelectedItem() );
      myButtonScan.setText( TEXT_BUTTON_CANCEL );
      myButtonScan.setToolTipText( TOOLTIP_BUTTON_CANCEL );
      myProgressScan.setMinimum( 0 );
      myProgressScan.setMaximum( myCurrentBundle.getTotalNumTraces() );
    }
    else {
      if( fireEvent ) fireCancelScanEvent();
      myButtonScan.setText( TEXT_BUTTON_SCAN );
      myButtonScan.setToolTipText( TOOLTIP_BUTTON_SCAN );
      myProgressScan.setValue( -1 );
    }
  }
  private void updateEnsSlider( int numEns ) {
    mySliderFirstEns.setMinimum( 0 );
    mySliderFirstEns.setMaximum( numEns-1 );
    int increment = Math.max( 1, numEns-1 );
    Dictionary<Integer,JComponent> table = mySliderFirstEns.createStandardLabels( increment );
    Enumeration<Integer> ee = table.keys();
    while( ee.hasMoreElements() ) {
      table.remove( ee.nextElement() );
    }
    table.put( 0, new JLabel("" + myScanHeaderInfo.getEnsValue(0) + " (#1)" ) );
    int num = (numEns-1)/2;
    table.put( num, new JLabel("" + myScanHeaderInfo.getEnsValue(num) +
            " (#" + (num+1) + ")" ) );
    table.put( numEns-1, new JLabel("" + myScanHeaderInfo.getEnsValue(numEns-1) +
            " (#" + numEns + ")" ) );
    mySliderFirstEns.setLabelTable( table );
  }
  private void updateHdrSlider( int numHeaders ) {
    mySliderHeaderValue.setMinimum( 0 );
    mySliderHeaderValue.setMaximum( numHeaders-1 );
    int increment = Math.max( 1, numHeaders-1 );
    Dictionary<Integer,JComponent> table = mySliderHeaderValue.createStandardLabels( increment );
    Enumeration<Integer> ee = table.keys();
    while( ee.hasMoreElements() ) {
      table.remove( ee.nextElement() );
    }
    table.put( 0, new JLabel("" + myScanHeaderInfo.getSortedHeaderValue(0) + " (#1)" ) );
    int num = (numHeaders-1)/2;
    table.put( num, new JLabel("" + myScanHeaderInfo.getSortedHeaderValue(num) +
            " (#" + (num+1) + ")" ) );
    table.put( numHeaders-1, new JLabel("" + myScanHeaderInfo.getSortedHeaderValue(numHeaders-1) +
            " (#" + numHeaders + ")" ) );
    mySliderHeaderValue.setLabelTable( table );
  }
  public void updateTraceSlider( int numTraces ) {
    if( numTraces == 0 ) numTraces = 1;
    mySliderFirstTrace.setMinimum( 1 );
    mySliderFirstTrace.setMaximum( numTraces );
    mySliderFirstTrace.setMajorTickSpacing( numTraces / 10 );
    if( numTraces > 1 ) {
      Dictionary<Integer,JComponent> table = mySliderFirstTrace.createStandardLabels(numTraces-1);
      table.put( numTraces/2, new JLabel(""+numTraces/2) );
      mySliderFirstTrace.setLabelTable( table );
    }
  }
  private void enableGUIElements() {
    boolean isTraceSelection   = myButtonDisplayTraces.isSelected();
    boolean isEnsSelection     = myButtonDisplayEns.isSelected();
    boolean isHeaderSelection = myButtonDisplayHeader.isSelected();
    boolean isScanned = ( myCurrentBundle != null && myScanHeaderInfo != null &&
            myCurrentBundle.isHeaderScanned( myScanHeaderInfo.getHeaderName() ) );

    if( myCurrentBundle != null ) {
      isTraceSelection  &= myCurrentBundle.hasRandomFileAccess();
      isEnsSelection    &= myCurrentBundle.hasRandomFileAccess();
      isHeaderSelection &= myCurrentBundle.hasRandomFileAccess();
    }

    myTextFirstTrace.setEnabled( isTraceSelection );
    myTextNumTraces.setEnabled( isTraceSelection );
    myTextStepTrace.setEnabled( isTraceSelection );
    mySliderFirstTrace.setEnabled( isTraceSelection );

    myTextFirstEns.setEnabled( isEnsSelection && isScanned );
    myTextNumEns.setEnabled( isEnsSelection );
    mySliderFirstEns.setEnabled( isEnsSelection && isScanned );

    myTextHeaderValue.setEnabled( isHeaderSelection );
    mySliderHeaderValue.setEnabled( isHeaderSelection && isScanned );
    myTextHeaderNumTraces.setEnabled( isHeaderSelection );

//    myButtonApply.setEnabled( myCurrentBundle != null && myCurrentBundle.hasRandomFileAccess() );
  }
  //----------------------------------------------------------------------------------------------
  //
  public void updateFirstTraceIndex( int firstTraceIndex ) {
    myTextFirstTrace.setText( "" + (firstTraceIndex+1) );
    mySliderFirstTrace.setValue( firstTraceIndex+1 );
  }
  public void updateFirstEnsIndex( int firstEnsIndex ) {
    myTextFirstEns.setText( "" + (firstEnsIndex+1) );
    myTextFirstEns.repaint();
    mySliderFirstEns.setValue( firstEnsIndex );
    mySliderFirstEns.repaint();
  }
  public void updateBundle( csSeisPaneBundle bundle ) {
    updateBundle( bundle, false );
  }
  public synchronized void updateBundle( csSeisPaneBundle bundle, boolean forceRefresh ) {
    if( bundle == null ) {
      myCurrentBundle = bundle;
      return;
    }
    csTraceSelectionParam paramNew = bundle.getTraceSelectionParam();
    csTraceSelectionParam paramPrev = null;
    csHeaderDef[] headerDefPrev = null;
    forceRefresh = forceRefresh || ( isVisible() && (!bundle.isSameAs(myCurrentBundle) ) || paramPrev == null );
    if( myCurrentBundle != null ) {
      paramPrev = myCurrentBundle.getTraceSelectionParam();
      headerDefPrev = myCurrentBundle.getSortedTraceHeaderDef();
    }
    myCurrentBundle = bundle;
    if( forceRefresh || paramPrev.firstTraceIndex != paramNew.firstTraceIndex ) updateFirstTraceIndex( paramNew.firstTraceIndex );
    if( forceRefresh || paramPrev.firstEnsIndex != paramNew.firstEnsIndex ) updateFirstEnsIndex( paramNew.firstEnsIndex );
    if( forceRefresh ) {
      updateTraceSlider( myCurrentBundle.getTotalNumTraces() );
      myTextNumTraces.setText( "" + paramNew.numTraces );
      myTextStepTrace.setText( "" + paramNew.traceStep );

      myTextFirstEns.setText( "" + (paramNew.firstEnsIndex+1) );
      myTextNumEns.setText( "" + paramNew.numEns );

      myTextHeaderNumTraces.setText( "" + paramNew.numTracesHdr );
      
      csHeaderDef[] headerDef = myCurrentBundle.getSortedTraceHeaderDef();
      if( headerDef != null ) {
        boolean renewCombo = ( forceRefresh || headerDefPrev.length != headerDef.length );
        if( !renewCombo ) {
          for( int i = 0; i < headerDef.length; i++ ) {
            if( headerDef[i].compareTo(headerDefPrev[i]) != 0 ) {
              renewCombo = true;
              break;
            }
          }
        }
        if( renewCombo ) {
          myComboScanHdr.setModel( new DefaultComboBoxModel<csHeaderDef>(headerDef) );
        }
        myComboScanHdr.setSelectedIndex( paramNew.selectedHdrIndex );
        myTextScanHdrName.setText( paramNew.selectedHdrName );
      }
      else {
        myComboScanHdr.removeAllItems();
      }
      if( paramNew.selectOption == csTraceSelectionParam.SELECT_TRACE ) myButtonDisplayTraces.setSelected(true);
      else if( paramNew.selectOption == csTraceSelectionParam.SELECT_ENS ) myButtonDisplayEns.setSelected(true);
      else if( paramNew.selectOption == csTraceSelectionParam.SELECT_HEADER ) myButtonDisplayHeader.setSelected(true);

      if( myCurrentBundle.isScanning() ) {
        updateScanAction( false, false );
      }
      enableGUIElements();
      updateScanGUIComponents();
    }
  }
  private void cancel() {
    dispose();
  }
  private boolean apply() {
    csTraceSelectionParam param = extractGUIEntries();
    if( param.selectOption != csTraceSelectionParam.SELECT_TRACE &&
        param.selectedHdrIndex < 0 ) {
      String text = (param.selectOption == csTraceSelectionParam.SELECT_ENS) ? "ensemble" : "header";
      JOptionPane.showMessageDialog( this, "No trace header selected for " + text +
              " selection.", "Error", JOptionPane.ERROR_MESSAGE );
      return false;
    }
    fireResetSelectionEvent( param );
    return true;
  }
  /**
   * Extract entries made in GUI and update mySelectParam
   */
  private csTraceSelectionParam extractGUIEntries() {
    csTraceSelectionParam param = new csTraceSelectionParam();
    if( myButtonDisplayTraces.isSelected() ) {
      param.selectOption = csTraceSelectionParam.SELECT_TRACE;
    }
    else if( myButtonDisplayEns.isSelected() ) {
      param.selectOption = csTraceSelectionParam.SELECT_ENS;
    }
    else if( myButtonDisplayHeader.isSelected() ) {
      param.selectOption = csTraceSelectionParam.SELECT_HEADER;
    }
    try {
      param.numTraces = Integer.parseInt( myTextNumTraces.getText() );
    } catch( NumberFormatException exc ) {}
    try {
      param.firstTraceIndex = Integer.parseInt( myTextFirstTrace.getText() ) - 1; // Convert from trace number to trace index
      if( myCurrentBundle != null && param.firstTraceIndex >= myCurrentBundle.getTotalNumTraces() ) {
        param.firstTraceIndex = myCurrentBundle.getTotalNumTraces()-1;
      }
    } catch( NumberFormatException exc ) {}
    try {
      param.traceStep = Integer.parseInt( myTextStepTrace.getText() );
    } catch( NumberFormatException exc ) {}
    try {
      param.numEns = Integer.parseInt( myTextNumEns.getText() );
    } catch( NumberFormatException exc ) {}
    try {
      param.firstEnsIndex = Integer.parseInt( myTextFirstEns.getText() ) - 1; // Convert from trace number to trace index
    } catch( NumberFormatException exc ) {}
    try {
      param.selectedHdrValue = Double.parseDouble( myTextHeaderValue.getText() );
    } catch( NumberFormatException exc2 ) {}
    try {
      param.numTracesHdr = Integer.parseInt( myTextHeaderNumTraces.getText() );
    } catch( NumberFormatException exc ) {}

    Object obj = myComboScanHdr.getSelectedItem();
    if( obj != null ) param.selectedHdrName = obj.toString();

    param.selectedHdrIndex = myComboScanHdr.getSelectedIndex();
    return param;
  }
  private boolean performScanComboAction( String hdrName ) {
    for( int i = 0; i < myComboScanHdr.getItemCount(); i++ ) {
      csHeaderDef headerDef = (csHeaderDef)myComboScanHdr.getItemAt(i);
      if( hdrName.compareTo( headerDef.name ) == 0 ) {
        myComboScanHdr.setSelectedIndex(i);
        setSelectedTraceHeader( headerDef );
        return true;
      }
    }
    myComboScanHdr.setSelectedIndex(-1);
    myButtonScan.setEnabled( false );
    return false;
  }
  public void updateScanProgress( int currentTraceIndex ) {
    myProgressScan.setValue( currentTraceIndex );
  }
  public void addTraceSelectionListener( csITraceSelectionListener listener ) {
    myListeners.add(listener);
  }
  private void fireResetSelectionEvent( csTraceSelectionParam param ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).updateTraceSelection( param );
    }
  }
  private void fireStartScanEvent( csHeaderDef headerDef ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).startScanTraceSelection( headerDef );
    }
  }
  private void fireCancelScanEvent() {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).cancelScanTraceSelection();
    }
  }
  public static void main( String[] args ) {
    csTraceSelectionDialog dialog = new csTraceSelectionDialog(null,null);
    dialog.setVisible(true);
  }
}

